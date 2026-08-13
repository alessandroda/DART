+ /bin/bash -x /tmp/tmp.4C5kmSMsm8
+ SCRIPT_PID=2346573
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
2026-06-21 19:27:39 INFO 
███    ███ ██ ███    ███ ███████ ███████ ██
████  ████ ██ ████  ████ ██      ██      ██
██ ████ ██ ██ ██ ████ ██ █████   ███████ ██
██  ██  ██ ██ ██  ██  ██ ██           ██ ██
██      ██ ██ ██      ██ ███████ ███████ ██



2026-06-21 19:27:39 INFO [PIPELINE] =======================================
2026-06-21 19:27:39 INFO [PIPELINE] Starting chimere–DART orchestrator
2026-06-21 19:27:39 INFO [PIPELINE] Config file: config/config_irene_IM.yaml
2026-06-21 19:27:39 INFO [PIPELINE] Run dir: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart
2026-06-21 19:27:39 INFO [PIPELINE] Log file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/mimesi_orchestrator_logs/chimere_DART_20260621_192739.log
2026-06-21 19:27:39 INFO [PIPELINE] =======================================
2026-06-21 19:27:39 INFO Running assimilation with model_type=ModelType.CHIMERE
2026-06-21 19:27:39 INFO Using scheduler=Scheduler.SLURM, queue=rome
2026-06-21 19:27:39 INFO [STEP] ---- TIME LOOP START ----
2026-06-21 19:27:39 INFO [TIME] step_start current_time=2020-02-14 15:00:00 simulated_time=None dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 19:27:39 INFO [TIME] window start=2020-02-14 15:00:00 end=2020-02-15 00:00:00 run_hours=9 has_assimilation=False
2026-06-21 19:27:39 INFO Creating directories and links for ENS1 to run chimere's parallel part
2026-06-21 19:27:39 INFO Linking EMIS ...
2026-06-21 19:27:39 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens01.nc
2026-06-21 19:27:40 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens01.nc
2026-06-21 19:27:40 INFO >> Checking links...
2026-06-21 19:27:42 INFO >> All links are good for ENS1  ...
2026-06-21 19:27:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:27:50 INFO Hourly dataset computed and listing created
2026-06-21 19:27:53 INFO Hourly dataset computed
2026-06-21 19:27:53 INFO Creating directories and links for ENS2 to run chimere's parallel part
2026-06-21 19:27:53 INFO Linking EMIS ...
2026-06-21 19:27:53 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens02.nc
2026-06-21 19:27:53 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens02.nc
2026-06-21 19:27:53 INFO >> Checking links...
2026-06-21 19:27:56 INFO >> All links are good for ENS2  ...
2026-06-21 19:27:56 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:27:57 INFO Hourly dataset computed and listing created
2026-06-21 19:28:00 INFO Hourly dataset computed
2026-06-21 19:28:00 INFO Creating directories and links for ENS3 to run chimere's parallel part
2026-06-21 19:28:00 INFO Linking EMIS ...
2026-06-21 19:28:00 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens03.nc
2026-06-21 19:28:00 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens03.nc
2026-06-21 19:28:00 INFO >> Checking links...
2026-06-21 19:28:03 INFO >> All links are good for ENS3  ...
2026-06-21 19:28:03 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:28:04 INFO Hourly dataset computed and listing created
2026-06-21 19:28:06 INFO Hourly dataset computed
2026-06-21 19:28:06 INFO Creating directories and links for ENS4 to run chimere's parallel part
2026-06-21 19:28:06 INFO Linking EMIS ...
2026-06-21 19:28:07 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens04.nc
2026-06-21 19:28:07 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens04.nc
2026-06-21 19:28:07 INFO >> Checking links...
2026-06-21 19:28:09 INFO >> All links are good for ENS4  ...
2026-06-21 19:28:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:28:11 INFO Hourly dataset computed and listing created
2026-06-21 19:28:13 INFO Hourly dataset computed
2026-06-21 19:28:13 INFO Creating directories and links for ENS5 to run chimere's parallel part
2026-06-21 19:28:13 INFO Linking EMIS ...
2026-06-21 19:28:14 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens05.nc
2026-06-21 19:28:14 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens05.nc
2026-06-21 19:28:14 INFO >> Checking links...
2026-06-21 19:28:16 INFO >> All links are good for ENS5  ...
2026-06-21 19:28:16 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:28:17 INFO Hourly dataset computed and listing created
2026-06-21 19:28:19 INFO Hourly dataset computed
2026-06-21 19:28:19 INFO Creating directories and links for ENS6 to run chimere's parallel part
2026-06-21 19:28:19 INFO Linking EMIS ...
2026-06-21 19:28:20 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens06.nc
2026-06-21 19:28:20 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens06.nc
2026-06-21 19:28:20 INFO >> Checking links...
2026-06-21 19:28:23 INFO >> All links are good for ENS6  ...
2026-06-21 19:28:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:28:24 INFO Hourly dataset computed and listing created
2026-06-21 19:28:26 INFO Hourly dataset computed
2026-06-21 19:28:26 INFO Creating directories and links for ENS7 to run chimere's parallel part
2026-06-21 19:28:26 INFO Linking EMIS ...
2026-06-21 19:28:26 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens07.nc
2026-06-21 19:28:27 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens07.nc
2026-06-21 19:28:27 INFO >> Checking links...
2026-06-21 19:28:29 INFO >> All links are good for ENS7  ...
2026-06-21 19:28:29 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:28:30 INFO Hourly dataset computed and listing created
2026-06-21 19:28:33 INFO Hourly dataset computed
2026-06-21 19:28:33 INFO Creating directories and links for ENS8 to run chimere's parallel part
2026-06-21 19:28:33 INFO Linking EMIS ...
2026-06-21 19:28:33 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens08.nc
2026-06-21 19:28:33 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens08.nc
2026-06-21 19:28:33 INFO >> Checking links...
2026-06-21 19:28:36 INFO >> All links are good for ENS8  ...
2026-06-21 19:28:36 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:28:37 INFO Hourly dataset computed and listing created
2026-06-21 19:28:39 INFO Hourly dataset computed
2026-06-21 19:28:39 INFO Creating directories and links for ENS9 to run chimere's parallel part
2026-06-21 19:28:39 INFO Linking EMIS ...
2026-06-21 19:28:40 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens09.nc
2026-06-21 19:28:40 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens09.nc
2026-06-21 19:28:40 INFO >> Checking links...
2026-06-21 19:28:42 INFO >> All links are good for ENS9  ...
2026-06-21 19:28:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:28:43 INFO Hourly dataset computed and listing created
2026-06-21 19:28:46 INFO Hourly dataset computed
2026-06-21 19:28:46 INFO Creating directories and links for ENS10 to run chimere's parallel part
2026-06-21 19:28:46 INFO Linking EMIS ...
2026-06-21 19:28:46 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens10.nc
2026-06-21 19:28:47 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens10.nc
2026-06-21 19:28:47 INFO >> Checking links...
2026-06-21 19:28:49 INFO >> All links are good for ENS10  ...
2026-06-21 19:28:49 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:28:50 INFO Hourly dataset computed and listing created
2026-06-21 19:28:52 INFO Hourly dataset computed
2026-06-21 19:28:52 INFO Creating directories and links for ENS11 to run chimere's parallel part
2026-06-21 19:28:52 INFO Linking EMIS ...
2026-06-21 19:28:53 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens11.nc
2026-06-21 19:28:53 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens11.nc
2026-06-21 19:28:53 INFO >> Checking links...
2026-06-21 19:28:55 INFO >> All links are good for ENS11  ...
2026-06-21 19:28:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:28:56 INFO Hourly dataset computed and listing created
2026-06-21 19:28:58 INFO Hourly dataset computed
2026-06-21 19:28:58 INFO Creating directories and links for ENS12 to run chimere's parallel part
2026-06-21 19:28:58 INFO Linking EMIS ...
2026-06-21 19:28:59 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens12.nc
2026-06-21 19:28:59 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens12.nc
2026-06-21 19:28:59 INFO >> Checking links...
2026-06-21 19:29:01 INFO >> All links are good for ENS12  ...
2026-06-21 19:29:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:29:03 INFO Hourly dataset computed and listing created
2026-06-21 19:29:05 INFO Hourly dataset computed
2026-06-21 19:29:05 INFO Creating directories and links for ENS13 to run chimere's parallel part
2026-06-21 19:29:05 INFO Linking EMIS ...
2026-06-21 19:29:05 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens13.nc
2026-06-21 19:29:05 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens13.nc
2026-06-21 19:29:05 INFO >> Checking links...
2026-06-21 19:29:08 INFO >> All links are good for ENS13  ...
2026-06-21 19:29:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:29:09 INFO Hourly dataset computed and listing created
2026-06-21 19:29:11 INFO Hourly dataset computed
2026-06-21 19:29:11 INFO Creating directories and links for ENS14 to run chimere's parallel part
2026-06-21 19:29:11 INFO Linking EMIS ...
2026-06-21 19:29:12 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens14.nc
2026-06-21 19:29:12 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens14.nc
2026-06-21 19:29:12 INFO >> Checking links...
2026-06-21 19:29:14 INFO >> All links are good for ENS14  ...
2026-06-21 19:29:14 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:29:16 INFO Hourly dataset computed and listing created
2026-06-21 19:29:18 INFO Hourly dataset computed
2026-06-21 19:29:18 INFO Creating directories and links for ENS15 to run chimere's parallel part
2026-06-21 19:29:18 INFO Linking EMIS ...
2026-06-21 19:29:19 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens15.nc
2026-06-21 19:29:19 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens15.nc
2026-06-21 19:29:19 INFO >> Checking links...
2026-06-21 19:29:21 INFO >> All links are good for ENS15  ...
2026-06-21 19:29:21 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:29:22 INFO Hourly dataset computed and listing created
2026-06-21 19:29:24 INFO Hourly dataset computed
2026-06-21 19:29:24 INFO ---------->>> Running CHIMERE model from 2020-02-14 15:00:00 to 2020-02-15 00:00:00
2026-06-21 19:29:24 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:29:24 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1
2026-06-21 19:29:24 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/end.2020021413_2_ENS1.nc
2026-06-21 19:29:24 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-06-21 19:29:24 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:29:24 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-06-21 19:29:24 INFO Queuing job for member 1...
2026-06-21 19:29:24 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:29:24 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-06-21 19:29:25 INFO Found: ['4937751']
2026-06-21 19:29:30 INFO [TGCC-IRENE] Submitted job with ID:['4937751']
2026-06-21 19:29:30 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:29:30 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2
2026-06-21 19:29:30 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/end.2020021413_2_ENS2.nc
2026-06-21 19:29:30 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-06-21 19:29:30 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:29:30 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-06-21 19:29:30 INFO Queuing job for member 2...
2026-06-21 19:29:30 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:29:30 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-06-21 19:29:31 INFO Found: ['4937752']
2026-06-21 19:29:36 INFO [TGCC-IRENE] Submitted job with ID:['4937752']
2026-06-21 19:29:36 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:29:36 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3
2026-06-21 19:29:36 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/end.2020021413_2_ENS3.nc
2026-06-21 19:29:36 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-06-21 19:29:36 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:29:36 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-06-21 19:29:36 INFO Queuing job for member 3...
2026-06-21 19:29:36 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:29:36 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-06-21 19:29:37 INFO Found: ['4937753']
2026-06-21 19:29:42 INFO [TGCC-IRENE] Submitted job with ID:['4937753']
2026-06-21 19:29:42 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:29:42 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4
2026-06-21 19:29:42 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/end.2020021413_2_ENS4.nc
2026-06-21 19:29:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-06-21 19:29:42 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:29:42 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-06-21 19:29:42 INFO Queuing job for member 4...
2026-06-21 19:29:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:29:42 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-06-21 19:29:43 INFO Found: ['4937754']
2026-06-21 19:29:48 INFO [TGCC-IRENE] Submitted job with ID:['4937754']
2026-06-21 19:29:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:29:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5
2026-06-21 19:29:48 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/end.2020021413_2_ENS5.nc
2026-06-21 19:29:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-06-21 19:29:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:29:48 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-06-21 19:29:48 INFO Queuing job for member 5...
2026-06-21 19:29:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:29:48 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-06-21 19:29:48 INFO Found: ['4937755']
2026-06-21 19:29:53 INFO [TGCC-IRENE] Submitted job with ID:['4937755']
2026-06-21 19:29:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:29:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6
2026-06-21 19:29:53 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/end.2020021413_2_ENS6.nc
2026-06-21 19:29:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-06-21 19:29:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:29:53 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-06-21 19:29:53 INFO Queuing job for member 6...
2026-06-21 19:29:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:29:53 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-06-21 19:29:54 INFO Found: ['4937756']
2026-06-21 19:29:59 INFO [TGCC-IRENE] Submitted job with ID:['4937756']
2026-06-21 19:29:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:29:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7
2026-06-21 19:29:59 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/end.2020021413_2_ENS7.nc
2026-06-21 19:29:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-06-21 19:29:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:29:59 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-06-21 19:29:59 INFO Queuing job for member 7...
2026-06-21 19:29:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:29:59 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-06-21 19:30:00 INFO Found: ['4937757']
2026-06-21 19:30:05 INFO [TGCC-IRENE] Submitted job with ID:['4937757']
2026-06-21 19:30:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:30:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8
2026-06-21 19:30:05 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/end.2020021413_2_ENS8.nc
2026-06-21 19:30:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-06-21 19:30:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:30:05 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-06-21 19:30:05 INFO Queuing job for member 8...
2026-06-21 19:30:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:30:05 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-06-21 19:30:06 INFO Found: ['4937759']
2026-06-21 19:30:11 INFO [TGCC-IRENE] Submitted job with ID:['4937759']
2026-06-21 19:30:11 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:30:11 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9
2026-06-21 19:30:11 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/end.2020021413_2_ENS9.nc
2026-06-21 19:30:11 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-06-21 19:30:11 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:30:11 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-06-21 19:30:12 INFO Queuing job for member 9...
2026-06-21 19:30:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:30:12 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-06-21 19:30:12 INFO Found: ['4937760']
2026-06-21 19:30:17 INFO [TGCC-IRENE] Submitted job with ID:['4937760']
2026-06-21 19:30:17 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:30:17 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10
2026-06-21 19:30:17 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/end.2020021413_2_ENS10.nc
2026-06-21 19:30:17 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-06-21 19:30:17 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:30:17 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-06-21 19:30:17 INFO Queuing job for member 10...
2026-06-21 19:30:17 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:30:17 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-06-21 19:30:18 INFO Found: ['4937761']
2026-06-21 19:30:23 INFO [TGCC-IRENE] Submitted job with ID:['4937761']
2026-06-21 19:30:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:30:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11
2026-06-21 19:30:23 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/end.2020021413_2_ENS11.nc
2026-06-21 19:30:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-06-21 19:30:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:30:23 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-06-21 19:30:23 INFO Queuing job for member 11...
2026-06-21 19:30:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:30:23 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-06-21 19:30:24 INFO Found: ['4937762']
2026-06-21 19:30:29 INFO [TGCC-IRENE] Submitted job with ID:['4937762']
2026-06-21 19:30:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:30:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12
2026-06-21 19:30:29 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/end.2020021413_2_ENS12.nc
2026-06-21 19:30:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-06-21 19:30:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:30:29 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-06-21 19:30:29 INFO Queuing job for member 12...
2026-06-21 19:30:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:30:29 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-06-21 19:30:30 INFO Found: ['4937764']
2026-06-21 19:30:35 INFO [TGCC-IRENE] Submitted job with ID:['4937764']
2026-06-21 19:30:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:30:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13
2026-06-21 19:30:35 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/end.2020021413_2_ENS13.nc
2026-06-21 19:30:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-06-21 19:30:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:30:35 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-06-21 19:30:35 INFO Queuing job for member 13...
2026-06-21 19:30:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:30:35 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-06-21 19:30:36 INFO Found: ['4937765']
2026-06-21 19:30:41 INFO [TGCC-IRENE] Submitted job with ID:['4937765']
2026-06-21 19:30:41 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:30:41 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14
2026-06-21 19:30:41 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/end.2020021413_2_ENS14.nc
2026-06-21 19:30:41 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-06-21 19:30:41 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:30:41 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-06-21 19:30:41 INFO Queuing job for member 14...
2026-06-21 19:30:41 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:30:41 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-06-21 19:30:42 INFO Found: ['4937766']
2026-06-21 19:30:47 INFO [TGCC-IRENE] Submitted job with ID:['4937766']
2026-06-21 19:30:47 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:30:47 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15
2026-06-21 19:30:47 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/end.2020021413_2_ENS15.nc
2026-06-21 19:30:47 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-06-21 19:30:47 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:30:47 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-06-21 19:30:47 INFO Queuing job for member 15...
2026-06-21 19:30:47 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:30:47 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-06-21 19:30:47 INFO Found: ['4937767']
2026-06-21 19:30:52 INFO [TGCC-IRENE] Submitted job with ID:['4937767']
2026-06-21 19:30:52 INFO Checking job status ...
2026-06-21 19:30:52 INFO None 4937751: status RUNNING/PENDING
2026-06-21 19:30:52 INFO None 4937752: status RUNNING/PENDING
2026-06-21 19:30:52 INFO None 4937753: status RUNNING/PENDING
2026-06-21 19:30:52 INFO None 4937754: status RUNNING/PENDING
2026-06-21 19:30:52 INFO None 4937755: status RUNNING/PENDING
2026-06-21 19:30:52 INFO None 4937756: status RUNNING/PENDING
2026-06-21 19:30:52 INFO None 4937757: status RUNNING/PENDING
2026-06-21 19:30:52 INFO None 4937759: status RUNNING/PENDING
2026-06-21 19:30:52 INFO None 4937760: status RUNNING/PENDING
2026-06-21 19:30:53 INFO None 4937761: status RUNNING/PENDING
2026-06-21 19:30:53 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:30:53 INFO None 4937764: status RUNNING/PENDING
2026-06-21 19:30:53 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:30:53 INFO None 4937766: status RUNNING/PENDING
2026-06-21 19:30:53 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:30:53 INFO Jobs still running: ['4937751', '4937752', '4937753', '4937754', '4937755', '4937756', '4937757', '4937759', '4937760', '4937761', '4937762', '4937764', '4937765', '4937766', '4937767']. Waiting...
2026-06-21 19:31:08 INFO None 4937751: status RUNNING/PENDING
2026-06-21 19:31:08 INFO None 4937752: status RUNNING/PENDING
2026-06-21 19:31:08 INFO None 4937753: status RUNNING/PENDING
2026-06-21 19:31:08 INFO None 4937754: status RUNNING/PENDING
2026-06-21 19:31:08 INFO None 4937755: status RUNNING/PENDING
2026-06-21 19:31:08 INFO None 4937756: status RUNNING/PENDING
2026-06-21 19:31:08 INFO None 4937757: status RUNNING/PENDING
2026-06-21 19:31:08 INFO None 4937759: status RUNNING/PENDING
2026-06-21 19:31:08 INFO None 4937760: status RUNNING/PENDING
2026-06-21 19:31:08 INFO None 4937761: status RUNNING/PENDING
2026-06-21 19:31:08 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:31:08 INFO None 4937764: status RUNNING/PENDING
2026-06-21 19:31:08 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:31:08 INFO None 4937766: status RUNNING/PENDING
2026-06-21 19:31:08 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:31:08 INFO Jobs still running: ['4937751', '4937752', '4937753', '4937754', '4937755', '4937756', '4937757', '4937759', '4937760', '4937761', '4937762', '4937764', '4937765', '4937766', '4937767']. Waiting...
2026-06-21 19:31:23 INFO None 4937751: status RUNNING/PENDING
2026-06-21 19:31:23 INFO None 4937752: status RUNNING/PENDING
2026-06-21 19:31:23 INFO None 4937753: status RUNNING/PENDING
2026-06-21 19:31:23 INFO None 4937754: status RUNNING/PENDING
2026-06-21 19:31:23 INFO None 4937755: status RUNNING/PENDING
2026-06-21 19:31:23 INFO None 4937756: status RUNNING/PENDING
2026-06-21 19:31:23 INFO None 4937757: status RUNNING/PENDING
2026-06-21 19:31:23 INFO None 4937759: status RUNNING/PENDING
2026-06-21 19:31:23 INFO None 4937760: status RUNNING/PENDING
2026-06-21 19:31:23 INFO None 4937761: status RUNNING/PENDING
2026-06-21 19:31:23 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:31:23 INFO None 4937764: status RUNNING/PENDING
2026-06-21 19:31:23 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:31:23 INFO None 4937766: status RUNNING/PENDING
2026-06-21 19:31:23 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:31:23 INFO Jobs still running: ['4937751', '4937752', '4937753', '4937754', '4937755', '4937756', '4937757', '4937759', '4937760', '4937761', '4937762', '4937764', '4937765', '4937766', '4937767']. Waiting...
2026-06-21 19:31:38 INFO None 4937751: status RUNNING/PENDING
2026-06-21 19:31:38 INFO None 4937752: status RUNNING/PENDING
2026-06-21 19:31:38 INFO None 4937753: status RUNNING/PENDING
2026-06-21 19:31:38 INFO None 4937754: status RUNNING/PENDING
2026-06-21 19:31:38 INFO None 4937755: status RUNNING/PENDING
2026-06-21 19:31:38 INFO None 4937756: status RUNNING/PENDING
2026-06-21 19:31:38 INFO None 4937757: status RUNNING/PENDING
2026-06-21 19:31:38 INFO None 4937759: status RUNNING/PENDING
2026-06-21 19:31:38 INFO None 4937760: status RUNNING/PENDING
2026-06-21 19:31:38 INFO None 4937761: status RUNNING/PENDING
2026-06-21 19:31:38 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:31:38 INFO None 4937764: status RUNNING/PENDING
2026-06-21 19:31:39 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:31:39 INFO None 4937766: status RUNNING/PENDING
2026-06-21 19:31:39 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:31:39 INFO Jobs still running: ['4937751', '4937752', '4937753', '4937754', '4937755', '4937756', '4937757', '4937759', '4937760', '4937761', '4937762', '4937764', '4937765', '4937766', '4937767']. Waiting...
2026-06-21 19:31:54 INFO None 4937751: status RUNNING/PENDING
2026-06-21 19:31:54 INFO None 4937752: status RUNNING/PENDING
2026-06-21 19:31:54 INFO None 4937753: status RUNNING/PENDING
2026-06-21 19:31:54 INFO None 4937754: status RUNNING/PENDING
2026-06-21 19:31:54 INFO None 4937755: status RUNNING/PENDING
2026-06-21 19:31:54 INFO None 4937756: status RUNNING/PENDING
2026-06-21 19:31:54 INFO None 4937757: status RUNNING/PENDING
2026-06-21 19:31:54 INFO None 4937759: status RUNNING/PENDING
2026-06-21 19:31:54 INFO None 4937760: status RUNNING/PENDING
2026-06-21 19:31:54 INFO None 4937761: status RUNNING/PENDING
2026-06-21 19:31:54 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:31:54 INFO None 4937764: status RUNNING/PENDING
2026-06-21 19:31:54 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:31:54 INFO None 4937766: status RUNNING/PENDING
2026-06-21 19:31:54 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:31:54 INFO Jobs still running: ['4937751', '4937752', '4937753', '4937754', '4937755', '4937756', '4937757', '4937759', '4937760', '4937761', '4937762', '4937764', '4937765', '4937766', '4937767']. Waiting...
2026-06-21 19:32:09 INFO None 4937751: status RUNNING/PENDING
2026-06-21 19:32:09 INFO None 4937752: status RUNNING/PENDING
2026-06-21 19:32:09 INFO None 4937753: status RUNNING/PENDING
2026-06-21 19:32:09 INFO None 4937754: status RUNNING/PENDING
2026-06-21 19:32:09 INFO None 4937755: status RUNNING/PENDING
2026-06-21 19:32:09 INFO None 4937756: status RUNNING/PENDING
2026-06-21 19:32:09 INFO None 4937757: status RUNNING/PENDING
2026-06-21 19:32:09 INFO None 4937759: status RUNNING/PENDING
2026-06-21 19:32:09 INFO None 4937760: status RUNNING/PENDING
2026-06-21 19:32:09 INFO None 4937761: status RUNNING/PENDING
2026-06-21 19:32:09 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:32:09 INFO None 4937764: status RUNNING/PENDING
2026-06-21 19:32:09 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:32:09 INFO None 4937766: status RUNNING/PENDING
2026-06-21 19:32:09 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:32:09 INFO Jobs still running: ['4937751', '4937752', '4937753', '4937754', '4937755', '4937756', '4937757', '4937759', '4937760', '4937761', '4937762', '4937764', '4937765', '4937766', '4937767']. Waiting...
2026-06-21 19:32:24 INFO None 4937751: status RUNNING/PENDING
2026-06-21 19:32:24 INFO None 4937752: status RUNNING/PENDING
2026-06-21 19:32:24 INFO None 4937753: status RUNNING/PENDING
2026-06-21 19:32:24 INFO None 4937754: status RUNNING/PENDING
2026-06-21 19:32:24 INFO None 4937755: status RUNNING/PENDING
2026-06-21 19:32:24 INFO None 4937756: status RUNNING/PENDING
2026-06-21 19:32:24 INFO None 4937757: status RUNNING/PENDING
2026-06-21 19:32:24 INFO None 4937759: status RUNNING/PENDING
2026-06-21 19:32:24 INFO None 4937760: status RUNNING/PENDING
2026-06-21 19:32:24 INFO None 4937761: status RUNNING/PENDING
2026-06-21 19:32:24 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:32:24 INFO None 4937764: status RUNNING/PENDING
2026-06-21 19:32:24 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:32:24 INFO None 4937766: status RUNNING/PENDING
2026-06-21 19:32:24 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:32:24 INFO Jobs still running: ['4937751', '4937752', '4937753', '4937754', '4937755', '4937756', '4937757', '4937759', '4937760', '4937761', '4937762', '4937764', '4937765', '4937766', '4937767']. Waiting...
2026-06-21 19:32:40 INFO None 4937751: status RUNNING/PENDING
2026-06-21 19:32:40 INFO None 4937752: status RUNNING/PENDING
2026-06-21 19:32:40 INFO None 4937753: status RUNNING/PENDING
2026-06-21 19:32:40 INFO None 4937754: status RUNNING/PENDING
2026-06-21 19:32:40 INFO None 4937755: status RUNNING/PENDING
2026-06-21 19:32:40 INFO None 4937756: status RUNNING/PENDING
2026-06-21 19:32:40 INFO None 4937757: status RUNNING/PENDING
2026-06-21 19:32:40 INFO None 4937759: status RUNNING/PENDING
2026-06-21 19:32:40 INFO None 4937760: status RUNNING/PENDING
2026-06-21 19:32:40 INFO None 4937761: status RUNNING/PENDING
2026-06-21 19:32:40 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:32:40 INFO None 4937764: status RUNNING/PENDING
2026-06-21 19:32:40 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:32:40 INFO None 4937766: status RUNNING/PENDING
2026-06-21 19:32:40 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:32:40 INFO Jobs still running: ['4937751', '4937752', '4937753', '4937754', '4937755', '4937756', '4937757', '4937759', '4937760', '4937761', '4937762', '4937764', '4937765', '4937766', '4937767']. Waiting...
2026-06-21 19:32:55 INFO None 4937751: status RUNNING/PENDING
2026-06-21 19:32:55 INFO None 4937752: status RUNNING/PENDING
2026-06-21 19:32:55 INFO None 4937753: status RUNNING/PENDING
2026-06-21 19:32:55 INFO None 4937754: status RUNNING/PENDING
2026-06-21 19:32:55 INFO None 4937755: status RUNNING/PENDING
2026-06-21 19:32:55 INFO None 4937756: status RUNNING/PENDING
2026-06-21 19:32:55 INFO None 4937757: status RUNNING/PENDING
2026-06-21 19:32:55 INFO None 4937759: status RUNNING/PENDING
2026-06-21 19:32:55 INFO None 4937760: status RUNNING/PENDING
2026-06-21 19:32:55 INFO None 4937761: status RUNNING/PENDING
2026-06-21 19:32:55 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:32:55 INFO None 4937764: status RUNNING/PENDING
2026-06-21 19:32:55 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:32:55 INFO None 4937766: status RUNNING/PENDING
2026-06-21 19:32:55 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:32:55 INFO Jobs still running: ['4937751', '4937752', '4937753', '4937754', '4937755', '4937756', '4937757', '4937759', '4937760', '4937761', '4937762', '4937764', '4937765', '4937766', '4937767']. Waiting...
2026-06-21 19:33:10 INFO None 4937751: status RUNNING/PENDING
2026-06-21 19:33:10 INFO None 4937752: status RUNNING/PENDING
2026-06-21 19:33:10 INFO None 4937753: status RUNNING/PENDING
2026-06-21 19:33:10 INFO None 4937754: status RUNNING/PENDING
2026-06-21 19:33:10 INFO None 4937755: status RUNNING/PENDING
2026-06-21 19:33:10 INFO None 4937756: status RUNNING/PENDING
2026-06-21 19:33:10 INFO None 4937757: status RUNNING/PENDING
2026-06-21 19:33:10 INFO None 4937759: status RUNNING/PENDING
2026-06-21 19:33:10 INFO None 4937760: status RUNNING/PENDING
2026-06-21 19:33:10 INFO None 4937761: status RUNNING/PENDING
2026-06-21 19:33:10 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:33:10 INFO None 4937764: status RUNNING/PENDING
2026-06-21 19:33:10 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:33:10 INFO None 4937766: status RUNNING/PENDING
2026-06-21 19:33:10 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:33:10 INFO Jobs still running: ['4937751', '4937752', '4937753', '4937754', '4937755', '4937756', '4937757', '4937759', '4937760', '4937761', '4937762', '4937764', '4937765', '4937766', '4937767']. Waiting...
2026-06-21 19:33:25 INFO None 4937751: status RUNNING/PENDING
2026-06-21 19:33:25 INFO None 4937752: status RUNNING/PENDING
2026-06-21 19:33:26 INFO None 4937753: status RUNNING/PENDING
2026-06-21 19:33:26 INFO None 4937754: status RUNNING/PENDING
2026-06-21 19:33:26 INFO None 4937755: status RUNNING/PENDING
2026-06-21 19:33:26 INFO None 4937756: status RUNNING/PENDING
2026-06-21 19:33:26 INFO None 4937757: status RUNNING/PENDING
2026-06-21 19:33:26 INFO None 4937759: status RUNNING/PENDING
2026-06-21 19:33:26 INFO None 4937760: status RUNNING/PENDING
2026-06-21 19:33:26 INFO None 4937761: status RUNNING/PENDING
2026-06-21 19:33:26 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:33:26 INFO None 4937764: status RUNNING/PENDING
2026-06-21 19:33:26 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:33:26 INFO None 4937766: status RUNNING/PENDING
2026-06-21 19:33:26 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:33:26 INFO Jobs still running: ['4937751', '4937752', '4937753', '4937754', '4937755', '4937756', '4937757', '4937759', '4937760', '4937761', '4937762', '4937764', '4937765', '4937766', '4937767']. Waiting...
2026-06-21 19:33:41 INFO None 4937751: status RUNNING/PENDING
2026-06-21 19:33:41 INFO None 4937752: status RUNNING/PENDING
2026-06-21 19:33:41 INFO None 4937753: status RUNNING/PENDING
2026-06-21 19:33:41 INFO None 4937754: status RUNNING/PENDING
2026-06-21 19:33:41 INFO None 4937755: status RUNNING/PENDING
2026-06-21 19:33:41 INFO None 4937756: status RUNNING/PENDING
2026-06-21 19:33:41 INFO None 4937757: status RUNNING/PENDING
2026-06-21 19:33:41 INFO None 4937759: status RUNNING/PENDING
2026-06-21 19:33:41 INFO None 4937760: status RUNNING/PENDING
2026-06-21 19:33:41 INFO None 4937761: status RUNNING/PENDING
2026-06-21 19:33:41 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:33:41 INFO None 4937764: status RUNNING/PENDING
2026-06-21 19:33:41 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:33:41 INFO None 4937766: status RUNNING/PENDING
2026-06-21 19:33:41 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:33:41 INFO Jobs still running: ['4937751', '4937752', '4937753', '4937754', '4937755', '4937756', '4937757', '4937759', '4937760', '4937761', '4937762', '4937764', '4937765', '4937766', '4937767']. Waiting...
2026-06-21 19:33:56 INFO None 4937751: status RUNNING/PENDING
2026-06-21 19:33:56 INFO None 4937752: status RUNNING/PENDING
2026-06-21 19:33:56 INFO None 4937753: status RUNNING/PENDING
2026-06-21 19:33:56 INFO None 4937754: status RUNNING/PENDING
2026-06-21 19:33:56 INFO None 4937755: status RUNNING/PENDING
2026-06-21 19:33:56 INFO None 4937756: status RUNNING/PENDING
2026-06-21 19:33:56 INFO None 4937757: status RUNNING/PENDING
2026-06-21 19:33:56 INFO None 4937759: status RUNNING/PENDING
2026-06-21 19:33:56 INFO None 4937760: status RUNNING/PENDING
2026-06-21 19:33:56 INFO None 4937761: status RUNNING/PENDING
2026-06-21 19:33:56 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:33:56 INFO None 4937764: status RUNNING/PENDING
2026-06-21 19:33:56 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:33:56 INFO None 4937766: status RUNNING/PENDING
2026-06-21 19:33:56 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:33:56 INFO Jobs still running: ['4937751', '4937752', '4937753', '4937754', '4937755', '4937756', '4937757', '4937759', '4937760', '4937761', '4937762', '4937764', '4937765', '4937766', '4937767']. Waiting...
2026-06-21 19:34:11 INFO None 4937751: status RUNNING/PENDING
2026-06-21 19:34:11 INFO None 4937752: status RUNNING/PENDING
2026-06-21 19:34:11 INFO None 4937753: status RUNNING/PENDING
2026-06-21 19:34:11 INFO None 4937754: status RUNNING/PENDING
2026-06-21 19:34:11 INFO None 4937755: status RUNNING/PENDING
2026-06-21 19:34:11 INFO None 4937756: status RUNNING/PENDING
2026-06-21 19:34:12 INFO None 4937757: status RUNNING/PENDING
2026-06-21 19:34:12 INFO None 4937759: status RUNNING/PENDING
2026-06-21 19:34:12 INFO None 4937760: status RUNNING/PENDING
2026-06-21 19:34:12 INFO None 4937761: status RUNNING/PENDING
2026-06-21 19:34:12 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:34:12 INFO None 4937764: status RUNNING/PENDING
2026-06-21 19:34:12 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:34:12 INFO None 4937766: status RUNNING/PENDING
2026-06-21 19:34:12 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:34:12 INFO Jobs still running: ['4937751', '4937752', '4937753', '4937754', '4937755', '4937756', '4937757', '4937759', '4937760', '4937761', '4937762', '4937764', '4937765', '4937766', '4937767']. Waiting...
2026-06-21 19:34:27 INFO None 4937751: status RUNNING/PENDING
2026-06-21 19:34:27 INFO None 4937752: status RUNNING/PENDING
2026-06-21 19:34:27 INFO None 4937753: status RUNNING/PENDING
2026-06-21 19:34:27 INFO None 4937754: status RUNNING/PENDING
2026-06-21 19:34:27 INFO None 4937755: status RUNNING/PENDING
2026-06-21 19:34:27 INFO None 4937756: status RUNNING/PENDING
2026-06-21 19:34:27 INFO None 4937757: status RUNNING/PENDING
2026-06-21 19:34:27 INFO None 4937759: status RUNNING/PENDING
2026-06-21 19:34:27 INFO None 4937760: status RUNNING/PENDING
2026-06-21 19:34:27 INFO None 4937761: status RUNNING/PENDING
2026-06-21 19:34:27 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:34:27 INFO None 4937764: status RUNNING/PENDING
2026-06-21 19:34:27 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:34:27 INFO None 4937766: status RUNNING/PENDING
2026-06-21 19:34:27 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:34:27 INFO Jobs still running: ['4937751', '4937752', '4937753', '4937754', '4937755', '4937756', '4937757', '4937759', '4937760', '4937761', '4937762', '4937764', '4937765', '4937766', '4937767']. Waiting...
2026-06-21 19:34:42 INFO None 4937751: status RUNNING/PENDING
2026-06-21 19:34:42 INFO None 4937752: status RUNNING/PENDING
2026-06-21 19:34:42 INFO None 4937753: status RUNNING/PENDING
2026-06-21 19:34:42 INFO None 4937754: status RUNNING/PENDING
2026-06-21 19:34:42 INFO None 4937755: status RUNNING/PENDING
2026-06-21 19:34:42 INFO None 4937756: status RUNNING/PENDING
2026-06-21 19:34:42 INFO None 4937757: status RUNNING/PENDING
2026-06-21 19:34:42 INFO None 4937759: status RUNNING/PENDING
2026-06-21 19:34:42 INFO None 4937760: status RUNNING/PENDING
2026-06-21 19:34:42 INFO None 4937761: status RUNNING/PENDING
2026-06-21 19:34:42 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:34:42 INFO None 4937764: status RUNNING/PENDING
2026-06-21 19:34:42 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:34:42 INFO None 4937766: status RUNNING/PENDING
2026-06-21 19:34:42 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:34:42 INFO Jobs still running: ['4937751', '4937752', '4937753', '4937754', '4937755', '4937756', '4937757', '4937759', '4937760', '4937761', '4937762', '4937764', '4937765', '4937766', '4937767']. Waiting...
2026-06-21 19:34:57 INFO None 4937751: status RUNNING/PENDING
2026-06-21 19:34:57 INFO None 4937752: status RUNNING/PENDING
2026-06-21 19:34:57 INFO None 4937753: status RUNNING/PENDING
2026-06-21 19:34:57 INFO None 4937754: status RUNNING/PENDING
2026-06-21 19:34:57 INFO None 4937755: status RUNNING/PENDING
2026-06-21 19:34:57 INFO None 4937756: status RUNNING/PENDING
2026-06-21 19:34:57 INFO None 4937757: status RUNNING/PENDING
2026-06-21 19:34:57 INFO None 4937759: status RUNNING/PENDING
2026-06-21 19:34:57 INFO None 4937760: status RUNNING/PENDING
2026-06-21 19:34:58 INFO None 4937761: status RUNNING/PENDING
2026-06-21 19:34:58 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:34:58 INFO None 4937764: status RUNNING/PENDING
2026-06-21 19:34:58 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:34:58 INFO None 4937766: status RUNNING/PENDING
2026-06-21 19:34:58 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:34:58 INFO Jobs still running: ['4937751', '4937752', '4937753', '4937754', '4937755', '4937756', '4937757', '4937759', '4937760', '4937761', '4937762', '4937764', '4937765', '4937766', '4937767']. Waiting...
2026-06-21 19:35:13 INFO None 4937751: status RUNNING/PENDING
2026-06-21 19:35:13 INFO None 4937752: status RUNNING/PENDING
2026-06-21 19:35:13 INFO None 4937753: status RUNNING/PENDING
2026-06-21 19:35:13 INFO None 4937754: status RUNNING/PENDING
2026-06-21 19:35:13 INFO None 4937755: status RUNNING/PENDING
2026-06-21 19:35:13 INFO None 4937756: status RUNNING/PENDING
2026-06-21 19:35:13 INFO None 4937757: status RUNNING/PENDING
2026-06-21 19:35:13 INFO None 4937759: status RUNNING/PENDING
2026-06-21 19:35:13 INFO None 4937760: status RUNNING/PENDING
2026-06-21 19:35:13 INFO None 4937761: status RUNNING/PENDING
2026-06-21 19:35:13 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:35:13 INFO None 4937764: status RUNNING/PENDING
2026-06-21 19:35:13 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:35:13 INFO None 4937766: status RUNNING/PENDING
2026-06-21 19:35:13 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:35:13 INFO Jobs still running: ['4937751', '4937752', '4937753', '4937754', '4937755', '4937756', '4937757', '4937759', '4937760', '4937761', '4937762', '4937764', '4937765', '4937766', '4937767']. Waiting...
2026-06-21 19:35:28 INFO None 4937751: status RUNNING/PENDING
2026-06-21 19:35:28 INFO None 4937752: status RUNNING/PENDING
2026-06-21 19:35:28 INFO None 4937753: status RUNNING/PENDING
2026-06-21 19:35:28 INFO None 4937754: status RUNNING/PENDING
2026-06-21 19:35:28 INFO None 4937755: status RUNNING/PENDING
2026-06-21 19:35:28 INFO None 4937756: status RUNNING/PENDING
2026-06-21 19:35:28 INFO None 4937757: status RUNNING/PENDING
2026-06-21 19:35:28 INFO None 4937759: status RUNNING/PENDING
2026-06-21 19:35:28 INFO None 4937760: status RUNNING/PENDING
2026-06-21 19:35:28 INFO None 4937761: status RUNNING/PENDING
2026-06-21 19:35:28 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:35:28 INFO None 4937764: status RUNNING/PENDING
2026-06-21 19:35:28 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:35:28 INFO None 4937766: status RUNNING/PENDING
2026-06-21 19:35:28 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:35:28 INFO Jobs still running: ['4937751', '4937752', '4937753', '4937754', '4937755', '4937756', '4937757', '4937759', '4937760', '4937761', '4937762', '4937764', '4937765', '4937766', '4937767']. Waiting...
2026-06-21 19:35:43 INFO None 4937751: status RUNNING/PENDING
2026-06-21 19:35:43 INFO None 4937752: status RUNNING/PENDING
2026-06-21 19:35:43 INFO None 4937753: status RUNNING/PENDING
2026-06-21 19:35:43 INFO None 4937754: status RUNNING/PENDING
2026-06-21 19:35:43 INFO None 4937755: status RUNNING/PENDING
2026-06-21 19:35:43 INFO None 4937756: status RUNNING/PENDING
2026-06-21 19:35:43 INFO None 4937757: status RUNNING/PENDING
2026-06-21 19:35:43 INFO None 4937759: status RUNNING/PENDING
2026-06-21 19:35:43 INFO None 4937760: status RUNNING/PENDING
2026-06-21 19:35:43 INFO None 4937761: status RUNNING/PENDING
2026-06-21 19:35:43 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:35:44 INFO None 4937764: status RUNNING/PENDING
2026-06-21 19:35:44 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:35:44 INFO None 4937766: status RUNNING/PENDING
2026-06-21 19:35:44 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:35:44 INFO Jobs still running: ['4937751', '4937752', '4937753', '4937754', '4937755', '4937756', '4937757', '4937759', '4937760', '4937761', '4937762', '4937764', '4937765', '4937766', '4937767']. Waiting...
2026-06-21 19:35:59 INFO None 4937751: status RUNNING/PENDING
2026-06-21 19:35:59 INFO None 4937752: status RUNNING/PENDING
2026-06-21 19:35:59 INFO None 4937753: status RUNNING/PENDING
2026-06-21 19:35:59 INFO None 4937754: status RUNNING/PENDING
2026-06-21 19:35:59 INFO None 4937755: status RUNNING/PENDING
2026-06-21 19:35:59 INFO None 4937756: status RUNNING/PENDING
2026-06-21 19:35:59 INFO None 4937757: status RUNNING/PENDING
2026-06-21 19:35:59 INFO None 4937759: status RUNNING/PENDING
2026-06-21 19:35:59 INFO None 4937760: status RUNNING/PENDING
2026-06-21 19:35:59 INFO None 4937761: status RUNNING/PENDING
2026-06-21 19:35:59 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:35:59 INFO None 4937764: status RUNNING/PENDING
2026-06-21 19:35:59 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:35:59 INFO None 4937766: status RUNNING/PENDING
2026-06-21 19:35:59 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:35:59 INFO Jobs still running: ['4937751', '4937752', '4937753', '4937754', '4937755', '4937756', '4937757', '4937759', '4937760', '4937761', '4937762', '4937764', '4937765', '4937766', '4937767']. Waiting...
2026-06-21 19:36:14 INFO None 4937751: status RUNNING/PENDING
2026-06-21 19:36:14 INFO None 4937752: status RUNNING/PENDING
2026-06-21 19:36:14 INFO None 4937753: status RUNNING/PENDING
2026-06-21 19:36:14 INFO None 4937754: status RUNNING/PENDING
2026-06-21 19:36:14 INFO None 4937755: status RUNNING/PENDING
2026-06-21 19:36:14 INFO None 4937756: status RUNNING/PENDING
2026-06-21 19:36:14 INFO None 4937757: status RUNNING/PENDING
2026-06-21 19:36:14 INFO None 4937759: status RUNNING/PENDING
2026-06-21 19:36:14 INFO None 4937760: status RUNNING/PENDING
2026-06-21 19:36:14 INFO None 4937761: status RUNNING/PENDING
2026-06-21 19:36:14 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:36:14 INFO None 4937764: status RUNNING/PENDING
2026-06-21 19:36:14 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:36:14 INFO None 4937766: status RUNNING/PENDING
2026-06-21 19:36:14 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:36:14 INFO Jobs still running: ['4937751', '4937752', '4937753', '4937754', '4937755', '4937756', '4937757', '4937759', '4937760', '4937761', '4937762', '4937764', '4937765', '4937766', '4937767']. Waiting...
2026-06-21 19:36:29 INFO None 4937751: status RUNNING/PENDING
2026-06-21 19:36:29 INFO None 4937752: status RUNNING/PENDING
2026-06-21 19:36:29 INFO None 4937753: status RUNNING/PENDING
2026-06-21 19:36:29 INFO None 4937754: status RUNNING/PENDING
2026-06-21 19:36:29 INFO None 4937755: status RUNNING/PENDING
2026-06-21 19:36:29 INFO None 4937756: status RUNNING/PENDING
2026-06-21 19:36:29 INFO None 4937757: status RUNNING/PENDING
2026-06-21 19:36:29 INFO None 4937759: status RUNNING/PENDING
2026-06-21 19:36:29 INFO None 4937760: status RUNNING/PENDING
2026-06-21 19:36:29 INFO None 4937761: status RUNNING/PENDING
2026-06-21 19:36:29 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:36:29 INFO None 4937764: status RUNNING/PENDING
2026-06-21 19:36:29 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:36:29 INFO None 4937766: status RUNNING/PENDING
2026-06-21 19:36:30 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:36:30 INFO Jobs still running: ['4937751', '4937752', '4937753', '4937754', '4937755', '4937756', '4937757', '4937759', '4937760', '4937761', '4937762', '4937764', '4937765', '4937766', '4937767']. Waiting...
2026-06-21 19:36:45 INFO None 4937751: status RUNNING/PENDING
2026-06-21 19:36:45 INFO None 4937752: status RUNNING/PENDING
2026-06-21 19:36:45 INFO None 4937753: status RUNNING/PENDING
2026-06-21 19:36:45 INFO None 4937754: status RUNNING/PENDING
2026-06-21 19:36:45 INFO None 4937755: status RUNNING/PENDING
2026-06-21 19:36:45 INFO None 4937756: status RUNNING/PENDING
2026-06-21 19:36:45 INFO None 4937757: status RUNNING/PENDING
2026-06-21 19:36:45 INFO None 4937759: status RUNNING/PENDING
2026-06-21 19:36:45 INFO None 4937760: status RUNNING/PENDING
2026-06-21 19:36:45 INFO None 4937761: status RUNNING/PENDING
2026-06-21 19:36:45 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:36:45 INFO None 4937764: status RUNNING/PENDING
2026-06-21 19:36:45 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:36:45 INFO None 4937766: status RUNNING/PENDING
2026-06-21 19:36:45 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:36:45 INFO Jobs still running: ['4937751', '4937752', '4937753', '4937754', '4937755', '4937756', '4937757', '4937759', '4937760', '4937761', '4937762', '4937764', '4937765', '4937766', '4937767']. Waiting...
2026-06-21 19:37:00 INFO None 4937751: status RUNNING/PENDING
2026-06-21 19:37:00 INFO None 4937752: status RUNNING/PENDING
2026-06-21 19:37:00 INFO None 4937753: status RUNNING/PENDING
2026-06-21 19:37:00 INFO None 4937754: status RUNNING/PENDING
2026-06-21 19:37:00 INFO None 4937755: status RUNNING/PENDING
2026-06-21 19:37:00 INFO None 4937756: status RUNNING/PENDING
2026-06-21 19:37:00 INFO None 4937757: status RUNNING/PENDING
2026-06-21 19:37:00 INFO None 4937759: status RUNNING/PENDING
2026-06-21 19:37:00 INFO None 4937760: status RUNNING/PENDING
2026-06-21 19:37:00 INFO None 4937761: status RUNNING/PENDING
2026-06-21 19:37:00 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:37:00 INFO None 4937764: status RUNNING/PENDING
2026-06-21 19:37:00 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:37:00 INFO None 4937766: status RUNNING/PENDING
2026-06-21 19:37:00 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:37:00 INFO Jobs still running: ['4937751', '4937752', '4937753', '4937754', '4937755', '4937756', '4937757', '4937759', '4937760', '4937761', '4937762', '4937764', '4937765', '4937766', '4937767']. Waiting...
2026-06-21 19:37:15 INFO None 4937751: status RUNNING/PENDING
2026-06-21 19:37:15 INFO None 4937752: status RUNNING/PENDING
2026-06-21 19:37:15 INFO None 4937753: status FINISHED
2026-06-21 19:37:15 INFO None 4937754: status RUNNING/PENDING
2026-06-21 19:37:15 INFO None 4937755: status RUNNING/PENDING
2026-06-21 19:37:15 INFO None 4937756: status RUNNING/PENDING
2026-06-21 19:37:15 INFO None 4937757: status RUNNING/PENDING
2026-06-21 19:37:15 INFO None 4937759: status RUNNING/PENDING
2026-06-21 19:37:15 INFO None 4937760: status RUNNING/PENDING
2026-06-21 19:37:15 INFO None 4937761: status RUNNING/PENDING
2026-06-21 19:37:15 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:37:15 INFO None 4937764: status RUNNING/PENDING
2026-06-21 19:37:15 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:37:15 INFO None 4937766: status RUNNING/PENDING
2026-06-21 19:37:15 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:37:15 INFO Jobs still running: ['4937751', '4937752', '4937754', '4937755', '4937756', '4937757', '4937759', '4937760', '4937761', '4937762', '4937764', '4937765', '4937766', '4937767']. Waiting...
2026-06-21 19:37:30 INFO None 4937751: status RUNNING/PENDING
2026-06-21 19:37:30 INFO None 4937752: status RUNNING/PENDING
2026-06-21 19:37:31 INFO None 4937753: status FINISHED
2026-06-21 19:37:31 INFO None 4937754: status FINISHED
2026-06-21 19:37:31 INFO None 4937755: status RUNNING/PENDING
2026-06-21 19:37:31 INFO None 4937756: status RUNNING/PENDING
2026-06-21 19:37:31 INFO None 4937757: status RUNNING/PENDING
2026-06-21 19:37:31 INFO None 4937759: status RUNNING/PENDING
2026-06-21 19:37:31 INFO None 4937760: status RUNNING/PENDING
2026-06-21 19:37:31 INFO None 4937761: status RUNNING/PENDING
2026-06-21 19:37:31 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:37:31 INFO None 4937764: status RUNNING/PENDING
2026-06-21 19:37:31 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:37:31 INFO None 4937766: status RUNNING/PENDING
2026-06-21 19:37:31 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:37:31 INFO Jobs still running: ['4937751', '4937752', '4937755', '4937756', '4937757', '4937759', '4937760', '4937761', '4937762', '4937764', '4937765', '4937766', '4937767']. Waiting...
2026-06-21 19:37:46 INFO None 4937751: status RUNNING/PENDING
2026-06-21 19:37:46 INFO None 4937752: status FINISHED
2026-06-21 19:37:46 INFO None 4937753: status FINISHED
2026-06-21 19:37:46 INFO None 4937754: status FINISHED
2026-06-21 19:37:46 INFO None 4937755: status RUNNING/PENDING
2026-06-21 19:37:46 INFO None 4937756: status RUNNING/PENDING
2026-06-21 19:37:46 INFO None 4937757: status RUNNING/PENDING
2026-06-21 19:37:46 INFO None 4937759: status RUNNING/PENDING
2026-06-21 19:37:46 INFO None 4937760: status RUNNING/PENDING
2026-06-21 19:37:46 INFO None 4937761: status RUNNING/PENDING
2026-06-21 19:37:46 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:37:46 INFO None 4937764: status RUNNING/PENDING
2026-06-21 19:37:46 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:37:46 INFO None 4937766: status RUNNING/PENDING
2026-06-21 19:37:46 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:37:46 INFO Jobs still running: ['4937751', '4937755', '4937756', '4937757', '4937759', '4937760', '4937761', '4937762', '4937764', '4937765', '4937766', '4937767']. Waiting...
2026-06-21 19:38:01 INFO None 4937751: status RUNNING/PENDING
2026-06-21 19:38:01 INFO None 4937752: status FINISHED
2026-06-21 19:38:01 INFO None 4937753: status FINISHED
2026-06-21 19:38:01 INFO None 4937754: status FINISHED
2026-06-21 19:38:01 INFO None 4937755: status RUNNING/PENDING
2026-06-21 19:38:01 INFO None 4937756: status FINISHED
2026-06-21 19:38:01 INFO None 4937757: status FINISHED
2026-06-21 19:38:01 INFO None 4937759: status RUNNING/PENDING
2026-06-21 19:38:01 INFO None 4937760: status RUNNING/PENDING
2026-06-21 19:38:01 INFO None 4937761: status RUNNING/PENDING
2026-06-21 19:38:01 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:38:01 INFO None 4937764: status RUNNING/PENDING
2026-06-21 19:38:01 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:38:01 INFO None 4937766: status RUNNING/PENDING
2026-06-21 19:38:01 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:38:01 INFO Jobs still running: ['4937751', '4937755', '4937759', '4937760', '4937761', '4937762', '4937764', '4937765', '4937766', '4937767']. Waiting...
2026-06-21 19:38:16 INFO None 4937751: status FINISHED
2026-06-21 19:38:16 INFO None 4937752: status FINISHED
2026-06-21 19:38:16 INFO None 4937753: status FINISHED
2026-06-21 19:38:16 INFO None 4937754: status FINISHED
2026-06-21 19:38:16 INFO None 4937755: status FINISHED
2026-06-21 19:38:16 INFO None 4937756: status FINISHED
2026-06-21 19:38:16 INFO None 4937757: status FINISHED
2026-06-21 19:38:16 INFO None 4937759: status FINISHED
2026-06-21 19:38:17 INFO None 4937760: status FINISHED
2026-06-21 19:38:17 INFO None 4937761: status RUNNING/PENDING
2026-06-21 19:38:17 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:38:17 INFO None 4937764: status FINISHED
2026-06-21 19:38:17 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:38:17 INFO None 4937766: status RUNNING/PENDING
2026-06-21 19:38:17 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:38:17 INFO Jobs still running: ['4937761', '4937762', '4937765', '4937766', '4937767']. Waiting...
2026-06-21 19:38:32 INFO None 4937751: status FINISHED
2026-06-21 19:38:32 INFO None 4937752: status FINISHED
2026-06-21 19:38:32 INFO None 4937753: status FINISHED
2026-06-21 19:38:32 INFO None 4937754: status FINISHED
2026-06-21 19:38:32 INFO None 4937755: status FINISHED
2026-06-21 19:38:32 INFO None 4937756: status FINISHED
2026-06-21 19:38:32 INFO None 4937757: status FINISHED
2026-06-21 19:38:32 INFO None 4937759: status FINISHED
2026-06-21 19:38:32 INFO None 4937760: status FINISHED
2026-06-21 19:38:32 INFO None 4937761: status FINISHED
2026-06-21 19:38:32 INFO None 4937762: status RUNNING/PENDING
2026-06-21 19:38:32 INFO None 4937764: status FINISHED
2026-06-21 19:38:32 INFO None 4937765: status RUNNING/PENDING
2026-06-21 19:38:32 INFO None 4937766: status FINISHED
2026-06-21 19:38:32 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:38:32 INFO Jobs still running: ['4937762', '4937765', '4937767']. Waiting...
2026-06-21 19:38:47 INFO None 4937751: status FINISHED
2026-06-21 19:38:47 INFO None 4937752: status FINISHED
2026-06-21 19:38:47 INFO None 4937753: status FINISHED
2026-06-21 19:38:47 INFO None 4937754: status FINISHED
2026-06-21 19:38:47 INFO None 4937755: status FINISHED
2026-06-21 19:38:47 INFO None 4937756: status FINISHED
2026-06-21 19:38:47 INFO None 4937757: status FINISHED
2026-06-21 19:38:47 INFO None 4937759: status FINISHED
2026-06-21 19:38:47 INFO None 4937760: status FINISHED
2026-06-21 19:38:47 INFO None 4937761: status FINISHED
2026-06-21 19:38:47 INFO None 4937762: status FINISHED
2026-06-21 19:38:47 INFO None 4937764: status FINISHED
2026-06-21 19:38:47 INFO None 4937765: status FINISHED
2026-06-21 19:38:47 INFO None 4937766: status FINISHED
2026-06-21 19:38:47 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:38:47 INFO Jobs still running: ['4937767']. Waiting...
2026-06-21 19:39:02 INFO None 4937751: status FINISHED
2026-06-21 19:39:02 INFO None 4937752: status FINISHED
2026-06-21 19:39:02 INFO None 4937753: status FINISHED
2026-06-21 19:39:02 INFO None 4937754: status FINISHED
2026-06-21 19:39:02 INFO None 4937755: status FINISHED
2026-06-21 19:39:02 INFO None 4937756: status FINISHED
2026-06-21 19:39:02 INFO None 4937757: status FINISHED
2026-06-21 19:39:02 INFO None 4937759: status FINISHED
2026-06-21 19:39:02 INFO None 4937760: status FINISHED
2026-06-21 19:39:02 INFO None 4937761: status FINISHED
2026-06-21 19:39:02 INFO None 4937762: status FINISHED
2026-06-21 19:39:02 INFO None 4937764: status FINISHED
2026-06-21 19:39:02 INFO None 4937765: status FINISHED
2026-06-21 19:39:02 INFO None 4937766: status FINISHED
2026-06-21 19:39:02 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:39:02 INFO Jobs still running: ['4937767']. Waiting...
2026-06-21 19:39:17 INFO None 4937751: status FINISHED
2026-06-21 19:39:18 INFO None 4937752: status FINISHED
2026-06-21 19:39:18 INFO None 4937753: status FINISHED
2026-06-21 19:39:18 INFO None 4937754: status FINISHED
2026-06-21 19:39:18 INFO None 4937755: status FINISHED
2026-06-21 19:39:18 INFO None 4937756: status FINISHED
2026-06-21 19:39:18 INFO None 4937757: status FINISHED
2026-06-21 19:39:18 INFO None 4937759: status FINISHED
2026-06-21 19:39:18 INFO None 4937760: status FINISHED
2026-06-21 19:39:18 INFO None 4937761: status FINISHED
2026-06-21 19:39:18 INFO None 4937762: status FINISHED
2026-06-21 19:39:18 INFO None 4937764: status FINISHED
2026-06-21 19:39:18 INFO None 4937765: status FINISHED
2026-06-21 19:39:18 INFO None 4937766: status FINISHED
2026-06-21 19:39:18 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:39:18 INFO Jobs still running: ['4937767']. Waiting...
2026-06-21 19:39:33 INFO None 4937751: status FINISHED
2026-06-21 19:39:33 INFO None 4937752: status FINISHED
2026-06-21 19:39:33 INFO None 4937753: status FINISHED
2026-06-21 19:39:33 INFO None 4937754: status FINISHED
2026-06-21 19:39:33 INFO None 4937755: status FINISHED
2026-06-21 19:39:33 INFO None 4937756: status FINISHED
2026-06-21 19:39:33 INFO None 4937757: status FINISHED
2026-06-21 19:39:33 INFO None 4937759: status FINISHED
2026-06-21 19:39:33 INFO None 4937760: status FINISHED
2026-06-21 19:39:33 INFO None 4937761: status FINISHED
2026-06-21 19:39:33 INFO None 4937762: status FINISHED
2026-06-21 19:39:33 INFO None 4937764: status FINISHED
2026-06-21 19:39:33 INFO None 4937765: status FINISHED
2026-06-21 19:39:33 INFO None 4937766: status FINISHED
2026-06-21 19:39:33 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:39:33 INFO Jobs still running: ['4937767']. Waiting...
2026-06-21 19:39:48 INFO None 4937751: status FINISHED
2026-06-21 19:39:48 INFO None 4937752: status FINISHED
2026-06-21 19:39:48 INFO None 4937753: status FINISHED
2026-06-21 19:39:48 INFO None 4937754: status FINISHED
2026-06-21 19:39:48 INFO None 4937755: status FINISHED
2026-06-21 19:39:48 INFO None 4937756: status FINISHED
2026-06-21 19:39:48 INFO None 4937757: status FINISHED
2026-06-21 19:39:48 INFO None 4937759: status FINISHED
2026-06-21 19:39:48 INFO None 4937760: status FINISHED
2026-06-21 19:39:48 INFO None 4937761: status FINISHED
2026-06-21 19:39:48 INFO None 4937762: status FINISHED
2026-06-21 19:39:48 INFO None 4937764: status FINISHED
2026-06-21 19:39:48 INFO None 4937765: status FINISHED
2026-06-21 19:39:48 INFO None 4937766: status FINISHED
2026-06-21 19:39:48 INFO None 4937767: status RUNNING/PENDING
2026-06-21 19:39:48 INFO Jobs still running: ['4937767']. Waiting...
2026-06-21 19:40:03 INFO None 4937751: status FINISHED
2026-06-21 19:40:03 INFO None 4937752: status FINISHED
2026-06-21 19:40:03 INFO None 4937753: status FINISHED
2026-06-21 19:40:03 INFO None 4937754: status FINISHED
2026-06-21 19:40:03 INFO None 4937755: status FINISHED
2026-06-21 19:40:03 INFO None 4937756: status FINISHED
2026-06-21 19:40:03 INFO None 4937757: status FINISHED
2026-06-21 19:40:03 INFO None 4937759: status FINISHED
2026-06-21 19:40:03 INFO None 4937760: status FINISHED
2026-06-21 19:40:03 INFO None 4937761: status FINISHED
2026-06-21 19:40:03 INFO None 4937762: status FINISHED
2026-06-21 19:40:04 INFO None 4937764: status FINISHED
2026-06-21 19:40:04 INFO None 4937765: status FINISHED
2026-06-21 19:40:04 INFO None 4937766: status FINISHED
2026-06-21 19:40:04 INFO None 4937767: status FINISHED
2026-06-21 19:40:04 INFO Jobs ['4937751', '4937752', '4937753', '4937754', '4937755', '4937756', '4937757', '4937759', '4937760', '4937761', '4937762', '4937764', '4937765', '4937766', '4937767'] have finished
2026-06-21 19:40:04 INFO Checking restart files were created ...
2026-06-21 19:40:04 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/end.2020021415_9_ENS1.nc(3339660275 bytes)
2026-06-21 19:40:04 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/end.2020021415_9_ENS2.nc(3339660275 bytes)
2026-06-21 19:40:04 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/end.2020021415_9_ENS3.nc(3339660275 bytes)
2026-06-21 19:40:04 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/end.2020021415_9_ENS4.nc(3339660275 bytes)
2026-06-21 19:40:04 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/end.2020021415_9_ENS5.nc(3339660275 bytes)
2026-06-21 19:40:04 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/end.2020021415_9_ENS6.nc(3339660275 bytes)
2026-06-21 19:40:04 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/end.2020021415_9_ENS7.nc(3339660275 bytes)
2026-06-21 19:40:04 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/end.2020021415_9_ENS8.nc(3339660275 bytes)
2026-06-21 19:40:04 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/end.2020021415_9_ENS9.nc(3339660275 bytes)
2026-06-21 19:40:04 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/end.2020021415_9_ENS10.nc(3339660275 bytes)
2026-06-21 19:40:04 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/end.2020021415_9_ENS11.nc(3339660275 bytes)
2026-06-21 19:40:04 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/end.2020021415_9_ENS12.nc(3339660275 bytes)
2026-06-21 19:40:04 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/end.2020021415_9_ENS13.nc(3339660275 bytes)
2026-06-21 19:40:04 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/end.2020021415_9_ENS14.nc(3339660275 bytes)
2026-06-21 19:40:04 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/end.2020021415_9_ENS15.nc(3339660275 bytes)
2026-06-21 19:40:04 INFO  Run_model() completed successfully.
2026-06-21 19:40:04 INFO [TIME] after_model_set_simulated_time current_time=2020-02-14 15:00:00 simulated_time=2020-02-15 00:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 19:40:04 INFO [TIME] gregorian_conversion simulated_time=2020-02-15 00:00:00 days=153081 seconds=0
2026-06-21 19:40:04 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-06-21 19:40:04 INFO [TIME] increment current_time 2020-02-14 15:00:00 -> 2020-02-15 00:00:00
2026-06-21 19:40:04 INFO [TIME] after_increment_before_assimilation current_time=2020-02-15 00:00:00 simulated_time=2020-02-15 00:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 19:40:04 INFO ---------->>> Running process_satellite_data()
2026-06-21 19:40:04 INFO [DART] No satellite data found, skipping assimilation
2026-06-21 19:40:04 INFO after_assimilation() skipped
2026-06-21 19:40:04 INFO Next run starts from 2020-02-15 00:00:00
2026-06-21 19:40:04 INFO Cycle is DONE; starting a new loop!
2026-06-21 19:40:04 INFO [TIME] step_end current_time=2020-02-15 00:00:00 simulated_time=2020-02-15 00:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 19:40:04 INFO [TIME] step_start current_time=2020-02-15 00:00:00 simulated_time=2020-02-15 00:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 19:40:04 INFO [TIME] window start=2020-02-15 00:00:00 end=2020-02-15 01:00:00 run_hours=1 has_assimilation=False
2026-06-21 19:40:04 INFO Creating directories and links for ENS1 to run chimere's parallel part
2026-06-21 19:40:04 INFO Linking EMIS ...
2026-06-21 19:40:04 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens01.nc
2026-06-21 19:40:04 INFO >> Checking links...
2026-06-21 19:40:07 INFO >> All links are good for ENS1  ...
2026-06-21 19:40:07 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:40:08 INFO Hourly dataset computed and listing created
2026-06-21 19:40:17 INFO Hourly dataset computed
2026-06-21 19:40:17 INFO Creating directories and links for ENS2 to run chimere's parallel part
2026-06-21 19:40:17 INFO Linking EMIS ...
2026-06-21 19:40:18 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens02.nc
2026-06-21 19:40:18 INFO >> Checking links...
2026-06-21 19:40:20 INFO >> All links are good for ENS2  ...
2026-06-21 19:40:20 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:40:21 INFO Hourly dataset computed and listing created
2026-06-21 19:40:22 INFO Hourly dataset computed
2026-06-21 19:40:22 INFO Creating directories and links for ENS3 to run chimere's parallel part
2026-06-21 19:40:22 INFO Linking EMIS ...
2026-06-21 19:40:22 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens03.nc
2026-06-21 19:40:22 INFO >> Checking links...
2026-06-21 19:40:24 INFO >> All links are good for ENS3  ...
2026-06-21 19:40:24 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:40:25 INFO Hourly dataset computed and listing created
2026-06-21 19:40:32 INFO Hourly dataset computed
2026-06-21 19:40:32 INFO Creating directories and links for ENS4 to run chimere's parallel part
2026-06-21 19:40:32 INFO Linking EMIS ...
2026-06-21 19:40:32 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens04.nc
2026-06-21 19:40:32 INFO >> Checking links...
2026-06-21 19:40:34 INFO >> All links are good for ENS4  ...
2026-06-21 19:40:34 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:40:35 INFO Hourly dataset computed and listing created
2026-06-21 19:40:36 INFO Hourly dataset computed
2026-06-21 19:40:36 INFO Creating directories and links for ENS5 to run chimere's parallel part
2026-06-21 19:40:36 INFO Linking EMIS ...
2026-06-21 19:40:36 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens05.nc
2026-06-21 19:40:36 INFO >> Checking links...
2026-06-21 19:40:39 INFO >> All links are good for ENS5  ...
2026-06-21 19:40:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:40:40 INFO Hourly dataset computed and listing created
2026-06-21 19:40:40 INFO Hourly dataset computed
2026-06-21 19:40:40 INFO Creating directories and links for ENS6 to run chimere's parallel part
2026-06-21 19:40:40 INFO Linking EMIS ...
2026-06-21 19:40:40 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens06.nc
2026-06-21 19:40:40 INFO >> Checking links...
2026-06-21 19:40:43 INFO >> All links are good for ENS6  ...
2026-06-21 19:40:43 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:40:44 INFO Hourly dataset computed and listing created
2026-06-21 19:40:44 INFO Hourly dataset computed
2026-06-21 19:40:44 INFO Creating directories and links for ENS7 to run chimere's parallel part
2026-06-21 19:40:44 INFO Linking EMIS ...
2026-06-21 19:40:44 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens07.nc
2026-06-21 19:40:44 INFO >> Checking links...
2026-06-21 19:40:47 INFO >> All links are good for ENS7  ...
2026-06-21 19:40:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:40:48 INFO Hourly dataset computed and listing created
2026-06-21 19:40:48 INFO Hourly dataset computed
2026-06-21 19:40:48 INFO Creating directories and links for ENS8 to run chimere's parallel part
2026-06-21 19:40:48 INFO Linking EMIS ...
2026-06-21 19:40:49 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens08.nc
2026-06-21 19:40:49 INFO >> Checking links...
2026-06-21 19:40:51 INFO >> All links are good for ENS8  ...
2026-06-21 19:40:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:40:52 INFO Hourly dataset computed and listing created
2026-06-21 19:40:52 INFO Hourly dataset computed
2026-06-21 19:40:52 INFO Creating directories and links for ENS9 to run chimere's parallel part
2026-06-21 19:40:52 INFO Linking EMIS ...
2026-06-21 19:40:53 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens09.nc
2026-06-21 19:40:53 INFO >> Checking links...
2026-06-21 19:40:55 INFO >> All links are good for ENS9  ...
2026-06-21 19:40:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:40:56 INFO Hourly dataset computed and listing created
2026-06-21 19:40:57 INFO Hourly dataset computed
2026-06-21 19:40:57 INFO Creating directories and links for ENS10 to run chimere's parallel part
2026-06-21 19:40:57 INFO Linking EMIS ...
2026-06-21 19:40:57 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens10.nc
2026-06-21 19:40:57 INFO >> Checking links...
2026-06-21 19:40:59 INFO >> All links are good for ENS10  ...
2026-06-21 19:40:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:41:00 INFO Hourly dataset computed and listing created
2026-06-21 19:41:01 INFO Hourly dataset computed
2026-06-21 19:41:01 INFO Creating directories and links for ENS11 to run chimere's parallel part
2026-06-21 19:41:01 INFO Linking EMIS ...
2026-06-21 19:41:01 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens11.nc
2026-06-21 19:41:01 INFO >> Checking links...
2026-06-21 19:41:04 INFO >> All links are good for ENS11  ...
2026-06-21 19:41:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:41:04 INFO Hourly dataset computed and listing created
2026-06-21 19:41:05 INFO Hourly dataset computed
2026-06-21 19:41:05 INFO Creating directories and links for ENS12 to run chimere's parallel part
2026-06-21 19:41:05 INFO Linking EMIS ...
2026-06-21 19:41:05 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens12.nc
2026-06-21 19:41:05 INFO >> Checking links...
2026-06-21 19:41:08 INFO >> All links are good for ENS12  ...
2026-06-21 19:41:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:41:09 INFO Hourly dataset computed and listing created
2026-06-21 19:41:12 INFO Hourly dataset computed
2026-06-21 19:41:12 INFO Creating directories and links for ENS13 to run chimere's parallel part
2026-06-21 19:41:12 INFO Linking EMIS ...
2026-06-21 19:41:12 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens13.nc
2026-06-21 19:41:12 INFO >> Checking links...
2026-06-21 19:41:15 INFO >> All links are good for ENS13  ...
2026-06-21 19:41:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:41:16 INFO Hourly dataset computed and listing created
2026-06-21 19:41:18 INFO Hourly dataset computed
2026-06-21 19:41:18 INFO Creating directories and links for ENS14 to run chimere's parallel part
2026-06-21 19:41:18 INFO Linking EMIS ...
2026-06-21 19:41:19 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens14.nc
2026-06-21 19:41:19 INFO >> Checking links...
2026-06-21 19:41:21 INFO >> All links are good for ENS14  ...
2026-06-21 19:41:21 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:41:22 INFO Hourly dataset computed and listing created
2026-06-21 19:41:25 INFO Hourly dataset computed
2026-06-21 19:41:25 INFO Creating directories and links for ENS15 to run chimere's parallel part
2026-06-21 19:41:25 INFO Linking EMIS ...
2026-06-21 19:41:25 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens15.nc
2026-06-21 19:41:25 INFO >> Checking links...
2026-06-21 19:41:28 INFO >> All links are good for ENS15  ...
2026-06-21 19:41:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:41:28 INFO Hourly dataset computed and listing created
2026-06-21 19:41:31 INFO Hourly dataset computed
2026-06-21 19:41:31 INFO ---------->>> Running CHIMERE model from 2020-02-15 00:00:00 to 2020-02-15 01:00:00
2026-06-21 19:41:31 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:41:31 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1
2026-06-21 19:41:31 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/end.2020021415_9_ENS1.nc
2026-06-21 19:41:31 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-06-21 19:41:31 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:41:31 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-06-21 19:41:31 INFO Queuing job for member 1...
2026-06-21 19:41:31 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:41:31 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-06-21 19:41:31 INFO Found: ['4937794']
2026-06-21 19:41:36 INFO [TGCC-IRENE] Submitted job with ID:['4937794']
2026-06-21 19:41:36 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:41:36 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2
2026-06-21 19:41:36 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/end.2020021415_9_ENS2.nc
2026-06-21 19:41:36 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-06-21 19:41:36 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:41:36 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-06-21 19:41:36 INFO Queuing job for member 2...
2026-06-21 19:41:36 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:41:36 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-06-21 19:41:37 INFO Found: ['4937795']
2026-06-21 19:41:42 INFO [TGCC-IRENE] Submitted job with ID:['4937795']
2026-06-21 19:41:42 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:41:42 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3
2026-06-21 19:41:42 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/end.2020021415_9_ENS3.nc
2026-06-21 19:41:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-06-21 19:41:42 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:41:42 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-06-21 19:41:42 INFO Queuing job for member 3...
2026-06-21 19:41:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:41:42 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-06-21 19:41:43 INFO Found: ['4937796']
2026-06-21 19:41:48 INFO [TGCC-IRENE] Submitted job with ID:['4937796']
2026-06-21 19:41:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:41:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4
2026-06-21 19:41:48 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/end.2020021415_9_ENS4.nc
2026-06-21 19:41:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-06-21 19:41:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:41:48 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-06-21 19:41:48 INFO Queuing job for member 4...
2026-06-21 19:41:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:41:48 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-06-21 19:41:49 INFO Found: ['4937797']
2026-06-21 19:41:54 INFO [TGCC-IRENE] Submitted job with ID:['4937797']
2026-06-21 19:41:54 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:41:54 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5
2026-06-21 19:41:54 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/end.2020021415_9_ENS5.nc
2026-06-21 19:41:54 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-06-21 19:41:54 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:41:54 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-06-21 19:41:54 INFO Queuing job for member 5...
2026-06-21 19:41:54 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:41:54 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-06-21 19:41:54 INFO Found: ['4937798']
2026-06-21 19:41:59 INFO [TGCC-IRENE] Submitted job with ID:['4937798']
2026-06-21 19:41:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:41:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6
2026-06-21 19:41:59 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/end.2020021415_9_ENS6.nc
2026-06-21 19:41:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-06-21 19:41:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:41:59 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-06-21 19:41:59 INFO Queuing job for member 6...
2026-06-21 19:41:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:41:59 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-06-21 19:42:00 INFO Found: ['4937800']
2026-06-21 19:42:05 INFO [TGCC-IRENE] Submitted job with ID:['4937800']
2026-06-21 19:42:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:42:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7
2026-06-21 19:42:05 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/end.2020021415_9_ENS7.nc
2026-06-21 19:42:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-06-21 19:42:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:42:05 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-06-21 19:42:05 INFO Queuing job for member 7...
2026-06-21 19:42:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:42:05 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-06-21 19:42:06 INFO Found: ['4937802']
2026-06-21 19:42:11 INFO [TGCC-IRENE] Submitted job with ID:['4937802']
2026-06-21 19:42:11 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:42:11 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8
2026-06-21 19:42:11 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/end.2020021415_9_ENS8.nc
2026-06-21 19:42:11 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-06-21 19:42:11 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:42:11 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-06-21 19:42:11 INFO Queuing job for member 8...
2026-06-21 19:42:11 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:42:11 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-06-21 19:42:12 INFO Found: ['4937804']
2026-06-21 19:42:17 INFO [TGCC-IRENE] Submitted job with ID:['4937804']
2026-06-21 19:42:17 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:42:17 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9
2026-06-21 19:42:17 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/end.2020021415_9_ENS9.nc
2026-06-21 19:42:17 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-06-21 19:42:17 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:42:17 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-06-21 19:42:17 INFO Queuing job for member 9...
2026-06-21 19:42:17 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:42:17 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-06-21 19:42:18 INFO Found: ['4937805']
2026-06-21 19:42:23 INFO [TGCC-IRENE] Submitted job with ID:['4937805']
2026-06-21 19:42:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:42:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10
2026-06-21 19:42:23 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/end.2020021415_9_ENS10.nc
2026-06-21 19:42:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-06-21 19:42:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:42:23 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-06-21 19:42:23 INFO Queuing job for member 10...
2026-06-21 19:42:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:42:23 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-06-21 19:42:24 INFO Found: ['4937807']
2026-06-21 19:42:29 INFO [TGCC-IRENE] Submitted job with ID:['4937807']
2026-06-21 19:42:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:42:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11
2026-06-21 19:42:29 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/end.2020021415_9_ENS11.nc
2026-06-21 19:42:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-06-21 19:42:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:42:29 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-06-21 19:42:29 INFO Queuing job for member 11...
2026-06-21 19:42:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:42:29 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-06-21 19:42:30 INFO Found: ['4937808']
2026-06-21 19:42:35 INFO [TGCC-IRENE] Submitted job with ID:['4937808']
2026-06-21 19:42:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:42:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12
2026-06-21 19:42:35 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/end.2020021415_9_ENS12.nc
2026-06-21 19:42:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-06-21 19:42:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:42:35 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-06-21 19:42:35 INFO Queuing job for member 12...
2026-06-21 19:42:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:42:35 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-06-21 19:42:35 INFO Found: ['4937810']
2026-06-21 19:42:40 INFO [TGCC-IRENE] Submitted job with ID:['4937810']
2026-06-21 19:42:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:42:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13
2026-06-21 19:42:40 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/end.2020021415_9_ENS13.nc
2026-06-21 19:42:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-06-21 19:42:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:42:40 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-06-21 19:42:40 INFO Queuing job for member 13...
2026-06-21 19:42:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:42:40 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-06-21 19:42:41 INFO Found: ['4937811']
2026-06-21 19:42:46 INFO [TGCC-IRENE] Submitted job with ID:['4937811']
2026-06-21 19:42:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:42:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14
2026-06-21 19:42:46 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/end.2020021415_9_ENS14.nc
2026-06-21 19:42:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-06-21 19:42:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:42:46 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-06-21 19:42:46 INFO Queuing job for member 14...
2026-06-21 19:42:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:42:46 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-06-21 19:42:47 INFO Found: ['4937813']
2026-06-21 19:42:52 INFO [TGCC-IRENE] Submitted job with ID:['4937813']
2026-06-21 19:42:52 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:42:52 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15
2026-06-21 19:42:52 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/end.2020021415_9_ENS15.nc
2026-06-21 19:42:52 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-06-21 19:42:52 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:42:52 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-06-21 19:42:52 INFO Queuing job for member 15...
2026-06-21 19:42:52 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:42:52 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-06-21 19:42:53 INFO Found: ['4937814']
2026-06-21 19:42:58 INFO [TGCC-IRENE] Submitted job with ID:['4937814']
2026-06-21 19:42:58 INFO Checking job status ...
2026-06-21 19:42:58 INFO None 4937794: status RUNNING/PENDING
2026-06-21 19:42:58 INFO None 4937795: status RUNNING/PENDING
2026-06-21 19:42:58 INFO None 4937796: status RUNNING/PENDING
2026-06-21 19:42:58 INFO None 4937797: status RUNNING/PENDING
2026-06-21 19:42:58 INFO None 4937798: status RUNNING/PENDING
2026-06-21 19:42:58 INFO None 4937800: status RUNNING/PENDING
2026-06-21 19:42:58 INFO None 4937802: status RUNNING/PENDING
2026-06-21 19:42:58 INFO None 4937804: status RUNNING/PENDING
2026-06-21 19:42:58 INFO None 4937805: status RUNNING/PENDING
2026-06-21 19:42:58 INFO None 4937807: status RUNNING/PENDING
2026-06-21 19:42:58 INFO None 4937808: status RUNNING/PENDING
2026-06-21 19:42:58 INFO None 4937810: status RUNNING/PENDING
2026-06-21 19:42:58 INFO None 4937811: status RUNNING/PENDING
2026-06-21 19:42:58 INFO None 4937813: status RUNNING/PENDING
2026-06-21 19:42:58 INFO None 4937814: status RUNNING/PENDING
2026-06-21 19:42:58 INFO Jobs still running: ['4937794', '4937795', '4937796', '4937797', '4937798', '4937800', '4937802', '4937804', '4937805', '4937807', '4937808', '4937810', '4937811', '4937813', '4937814']. Waiting...
2026-06-21 19:43:13 INFO None 4937794: status RUNNING/PENDING
2026-06-21 19:43:13 INFO None 4937795: status RUNNING/PENDING
2026-06-21 19:43:13 INFO None 4937796: status FINISHED
2026-06-21 19:43:13 INFO None 4937797: status RUNNING/PENDING
2026-06-21 19:43:13 INFO None 4937798: status RUNNING/PENDING
2026-06-21 19:43:13 INFO None 4937800: status RUNNING/PENDING
2026-06-21 19:43:13 INFO None 4937802: status RUNNING/PENDING
2026-06-21 19:43:13 INFO None 4937804: status RUNNING/PENDING
2026-06-21 19:43:13 INFO None 4937805: status RUNNING/PENDING
2026-06-21 19:43:13 INFO None 4937807: status RUNNING/PENDING
2026-06-21 19:43:13 INFO None 4937808: status RUNNING/PENDING
2026-06-21 19:43:13 INFO None 4937810: status RUNNING/PENDING
2026-06-21 19:43:13 INFO None 4937811: status RUNNING/PENDING
2026-06-21 19:43:13 INFO None 4937813: status RUNNING/PENDING
2026-06-21 19:43:13 INFO None 4937814: status RUNNING/PENDING
2026-06-21 19:43:13 INFO Jobs still running: ['4937794', '4937795', '4937797', '4937798', '4937800', '4937802', '4937804', '4937805', '4937807', '4937808', '4937810', '4937811', '4937813', '4937814']. Waiting...
2026-06-21 19:43:28 INFO None 4937794: status RUNNING/PENDING
2026-06-21 19:43:28 INFO None 4937795: status RUNNING/PENDING
2026-06-21 19:43:28 INFO None 4937796: status FINISHED
2026-06-21 19:43:28 INFO None 4937797: status FINISHED
2026-06-21 19:43:28 INFO None 4937798: status RUNNING/PENDING
2026-06-21 19:43:29 INFO None 4937800: status RUNNING/PENDING
2026-06-21 19:43:29 INFO None 4937802: status RUNNING/PENDING
2026-06-21 19:43:29 INFO None 4937804: status RUNNING/PENDING
2026-06-21 19:43:29 INFO None 4937805: status RUNNING/PENDING
2026-06-21 19:43:29 INFO None 4937807: status RUNNING/PENDING
2026-06-21 19:43:29 INFO None 4937808: status RUNNING/PENDING
2026-06-21 19:43:29 INFO None 4937810: status RUNNING/PENDING
2026-06-21 19:43:29 INFO None 4937811: status RUNNING/PENDING
2026-06-21 19:43:29 INFO None 4937813: status RUNNING/PENDING
2026-06-21 19:43:29 INFO None 4937814: status RUNNING/PENDING
2026-06-21 19:43:29 INFO Jobs still running: ['4937794', '4937795', '4937798', '4937800', '4937802', '4937804', '4937805', '4937807', '4937808', '4937810', '4937811', '4937813', '4937814']. Waiting...
2026-06-21 19:43:44 INFO None 4937794: status RUNNING/PENDING
2026-06-21 19:43:44 INFO None 4937795: status RUNNING/PENDING
2026-06-21 19:43:44 INFO None 4937796: status FINISHED
2026-06-21 19:43:44 INFO None 4937797: status FINISHED
2026-06-21 19:43:44 INFO None 4937798: status RUNNING/PENDING
2026-06-21 19:43:44 INFO None 4937800: status FINISHED
2026-06-21 19:43:44 INFO None 4937802: status FINISHED
2026-06-21 19:43:44 INFO None 4937804: status FINISHED
2026-06-21 19:43:44 INFO None 4937805: status FINISHED
2026-06-21 19:43:44 INFO None 4937807: status RUNNING/PENDING
2026-06-21 19:43:44 INFO None 4937808: status RUNNING/PENDING
2026-06-21 19:43:44 INFO None 4937810: status RUNNING/PENDING
2026-06-21 19:43:44 INFO None 4937811: status RUNNING/PENDING
2026-06-21 19:43:44 INFO None 4937813: status RUNNING/PENDING
2026-06-21 19:43:44 INFO None 4937814: status RUNNING/PENDING
2026-06-21 19:43:44 INFO Jobs still running: ['4937794', '4937795', '4937798', '4937807', '4937808', '4937810', '4937811', '4937813', '4937814']. Waiting...
2026-06-21 19:43:59 INFO None 4937794: status RUNNING/PENDING
2026-06-21 19:43:59 INFO None 4937795: status FINISHED
2026-06-21 19:43:59 INFO None 4937796: status FINISHED
2026-06-21 19:43:59 INFO None 4937797: status FINISHED
2026-06-21 19:43:59 INFO None 4937798: status RUNNING/PENDING
2026-06-21 19:43:59 INFO None 4937800: status FINISHED
2026-06-21 19:43:59 INFO None 4937802: status FINISHED
2026-06-21 19:43:59 INFO None 4937804: status FINISHED
2026-06-21 19:43:59 INFO None 4937805: status FINISHED
2026-06-21 19:43:59 INFO None 4937807: status FINISHED
2026-06-21 19:43:59 INFO None 4937808: status RUNNING/PENDING
2026-06-21 19:43:59 INFO None 4937810: status FINISHED
2026-06-21 19:43:59 INFO None 4937811: status RUNNING/PENDING
2026-06-21 19:43:59 INFO None 4937813: status RUNNING/PENDING
2026-06-21 19:43:59 INFO None 4937814: status RUNNING/PENDING
2026-06-21 19:43:59 INFO Jobs still running: ['4937794', '4937798', '4937808', '4937811', '4937813', '4937814']. Waiting...
2026-06-21 19:44:14 INFO None 4937794: status FINISHED
2026-06-21 19:44:14 INFO None 4937795: status FINISHED
2026-06-21 19:44:14 INFO None 4937796: status FINISHED
2026-06-21 19:44:14 INFO None 4937797: status FINISHED
2026-06-21 19:44:14 INFO None 4937798: status FINISHED
2026-06-21 19:44:14 INFO None 4937800: status FINISHED
2026-06-21 19:44:14 INFO None 4937802: status FINISHED
2026-06-21 19:44:14 INFO None 4937804: status FINISHED
2026-06-21 19:44:14 INFO None 4937805: status FINISHED
2026-06-21 19:44:14 INFO None 4937807: status FINISHED
2026-06-21 19:44:14 INFO None 4937808: status FINISHED
2026-06-21 19:44:15 INFO None 4937810: status FINISHED
2026-06-21 19:44:15 INFO None 4937811: status RUNNING/PENDING
2026-06-21 19:44:15 INFO None 4937813: status FINISHED
2026-06-21 19:44:15 INFO None 4937814: status RUNNING/PENDING
2026-06-21 19:44:15 INFO Jobs still running: ['4937811', '4937814']. Waiting...
2026-06-21 19:44:30 INFO None 4937794: status FINISHED
2026-06-21 19:44:30 INFO None 4937795: status FINISHED
2026-06-21 19:44:30 INFO None 4937796: status FINISHED
2026-06-21 19:44:30 INFO None 4937797: status FINISHED
2026-06-21 19:44:30 INFO None 4937798: status FINISHED
2026-06-21 19:44:30 INFO None 4937800: status FINISHED
2026-06-21 19:44:30 INFO None 4937802: status FINISHED
2026-06-21 19:44:30 INFO None 4937804: status FINISHED
2026-06-21 19:44:30 INFO None 4937805: status FINISHED
2026-06-21 19:44:30 INFO None 4937807: status FINISHED
2026-06-21 19:44:30 INFO None 4937808: status FINISHED
2026-06-21 19:44:30 INFO None 4937810: status FINISHED
2026-06-21 19:44:30 INFO None 4937811: status RUNNING/PENDING
2026-06-21 19:44:30 INFO None 4937813: status FINISHED
2026-06-21 19:44:30 INFO None 4937814: status RUNNING/PENDING
2026-06-21 19:44:30 INFO Jobs still running: ['4937811', '4937814']. Waiting...
2026-06-21 19:44:45 INFO None 4937794: status FINISHED
2026-06-21 19:44:45 INFO None 4937795: status FINISHED
2026-06-21 19:44:45 INFO None 4937796: status FINISHED
2026-06-21 19:44:45 INFO None 4937797: status FINISHED
2026-06-21 19:44:45 INFO None 4937798: status FINISHED
2026-06-21 19:44:45 INFO None 4937800: status FINISHED
2026-06-21 19:44:45 INFO None 4937802: status FINISHED
2026-06-21 19:44:45 INFO None 4937804: status FINISHED
2026-06-21 19:44:45 INFO None 4937805: status FINISHED
2026-06-21 19:44:45 INFO None 4937807: status FINISHED
2026-06-21 19:44:45 INFO None 4937808: status FINISHED
2026-06-21 19:44:45 INFO None 4937810: status FINISHED
2026-06-21 19:44:45 INFO None 4937811: status RUNNING/PENDING
2026-06-21 19:44:45 INFO None 4937813: status FINISHED
2026-06-21 19:44:45 INFO None 4937814: status FINISHED
2026-06-21 19:44:45 INFO Jobs still running: ['4937811']. Waiting...
2026-06-21 19:45:00 INFO None 4937794: status FINISHED
2026-06-21 19:45:00 INFO None 4937795: status FINISHED
2026-06-21 19:45:00 INFO None 4937796: status FINISHED
2026-06-21 19:45:00 INFO None 4937797: status FINISHED
2026-06-21 19:45:00 INFO None 4937798: status FINISHED
2026-06-21 19:45:00 INFO None 4937800: status FINISHED
2026-06-21 19:45:00 INFO None 4937802: status FINISHED
2026-06-21 19:45:00 INFO None 4937804: status FINISHED
2026-06-21 19:45:00 INFO None 4937805: status FINISHED
2026-06-21 19:45:00 INFO None 4937807: status FINISHED
2026-06-21 19:45:00 INFO None 4937808: status FINISHED
2026-06-21 19:45:00 INFO None 4937810: status FINISHED
2026-06-21 19:45:00 INFO None 4937811: status RUNNING/PENDING
2026-06-21 19:45:00 INFO None 4937813: status FINISHED
2026-06-21 19:45:00 INFO None 4937814: status FINISHED
2026-06-21 19:45:00 INFO Jobs still running: ['4937811']. Waiting...
2026-06-21 19:45:15 INFO None 4937794: status FINISHED
2026-06-21 19:45:15 INFO None 4937795: status FINISHED
2026-06-21 19:45:15 INFO None 4937796: status FINISHED
2026-06-21 19:45:15 INFO None 4937797: status FINISHED
2026-06-21 19:45:15 INFO None 4937798: status FINISHED
2026-06-21 19:45:16 INFO None 4937800: status FINISHED
2026-06-21 19:45:16 INFO None 4937802: status FINISHED
2026-06-21 19:45:16 INFO None 4937804: status FINISHED
2026-06-21 19:45:16 INFO None 4937805: status FINISHED
2026-06-21 19:45:16 INFO None 4937807: status FINISHED
2026-06-21 19:45:16 INFO None 4937808: status FINISHED
2026-06-21 19:45:16 INFO None 4937810: status FINISHED
2026-06-21 19:45:16 INFO None 4937811: status RUNNING/PENDING
2026-06-21 19:45:16 INFO None 4937813: status FINISHED
2026-06-21 19:45:16 INFO None 4937814: status FINISHED
2026-06-21 19:45:16 INFO Jobs still running: ['4937811']. Waiting...
2026-06-21 19:45:31 INFO None 4937794: status FINISHED
2026-06-21 19:45:31 INFO None 4937795: status FINISHED
2026-06-21 19:45:31 INFO None 4937796: status FINISHED
2026-06-21 19:45:31 INFO None 4937797: status FINISHED
2026-06-21 19:45:31 INFO None 4937798: status FINISHED
2026-06-21 19:45:31 INFO None 4937800: status FINISHED
2026-06-21 19:45:31 INFO None 4937802: status FINISHED
2026-06-21 19:45:31 INFO None 4937804: status FINISHED
2026-06-21 19:45:31 INFO None 4937805: status FINISHED
2026-06-21 19:45:31 INFO None 4937807: status FINISHED
2026-06-21 19:45:31 INFO None 4937808: status FINISHED
2026-06-21 19:45:31 INFO None 4937810: status FINISHED
2026-06-21 19:45:31 INFO None 4937811: status RUNNING/PENDING
2026-06-21 19:45:31 INFO None 4937813: status FINISHED
2026-06-21 19:45:31 INFO None 4937814: status FINISHED
2026-06-21 19:45:31 INFO Jobs still running: ['4937811']. Waiting...
2026-06-21 19:45:46 INFO None 4937794: status FINISHED
2026-06-21 19:45:46 INFO None 4937795: status FINISHED
2026-06-21 19:45:46 INFO None 4937796: status FINISHED
2026-06-21 19:45:46 INFO None 4937797: status FINISHED
2026-06-21 19:45:46 INFO None 4937798: status FINISHED
2026-06-21 19:45:46 INFO None 4937800: status FINISHED
2026-06-21 19:45:46 INFO None 4937802: status FINISHED
2026-06-21 19:45:46 INFO None 4937804: status FINISHED
2026-06-21 19:45:46 INFO None 4937805: status FINISHED
2026-06-21 19:45:46 INFO None 4937807: status FINISHED
2026-06-21 19:45:46 INFO None 4937808: status FINISHED
2026-06-21 19:45:46 INFO None 4937810: status FINISHED
2026-06-21 19:45:46 INFO None 4937811: status FINISHED
2026-06-21 19:45:46 INFO None 4937813: status FINISHED
2026-06-21 19:45:46 INFO None 4937814: status FINISHED
2026-06-21 19:45:46 INFO Jobs ['4937794', '4937795', '4937796', '4937797', '4937798', '4937800', '4937802', '4937804', '4937805', '4937807', '4937808', '4937810', '4937811', '4937813', '4937814'] have finished
2026-06-21 19:45:46 INFO Checking restart files were created ...
2026-06-21 19:45:46 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/end.2020021500_1_ENS1.nc(668832435 bytes)
2026-06-21 19:45:46 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/end.2020021500_1_ENS2.nc(668832435 bytes)
2026-06-21 19:45:46 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/end.2020021500_1_ENS3.nc(668832435 bytes)
2026-06-21 19:45:46 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/end.2020021500_1_ENS4.nc(668832435 bytes)
2026-06-21 19:45:46 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/end.2020021500_1_ENS5.nc(668832435 bytes)
2026-06-21 19:45:46 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/end.2020021500_1_ENS6.nc(668832435 bytes)
2026-06-21 19:45:46 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/end.2020021500_1_ENS7.nc(668832435 bytes)
2026-06-21 19:45:46 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/end.2020021500_1_ENS8.nc(668832435 bytes)
2026-06-21 19:45:46 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/end.2020021500_1_ENS9.nc(668832435 bytes)
2026-06-21 19:45:46 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/end.2020021500_1_ENS10.nc(668832435 bytes)
2026-06-21 19:45:46 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/end.2020021500_1_ENS11.nc(668832435 bytes)
2026-06-21 19:45:46 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/end.2020021500_1_ENS12.nc(668832435 bytes)
2026-06-21 19:45:46 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/end.2020021500_1_ENS13.nc(668832435 bytes)
2026-06-21 19:45:46 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/end.2020021500_1_ENS14.nc(668832435 bytes)
2026-06-21 19:45:46 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/end.2020021500_1_ENS15.nc(668832435 bytes)
2026-06-21 19:45:46 INFO  Run_model() completed successfully.
2026-06-21 19:45:46 INFO [TIME] after_model_set_simulated_time current_time=2020-02-15 00:00:00 simulated_time=2020-02-15 01:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 19:45:46 INFO [TIME] gregorian_conversion simulated_time=2020-02-15 01:00:00 days=153081 seconds=3600
2026-06-21 19:45:46 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-06-21 19:45:46 INFO [TIME] increment current_time 2020-02-15 00:00:00 -> 2020-02-15 01:00:00
2026-06-21 19:45:46 INFO [TIME] after_increment_before_assimilation current_time=2020-02-15 01:00:00 simulated_time=2020-02-15 01:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 19:45:46 INFO ---------->>> Running process_satellite_data()
2026-06-21 19:45:46 INFO [DART] No satellite data found, skipping assimilation
2026-06-21 19:45:46 INFO after_assimilation() skipped
2026-06-21 19:45:46 INFO Next run starts from 2020-02-15 01:00:00
2026-06-21 19:45:46 INFO Cycle is DONE; starting a new loop!
2026-06-21 19:45:46 INFO [TIME] step_end current_time=2020-02-15 01:00:00 simulated_time=2020-02-15 01:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 19:45:46 INFO [TIME] step_start current_time=2020-02-15 01:00:00 simulated_time=2020-02-15 01:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 19:45:46 INFO [TIME] window start=2020-02-15 01:00:00 end=2020-02-15 08:00:00 run_hours=7 has_assimilation=True
2026-06-21 19:45:46 INFO Creating directories and links for ENS1 to run chimere's parallel part
2026-06-21 19:45:46 INFO Linking EMIS ...
2026-06-21 19:45:47 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens01.nc
2026-06-21 19:45:47 INFO >> Checking links...
2026-06-21 19:45:49 INFO >> All links are good for ENS1  ...
2026-06-21 19:45:49 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:45:53 INFO Hourly dataset computed and listing created
2026-06-21 19:46:02 INFO Hourly dataset computed
2026-06-21 19:46:02 INFO Creating directories and links for ENS2 to run chimere's parallel part
2026-06-21 19:46:02 INFO Linking EMIS ...
2026-06-21 19:46:02 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens02.nc
2026-06-21 19:46:02 INFO >> Checking links...
2026-06-21 19:46:05 INFO >> All links are good for ENS2  ...
2026-06-21 19:46:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:46:06 INFO Hourly dataset computed and listing created
2026-06-21 19:46:08 INFO Hourly dataset computed
2026-06-21 19:46:08 INFO Creating directories and links for ENS3 to run chimere's parallel part
2026-06-21 19:46:08 INFO Linking EMIS ...
2026-06-21 19:46:08 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens03.nc
2026-06-21 19:46:08 INFO >> Checking links...
2026-06-21 19:46:10 INFO >> All links are good for ENS3  ...
2026-06-21 19:46:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:46:12 INFO Hourly dataset computed and listing created
2026-06-21 19:46:20 INFO Hourly dataset computed
2026-06-21 19:46:20 INFO Creating directories and links for ENS4 to run chimere's parallel part
2026-06-21 19:46:20 INFO Linking EMIS ...
2026-06-21 19:46:21 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens04.nc
2026-06-21 19:46:21 INFO >> Checking links...
2026-06-21 19:46:23 INFO >> All links are good for ENS4  ...
2026-06-21 19:46:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:46:24 INFO Hourly dataset computed and listing created
2026-06-21 19:46:34 INFO Hourly dataset computed
2026-06-21 19:46:34 INFO Creating directories and links for ENS5 to run chimere's parallel part
2026-06-21 19:46:34 INFO Linking EMIS ...
2026-06-21 19:46:34 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens05.nc
2026-06-21 19:46:34 INFO >> Checking links...
2026-06-21 19:46:37 INFO >> All links are good for ENS5  ...
2026-06-21 19:46:37 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:46:37 INFO Hourly dataset computed and listing created
2026-06-21 19:46:46 INFO Hourly dataset computed
2026-06-21 19:46:46 INFO Creating directories and links for ENS6 to run chimere's parallel part
2026-06-21 19:46:46 INFO Linking EMIS ...
2026-06-21 19:46:46 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens06.nc
2026-06-21 19:46:46 INFO >> Checking links...
2026-06-21 19:46:49 INFO >> All links are good for ENS6  ...
2026-06-21 19:46:49 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:46:50 INFO Hourly dataset computed and listing created
2026-06-21 19:46:58 INFO Hourly dataset computed
2026-06-21 19:46:58 INFO Creating directories and links for ENS7 to run chimere's parallel part
2026-06-21 19:46:58 INFO Linking EMIS ...
2026-06-21 19:46:59 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens07.nc
2026-06-21 19:46:59 INFO >> Checking links...
2026-06-21 19:47:01 INFO >> All links are good for ENS7  ...
2026-06-21 19:47:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:47:02 INFO Hourly dataset computed and listing created
2026-06-21 19:47:10 INFO Hourly dataset computed
2026-06-21 19:47:10 INFO Creating directories and links for ENS8 to run chimere's parallel part
2026-06-21 19:47:10 INFO Linking EMIS ...
2026-06-21 19:47:11 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens08.nc
2026-06-21 19:47:11 INFO >> Checking links...
2026-06-21 19:47:13 INFO >> All links are good for ENS8  ...
2026-06-21 19:47:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:47:14 INFO Hourly dataset computed and listing created
2026-06-21 19:47:23 INFO Hourly dataset computed
2026-06-21 19:47:23 INFO Creating directories and links for ENS9 to run chimere's parallel part
2026-06-21 19:47:23 INFO Linking EMIS ...
2026-06-21 19:47:23 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens09.nc
2026-06-21 19:47:23 INFO >> Checking links...
2026-06-21 19:47:25 INFO >> All links are good for ENS9  ...
2026-06-21 19:47:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:47:26 INFO Hourly dataset computed and listing created
2026-06-21 19:47:35 INFO Hourly dataset computed
2026-06-21 19:47:35 INFO Creating directories and links for ENS10 to run chimere's parallel part
2026-06-21 19:47:35 INFO Linking EMIS ...
2026-06-21 19:47:36 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens10.nc
2026-06-21 19:47:36 INFO >> Checking links...
2026-06-21 19:47:38 INFO >> All links are good for ENS10  ...
2026-06-21 19:47:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:47:39 INFO Hourly dataset computed and listing created
2026-06-21 19:47:47 INFO Hourly dataset computed
2026-06-21 19:47:47 INFO Creating directories and links for ENS11 to run chimere's parallel part
2026-06-21 19:47:47 INFO Linking EMIS ...
2026-06-21 19:47:48 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens11.nc
2026-06-21 19:47:48 INFO >> Checking links...
2026-06-21 19:47:50 INFO >> All links are good for ENS11  ...
2026-06-21 19:47:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:47:51 INFO Hourly dataset computed and listing created
2026-06-21 19:48:00 INFO Hourly dataset computed
2026-06-21 19:48:00 INFO Creating directories and links for ENS12 to run chimere's parallel part
2026-06-21 19:48:00 INFO Linking EMIS ...
2026-06-21 19:48:01 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens12.nc
2026-06-21 19:48:01 INFO >> Checking links...
2026-06-21 19:48:03 INFO >> All links are good for ENS12  ...
2026-06-21 19:48:03 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:48:04 INFO Hourly dataset computed and listing created
2026-06-21 19:48:13 INFO Hourly dataset computed
2026-06-21 19:48:13 INFO Creating directories and links for ENS13 to run chimere's parallel part
2026-06-21 19:48:13 INFO Linking EMIS ...
2026-06-21 19:48:13 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens13.nc
2026-06-21 19:48:13 INFO >> Checking links...
2026-06-21 19:48:16 INFO >> All links are good for ENS13  ...
2026-06-21 19:48:16 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:48:17 INFO Hourly dataset computed and listing created
2026-06-21 19:48:25 INFO Hourly dataset computed
2026-06-21 19:48:25 INFO Creating directories and links for ENS14 to run chimere's parallel part
2026-06-21 19:48:25 INFO Linking EMIS ...
2026-06-21 19:48:25 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens14.nc
2026-06-21 19:48:25 INFO >> Checking links...
2026-06-21 19:48:28 INFO >> All links are good for ENS14  ...
2026-06-21 19:48:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:48:29 INFO Hourly dataset computed and listing created
2026-06-21 19:48:37 INFO Hourly dataset computed
2026-06-21 19:48:37 INFO Creating directories and links for ENS15 to run chimere's parallel part
2026-06-21 19:48:37 INFO Linking EMIS ...
2026-06-21 19:48:38 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens15.nc
2026-06-21 19:48:38 INFO >> Checking links...
2026-06-21 19:48:40 INFO >> All links are good for ENS15  ...
2026-06-21 19:48:40 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:48:41 INFO Hourly dataset computed and listing created
2026-06-21 19:48:49 INFO Hourly dataset computed
2026-06-21 19:48:49 INFO ---------->>> Running CHIMERE model from 2020-02-15 01:00:00 to 2020-02-15 08:00:00
2026-06-21 19:48:49 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:48:49 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1
2026-06-21 19:48:49 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/end.2020021500_1_ENS1.nc
2026-06-21 19:48:49 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-06-21 19:48:49 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:48:49 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-06-21 19:48:49 INFO Queuing job for member 1...
2026-06-21 19:48:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:48:49 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-06-21 19:48:50 INFO Found: ['4937865']
2026-06-21 19:48:55 INFO [TGCC-IRENE] Submitted job with ID:['4937865']
2026-06-21 19:48:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:48:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2
2026-06-21 19:48:55 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/end.2020021500_1_ENS2.nc
2026-06-21 19:48:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-06-21 19:48:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:48:55 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-06-21 19:48:55 INFO Queuing job for member 2...
2026-06-21 19:48:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:48:55 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-06-21 19:48:56 INFO Found: ['4937868']
2026-06-21 19:49:01 INFO [TGCC-IRENE] Submitted job with ID:['4937868']
2026-06-21 19:49:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:49:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3
2026-06-21 19:49:01 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/end.2020021500_1_ENS3.nc
2026-06-21 19:49:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-06-21 19:49:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:49:01 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-06-21 19:49:01 INFO Queuing job for member 3...
2026-06-21 19:49:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:49:01 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-06-21 19:49:01 INFO Found: ['4937870']
2026-06-21 19:49:06 INFO [TGCC-IRENE] Submitted job with ID:['4937870']
2026-06-21 19:49:06 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:49:06 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4
2026-06-21 19:49:06 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/end.2020021500_1_ENS4.nc
2026-06-21 19:49:06 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-06-21 19:49:06 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:49:06 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-06-21 19:49:06 INFO Queuing job for member 4...
2026-06-21 19:49:06 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:49:06 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-06-21 19:49:07 INFO Found: ['4937872']
2026-06-21 19:49:12 INFO [TGCC-IRENE] Submitted job with ID:['4937872']
2026-06-21 19:49:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:49:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5
2026-06-21 19:49:12 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/end.2020021500_1_ENS5.nc
2026-06-21 19:49:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-06-21 19:49:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:49:12 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-06-21 19:49:12 INFO Queuing job for member 5...
2026-06-21 19:49:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:49:12 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-06-21 19:49:13 INFO Found: ['4937873']
2026-06-21 19:49:18 INFO [TGCC-IRENE] Submitted job with ID:['4937873']
2026-06-21 19:49:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:49:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6
2026-06-21 19:49:18 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/end.2020021500_1_ENS6.nc
2026-06-21 19:49:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-06-21 19:49:18 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:49:18 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-06-21 19:49:18 INFO Queuing job for member 6...
2026-06-21 19:49:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:49:18 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-06-21 19:49:19 INFO Found: ['4937875']
2026-06-21 19:49:24 INFO [TGCC-IRENE] Submitted job with ID:['4937875']
2026-06-21 19:49:24 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:49:24 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7
2026-06-21 19:49:24 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/end.2020021500_1_ENS7.nc
2026-06-21 19:49:24 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-06-21 19:49:24 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:49:24 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-06-21 19:49:24 INFO Queuing job for member 7...
2026-06-21 19:49:24 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:49:24 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-06-21 19:49:24 INFO Found: ['4937876']
2026-06-21 19:49:29 INFO [TGCC-IRENE] Submitted job with ID:['4937876']
2026-06-21 19:49:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:49:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8
2026-06-21 19:49:29 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/end.2020021500_1_ENS8.nc
2026-06-21 19:49:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-06-21 19:49:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:49:29 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-06-21 19:49:29 INFO Queuing job for member 8...
2026-06-21 19:49:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:49:29 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-06-21 19:49:30 INFO Found: ['4937878']
2026-06-21 19:49:35 INFO [TGCC-IRENE] Submitted job with ID:['4937878']
2026-06-21 19:49:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:49:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9
2026-06-21 19:49:35 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/end.2020021500_1_ENS9.nc
2026-06-21 19:49:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-06-21 19:49:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:49:35 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-06-21 19:49:35 INFO Queuing job for member 9...
2026-06-21 19:49:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:49:35 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-06-21 19:49:36 INFO Found: ['4937880']
2026-06-21 19:49:41 INFO [TGCC-IRENE] Submitted job with ID:['4937880']
2026-06-21 19:49:41 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:49:41 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10
2026-06-21 19:49:41 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/end.2020021500_1_ENS10.nc
2026-06-21 19:49:41 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-06-21 19:49:41 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:49:41 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-06-21 19:49:41 INFO Queuing job for member 10...
2026-06-21 19:49:41 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:49:41 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-06-21 19:49:42 INFO Found: ['4937881']
2026-06-21 19:49:47 INFO [TGCC-IRENE] Submitted job with ID:['4937881']
2026-06-21 19:49:47 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:49:47 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11
2026-06-21 19:49:47 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/end.2020021500_1_ENS11.nc
2026-06-21 19:49:47 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-06-21 19:49:47 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:49:47 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-06-21 19:49:47 INFO Queuing job for member 11...
2026-06-21 19:49:47 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:49:47 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-06-21 19:49:48 INFO Found: ['4937883']
2026-06-21 19:49:53 INFO [TGCC-IRENE] Submitted job with ID:['4937883']
2026-06-21 19:49:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:49:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12
2026-06-21 19:49:53 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/end.2020021500_1_ENS12.nc
2026-06-21 19:49:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-06-21 19:49:54 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:49:54 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-06-21 19:49:54 INFO Queuing job for member 12...
2026-06-21 19:49:54 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:49:54 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-06-21 19:49:54 INFO Found: ['4937884']
2026-06-21 19:49:59 INFO [TGCC-IRENE] Submitted job with ID:['4937884']
2026-06-21 19:49:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:49:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13
2026-06-21 19:49:59 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/end.2020021500_1_ENS13.nc
2026-06-21 19:49:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-06-21 19:49:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:49:59 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-06-21 19:49:59 INFO Queuing job for member 13...
2026-06-21 19:49:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:49:59 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-06-21 19:50:00 INFO Found: ['4937887']
2026-06-21 19:50:05 INFO [TGCC-IRENE] Submitted job with ID:['4937887']
2026-06-21 19:50:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:50:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14
2026-06-21 19:50:05 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/end.2020021500_1_ENS14.nc
2026-06-21 19:50:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-06-21 19:50:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:50:05 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-06-21 19:50:05 INFO Queuing job for member 14...
2026-06-21 19:50:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:50:05 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-06-21 19:50:06 INFO Found: ['4937889']
2026-06-21 19:50:11 INFO [TGCC-IRENE] Submitted job with ID:['4937889']
2026-06-21 19:50:11 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:50:11 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15
2026-06-21 19:50:11 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/end.2020021500_1_ENS15.nc
2026-06-21 19:50:11 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-06-21 19:50:11 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:50:11 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-06-21 19:50:11 INFO Queuing job for member 15...
2026-06-21 19:50:11 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:50:11 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-06-21 19:50:12 INFO Found: ['4937891']
2026-06-21 19:50:17 INFO [TGCC-IRENE] Submitted job with ID:['4937891']
2026-06-21 19:50:17 INFO Checking job status ...
2026-06-21 19:50:17 INFO None 4937865: status RUNNING/PENDING
2026-06-21 19:50:17 INFO None 4937868: status RUNNING/PENDING
2026-06-21 19:50:17 INFO None 4937870: status RUNNING/PENDING
2026-06-21 19:50:17 INFO None 4937872: status RUNNING/PENDING
2026-06-21 19:50:17 INFO None 4937873: status RUNNING/PENDING
2026-06-21 19:50:17 INFO None 4937875: status RUNNING/PENDING
2026-06-21 19:50:17 INFO None 4937876: status RUNNING/PENDING
2026-06-21 19:50:17 INFO None 4937878: status RUNNING/PENDING
2026-06-21 19:50:17 INFO None 4937880: status RUNNING/PENDING
2026-06-21 19:50:17 INFO None 4937881: status RUNNING/PENDING
2026-06-21 19:50:17 INFO None 4937883: status RUNNING/PENDING
2026-06-21 19:50:17 INFO None 4937884: status RUNNING/PENDING
2026-06-21 19:50:17 INFO None 4937887: status RUNNING/PENDING
2026-06-21 19:50:17 INFO None 4937889: status RUNNING/PENDING
2026-06-21 19:50:17 INFO None 4937891: status RUNNING/PENDING
2026-06-21 19:50:17 INFO Jobs still running: ['4937865', '4937868', '4937870', '4937872', '4937873', '4937875', '4937876', '4937878', '4937880', '4937881', '4937883', '4937884', '4937887', '4937889', '4937891']. Waiting...
2026-06-21 19:50:32 INFO None 4937865: status RUNNING/PENDING
2026-06-21 19:50:32 INFO None 4937868: status RUNNING/PENDING
2026-06-21 19:50:32 INFO None 4937870: status RUNNING/PENDING
2026-06-21 19:50:32 INFO None 4937872: status RUNNING/PENDING
2026-06-21 19:50:32 INFO None 4937873: status RUNNING/PENDING
2026-06-21 19:50:32 INFO None 4937875: status RUNNING/PENDING
2026-06-21 19:50:32 INFO None 4937876: status RUNNING/PENDING
2026-06-21 19:50:32 INFO None 4937878: status RUNNING/PENDING
2026-06-21 19:50:32 INFO None 4937880: status RUNNING/PENDING
2026-06-21 19:50:32 INFO None 4937881: status RUNNING/PENDING
2026-06-21 19:50:32 INFO None 4937883: status RUNNING/PENDING
2026-06-21 19:50:32 INFO None 4937884: status RUNNING/PENDING
2026-06-21 19:50:32 INFO None 4937887: status RUNNING/PENDING
2026-06-21 19:50:32 INFO None 4937889: status RUNNING/PENDING
2026-06-21 19:50:32 INFO None 4937891: status RUNNING/PENDING
2026-06-21 19:50:33 INFO Jobs still running: ['4937865', '4937868', '4937870', '4937872', '4937873', '4937875', '4937876', '4937878', '4937880', '4937881', '4937883', '4937884', '4937887', '4937889', '4937891']. Waiting...
2026-06-21 19:50:48 INFO None 4937865: status RUNNING/PENDING
2026-06-21 19:50:48 INFO None 4937868: status RUNNING/PENDING
2026-06-21 19:50:48 INFO None 4937870: status RUNNING/PENDING
2026-06-21 19:50:48 INFO None 4937872: status RUNNING/PENDING
2026-06-21 19:50:48 INFO None 4937873: status RUNNING/PENDING
2026-06-21 19:50:48 INFO None 4937875: status RUNNING/PENDING
2026-06-21 19:50:48 INFO None 4937876: status RUNNING/PENDING
2026-06-21 19:50:48 INFO None 4937878: status RUNNING/PENDING
2026-06-21 19:50:48 INFO None 4937880: status RUNNING/PENDING
2026-06-21 19:50:48 INFO None 4937881: status RUNNING/PENDING
2026-06-21 19:50:48 INFO None 4937883: status RUNNING/PENDING
2026-06-21 19:50:48 INFO None 4937884: status RUNNING/PENDING
2026-06-21 19:50:48 INFO None 4937887: status RUNNING/PENDING
2026-06-21 19:50:48 INFO None 4937889: status RUNNING/PENDING
2026-06-21 19:50:48 INFO None 4937891: status RUNNING/PENDING
2026-06-21 19:50:48 INFO Jobs still running: ['4937865', '4937868', '4937870', '4937872', '4937873', '4937875', '4937876', '4937878', '4937880', '4937881', '4937883', '4937884', '4937887', '4937889', '4937891']. Waiting...
2026-06-21 19:51:03 INFO None 4937865: status RUNNING/PENDING
2026-06-21 19:51:03 INFO None 4937868: status RUNNING/PENDING
2026-06-21 19:51:03 INFO None 4937870: status RUNNING/PENDING
2026-06-21 19:51:03 INFO None 4937872: status RUNNING/PENDING
2026-06-21 19:51:03 INFO None 4937873: status RUNNING/PENDING
2026-06-21 19:51:03 INFO None 4937875: status RUNNING/PENDING
2026-06-21 19:51:03 INFO None 4937876: status RUNNING/PENDING
2026-06-21 19:51:03 INFO None 4937878: status RUNNING/PENDING
2026-06-21 19:51:03 INFO None 4937880: status RUNNING/PENDING
2026-06-21 19:51:03 INFO None 4937881: status RUNNING/PENDING
2026-06-21 19:51:03 INFO None 4937883: status RUNNING/PENDING
2026-06-21 19:51:03 INFO None 4937884: status RUNNING/PENDING
2026-06-21 19:51:03 INFO None 4937887: status RUNNING/PENDING
2026-06-21 19:51:03 INFO None 4937889: status RUNNING/PENDING
2026-06-21 19:51:03 INFO None 4937891: status RUNNING/PENDING
2026-06-21 19:51:03 INFO Jobs still running: ['4937865', '4937868', '4937870', '4937872', '4937873', '4937875', '4937876', '4937878', '4937880', '4937881', '4937883', '4937884', '4937887', '4937889', '4937891']. Waiting...
2026-06-21 19:51:18 INFO None 4937865: status RUNNING/PENDING
2026-06-21 19:51:18 INFO None 4937868: status RUNNING/PENDING
2026-06-21 19:51:18 INFO None 4937870: status RUNNING/PENDING
2026-06-21 19:51:18 INFO None 4937872: status RUNNING/PENDING
2026-06-21 19:51:18 INFO None 4937873: status RUNNING/PENDING
2026-06-21 19:51:18 INFO None 4937875: status RUNNING/PENDING
2026-06-21 19:51:18 INFO None 4937876: status RUNNING/PENDING
2026-06-21 19:51:18 INFO None 4937878: status RUNNING/PENDING
2026-06-21 19:51:18 INFO None 4937880: status RUNNING/PENDING
2026-06-21 19:51:18 INFO None 4937881: status RUNNING/PENDING
2026-06-21 19:51:18 INFO None 4937883: status RUNNING/PENDING
2026-06-21 19:51:18 INFO None 4937884: status RUNNING/PENDING
2026-06-21 19:51:18 INFO None 4937887: status RUNNING/PENDING
2026-06-21 19:51:18 INFO None 4937889: status RUNNING/PENDING
2026-06-21 19:51:18 INFO None 4937891: status RUNNING/PENDING
2026-06-21 19:51:18 INFO Jobs still running: ['4937865', '4937868', '4937870', '4937872', '4937873', '4937875', '4937876', '4937878', '4937880', '4937881', '4937883', '4937884', '4937887', '4937889', '4937891']. Waiting...
2026-06-21 19:51:33 INFO None 4937865: status RUNNING/PENDING
2026-06-21 19:51:33 INFO None 4937868: status RUNNING/PENDING
2026-06-21 19:51:34 INFO None 4937870: status RUNNING/PENDING
2026-06-21 19:51:34 INFO None 4937872: status RUNNING/PENDING
2026-06-21 19:51:34 INFO None 4937873: status RUNNING/PENDING
2026-06-21 19:51:34 INFO None 4937875: status RUNNING/PENDING
2026-06-21 19:51:34 INFO None 4937876: status RUNNING/PENDING
2026-06-21 19:51:34 INFO None 4937878: status RUNNING/PENDING
2026-06-21 19:51:34 INFO None 4937880: status RUNNING/PENDING
2026-06-21 19:51:34 INFO None 4937881: status RUNNING/PENDING
2026-06-21 19:51:34 INFO None 4937883: status RUNNING/PENDING
2026-06-21 19:51:34 INFO None 4937884: status RUNNING/PENDING
2026-06-21 19:51:34 INFO None 4937887: status RUNNING/PENDING
2026-06-21 19:51:34 INFO None 4937889: status RUNNING/PENDING
2026-06-21 19:51:34 INFO None 4937891: status RUNNING/PENDING
2026-06-21 19:51:34 INFO Jobs still running: ['4937865', '4937868', '4937870', '4937872', '4937873', '4937875', '4937876', '4937878', '4937880', '4937881', '4937883', '4937884', '4937887', '4937889', '4937891']. Waiting...
2026-06-21 19:51:49 INFO None 4937865: status RUNNING/PENDING
2026-06-21 19:51:49 INFO None 4937868: status RUNNING/PENDING
2026-06-21 19:51:49 INFO None 4937870: status RUNNING/PENDING
2026-06-21 19:51:49 INFO None 4937872: status RUNNING/PENDING
2026-06-21 19:51:49 INFO None 4937873: status RUNNING/PENDING
2026-06-21 19:51:49 INFO None 4937875: status RUNNING/PENDING
2026-06-21 19:51:49 INFO None 4937876: status RUNNING/PENDING
2026-06-21 19:51:49 INFO None 4937878: status RUNNING/PENDING
2026-06-21 19:51:49 INFO None 4937880: status RUNNING/PENDING
2026-06-21 19:51:49 INFO None 4937881: status RUNNING/PENDING
2026-06-21 19:51:49 INFO None 4937883: status RUNNING/PENDING
2026-06-21 19:51:49 INFO None 4937884: status RUNNING/PENDING
2026-06-21 19:51:49 INFO None 4937887: status RUNNING/PENDING
2026-06-21 19:51:49 INFO None 4937889: status RUNNING/PENDING
2026-06-21 19:51:49 INFO None 4937891: status RUNNING/PENDING
2026-06-21 19:51:49 INFO Jobs still running: ['4937865', '4937868', '4937870', '4937872', '4937873', '4937875', '4937876', '4937878', '4937880', '4937881', '4937883', '4937884', '4937887', '4937889', '4937891']. Waiting...
2026-06-21 19:52:04 INFO None 4937865: status RUNNING/PENDING
2026-06-21 19:52:04 INFO None 4937868: status RUNNING/PENDING
2026-06-21 19:52:04 INFO None 4937870: status RUNNING/PENDING
2026-06-21 19:52:04 INFO None 4937872: status RUNNING/PENDING
2026-06-21 19:52:04 INFO None 4937873: status RUNNING/PENDING
2026-06-21 19:52:04 INFO None 4937875: status RUNNING/PENDING
2026-06-21 19:52:04 INFO None 4937876: status RUNNING/PENDING
2026-06-21 19:52:04 INFO None 4937878: status RUNNING/PENDING
2026-06-21 19:52:04 INFO None 4937880: status RUNNING/PENDING
2026-06-21 19:52:04 INFO None 4937881: status RUNNING/PENDING
2026-06-21 19:52:04 INFO None 4937883: status RUNNING/PENDING
2026-06-21 19:52:04 INFO None 4937884: status RUNNING/PENDING
2026-06-21 19:52:04 INFO None 4937887: status RUNNING/PENDING
2026-06-21 19:52:04 INFO None 4937889: status RUNNING/PENDING
2026-06-21 19:52:04 INFO None 4937891: status RUNNING/PENDING
2026-06-21 19:52:04 INFO Jobs still running: ['4937865', '4937868', '4937870', '4937872', '4937873', '4937875', '4937876', '4937878', '4937880', '4937881', '4937883', '4937884', '4937887', '4937889', '4937891']. Waiting...
2026-06-21 19:52:19 INFO None 4937865: status RUNNING/PENDING
2026-06-21 19:52:19 INFO None 4937868: status RUNNING/PENDING
2026-06-21 19:52:19 INFO None 4937870: status RUNNING/PENDING
2026-06-21 19:52:19 INFO None 4937872: status RUNNING/PENDING
2026-06-21 19:52:19 INFO None 4937873: status RUNNING/PENDING
2026-06-21 19:52:19 INFO None 4937875: status RUNNING/PENDING
2026-06-21 19:52:19 INFO None 4937876: status RUNNING/PENDING
2026-06-21 19:52:20 INFO None 4937878: status RUNNING/PENDING
2026-06-21 19:52:20 INFO None 4937880: status RUNNING/PENDING
2026-06-21 19:52:20 INFO None 4937881: status RUNNING/PENDING
2026-06-21 19:52:20 INFO None 4937883: status RUNNING/PENDING
2026-06-21 19:52:20 INFO None 4937884: status RUNNING/PENDING
2026-06-21 19:52:20 INFO None 4937887: status RUNNING/PENDING
2026-06-21 19:52:20 INFO None 4937889: status RUNNING/PENDING
2026-06-21 19:52:20 INFO None 4937891: status RUNNING/PENDING
2026-06-21 19:52:20 INFO Jobs still running: ['4937865', '4937868', '4937870', '4937872', '4937873', '4937875', '4937876', '4937878', '4937880', '4937881', '4937883', '4937884', '4937887', '4937889', '4937891']. Waiting...
2026-06-21 19:52:35 INFO None 4937865: status RUNNING/PENDING
2026-06-21 19:52:35 INFO None 4937868: status RUNNING/PENDING
2026-06-21 19:52:35 INFO None 4937870: status RUNNING/PENDING
2026-06-21 19:52:35 INFO None 4937872: status RUNNING/PENDING
2026-06-21 19:52:35 INFO None 4937873: status RUNNING/PENDING
2026-06-21 19:52:35 INFO None 4937875: status RUNNING/PENDING
2026-06-21 19:52:35 INFO None 4937876: status RUNNING/PENDING
2026-06-21 19:52:35 INFO None 4937878: status RUNNING/PENDING
2026-06-21 19:52:35 INFO None 4937880: status RUNNING/PENDING
2026-06-21 19:52:35 INFO None 4937881: status RUNNING/PENDING
2026-06-21 19:52:35 INFO None 4937883: status RUNNING/PENDING
2026-06-21 19:52:35 INFO None 4937884: status RUNNING/PENDING
2026-06-21 19:52:35 INFO None 4937887: status RUNNING/PENDING
2026-06-21 19:52:35 INFO None 4937889: status RUNNING/PENDING
2026-06-21 19:52:35 INFO None 4937891: status RUNNING/PENDING
2026-06-21 19:52:35 INFO Jobs still running: ['4937865', '4937868', '4937870', '4937872', '4937873', '4937875', '4937876', '4937878', '4937880', '4937881', '4937883', '4937884', '4937887', '4937889', '4937891']. Waiting...
2026-06-21 19:52:50 INFO None 4937865: status RUNNING/PENDING
2026-06-21 19:52:50 INFO None 4937868: status RUNNING/PENDING
2026-06-21 19:52:50 INFO None 4937870: status RUNNING/PENDING
2026-06-21 19:52:50 INFO None 4937872: status RUNNING/PENDING
2026-06-21 19:52:50 INFO None 4937873: status RUNNING/PENDING
2026-06-21 19:52:50 INFO None 4937875: status RUNNING/PENDING
2026-06-21 19:52:50 INFO None 4937876: status RUNNING/PENDING
2026-06-21 19:52:50 INFO None 4937878: status RUNNING/PENDING
2026-06-21 19:52:50 INFO None 4937880: status RUNNING/PENDING
2026-06-21 19:52:50 INFO None 4937881: status RUNNING/PENDING
2026-06-21 19:52:50 INFO None 4937883: status RUNNING/PENDING
2026-06-21 19:52:50 INFO None 4937884: status RUNNING/PENDING
2026-06-21 19:52:50 INFO None 4937887: status RUNNING/PENDING
2026-06-21 19:52:50 INFO None 4937889: status RUNNING/PENDING
2026-06-21 19:52:50 INFO None 4937891: status RUNNING/PENDING
2026-06-21 19:52:50 INFO Jobs still running: ['4937865', '4937868', '4937870', '4937872', '4937873', '4937875', '4937876', '4937878', '4937880', '4937881', '4937883', '4937884', '4937887', '4937889', '4937891']. Waiting...
2026-06-21 19:53:05 INFO None 4937865: status RUNNING/PENDING
2026-06-21 19:53:05 INFO None 4937868: status RUNNING/PENDING
2026-06-21 19:53:05 INFO None 4937870: status RUNNING/PENDING
2026-06-21 19:53:05 INFO None 4937872: status RUNNING/PENDING
2026-06-21 19:53:05 INFO None 4937873: status RUNNING/PENDING
2026-06-21 19:53:05 INFO None 4937875: status RUNNING/PENDING
2026-06-21 19:53:05 INFO None 4937876: status RUNNING/PENDING
2026-06-21 19:53:05 INFO None 4937878: status RUNNING/PENDING
2026-06-21 19:53:05 INFO None 4937880: status RUNNING/PENDING
2026-06-21 19:53:05 INFO None 4937881: status RUNNING/PENDING
2026-06-21 19:53:05 INFO None 4937883: status RUNNING/PENDING
2026-06-21 19:53:05 INFO None 4937884: status RUNNING/PENDING
2026-06-21 19:53:05 INFO None 4937887: status RUNNING/PENDING
2026-06-21 19:53:06 INFO None 4937889: status RUNNING/PENDING
2026-06-21 19:53:06 INFO None 4937891: status RUNNING/PENDING
2026-06-21 19:53:06 INFO Jobs still running: ['4937865', '4937868', '4937870', '4937872', '4937873', '4937875', '4937876', '4937878', '4937880', '4937881', '4937883', '4937884', '4937887', '4937889', '4937891']. Waiting...
2026-06-21 19:53:21 INFO None 4937865: status RUNNING/PENDING
2026-06-21 19:53:21 INFO None 4937868: status RUNNING/PENDING
2026-06-21 19:53:21 INFO None 4937870: status RUNNING/PENDING
2026-06-21 19:53:21 INFO None 4937872: status RUNNING/PENDING
2026-06-21 19:53:21 INFO None 4937873: status RUNNING/PENDING
2026-06-21 19:53:21 INFO None 4937875: status RUNNING/PENDING
2026-06-21 19:53:21 INFO None 4937876: status RUNNING/PENDING
2026-06-21 19:53:21 INFO None 4937878: status RUNNING/PENDING
2026-06-21 19:53:21 INFO None 4937880: status RUNNING/PENDING
2026-06-21 19:53:21 INFO None 4937881: status RUNNING/PENDING
2026-06-21 19:53:21 INFO None 4937883: status RUNNING/PENDING
2026-06-21 19:53:21 INFO None 4937884: status RUNNING/PENDING
2026-06-21 19:53:21 INFO None 4937887: status RUNNING/PENDING
2026-06-21 19:53:21 INFO None 4937889: status RUNNING/PENDING
2026-06-21 19:53:21 INFO None 4937891: status RUNNING/PENDING
2026-06-21 19:53:21 INFO Jobs still running: ['4937865', '4937868', '4937870', '4937872', '4937873', '4937875', '4937876', '4937878', '4937880', '4937881', '4937883', '4937884', '4937887', '4937889', '4937891']. Waiting...
2026-06-21 19:53:36 INFO None 4937865: status RUNNING/PENDING
2026-06-21 19:53:36 INFO None 4937868: status RUNNING/PENDING
2026-06-21 19:53:36 INFO None 4937870: status RUNNING/PENDING
2026-06-21 19:53:36 INFO None 4937872: status RUNNING/PENDING
2026-06-21 19:53:36 INFO None 4937873: status RUNNING/PENDING
2026-06-21 19:53:36 INFO None 4937875: status RUNNING/PENDING
2026-06-21 19:53:36 INFO None 4937876: status RUNNING/PENDING
2026-06-21 19:53:36 INFO None 4937878: status RUNNING/PENDING
2026-06-21 19:53:36 INFO None 4937880: status RUNNING/PENDING
2026-06-21 19:53:36 INFO None 4937881: status RUNNING/PENDING
2026-06-21 19:53:36 INFO None 4937883: status RUNNING/PENDING
2026-06-21 19:53:36 INFO None 4937884: status RUNNING/PENDING
2026-06-21 19:53:36 INFO None 4937887: status RUNNING/PENDING
2026-06-21 19:53:36 INFO None 4937889: status RUNNING/PENDING
2026-06-21 19:53:36 INFO None 4937891: status RUNNING/PENDING
2026-06-21 19:53:36 INFO Jobs still running: ['4937865', '4937868', '4937870', '4937872', '4937873', '4937875', '4937876', '4937878', '4937880', '4937881', '4937883', '4937884', '4937887', '4937889', '4937891']. Waiting...
2026-06-21 19:53:51 INFO None 4937865: status RUNNING/PENDING
2026-06-21 19:53:51 INFO None 4937868: status RUNNING/PENDING
2026-06-21 19:53:51 INFO None 4937870: status RUNNING/PENDING
2026-06-21 19:53:51 INFO None 4937872: status RUNNING/PENDING
2026-06-21 19:53:51 INFO None 4937873: status RUNNING/PENDING
2026-06-21 19:53:51 INFO None 4937875: status RUNNING/PENDING
2026-06-21 19:53:51 INFO None 4937876: status RUNNING/PENDING
2026-06-21 19:53:51 INFO None 4937878: status RUNNING/PENDING
2026-06-21 19:53:51 INFO None 4937880: status RUNNING/PENDING
2026-06-21 19:53:51 INFO None 4937881: status RUNNING/PENDING
2026-06-21 19:53:51 INFO None 4937883: status RUNNING/PENDING
2026-06-21 19:53:51 INFO None 4937884: status RUNNING/PENDING
2026-06-21 19:53:51 INFO None 4937887: status RUNNING/PENDING
2026-06-21 19:53:52 INFO None 4937889: status RUNNING/PENDING
2026-06-21 19:53:52 INFO None 4937891: status RUNNING/PENDING
2026-06-21 19:53:52 INFO Jobs still running: ['4937865', '4937868', '4937870', '4937872', '4937873', '4937875', '4937876', '4937878', '4937880', '4937881', '4937883', '4937884', '4937887', '4937889', '4937891']. Waiting...
2026-06-21 19:54:07 INFO None 4937865: status RUNNING/PENDING
2026-06-21 19:54:07 INFO None 4937868: status RUNNING/PENDING
2026-06-21 19:54:07 INFO None 4937870: status RUNNING/PENDING
2026-06-21 19:54:07 INFO None 4937872: status RUNNING/PENDING
2026-06-21 19:54:07 INFO None 4937873: status RUNNING/PENDING
2026-06-21 19:54:07 INFO None 4937875: status RUNNING/PENDING
2026-06-21 19:54:07 INFO None 4937876: status RUNNING/PENDING
2026-06-21 19:54:07 INFO None 4937878: status RUNNING/PENDING
2026-06-21 19:54:07 INFO None 4937880: status RUNNING/PENDING
2026-06-21 19:54:07 INFO None 4937881: status RUNNING/PENDING
2026-06-21 19:54:07 INFO None 4937883: status RUNNING/PENDING
2026-06-21 19:54:07 INFO None 4937884: status RUNNING/PENDING
2026-06-21 19:54:07 INFO None 4937887: status RUNNING/PENDING
2026-06-21 19:54:07 INFO None 4937889: status RUNNING/PENDING
2026-06-21 19:54:07 INFO None 4937891: status RUNNING/PENDING
2026-06-21 19:54:07 INFO Jobs still running: ['4937865', '4937868', '4937870', '4937872', '4937873', '4937875', '4937876', '4937878', '4937880', '4937881', '4937883', '4937884', '4937887', '4937889', '4937891']. Waiting...
2026-06-21 19:54:22 INFO None 4937865: status RUNNING/PENDING
2026-06-21 19:54:22 INFO None 4937868: status RUNNING/PENDING
2026-06-21 19:54:22 INFO None 4937870: status RUNNING/PENDING
2026-06-21 19:54:22 INFO None 4937872: status RUNNING/PENDING
2026-06-21 19:54:22 INFO None 4937873: status RUNNING/PENDING
2026-06-21 19:54:22 INFO None 4937875: status RUNNING/PENDING
2026-06-21 19:54:22 INFO None 4937876: status RUNNING/PENDING
2026-06-21 19:54:22 INFO None 4937878: status RUNNING/PENDING
2026-06-21 19:54:22 INFO None 4937880: status RUNNING/PENDING
2026-06-21 19:54:22 INFO None 4937881: status RUNNING/PENDING
2026-06-21 19:54:22 INFO None 4937883: status RUNNING/PENDING
2026-06-21 19:54:22 INFO None 4937884: status RUNNING/PENDING
2026-06-21 19:54:22 INFO None 4937887: status RUNNING/PENDING
2026-06-21 19:54:22 INFO None 4937889: status RUNNING/PENDING
2026-06-21 19:54:22 INFO None 4937891: status RUNNING/PENDING
2026-06-21 19:54:22 INFO Jobs still running: ['4937865', '4937868', '4937870', '4937872', '4937873', '4937875', '4937876', '4937878', '4937880', '4937881', '4937883', '4937884', '4937887', '4937889', '4937891']. Waiting...
2026-06-21 19:54:37 INFO None 4937865: status RUNNING/PENDING
2026-06-21 19:54:37 INFO None 4937868: status RUNNING/PENDING
2026-06-21 19:54:37 INFO None 4937870: status RUNNING/PENDING
2026-06-21 19:54:37 INFO None 4937872: status RUNNING/PENDING
2026-06-21 19:54:37 INFO None 4937873: status RUNNING/PENDING
2026-06-21 19:54:37 INFO None 4937875: status RUNNING/PENDING
2026-06-21 19:54:37 INFO None 4937876: status RUNNING/PENDING
2026-06-21 19:54:37 INFO None 4937878: status RUNNING/PENDING
2026-06-21 19:54:37 INFO None 4937880: status RUNNING/PENDING
2026-06-21 19:54:37 INFO None 4937881: status RUNNING/PENDING
2026-06-21 19:54:37 INFO None 4937883: status RUNNING/PENDING
2026-06-21 19:54:37 INFO None 4937884: status RUNNING/PENDING
2026-06-21 19:54:37 INFO None 4937887: status RUNNING/PENDING
2026-06-21 19:54:37 INFO None 4937889: status RUNNING/PENDING
2026-06-21 19:54:37 INFO None 4937891: status RUNNING/PENDING
2026-06-21 19:54:37 INFO Jobs still running: ['4937865', '4937868', '4937870', '4937872', '4937873', '4937875', '4937876', '4937878', '4937880', '4937881', '4937883', '4937884', '4937887', '4937889', '4937891']. Waiting...
2026-06-21 19:54:52 INFO None 4937865: status RUNNING/PENDING
2026-06-21 19:54:53 INFO None 4937868: status FINISHED
2026-06-21 19:54:53 INFO None 4937870: status FINISHED
2026-06-21 19:54:53 INFO None 4937872: status RUNNING/PENDING
2026-06-21 19:54:53 INFO None 4937873: status RUNNING/PENDING
2026-06-21 19:54:53 INFO None 4937875: status FINISHED
2026-06-21 19:54:53 INFO None 4937876: status RUNNING/PENDING
2026-06-21 19:54:53 INFO None 4937878: status RUNNING/PENDING
2026-06-21 19:54:53 INFO None 4937880: status RUNNING/PENDING
2026-06-21 19:54:53 INFO None 4937881: status RUNNING/PENDING
2026-06-21 19:54:53 INFO None 4937883: status RUNNING/PENDING
2026-06-21 19:54:53 INFO None 4937884: status RUNNING/PENDING
2026-06-21 19:54:53 INFO None 4937887: status RUNNING/PENDING
2026-06-21 19:54:53 INFO None 4937889: status RUNNING/PENDING
2026-06-21 19:54:53 INFO None 4937891: status RUNNING/PENDING
2026-06-21 19:54:53 INFO Jobs still running: ['4937865', '4937872', '4937873', '4937876', '4937878', '4937880', '4937881', '4937883', '4937884', '4937887', '4937889', '4937891']. Waiting...
2026-06-21 19:55:08 INFO None 4937865: status FINISHED
2026-06-21 19:55:08 INFO None 4937868: status FINISHED
2026-06-21 19:55:08 INFO None 4937870: status FINISHED
2026-06-21 19:55:08 INFO None 4937872: status FINISHED
2026-06-21 19:55:08 INFO None 4937873: status RUNNING/PENDING
2026-06-21 19:55:08 INFO None 4937875: status FINISHED
2026-06-21 19:55:08 INFO None 4937876: status FINISHED
2026-06-21 19:55:08 INFO None 4937878: status RUNNING/PENDING
2026-06-21 19:55:08 INFO None 4937880: status RUNNING/PENDING
2026-06-21 19:55:08 INFO None 4937881: status RUNNING/PENDING
2026-06-21 19:55:08 INFO None 4937883: status FINISHED
2026-06-21 19:55:08 INFO None 4937884: status RUNNING/PENDING
2026-06-21 19:55:08 INFO None 4937887: status RUNNING/PENDING
2026-06-21 19:55:08 INFO None 4937889: status RUNNING/PENDING
2026-06-21 19:55:08 INFO None 4937891: status RUNNING/PENDING
2026-06-21 19:55:08 INFO Jobs still running: ['4937873', '4937878', '4937880', '4937881', '4937884', '4937887', '4937889', '4937891']. Waiting...
2026-06-21 19:55:23 INFO None 4937865: status FINISHED
2026-06-21 19:55:23 INFO None 4937868: status FINISHED
2026-06-21 19:55:23 INFO None 4937870: status FINISHED
2026-06-21 19:55:23 INFO None 4937872: status FINISHED
2026-06-21 19:55:23 INFO None 4937873: status RUNNING/PENDING
2026-06-21 19:55:23 INFO None 4937875: status FINISHED
2026-06-21 19:55:23 INFO None 4937876: status FINISHED
2026-06-21 19:55:23 INFO None 4937878: status RUNNING/PENDING
2026-06-21 19:55:23 INFO None 4937880: status FINISHED
2026-06-21 19:55:23 INFO None 4937881: status FINISHED
2026-06-21 19:55:23 INFO None 4937883: status FINISHED
2026-06-21 19:55:23 INFO None 4937884: status FINISHED
2026-06-21 19:55:23 INFO None 4937887: status FINISHED
2026-06-21 19:55:23 INFO None 4937889: status RUNNING/PENDING
2026-06-21 19:55:23 INFO None 4937891: status RUNNING/PENDING
2026-06-21 19:55:23 INFO Jobs still running: ['4937873', '4937878', '4937889', '4937891']. Waiting...
2026-06-21 19:55:38 INFO None 4937865: status FINISHED
2026-06-21 19:55:38 INFO None 4937868: status FINISHED
2026-06-21 19:55:38 INFO None 4937870: status FINISHED
2026-06-21 19:55:38 INFO None 4937872: status FINISHED
2026-06-21 19:55:38 INFO None 4937873: status RUNNING/PENDING
2026-06-21 19:55:38 INFO None 4937875: status FINISHED
2026-06-21 19:55:38 INFO None 4937876: status FINISHED
2026-06-21 19:55:39 INFO None 4937878: status RUNNING/PENDING
2026-06-21 19:55:39 INFO None 4937880: status FINISHED
2026-06-21 19:55:39 INFO None 4937881: status FINISHED
2026-06-21 19:55:39 INFO None 4937883: status FINISHED
2026-06-21 19:55:39 INFO None 4937884: status FINISHED
2026-06-21 19:55:39 INFO None 4937887: status FINISHED
2026-06-21 19:55:39 INFO None 4937889: status RUNNING/PENDING
2026-06-21 19:55:39 INFO None 4937891: status RUNNING/PENDING
2026-06-21 19:55:39 INFO Jobs still running: ['4937873', '4937878', '4937889', '4937891']. Waiting...
2026-06-21 19:55:54 INFO None 4937865: status FINISHED
2026-06-21 19:55:54 INFO None 4937868: status FINISHED
2026-06-21 19:55:54 INFO None 4937870: status FINISHED
2026-06-21 19:55:54 INFO None 4937872: status FINISHED
2026-06-21 19:55:54 INFO None 4937873: status RUNNING/PENDING
2026-06-21 19:55:54 INFO None 4937875: status FINISHED
2026-06-21 19:55:54 INFO None 4937876: status FINISHED
2026-06-21 19:55:54 INFO None 4937878: status FINISHED
2026-06-21 19:55:54 INFO None 4937880: status FINISHED
2026-06-21 19:55:54 INFO None 4937881: status FINISHED
2026-06-21 19:55:54 INFO None 4937883: status FINISHED
2026-06-21 19:55:54 INFO None 4937884: status FINISHED
2026-06-21 19:55:54 INFO None 4937887: status FINISHED
2026-06-21 19:55:54 INFO None 4937889: status RUNNING/PENDING
2026-06-21 19:55:54 INFO None 4937891: status FINISHED
2026-06-21 19:55:54 INFO Jobs still running: ['4937873', '4937889']. Waiting...
2026-06-21 19:56:09 INFO None 4937865: status FINISHED
2026-06-21 19:56:09 INFO None 4937868: status FINISHED
2026-06-21 19:56:09 INFO None 4937870: status FINISHED
2026-06-21 19:56:09 INFO None 4937872: status FINISHED
2026-06-21 19:56:09 INFO None 4937873: status FINISHED
2026-06-21 19:56:09 INFO None 4937875: status FINISHED
2026-06-21 19:56:09 INFO None 4937876: status FINISHED
2026-06-21 19:56:09 INFO None 4937878: status FINISHED
2026-06-21 19:56:09 INFO None 4937880: status FINISHED
2026-06-21 19:56:09 INFO None 4937881: status FINISHED
2026-06-21 19:56:09 INFO None 4937883: status FINISHED
2026-06-21 19:56:09 INFO None 4937884: status FINISHED
2026-06-21 19:56:09 INFO None 4937887: status FINISHED
2026-06-21 19:56:09 INFO None 4937889: status FINISHED
2026-06-21 19:56:09 INFO None 4937891: status FINISHED
2026-06-21 19:56:09 INFO Jobs ['4937865', '4937868', '4937870', '4937872', '4937873', '4937875', '4937876', '4937878', '4937880', '4937881', '4937883', '4937884', '4937887', '4937889', '4937891'] have finished
2026-06-21 19:56:09 INFO Checking restart files were created ...
2026-06-21 19:56:09 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/end.2020021501_7_ENS1.nc(2671953315 bytes)
2026-06-21 19:56:09 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/end.2020021501_7_ENS2.nc(2671953315 bytes)
2026-06-21 19:56:09 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/end.2020021501_7_ENS3.nc(2671953315 bytes)
2026-06-21 19:56:09 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/end.2020021501_7_ENS4.nc(2671953315 bytes)
2026-06-21 19:56:09 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/end.2020021501_7_ENS5.nc(2671953315 bytes)
2026-06-21 19:56:09 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/end.2020021501_7_ENS6.nc(2671953315 bytes)
2026-06-21 19:56:09 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/end.2020021501_7_ENS7.nc(2671953315 bytes)
2026-06-21 19:56:09 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/end.2020021501_7_ENS8.nc(2671953315 bytes)
2026-06-21 19:56:09 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/end.2020021501_7_ENS9.nc(2671953315 bytes)
2026-06-21 19:56:09 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/end.2020021501_7_ENS10.nc(2671953315 bytes)
2026-06-21 19:56:09 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/end.2020021501_7_ENS11.nc(2671953315 bytes)
2026-06-21 19:56:09 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/end.2020021501_7_ENS12.nc(2671953315 bytes)
2026-06-21 19:56:09 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/end.2020021501_7_ENS13.nc(2671953315 bytes)
2026-06-21 19:56:09 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/end.2020021501_7_ENS14.nc(2671953315 bytes)
2026-06-21 19:56:09 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/end.2020021501_7_ENS15.nc(2671953315 bytes)
2026-06-21 19:56:09 INFO  Run_model() completed successfully.
2026-06-21 19:56:09 INFO [TIME] after_model_set_simulated_time current_time=2020-02-15 01:00:00 simulated_time=2020-02-15 08:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 19:56:09 INFO [TIME] gregorian_conversion simulated_time=2020-02-15 08:00:00 days=153081 seconds=28800
2026-06-21 19:56:09 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-06-21 19:56:09 INFO [TIME] increment current_time 2020-02-15 01:00:00 -> 2020-02-15 08:00:00
2026-06-21 19:56:09 INFO [TIME] after_increment_before_assimilation current_time=2020-02-15 08:00:00 simulated_time=2020-02-15 08:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 19:56:09 INFO ---------->>> Running process_satellite_data()
2026-06-21 19:56:09 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12129.nc
2026-06-21 19:56:09 INFO ---------->>> Running run_obs_converter()
2026-06-21 19:56:09 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_29326_153081.out
2026-06-21 19:56:09 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_29326_153081.out
2026-06-21 19:56:09 INFO ---------->>> Running DART
2026-06-21 19:56:09 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-06-21 19:56:09 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/chim_ENS1_2020021501_7_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/chim_ENS1_2020021508_1_out_toDART.nc
2026-06-21 19:56:10 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/chim_ENS2_2020021501_7_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/chim_ENS2_2020021508_1_out_toDART.nc
2026-06-21 19:56:10 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/chim_ENS3_2020021501_7_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/chim_ENS3_2020021508_1_out_toDART.nc
2026-06-21 19:56:11 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/chim_ENS4_2020021501_7_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/chim_ENS4_2020021508_1_out_toDART.nc
2026-06-21 19:56:11 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/chim_ENS5_2020021501_7_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/chim_ENS5_2020021508_1_out_toDART.nc
2026-06-21 19:56:11 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/chim_ENS6_2020021501_7_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/chim_ENS6_2020021508_1_out_toDART.nc
2026-06-21 19:56:12 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/chim_ENS7_2020021501_7_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/chim_ENS7_2020021508_1_out_toDART.nc
2026-06-21 19:56:12 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/chim_ENS8_2020021501_7_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/chim_ENS8_2020021508_1_out_toDART.nc
2026-06-21 19:56:13 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/chim_ENS9_2020021501_7_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/chim_ENS9_2020021508_1_out_toDART.nc
2026-06-21 19:56:13 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/chim_ENS10_2020021501_7_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/chim_ENS10_2020021508_1_out_toDART.nc
2026-06-21 19:56:13 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/chim_ENS11_2020021501_7_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/chim_ENS11_2020021508_1_out_toDART.nc
2026-06-21 19:56:14 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/chim_ENS12_2020021501_7_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/chim_ENS12_2020021508_1_out_toDART.nc
2026-06-21 19:56:14 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/chim_ENS13_2020021501_7_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/chim_ENS13_2020021508_1_out_toDART.nc
2026-06-21 19:56:15 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/chim_ENS14_2020021501_7_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/chim_ENS14_2020021508_1_out_toDART.nc
2026-06-21 19:56:15 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/chim_ENS15_2020021501_7_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/chim_ENS15_2020021508_1_out_toDART.nc
2026-06-21 19:56:15 INFO Replacement input_template.nml → input.nml completed successfully.
2026-06-21 19:56:15 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-06-21 19:56:15 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-06-21 19:56:15 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-06-21 19:56:15 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-06-21 19:56:16 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-06-21 19:56:22 INFO Found: []
2026-06-21 19:56:22 INFO No job id returned by command ./run_filter.bsh
2026-06-21 19:56:22 INFO No monitoring will be performed
2026-06-21 19:56:22 INFO Moving DART output files to analysis and preassim directories for date 2020021508 if present ...
2026-06-21 19:56:22 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021508'
2026-06-21 19:56:22 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021508'
2026-06-21 19:56:22 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021508'
2026-06-21 19:56:22 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021508'
2026-06-21 19:56:22 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021508'
2026-06-21 19:56:22 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021508'
2026-06-21 19:56:22 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021508'
2026-06-21 19:56:22 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021508'
2026-06-21 19:56:22 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021508'
2026-06-21 19:56:22 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021508'
2026-06-21 19:56:22 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021508'
2026-06-21 19:56:22 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021508'
2026-06-21 19:56:22 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021508'
2026-06-21 19:56:22 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021508'
2026-06-21 19:56:22 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021508'
2026-06-21 19:56:22 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021508'
2026-06-21 19:56:22 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021508'
2026-06-21 19:56:22 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021508'
2026-06-21 19:56:22 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021508'
2026-06-21 19:56:22 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021508'
2026-06-21 19:56:22 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021508'
2026-06-21 19:56:22 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021508'
2026-06-21 19:56:22 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021508'
2026-06-21 19:56:22 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021508'
2026-06-21 19:56:22 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021508'
2026-06-21 19:56:22 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021508'
2026-06-21 19:56:22 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021508'
2026-06-21 19:56:22 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021508'
2026-06-21 19:56:22 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021508'
2026-06-21 19:56:22 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021508'
2026-06-21 19:56:22 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021508'
2026-06-21 19:56:22 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021508'
2026-06-21 19:56:22 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021508'
2026-06-21 19:56:22 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021508'
2026-06-21 19:56:22 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-06-21 19:56:25 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-06-21 19:56:25 INFO run_dart() is DONE.
2026-06-21 19:56:25 INFO ---------->>> Running update_pollutant_in_end()
2026-06-21 19:56:25 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/end.2020021501_7_ENS1.nc
2026-06-21 19:56:37 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 19:56:38 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 19:56:38 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 19:56:38 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 19:56:39 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 19:56:39 INFO [Posterior Diff ENS1] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021508/diff_posterior_ENS1_2020021508.nc
2026-06-21 19:56:39 INFO [Posterior Diff ENS1] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021508/diff_posterior_ENS1_2020021508.relative.nc
2026-06-21 19:56:40 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/end.2020021501_7_ENS2.nc
2026-06-21 19:56:52 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 19:56:52 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 19:56:53 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 19:56:53 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 19:56:54 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 19:56:54 INFO [Posterior Diff ENS2] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021508/diff_posterior_ENS2_2020021508.nc
2026-06-21 19:56:54 INFO [Posterior Diff ENS2] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021508/diff_posterior_ENS2_2020021508.relative.nc
2026-06-21 19:56:54 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/end.2020021501_7_ENS3.nc
2026-06-21 19:57:07 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 19:57:07 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 19:57:08 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 19:57:08 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 19:57:08 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 19:57:09 INFO [Posterior Diff ENS3] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021508/diff_posterior_ENS3_2020021508.nc
2026-06-21 19:57:09 INFO [Posterior Diff ENS3] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021508/diff_posterior_ENS3_2020021508.relative.nc
2026-06-21 19:57:09 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/end.2020021501_7_ENS4.nc
2026-06-21 19:57:22 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 19:57:22 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 19:57:23 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 19:57:23 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 19:57:23 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 19:57:23 INFO [Posterior Diff ENS4] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021508/diff_posterior_ENS4_2020021508.nc
2026-06-21 19:57:23 INFO [Posterior Diff ENS4] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021508/diff_posterior_ENS4_2020021508.relative.nc
2026-06-21 19:57:24 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/end.2020021501_7_ENS5.nc
2026-06-21 19:57:36 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 19:57:36 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 19:57:37 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 19:57:37 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 19:57:38 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 19:57:38 INFO [Posterior Diff ENS5] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021508/diff_posterior_ENS5_2020021508.nc
2026-06-21 19:57:38 INFO [Posterior Diff ENS5] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021508/diff_posterior_ENS5_2020021508.relative.nc
2026-06-21 19:57:38 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/end.2020021501_7_ENS6.nc
2026-06-21 19:57:51 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 19:57:51 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 19:57:52 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 19:57:52 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 19:57:52 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 19:57:53 INFO [Posterior Diff ENS6] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021508/diff_posterior_ENS6_2020021508.nc
2026-06-21 19:57:53 INFO [Posterior Diff ENS6] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021508/diff_posterior_ENS6_2020021508.relative.nc
2026-06-21 19:57:53 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/end.2020021501_7_ENS7.nc
2026-06-21 19:58:06 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 19:58:06 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 19:58:07 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 19:58:07 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 19:58:07 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 19:58:07 INFO [Posterior Diff ENS7] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021508/diff_posterior_ENS7_2020021508.nc
2026-06-21 19:58:07 INFO [Posterior Diff ENS7] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021508/diff_posterior_ENS7_2020021508.relative.nc
2026-06-21 19:58:08 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/end.2020021501_7_ENS8.nc
2026-06-21 19:58:20 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 19:58:20 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 19:58:21 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 19:58:21 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 19:58:21 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 19:58:22 INFO [Posterior Diff ENS8] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021508/diff_posterior_ENS8_2020021508.nc
2026-06-21 19:58:22 INFO [Posterior Diff ENS8] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021508/diff_posterior_ENS8_2020021508.relative.nc
2026-06-21 19:58:22 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/end.2020021501_7_ENS9.nc
2026-06-21 19:58:34 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 19:58:34 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 19:58:35 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 19:58:35 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 19:58:36 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 19:58:36 INFO [Posterior Diff ENS9] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021508/diff_posterior_ENS9_2020021508.nc
2026-06-21 19:58:36 INFO [Posterior Diff ENS9] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021508/diff_posterior_ENS9_2020021508.relative.nc
2026-06-21 19:58:36 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/end.2020021501_7_ENS10.nc
2026-06-21 19:58:49 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 19:58:49 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 19:58:49 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 19:58:49 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 19:58:50 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 19:58:51 INFO [Posterior Diff ENS10] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021508/diff_posterior_ENS10_2020021508.nc
2026-06-21 19:58:51 INFO [Posterior Diff ENS10] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021508/diff_posterior_ENS10_2020021508.relative.nc
2026-06-21 19:58:51 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/end.2020021501_7_ENS11.nc
2026-06-21 19:59:03 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 19:59:03 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 19:59:04 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 19:59:04 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 19:59:04 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 19:59:05 INFO [Posterior Diff ENS11] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021508/diff_posterior_ENS11_2020021508.nc
2026-06-21 19:59:05 INFO [Posterior Diff ENS11] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021508/diff_posterior_ENS11_2020021508.relative.nc
2026-06-21 19:59:05 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/end.2020021501_7_ENS12.nc
2026-06-21 19:59:18 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 19:59:18 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 19:59:19 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 19:59:19 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 19:59:19 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 19:59:20 INFO [Posterior Diff ENS12] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021508/diff_posterior_ENS12_2020021508.nc
2026-06-21 19:59:20 INFO [Posterior Diff ENS12] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021508/diff_posterior_ENS12_2020021508.relative.nc
2026-06-21 19:59:20 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/end.2020021501_7_ENS13.nc
2026-06-21 19:59:33 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 19:59:33 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 19:59:34 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 19:59:34 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 19:59:34 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 19:59:35 INFO [Posterior Diff ENS13] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021508/diff_posterior_ENS13_2020021508.nc
2026-06-21 19:59:35 INFO [Posterior Diff ENS13] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021508/diff_posterior_ENS13_2020021508.relative.nc
2026-06-21 19:59:35 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/end.2020021501_7_ENS14.nc
2026-06-21 19:59:47 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 19:59:47 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 19:59:48 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 19:59:48 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 19:59:49 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 19:59:49 INFO [Posterior Diff ENS14] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021508/diff_posterior_ENS14_2020021508.nc
2026-06-21 19:59:49 INFO [Posterior Diff ENS14] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021508/diff_posterior_ENS14_2020021508.relative.nc
2026-06-21 19:59:50 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/end.2020021501_7_ENS15.nc
2026-06-21 20:00:02 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:00:02 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:00:03 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:00:03 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:00:04 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:00:04 INFO [Posterior Diff ENS15] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021508/diff_posterior_ENS15_2020021508.nc
2026-06-21 20:00:04 INFO [Posterior Diff ENS15] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021508/diff_posterior_ENS15_2020021508.relative.nc
2026-06-21 20:00:04 INFO Next run starts from 2020-02-15 08:00:00
2026-06-21 20:00:04 INFO Cycle is DONE; starting a new loop!
2026-06-21 20:00:04 INFO [TIME] step_end current_time=2020-02-15 08:00:00 simulated_time=2020-02-15 08:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 20:00:04 INFO [TIME] step_start current_time=2020-02-15 08:00:00 simulated_time=2020-02-15 08:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 20:00:04 INFO [TIME] window start=2020-02-15 08:00:00 end=2020-02-15 10:00:00 run_hours=2 has_assimilation=True
2026-06-21 20:00:04 INFO Creating directories and links for ENS1 to run chimere's parallel part
2026-06-21 20:00:04 INFO Linking EMIS ...
2026-06-21 20:00:04 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens01.nc
2026-06-21 20:00:04 INFO >> Checking links...
2026-06-21 20:00:07 INFO >> All links are good for ENS1  ...
2026-06-21 20:00:07 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:00:08 INFO Hourly dataset computed and listing created
2026-06-21 20:00:12 INFO Hourly dataset computed
2026-06-21 20:00:12 INFO Creating directories and links for ENS2 to run chimere's parallel part
2026-06-21 20:00:12 INFO Linking EMIS ...
2026-06-21 20:00:12 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens02.nc
2026-06-21 20:00:12 INFO >> Checking links...
2026-06-21 20:00:15 INFO >> All links are good for ENS2  ...
2026-06-21 20:00:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:00:16 INFO Hourly dataset computed and listing created
2026-06-21 20:00:17 INFO Hourly dataset computed
2026-06-21 20:00:17 INFO Creating directories and links for ENS3 to run chimere's parallel part
2026-06-21 20:00:17 INFO Linking EMIS ...
2026-06-21 20:00:17 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens03.nc
2026-06-21 20:00:17 INFO >> Checking links...
2026-06-21 20:00:20 INFO >> All links are good for ENS3  ...
2026-06-21 20:00:20 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:00:21 INFO Hourly dataset computed and listing created
2026-06-21 20:00:21 INFO Hourly dataset computed
2026-06-21 20:00:21 INFO Creating directories and links for ENS4 to run chimere's parallel part
2026-06-21 20:00:21 INFO Linking EMIS ...
2026-06-21 20:00:22 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens04.nc
2026-06-21 20:00:22 INFO >> Checking links...
2026-06-21 20:00:24 INFO >> All links are good for ENS4  ...
2026-06-21 20:00:24 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:00:25 INFO Hourly dataset computed and listing created
2026-06-21 20:00:26 INFO Hourly dataset computed
2026-06-21 20:00:26 INFO Creating directories and links for ENS5 to run chimere's parallel part
2026-06-21 20:00:26 INFO Linking EMIS ...
2026-06-21 20:00:27 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens05.nc
2026-06-21 20:00:27 INFO >> Checking links...
2026-06-21 20:00:29 INFO >> All links are good for ENS5  ...
2026-06-21 20:00:29 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:00:30 INFO Hourly dataset computed and listing created
2026-06-21 20:00:31 INFO Hourly dataset computed
2026-06-21 20:00:31 INFO Creating directories and links for ENS6 to run chimere's parallel part
2026-06-21 20:00:31 INFO Linking EMIS ...
2026-06-21 20:00:31 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens06.nc
2026-06-21 20:00:31 INFO >> Checking links...
2026-06-21 20:00:33 INFO >> All links are good for ENS6  ...
2026-06-21 20:00:33 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:00:34 INFO Hourly dataset computed and listing created
2026-06-21 20:00:35 INFO Hourly dataset computed
2026-06-21 20:00:35 INFO Creating directories and links for ENS7 to run chimere's parallel part
2026-06-21 20:00:35 INFO Linking EMIS ...
2026-06-21 20:00:36 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens07.nc
2026-06-21 20:00:36 INFO >> Checking links...
2026-06-21 20:00:38 INFO >> All links are good for ENS7  ...
2026-06-21 20:00:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:00:39 INFO Hourly dataset computed and listing created
2026-06-21 20:00:40 INFO Hourly dataset computed
2026-06-21 20:00:40 INFO Creating directories and links for ENS8 to run chimere's parallel part
2026-06-21 20:00:40 INFO Linking EMIS ...
2026-06-21 20:00:40 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens08.nc
2026-06-21 20:00:40 INFO >> Checking links...
2026-06-21 20:00:43 INFO >> All links are good for ENS8  ...
2026-06-21 20:00:43 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:00:44 INFO Hourly dataset computed and listing created
2026-06-21 20:00:44 INFO Hourly dataset computed
2026-06-21 20:00:44 INFO Creating directories and links for ENS9 to run chimere's parallel part
2026-06-21 20:00:44 INFO Linking EMIS ...
2026-06-21 20:00:45 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens09.nc
2026-06-21 20:00:45 INFO >> Checking links...
2026-06-21 20:00:47 INFO >> All links are good for ENS9  ...
2026-06-21 20:00:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:00:48 INFO Hourly dataset computed and listing created
2026-06-21 20:00:49 INFO Hourly dataset computed
2026-06-21 20:00:49 INFO Creating directories and links for ENS10 to run chimere's parallel part
2026-06-21 20:00:49 INFO Linking EMIS ...
2026-06-21 20:00:49 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens10.nc
2026-06-21 20:00:49 INFO >> Checking links...
2026-06-21 20:00:52 INFO >> All links are good for ENS10  ...
2026-06-21 20:00:52 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:00:53 INFO Hourly dataset computed and listing created
2026-06-21 20:00:54 INFO Hourly dataset computed
2026-06-21 20:00:54 INFO Creating directories and links for ENS11 to run chimere's parallel part
2026-06-21 20:00:54 INFO Linking EMIS ...
2026-06-21 20:00:54 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens11.nc
2026-06-21 20:00:54 INFO >> Checking links...
2026-06-21 20:00:56 INFO >> All links are good for ENS11  ...
2026-06-21 20:00:56 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:00:57 INFO Hourly dataset computed and listing created
2026-06-21 20:00:58 INFO Hourly dataset computed
2026-06-21 20:00:58 INFO Creating directories and links for ENS12 to run chimere's parallel part
2026-06-21 20:00:58 INFO Linking EMIS ...
2026-06-21 20:00:58 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens12.nc
2026-06-21 20:00:58 INFO >> Checking links...
2026-06-21 20:01:01 INFO >> All links are good for ENS12  ...
2026-06-21 20:01:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:01:02 INFO Hourly dataset computed and listing created
2026-06-21 20:01:06 INFO Hourly dataset computed
2026-06-21 20:01:06 INFO Creating directories and links for ENS13 to run chimere's parallel part
2026-06-21 20:01:06 INFO Linking EMIS ...
2026-06-21 20:01:06 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens13.nc
2026-06-21 20:01:06 INFO >> Checking links...
2026-06-21 20:01:09 INFO >> All links are good for ENS13  ...
2026-06-21 20:01:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:01:10 INFO Hourly dataset computed and listing created
2026-06-21 20:01:11 INFO Hourly dataset computed
2026-06-21 20:01:11 INFO Creating directories and links for ENS14 to run chimere's parallel part
2026-06-21 20:01:11 INFO Linking EMIS ...
2026-06-21 20:01:11 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens14.nc
2026-06-21 20:01:11 INFO >> Checking links...
2026-06-21 20:01:13 INFO >> All links are good for ENS14  ...
2026-06-21 20:01:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:01:14 INFO Hourly dataset computed and listing created
2026-06-21 20:01:15 INFO Hourly dataset computed
2026-06-21 20:01:15 INFO Creating directories and links for ENS15 to run chimere's parallel part
2026-06-21 20:01:15 INFO Linking EMIS ...
2026-06-21 20:01:16 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens15.nc
2026-06-21 20:01:16 INFO >> Checking links...
2026-06-21 20:01:18 INFO >> All links are good for ENS15  ...
2026-06-21 20:01:18 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:01:19 INFO Hourly dataset computed and listing created
2026-06-21 20:01:20 INFO Hourly dataset computed
2026-06-21 20:01:20 INFO ---------->>> Running CHIMERE model from 2020-02-15 08:00:00 to 2020-02-15 10:00:00
2026-06-21 20:01:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:01:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1
2026-06-21 20:01:20 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/end.2020021501_7_ENS1.nc
2026-06-21 20:01:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-06-21 20:01:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:01:20 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-06-21 20:01:20 INFO Queuing job for member 1...
2026-06-21 20:01:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:01:20 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-06-21 20:01:20 INFO Found: ['4937956']
2026-06-21 20:01:25 INFO [TGCC-IRENE] Submitted job with ID:['4937956']
2026-06-21 20:01:25 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:01:25 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2
2026-06-21 20:01:25 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/end.2020021501_7_ENS2.nc
2026-06-21 20:01:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-06-21 20:01:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:01:26 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-06-21 20:01:26 INFO Queuing job for member 2...
2026-06-21 20:01:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:01:26 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-06-21 20:01:26 INFO Found: ['4937957']
2026-06-21 20:01:31 INFO [TGCC-IRENE] Submitted job with ID:['4937957']
2026-06-21 20:01:31 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:01:31 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3
2026-06-21 20:01:31 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/end.2020021501_7_ENS3.nc
2026-06-21 20:01:31 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-06-21 20:01:31 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:01:31 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-06-21 20:01:31 INFO Queuing job for member 3...
2026-06-21 20:01:31 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:01:31 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-06-21 20:01:32 INFO Found: ['4937958']
2026-06-21 20:01:37 INFO [TGCC-IRENE] Submitted job with ID:['4937958']
2026-06-21 20:01:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:01:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4
2026-06-21 20:01:37 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/end.2020021501_7_ENS4.nc
2026-06-21 20:01:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-06-21 20:01:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:01:37 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-06-21 20:01:37 INFO Queuing job for member 4...
2026-06-21 20:01:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:01:37 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-06-21 20:01:38 INFO Found: ['4937959']
2026-06-21 20:01:43 INFO [TGCC-IRENE] Submitted job with ID:['4937959']
2026-06-21 20:01:43 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:01:43 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5
2026-06-21 20:01:43 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/end.2020021501_7_ENS5.nc
2026-06-21 20:01:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-06-21 20:01:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:01:43 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-06-21 20:01:43 INFO Queuing job for member 5...
2026-06-21 20:01:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:01:43 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-06-21 20:01:43 INFO Found: ['4937960']
2026-06-21 20:01:48 INFO [TGCC-IRENE] Submitted job with ID:['4937960']
2026-06-21 20:01:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:01:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6
2026-06-21 20:01:48 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/end.2020021501_7_ENS6.nc
2026-06-21 20:01:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-06-21 20:01:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:01:48 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-06-21 20:01:48 INFO Queuing job for member 6...
2026-06-21 20:01:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:01:48 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-06-21 20:01:49 INFO Found: ['4937961']
2026-06-21 20:01:54 INFO [TGCC-IRENE] Submitted job with ID:['4937961']
2026-06-21 20:01:54 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:01:54 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7
2026-06-21 20:01:54 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/end.2020021501_7_ENS7.nc
2026-06-21 20:01:54 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-06-21 20:01:54 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:01:54 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-06-21 20:01:54 INFO Queuing job for member 7...
2026-06-21 20:01:54 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:01:54 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-06-21 20:01:55 INFO Found: ['4937963']
2026-06-21 20:02:00 INFO [TGCC-IRENE] Submitted job with ID:['4937963']
2026-06-21 20:02:00 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:02:00 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8
2026-06-21 20:02:00 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/end.2020021501_7_ENS8.nc
2026-06-21 20:02:00 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-06-21 20:02:00 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:02:00 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-06-21 20:02:00 INFO Queuing job for member 8...
2026-06-21 20:02:00 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:02:00 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-06-21 20:02:01 INFO Found: ['4937964']
2026-06-21 20:02:06 INFO [TGCC-IRENE] Submitted job with ID:['4937964']
2026-06-21 20:02:06 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:02:06 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9
2026-06-21 20:02:06 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/end.2020021501_7_ENS9.nc
2026-06-21 20:02:06 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-06-21 20:02:06 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:02:06 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-06-21 20:02:06 INFO Queuing job for member 9...
2026-06-21 20:02:06 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:02:06 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-06-21 20:02:07 INFO Found: ['4937966']
2026-06-21 20:02:12 INFO [TGCC-IRENE] Submitted job with ID:['4937966']
2026-06-21 20:02:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:02:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10
2026-06-21 20:02:12 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/end.2020021501_7_ENS10.nc
2026-06-21 20:02:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-06-21 20:02:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:02:12 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-06-21 20:02:12 INFO Queuing job for member 10...
2026-06-21 20:02:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:02:12 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-06-21 20:02:12 INFO Found: ['4937967']
2026-06-21 20:02:17 INFO [TGCC-IRENE] Submitted job with ID:['4937967']
2026-06-21 20:02:17 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:02:17 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11
2026-06-21 20:02:17 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/end.2020021501_7_ENS11.nc
2026-06-21 20:02:17 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-06-21 20:02:17 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:02:17 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-06-21 20:02:17 INFO Queuing job for member 11...
2026-06-21 20:02:17 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:02:17 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-06-21 20:02:18 INFO Found: ['4937969']
2026-06-21 20:02:23 INFO [TGCC-IRENE] Submitted job with ID:['4937969']
2026-06-21 20:02:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:02:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12
2026-06-21 20:02:23 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/end.2020021501_7_ENS12.nc
2026-06-21 20:02:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-06-21 20:02:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:02:23 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-06-21 20:02:23 INFO Queuing job for member 12...
2026-06-21 20:02:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:02:23 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-06-21 20:02:24 INFO Found: ['4937970']
2026-06-21 20:02:29 INFO [TGCC-IRENE] Submitted job with ID:['4937970']
2026-06-21 20:02:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:02:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13
2026-06-21 20:02:29 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/end.2020021501_7_ENS13.nc
2026-06-21 20:02:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-06-21 20:02:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:02:29 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-06-21 20:02:29 INFO Queuing job for member 13...
2026-06-21 20:02:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:02:29 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-06-21 20:02:30 INFO Found: ['4937971']
2026-06-21 20:02:35 INFO [TGCC-IRENE] Submitted job with ID:['4937971']
2026-06-21 20:02:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:02:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14
2026-06-21 20:02:35 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/end.2020021501_7_ENS14.nc
2026-06-21 20:02:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-06-21 20:02:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:02:35 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-06-21 20:02:35 INFO Queuing job for member 14...
2026-06-21 20:02:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:02:35 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-06-21 20:02:35 INFO Found: ['4937973']
2026-06-21 20:02:40 INFO [TGCC-IRENE] Submitted job with ID:['4937973']
2026-06-21 20:02:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:02:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15
2026-06-21 20:02:40 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/end.2020021501_7_ENS15.nc
2026-06-21 20:02:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-06-21 20:02:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:02:40 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-06-21 20:02:41 INFO Queuing job for member 15...
2026-06-21 20:02:41 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:02:41 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-06-21 20:02:42 INFO Found: ['4937974']
2026-06-21 20:02:47 INFO [TGCC-IRENE] Submitted job with ID:['4937974']
2026-06-21 20:02:47 INFO Checking job status ...
2026-06-21 20:02:47 INFO None 4937956: status RUNNING/PENDING
2026-06-21 20:02:47 INFO None 4937957: status RUNNING/PENDING
2026-06-21 20:02:47 INFO None 4937958: status RUNNING/PENDING
2026-06-21 20:02:47 INFO None 4937959: status RUNNING/PENDING
2026-06-21 20:02:47 INFO None 4937960: status RUNNING/PENDING
2026-06-21 20:02:47 INFO None 4937961: status RUNNING/PENDING
2026-06-21 20:02:47 INFO None 4937963: status RUNNING/PENDING
2026-06-21 20:02:47 INFO None 4937964: status RUNNING/PENDING
2026-06-21 20:02:47 INFO None 4937966: status RUNNING/PENDING
2026-06-21 20:02:47 INFO None 4937967: status RUNNING/PENDING
2026-06-21 20:02:47 INFO None 4937969: status RUNNING/PENDING
2026-06-21 20:02:47 INFO None 4937970: status RUNNING/PENDING
2026-06-21 20:02:47 INFO None 4937971: status RUNNING/PENDING
2026-06-21 20:02:47 INFO None 4937973: status RUNNING/PENDING
2026-06-21 20:02:47 INFO None 4937974: status RUNNING/PENDING
2026-06-21 20:02:47 INFO Jobs still running: ['4937956', '4937957', '4937958', '4937959', '4937960', '4937961', '4937963', '4937964', '4937966', '4937967', '4937969', '4937970', '4937971', '4937973', '4937974']. Waiting...
2026-06-21 20:03:02 INFO None 4937956: status RUNNING/PENDING
2026-06-21 20:03:02 INFO None 4937957: status RUNNING/PENDING
2026-06-21 20:03:02 INFO None 4937958: status RUNNING/PENDING
2026-06-21 20:03:02 INFO None 4937959: status RUNNING/PENDING
2026-06-21 20:03:02 INFO None 4937960: status RUNNING/PENDING
2026-06-21 20:03:02 INFO None 4937961: status RUNNING/PENDING
2026-06-21 20:03:02 INFO None 4937963: status RUNNING/PENDING
2026-06-21 20:03:02 INFO None 4937964: status RUNNING/PENDING
2026-06-21 20:03:02 INFO None 4937966: status RUNNING/PENDING
2026-06-21 20:03:02 INFO None 4937967: status RUNNING/PENDING
2026-06-21 20:03:02 INFO None 4937969: status RUNNING/PENDING
2026-06-21 20:03:02 INFO None 4937970: status RUNNING/PENDING
2026-06-21 20:03:02 INFO None 4937971: status RUNNING/PENDING
2026-06-21 20:03:02 INFO None 4937973: status RUNNING/PENDING
2026-06-21 20:03:02 INFO None 4937974: status RUNNING/PENDING
2026-06-21 20:03:02 INFO Jobs still running: ['4937956', '4937957', '4937958', '4937959', '4937960', '4937961', '4937963', '4937964', '4937966', '4937967', '4937969', '4937970', '4937971', '4937973', '4937974']. Waiting...
2026-06-21 20:03:17 INFO None 4937956: status RUNNING/PENDING
2026-06-21 20:03:17 INFO None 4937957: status RUNNING/PENDING
2026-06-21 20:03:17 INFO None 4937958: status RUNNING/PENDING
2026-06-21 20:03:17 INFO None 4937959: status RUNNING/PENDING
2026-06-21 20:03:17 INFO None 4937960: status RUNNING/PENDING
2026-06-21 20:03:17 INFO None 4937961: status RUNNING/PENDING
2026-06-21 20:03:17 INFO None 4937963: status RUNNING/PENDING
2026-06-21 20:03:17 INFO None 4937964: status RUNNING/PENDING
2026-06-21 20:03:17 INFO None 4937966: status RUNNING/PENDING
2026-06-21 20:03:17 INFO None 4937967: status RUNNING/PENDING
2026-06-21 20:03:18 INFO None 4937969: status RUNNING/PENDING
2026-06-21 20:03:18 INFO None 4937970: status RUNNING/PENDING
2026-06-21 20:03:18 INFO None 4937971: status RUNNING/PENDING
2026-06-21 20:03:18 INFO None 4937973: status RUNNING/PENDING
2026-06-21 20:03:18 INFO None 4937974: status RUNNING/PENDING
2026-06-21 20:03:18 INFO Jobs still running: ['4937956', '4937957', '4937958', '4937959', '4937960', '4937961', '4937963', '4937964', '4937966', '4937967', '4937969', '4937970', '4937971', '4937973', '4937974']. Waiting...
2026-06-21 20:03:33 INFO None 4937956: status RUNNING/PENDING
2026-06-21 20:03:33 INFO None 4937957: status RUNNING/PENDING
2026-06-21 20:03:33 INFO None 4937958: status RUNNING/PENDING
2026-06-21 20:03:33 INFO None 4937959: status RUNNING/PENDING
2026-06-21 20:03:33 INFO None 4937960: status RUNNING/PENDING
2026-06-21 20:03:33 INFO None 4937961: status RUNNING/PENDING
2026-06-21 20:03:33 INFO None 4937963: status RUNNING/PENDING
2026-06-21 20:03:33 INFO None 4937964: status RUNNING/PENDING
2026-06-21 20:03:33 INFO None 4937966: status RUNNING/PENDING
2026-06-21 20:03:33 INFO None 4937967: status RUNNING/PENDING
2026-06-21 20:03:33 INFO None 4937969: status RUNNING/PENDING
2026-06-21 20:03:33 INFO None 4937970: status RUNNING/PENDING
2026-06-21 20:03:33 INFO None 4937971: status RUNNING/PENDING
2026-06-21 20:03:33 INFO None 4937973: status RUNNING/PENDING
2026-06-21 20:03:33 INFO None 4937974: status RUNNING/PENDING
2026-06-21 20:03:33 INFO Jobs still running: ['4937956', '4937957', '4937958', '4937959', '4937960', '4937961', '4937963', '4937964', '4937966', '4937967', '4937969', '4937970', '4937971', '4937973', '4937974']. Waiting...
2026-06-21 20:03:48 INFO None 4937956: status RUNNING/PENDING
2026-06-21 20:03:48 INFO None 4937957: status FINISHED
2026-06-21 20:03:48 INFO None 4937958: status FINISHED
2026-06-21 20:03:48 INFO None 4937959: status FINISHED
2026-06-21 20:03:48 INFO None 4937960: status RUNNING/PENDING
2026-06-21 20:03:48 INFO None 4937961: status RUNNING/PENDING
2026-06-21 20:03:48 INFO None 4937963: status RUNNING/PENDING
2026-06-21 20:03:48 INFO None 4937964: status RUNNING/PENDING
2026-06-21 20:03:48 INFO None 4937966: status RUNNING/PENDING
2026-06-21 20:03:48 INFO None 4937967: status RUNNING/PENDING
2026-06-21 20:03:48 INFO None 4937969: status RUNNING/PENDING
2026-06-21 20:03:48 INFO None 4937970: status RUNNING/PENDING
2026-06-21 20:03:48 INFO None 4937971: status RUNNING/PENDING
2026-06-21 20:03:48 INFO None 4937973: status RUNNING/PENDING
2026-06-21 20:03:48 INFO None 4937974: status RUNNING/PENDING
2026-06-21 20:03:48 INFO Jobs still running: ['4937956', '4937960', '4937961', '4937963', '4937964', '4937966', '4937967', '4937969', '4937970', '4937971', '4937973', '4937974']. Waiting...
2026-06-21 20:04:03 INFO None 4937956: status RUNNING/PENDING
2026-06-21 20:04:03 INFO None 4937957: status FINISHED
2026-06-21 20:04:03 INFO None 4937958: status FINISHED
2026-06-21 20:04:03 INFO None 4937959: status FINISHED
2026-06-21 20:04:03 INFO None 4937960: status FINISHED
2026-06-21 20:04:03 INFO None 4937961: status FINISHED
2026-06-21 20:04:03 INFO None 4937963: status RUNNING/PENDING
2026-06-21 20:04:03 INFO None 4937964: status RUNNING/PENDING
2026-06-21 20:04:03 INFO None 4937966: status RUNNING/PENDING
2026-06-21 20:04:03 INFO None 4937967: status RUNNING/PENDING
2026-06-21 20:04:03 INFO None 4937969: status RUNNING/PENDING
2026-06-21 20:04:04 INFO None 4937970: status RUNNING/PENDING
2026-06-21 20:04:04 INFO None 4937971: status RUNNING/PENDING
2026-06-21 20:04:04 INFO None 4937973: status RUNNING/PENDING
2026-06-21 20:04:04 INFO None 4937974: status RUNNING/PENDING
2026-06-21 20:04:04 INFO Jobs still running: ['4937956', '4937963', '4937964', '4937966', '4937967', '4937969', '4937970', '4937971', '4937973', '4937974']. Waiting...
2026-06-21 20:04:19 INFO None 4937956: status RUNNING/PENDING
2026-06-21 20:04:19 INFO None 4937957: status FINISHED
2026-06-21 20:04:19 INFO None 4937958: status FINISHED
2026-06-21 20:04:19 INFO None 4937959: status FINISHED
2026-06-21 20:04:19 INFO None 4937960: status FINISHED
2026-06-21 20:04:19 INFO None 4937961: status FINISHED
2026-06-21 20:04:19 INFO None 4937963: status FINISHED
2026-06-21 20:04:19 INFO None 4937964: status FINISHED
2026-06-21 20:04:19 INFO None 4937966: status RUNNING/PENDING
2026-06-21 20:04:19 INFO None 4937967: status RUNNING/PENDING
2026-06-21 20:04:19 INFO None 4937969: status RUNNING/PENDING
2026-06-21 20:04:19 INFO None 4937970: status RUNNING/PENDING
2026-06-21 20:04:19 INFO None 4937971: status RUNNING/PENDING
2026-06-21 20:04:19 INFO None 4937973: status RUNNING/PENDING
2026-06-21 20:04:19 INFO None 4937974: status RUNNING/PENDING
2026-06-21 20:04:19 INFO Jobs still running: ['4937956', '4937966', '4937967', '4937969', '4937970', '4937971', '4937973', '4937974']. Waiting...
2026-06-21 20:04:34 INFO None 4937956: status FINISHED
2026-06-21 20:04:34 INFO None 4937957: status FINISHED
2026-06-21 20:04:34 INFO None 4937958: status FINISHED
2026-06-21 20:04:34 INFO None 4937959: status FINISHED
2026-06-21 20:04:34 INFO None 4937960: status FINISHED
2026-06-21 20:04:34 INFO None 4937961: status FINISHED
2026-06-21 20:04:34 INFO None 4937963: status FINISHED
2026-06-21 20:04:34 INFO None 4937964: status FINISHED
2026-06-21 20:04:34 INFO None 4937966: status FINISHED
2026-06-21 20:04:34 INFO None 4937967: status FINISHED
2026-06-21 20:04:34 INFO None 4937969: status FINISHED
2026-06-21 20:04:34 INFO None 4937970: status RUNNING/PENDING
2026-06-21 20:04:34 INFO None 4937971: status RUNNING/PENDING
2026-06-21 20:04:34 INFO None 4937973: status RUNNING/PENDING
2026-06-21 20:04:34 INFO None 4937974: status RUNNING/PENDING
2026-06-21 20:04:34 INFO Jobs still running: ['4937970', '4937971', '4937973', '4937974']. Waiting...
2026-06-21 20:04:49 INFO None 4937956: status FINISHED
2026-06-21 20:04:49 INFO None 4937957: status FINISHED
2026-06-21 20:04:49 INFO None 4937958: status FINISHED
2026-06-21 20:04:49 INFO None 4937959: status FINISHED
2026-06-21 20:04:49 INFO None 4937960: status FINISHED
2026-06-21 20:04:49 INFO None 4937961: status FINISHED
2026-06-21 20:04:49 INFO None 4937963: status FINISHED
2026-06-21 20:04:49 INFO None 4937964: status FINISHED
2026-06-21 20:04:49 INFO None 4937966: status FINISHED
2026-06-21 20:04:49 INFO None 4937967: status FINISHED
2026-06-21 20:04:49 INFO None 4937969: status FINISHED
2026-06-21 20:04:49 INFO None 4937970: status FINISHED
2026-06-21 20:04:49 INFO None 4937971: status RUNNING/PENDING
2026-06-21 20:04:49 INFO None 4937973: status RUNNING/PENDING
2026-06-21 20:04:49 INFO None 4937974: status RUNNING/PENDING
2026-06-21 20:04:49 INFO Jobs still running: ['4937971', '4937973', '4937974']. Waiting...
2026-06-21 20:05:04 INFO None 4937956: status FINISHED
2026-06-21 20:05:05 INFO None 4937957: status FINISHED
2026-06-21 20:05:05 INFO None 4937958: status FINISHED
2026-06-21 20:05:05 INFO None 4937959: status FINISHED
2026-06-21 20:05:05 INFO None 4937960: status FINISHED
2026-06-21 20:05:05 INFO None 4937961: status FINISHED
2026-06-21 20:05:05 INFO None 4937963: status FINISHED
2026-06-21 20:05:05 INFO None 4937964: status FINISHED
2026-06-21 20:05:05 INFO None 4937966: status FINISHED
2026-06-21 20:05:05 INFO None 4937967: status FINISHED
2026-06-21 20:05:05 INFO None 4937969: status FINISHED
2026-06-21 20:05:05 INFO None 4937970: status FINISHED
2026-06-21 20:05:05 INFO None 4937971: status RUNNING/PENDING
2026-06-21 20:05:05 INFO None 4937973: status FINISHED
2026-06-21 20:05:05 INFO None 4937974: status FINISHED
2026-06-21 20:05:05 INFO Jobs still running: ['4937971']. Waiting...
2026-06-21 20:05:20 INFO None 4937956: status FINISHED
2026-06-21 20:05:20 INFO None 4937957: status FINISHED
2026-06-21 20:05:20 INFO None 4937958: status FINISHED
2026-06-21 20:05:20 INFO None 4937959: status FINISHED
2026-06-21 20:05:20 INFO None 4937960: status FINISHED
2026-06-21 20:05:20 INFO None 4937961: status FINISHED
2026-06-21 20:05:20 INFO None 4937963: status FINISHED
2026-06-21 20:05:20 INFO None 4937964: status FINISHED
2026-06-21 20:05:20 INFO None 4937966: status FINISHED
2026-06-21 20:05:20 INFO None 4937967: status FINISHED
2026-06-21 20:05:20 INFO None 4937969: status FINISHED
2026-06-21 20:05:20 INFO None 4937970: status FINISHED
2026-06-21 20:05:20 INFO None 4937971: status RUNNING/PENDING
2026-06-21 20:05:20 INFO None 4937973: status FINISHED
2026-06-21 20:05:20 INFO None 4937974: status FINISHED
2026-06-21 20:05:20 INFO Jobs still running: ['4937971']. Waiting...
2026-06-21 20:05:35 INFO None 4937956: status FINISHED
2026-06-21 20:05:35 INFO None 4937957: status FINISHED
2026-06-21 20:05:35 INFO None 4937958: status FINISHED
2026-06-21 20:05:35 INFO None 4937959: status FINISHED
2026-06-21 20:05:35 INFO None 4937960: status FINISHED
2026-06-21 20:05:35 INFO None 4937961: status FINISHED
2026-06-21 20:05:35 INFO None 4937963: status FINISHED
2026-06-21 20:05:35 INFO None 4937964: status FINISHED
2026-06-21 20:05:35 INFO None 4937966: status FINISHED
2026-06-21 20:05:35 INFO None 4937967: status FINISHED
2026-06-21 20:05:35 INFO None 4937969: status FINISHED
2026-06-21 20:05:35 INFO None 4937970: status FINISHED
2026-06-21 20:05:35 INFO None 4937971: status FINISHED
2026-06-21 20:05:35 INFO None 4937973: status FINISHED
2026-06-21 20:05:35 INFO None 4937974: status FINISHED
2026-06-21 20:05:35 INFO Jobs ['4937956', '4937957', '4937958', '4937959', '4937960', '4937961', '4937963', '4937964', '4937966', '4937967', '4937969', '4937970', '4937971', '4937973', '4937974'] have finished
2026-06-21 20:05:35 INFO Checking restart files were created ...
2026-06-21 20:05:35 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/end.2020021508_2_ENS1.nc(1002685915 bytes)
2026-06-21 20:05:35 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/end.2020021508_2_ENS2.nc(1002685915 bytes)
2026-06-21 20:05:35 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/end.2020021508_2_ENS3.nc(1002685915 bytes)
2026-06-21 20:05:35 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/end.2020021508_2_ENS4.nc(1002685915 bytes)
2026-06-21 20:05:35 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/end.2020021508_2_ENS5.nc(1002685915 bytes)
2026-06-21 20:05:35 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/end.2020021508_2_ENS6.nc(1002685915 bytes)
2026-06-21 20:05:35 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/end.2020021508_2_ENS7.nc(1002685915 bytes)
2026-06-21 20:05:35 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/end.2020021508_2_ENS8.nc(1002685915 bytes)
2026-06-21 20:05:35 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/end.2020021508_2_ENS9.nc(1002685915 bytes)
2026-06-21 20:05:35 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/end.2020021508_2_ENS10.nc(1002685915 bytes)
2026-06-21 20:05:35 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/end.2020021508_2_ENS11.nc(1002685915 bytes)
2026-06-21 20:05:35 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/end.2020021508_2_ENS12.nc(1002685915 bytes)
2026-06-21 20:05:35 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/end.2020021508_2_ENS13.nc(1002685915 bytes)
2026-06-21 20:05:35 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/end.2020021508_2_ENS14.nc(1002685915 bytes)
2026-06-21 20:05:35 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/end.2020021508_2_ENS15.nc(1002685915 bytes)
2026-06-21 20:05:35 INFO  Run_model() completed successfully.
2026-06-21 20:05:35 INFO [TIME] after_model_set_simulated_time current_time=2020-02-15 08:00:00 simulated_time=2020-02-15 10:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 20:05:35 INFO [TIME] gregorian_conversion simulated_time=2020-02-15 10:00:00 days=153081 seconds=36000
2026-06-21 20:05:35 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-06-21 20:05:35 INFO [TIME] increment current_time 2020-02-15 08:00:00 -> 2020-02-15 10:00:00
2026-06-21 20:05:35 INFO [TIME] after_increment_before_assimilation current_time=2020-02-15 10:00:00 simulated_time=2020-02-15 10:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 20:05:35 INFO ---------->>> Running process_satellite_data()
2026-06-21 20:05:35 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12130.nc
2026-06-21 20:05:35 INFO ---------->>> Running run_obs_converter()
2026-06-21 20:05:35 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_34908_153081.out
2026-06-21 20:05:35 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_34908_153081.out
2026-06-21 20:05:35 INFO ---------->>> Running DART
2026-06-21 20:05:35 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-06-21 20:05:35 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/chim_ENS1_2020021508_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/chim_ENS1_2020021510_1_out_toDART.nc
2026-06-21 20:05:36 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/chim_ENS2_2020021508_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/chim_ENS2_2020021510_1_out_toDART.nc
2026-06-21 20:05:36 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/chim_ENS3_2020021508_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/chim_ENS3_2020021510_1_out_toDART.nc
2026-06-21 20:05:37 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/chim_ENS4_2020021508_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/chim_ENS4_2020021510_1_out_toDART.nc
2026-06-21 20:05:37 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/chim_ENS5_2020021508_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/chim_ENS5_2020021510_1_out_toDART.nc
2026-06-21 20:05:37 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/chim_ENS6_2020021508_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/chim_ENS6_2020021510_1_out_toDART.nc
2026-06-21 20:05:38 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/chim_ENS7_2020021508_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/chim_ENS7_2020021510_1_out_toDART.nc
2026-06-21 20:05:38 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/chim_ENS8_2020021508_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/chim_ENS8_2020021510_1_out_toDART.nc
2026-06-21 20:05:38 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/chim_ENS9_2020021508_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/chim_ENS9_2020021510_1_out_toDART.nc
2026-06-21 20:05:39 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/chim_ENS10_2020021508_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/chim_ENS10_2020021510_1_out_toDART.nc
2026-06-21 20:05:39 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/chim_ENS11_2020021508_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/chim_ENS11_2020021510_1_out_toDART.nc
2026-06-21 20:05:39 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/chim_ENS12_2020021508_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/chim_ENS12_2020021510_1_out_toDART.nc
2026-06-21 20:05:40 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/chim_ENS13_2020021508_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/chim_ENS13_2020021510_1_out_toDART.nc
2026-06-21 20:05:40 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/chim_ENS14_2020021508_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/chim_ENS14_2020021510_1_out_toDART.nc
2026-06-21 20:05:40 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/chim_ENS15_2020021508_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/chim_ENS15_2020021510_1_out_toDART.nc
2026-06-21 20:05:41 INFO Replacement input_template.nml → input.nml completed successfully.
2026-06-21 20:05:41 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-06-21 20:05:41 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-06-21 20:05:41 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-06-21 20:05:41 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-06-21 20:05:41 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-06-21 20:05:48 INFO Found: []
2026-06-21 20:05:48 INFO No job id returned by command ./run_filter.bsh
2026-06-21 20:05:48 INFO No monitoring will be performed
2026-06-21 20:05:48 INFO Moving DART output files to analysis and preassim directories for date 2020021510 if present ...
2026-06-21 20:05:48 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021510'
2026-06-21 20:05:48 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021510'
2026-06-21 20:05:48 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021510'
2026-06-21 20:05:48 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021510'
2026-06-21 20:05:48 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021510'
2026-06-21 20:05:48 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021510'
2026-06-21 20:05:48 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021510'
2026-06-21 20:05:48 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021510'
2026-06-21 20:05:48 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021510'
2026-06-21 20:05:48 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021510'
2026-06-21 20:05:48 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021510'
2026-06-21 20:05:48 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021510'
2026-06-21 20:05:48 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021510'
2026-06-21 20:05:48 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021510'
2026-06-21 20:05:48 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021510'
2026-06-21 20:05:48 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021510'
2026-06-21 20:05:48 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021510'
2026-06-21 20:05:48 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021510'
2026-06-21 20:05:48 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021510'
2026-06-21 20:05:48 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021510'
2026-06-21 20:05:48 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021510'
2026-06-21 20:05:48 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021510'
2026-06-21 20:05:48 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021510'
2026-06-21 20:05:48 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021510'
2026-06-21 20:05:48 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021510'
2026-06-21 20:05:49 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021510'
2026-06-21 20:05:49 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021510'
2026-06-21 20:05:49 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021510'
2026-06-21 20:05:49 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021510'
2026-06-21 20:05:49 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021510'
2026-06-21 20:05:49 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021510'
2026-06-21 20:05:49 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021510'
2026-06-21 20:05:49 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021510'
2026-06-21 20:05:49 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021510'
2026-06-21 20:05:49 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-06-21 20:05:49 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-06-21 20:05:49 INFO run_dart() is DONE.
2026-06-21 20:05:49 INFO ---------->>> Running update_pollutant_in_end()
2026-06-21 20:05:49 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/end.2020021508_2_ENS1.nc
2026-06-21 20:05:54 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:05:54 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:05:55 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:05:55 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:05:55 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:05:56 INFO [Posterior Diff ENS1] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021510/diff_posterior_ENS1_2020021510.nc
2026-06-21 20:05:56 INFO [Posterior Diff ENS1] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021510/diff_posterior_ENS1_2020021510.relative.nc
2026-06-21 20:05:56 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/end.2020021508_2_ENS2.nc
2026-06-21 20:06:01 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:06:01 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:06:02 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:06:02 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:06:03 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:06:03 INFO [Posterior Diff ENS2] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021510/diff_posterior_ENS2_2020021510.nc
2026-06-21 20:06:03 INFO [Posterior Diff ENS2] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021510/diff_posterior_ENS2_2020021510.relative.nc
2026-06-21 20:06:04 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/end.2020021508_2_ENS3.nc
2026-06-21 20:06:09 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:06:09 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:06:09 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:06:09 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:06:10 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:06:10 INFO [Posterior Diff ENS3] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021510/diff_posterior_ENS3_2020021510.nc
2026-06-21 20:06:10 INFO [Posterior Diff ENS3] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021510/diff_posterior_ENS3_2020021510.relative.nc
2026-06-21 20:06:10 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/end.2020021508_2_ENS4.nc
2026-06-21 20:06:16 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:06:16 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:06:17 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:06:17 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:06:18 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:06:18 INFO [Posterior Diff ENS4] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021510/diff_posterior_ENS4_2020021510.nc
2026-06-21 20:06:18 INFO [Posterior Diff ENS4] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021510/diff_posterior_ENS4_2020021510.relative.nc
2026-06-21 20:06:18 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/end.2020021508_2_ENS5.nc
2026-06-21 20:06:23 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:06:23 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:06:24 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:06:24 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:06:25 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:06:25 INFO [Posterior Diff ENS5] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021510/diff_posterior_ENS5_2020021510.nc
2026-06-21 20:06:25 INFO [Posterior Diff ENS5] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021510/diff_posterior_ENS5_2020021510.relative.nc
2026-06-21 20:06:26 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/end.2020021508_2_ENS6.nc
2026-06-21 20:06:31 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:06:31 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:06:32 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:06:32 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:06:33 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:06:33 INFO [Posterior Diff ENS6] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021510/diff_posterior_ENS6_2020021510.nc
2026-06-21 20:06:33 INFO [Posterior Diff ENS6] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021510/diff_posterior_ENS6_2020021510.relative.nc
2026-06-21 20:06:33 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/end.2020021508_2_ENS7.nc
2026-06-21 20:06:39 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:06:39 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:06:40 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:06:40 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:06:40 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:06:41 INFO [Posterior Diff ENS7] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021510/diff_posterior_ENS7_2020021510.nc
2026-06-21 20:06:41 INFO [Posterior Diff ENS7] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021510/diff_posterior_ENS7_2020021510.relative.nc
2026-06-21 20:06:41 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/end.2020021508_2_ENS8.nc
2026-06-21 20:06:46 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:06:46 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:06:47 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:06:47 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:06:48 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:06:48 INFO [Posterior Diff ENS8] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021510/diff_posterior_ENS8_2020021510.nc
2026-06-21 20:06:48 INFO [Posterior Diff ENS8] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021510/diff_posterior_ENS8_2020021510.relative.nc
2026-06-21 20:06:48 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/end.2020021508_2_ENS9.nc
2026-06-21 20:06:54 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:06:54 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:06:55 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:06:55 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:06:56 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:06:56 INFO [Posterior Diff ENS9] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021510/diff_posterior_ENS9_2020021510.nc
2026-06-21 20:06:56 INFO [Posterior Diff ENS9] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021510/diff_posterior_ENS9_2020021510.relative.nc
2026-06-21 20:06:56 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/end.2020021508_2_ENS10.nc
2026-06-21 20:07:02 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:07:02 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:07:03 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:07:03 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:07:03 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:07:04 INFO [Posterior Diff ENS10] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021510/diff_posterior_ENS10_2020021510.nc
2026-06-21 20:07:04 INFO [Posterior Diff ENS10] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021510/diff_posterior_ENS10_2020021510.relative.nc
2026-06-21 20:07:04 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/end.2020021508_2_ENS11.nc
2026-06-21 20:07:09 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:07:10 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:07:10 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:07:10 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:07:11 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:07:11 INFO [Posterior Diff ENS11] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021510/diff_posterior_ENS11_2020021510.nc
2026-06-21 20:07:11 INFO [Posterior Diff ENS11] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021510/diff_posterior_ENS11_2020021510.relative.nc
2026-06-21 20:07:11 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/end.2020021508_2_ENS12.nc
2026-06-21 20:07:17 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:07:17 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:07:18 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:07:18 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:07:19 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:07:19 INFO [Posterior Diff ENS12] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021510/diff_posterior_ENS12_2020021510.nc
2026-06-21 20:07:19 INFO [Posterior Diff ENS12] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021510/diff_posterior_ENS12_2020021510.relative.nc
2026-06-21 20:07:19 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/end.2020021508_2_ENS13.nc
2026-06-21 20:07:24 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:07:25 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:07:25 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:07:25 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:07:26 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:07:26 INFO [Posterior Diff ENS13] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021510/diff_posterior_ENS13_2020021510.nc
2026-06-21 20:07:26 INFO [Posterior Diff ENS13] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021510/diff_posterior_ENS13_2020021510.relative.nc
2026-06-21 20:07:27 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/end.2020021508_2_ENS14.nc
2026-06-21 20:07:32 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:07:32 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:07:33 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:07:33 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:07:34 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:07:34 INFO [Posterior Diff ENS14] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021510/diff_posterior_ENS14_2020021510.nc
2026-06-21 20:07:34 INFO [Posterior Diff ENS14] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021510/diff_posterior_ENS14_2020021510.relative.nc
2026-06-21 20:07:34 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/end.2020021508_2_ENS15.nc
2026-06-21 20:07:40 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:07:40 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:07:40 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:07:40 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:07:41 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:07:41 INFO [Posterior Diff ENS15] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021510/diff_posterior_ENS15_2020021510.nc
2026-06-21 20:07:41 INFO [Posterior Diff ENS15] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021510/diff_posterior_ENS15_2020021510.relative.nc
2026-06-21 20:07:41 INFO Next run starts from 2020-02-15 10:00:00
2026-06-21 20:07:41 INFO Cycle is DONE; starting a new loop!
2026-06-21 20:07:41 INFO [TIME] step_end current_time=2020-02-15 10:00:00 simulated_time=2020-02-15 10:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 20:07:41 INFO [TIME] step_start current_time=2020-02-15 10:00:00 simulated_time=2020-02-15 10:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 20:07:41 INFO [TIME] window start=2020-02-15 10:00:00 end=2020-02-15 11:00:00 run_hours=1 has_assimilation=True
2026-06-21 20:07:41 INFO Creating directories and links for ENS1 to run chimere's parallel part
2026-06-21 20:07:41 INFO Linking EMIS ...
2026-06-21 20:07:42 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens01.nc
2026-06-21 20:07:42 INFO >> Checking links...
2026-06-21 20:07:44 INFO >> All links are good for ENS1  ...
2026-06-21 20:07:44 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:07:45 INFO Hourly dataset computed and listing created
2026-06-21 20:07:47 INFO Hourly dataset computed
2026-06-21 20:07:47 INFO Creating directories and links for ENS2 to run chimere's parallel part
2026-06-21 20:07:47 INFO Linking EMIS ...
2026-06-21 20:07:47 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens02.nc
2026-06-21 20:07:47 INFO >> Checking links...
2026-06-21 20:07:50 INFO >> All links are good for ENS2  ...
2026-06-21 20:07:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:07:51 INFO Hourly dataset computed and listing created
2026-06-21 20:07:51 INFO Hourly dataset computed
2026-06-21 20:07:51 INFO Creating directories and links for ENS3 to run chimere's parallel part
2026-06-21 20:07:51 INFO Linking EMIS ...
2026-06-21 20:07:52 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens03.nc
2026-06-21 20:07:52 INFO >> Checking links...
2026-06-21 20:07:54 INFO >> All links are good for ENS3  ...
2026-06-21 20:07:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:07:55 INFO Hourly dataset computed and listing created
2026-06-21 20:07:58 INFO Hourly dataset computed
2026-06-21 20:07:58 INFO Creating directories and links for ENS4 to run chimere's parallel part
2026-06-21 20:07:58 INFO Linking EMIS ...
2026-06-21 20:07:58 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens04.nc
2026-06-21 20:07:58 INFO >> Checking links...
2026-06-21 20:08:01 INFO >> All links are good for ENS4  ...
2026-06-21 20:08:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:08:02 INFO Hourly dataset computed and listing created
2026-06-21 20:08:02 INFO Hourly dataset computed
2026-06-21 20:08:02 INFO Creating directories and links for ENS5 to run chimere's parallel part
2026-06-21 20:08:02 INFO Linking EMIS ...
2026-06-21 20:08:03 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens05.nc
2026-06-21 20:08:03 INFO >> Checking links...
2026-06-21 20:08:05 INFO >> All links are good for ENS5  ...
2026-06-21 20:08:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:08:06 INFO Hourly dataset computed and listing created
2026-06-21 20:08:06 INFO Hourly dataset computed
2026-06-21 20:08:06 INFO Creating directories and links for ENS6 to run chimere's parallel part
2026-06-21 20:08:06 INFO Linking EMIS ...
2026-06-21 20:08:07 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens06.nc
2026-06-21 20:08:07 INFO >> Checking links...
2026-06-21 20:08:09 INFO >> All links are good for ENS6  ...
2026-06-21 20:08:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:08:10 INFO Hourly dataset computed and listing created
2026-06-21 20:08:11 INFO Hourly dataset computed
2026-06-21 20:08:11 INFO Creating directories and links for ENS7 to run chimere's parallel part
2026-06-21 20:08:11 INFO Linking EMIS ...
2026-06-21 20:08:11 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens07.nc
2026-06-21 20:08:11 INFO >> Checking links...
2026-06-21 20:08:14 INFO >> All links are good for ENS7  ...
2026-06-21 20:08:14 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:08:15 INFO Hourly dataset computed and listing created
2026-06-21 20:08:15 INFO Hourly dataset computed
2026-06-21 20:08:15 INFO Creating directories and links for ENS8 to run chimere's parallel part
2026-06-21 20:08:15 INFO Linking EMIS ...
2026-06-21 20:08:16 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens08.nc
2026-06-21 20:08:16 INFO >> Checking links...
2026-06-21 20:08:18 INFO >> All links are good for ENS8  ...
2026-06-21 20:08:18 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:08:19 INFO Hourly dataset computed and listing created
2026-06-21 20:08:20 INFO Hourly dataset computed
2026-06-21 20:08:20 INFO Creating directories and links for ENS9 to run chimere's parallel part
2026-06-21 20:08:20 INFO Linking EMIS ...
2026-06-21 20:08:20 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens09.nc
2026-06-21 20:08:20 INFO >> Checking links...
2026-06-21 20:08:22 INFO >> All links are good for ENS9  ...
2026-06-21 20:08:22 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:08:23 INFO Hourly dataset computed and listing created
2026-06-21 20:08:24 INFO Hourly dataset computed
2026-06-21 20:08:24 INFO Creating directories and links for ENS10 to run chimere's parallel part
2026-06-21 20:08:24 INFO Linking EMIS ...
2026-06-21 20:08:24 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens10.nc
2026-06-21 20:08:24 INFO >> Checking links...
2026-06-21 20:08:27 INFO >> All links are good for ENS10  ...
2026-06-21 20:08:27 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:08:27 INFO Hourly dataset computed and listing created
2026-06-21 20:08:28 INFO Hourly dataset computed
2026-06-21 20:08:28 INFO Creating directories and links for ENS11 to run chimere's parallel part
2026-06-21 20:08:28 INFO Linking EMIS ...
2026-06-21 20:08:28 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens11.nc
2026-06-21 20:08:28 INFO >> Checking links...
2026-06-21 20:08:31 INFO >> All links are good for ENS11  ...
2026-06-21 20:08:31 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:08:32 INFO Hourly dataset computed and listing created
2026-06-21 20:08:32 INFO Hourly dataset computed
2026-06-21 20:08:32 INFO Creating directories and links for ENS12 to run chimere's parallel part
2026-06-21 20:08:32 INFO Linking EMIS ...
2026-06-21 20:08:33 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens12.nc
2026-06-21 20:08:33 INFO >> Checking links...
2026-06-21 20:08:35 INFO >> All links are good for ENS12  ...
2026-06-21 20:08:35 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:08:36 INFO Hourly dataset computed and listing created
2026-06-21 20:08:39 INFO Hourly dataset computed
2026-06-21 20:08:39 INFO Creating directories and links for ENS13 to run chimere's parallel part
2026-06-21 20:08:39 INFO Linking EMIS ...
2026-06-21 20:08:40 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens13.nc
2026-06-21 20:08:40 INFO >> Checking links...
2026-06-21 20:08:42 INFO >> All links are good for ENS13  ...
2026-06-21 20:08:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:08:43 INFO Hourly dataset computed and listing created
2026-06-21 20:08:46 INFO Hourly dataset computed
2026-06-21 20:08:46 INFO Creating directories and links for ENS14 to run chimere's parallel part
2026-06-21 20:08:46 INFO Linking EMIS ...
2026-06-21 20:08:46 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens14.nc
2026-06-21 20:08:46 INFO >> Checking links...
2026-06-21 20:08:49 INFO >> All links are good for ENS14  ...
2026-06-21 20:08:49 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:08:49 INFO Hourly dataset computed and listing created
2026-06-21 20:08:52 INFO Hourly dataset computed
2026-06-21 20:08:52 INFO Creating directories and links for ENS15 to run chimere's parallel part
2026-06-21 20:08:52 INFO Linking EMIS ...
2026-06-21 20:08:52 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens15.nc
2026-06-21 20:08:52 INFO >> Checking links...
2026-06-21 20:08:55 INFO >> All links are good for ENS15  ...
2026-06-21 20:08:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:08:56 INFO Hourly dataset computed and listing created
2026-06-21 20:08:56 INFO Hourly dataset computed
2026-06-21 20:08:56 INFO ---------->>> Running CHIMERE model from 2020-02-15 10:00:00 to 2020-02-15 11:00:00
2026-06-21 20:08:56 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:08:56 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1
2026-06-21 20:08:56 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/end.2020021508_2_ENS1.nc
2026-06-21 20:08:56 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-06-21 20:08:56 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:08:56 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-06-21 20:08:56 INFO Queuing job for member 1...
2026-06-21 20:08:56 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:08:56 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-06-21 20:08:57 INFO Found: ['4938018']
2026-06-21 20:09:02 INFO [TGCC-IRENE] Submitted job with ID:['4938018']
2026-06-21 20:09:02 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:09:02 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2
2026-06-21 20:09:02 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/end.2020021508_2_ENS2.nc
2026-06-21 20:09:02 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-06-21 20:09:02 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:09:02 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-06-21 20:09:02 INFO Queuing job for member 2...
2026-06-21 20:09:02 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:09:02 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-06-21 20:09:03 INFO Found: ['4938020']
2026-06-21 20:09:08 INFO [TGCC-IRENE] Submitted job with ID:['4938020']
2026-06-21 20:09:08 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:09:08 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3
2026-06-21 20:09:08 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/end.2020021508_2_ENS3.nc
2026-06-21 20:09:08 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-06-21 20:09:08 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:09:08 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-06-21 20:09:08 INFO Queuing job for member 3...
2026-06-21 20:09:08 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:09:08 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-06-21 20:09:09 INFO Found: ['4938021']
2026-06-21 20:09:14 INFO [TGCC-IRENE] Submitted job with ID:['4938021']
2026-06-21 20:09:14 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:09:14 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4
2026-06-21 20:09:14 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/end.2020021508_2_ENS4.nc
2026-06-21 20:09:14 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-06-21 20:09:14 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:09:14 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-06-21 20:09:14 INFO Queuing job for member 4...
2026-06-21 20:09:14 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:09:14 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-06-21 20:09:14 INFO Found: ['4938022']
2026-06-21 20:09:19 INFO [TGCC-IRENE] Submitted job with ID:['4938022']
2026-06-21 20:09:19 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:09:19 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5
2026-06-21 20:09:19 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/end.2020021508_2_ENS5.nc
2026-06-21 20:09:19 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-06-21 20:09:19 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:09:19 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-06-21 20:09:19 INFO Queuing job for member 5...
2026-06-21 20:09:19 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:09:19 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-06-21 20:09:20 INFO Found: ['4938023']
2026-06-21 20:09:25 INFO [TGCC-IRENE] Submitted job with ID:['4938023']
2026-06-21 20:09:25 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:09:25 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6
2026-06-21 20:09:25 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/end.2020021508_2_ENS6.nc
2026-06-21 20:09:25 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-06-21 20:09:25 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:09:25 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-06-21 20:09:25 INFO Queuing job for member 6...
2026-06-21 20:09:25 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:09:25 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-06-21 20:09:26 INFO Found: ['4938024']
2026-06-21 20:09:31 INFO [TGCC-IRENE] Submitted job with ID:['4938024']
2026-06-21 20:09:31 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:09:31 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7
2026-06-21 20:09:31 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/end.2020021508_2_ENS7.nc
2026-06-21 20:09:31 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-06-21 20:09:31 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:09:31 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-06-21 20:09:31 INFO Queuing job for member 7...
2026-06-21 20:09:31 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:09:31 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-06-21 20:09:32 INFO Found: ['4938025']
2026-06-21 20:09:37 INFO [TGCC-IRENE] Submitted job with ID:['4938025']
2026-06-21 20:09:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:09:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8
2026-06-21 20:09:37 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/end.2020021508_2_ENS8.nc
2026-06-21 20:09:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-06-21 20:09:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:09:37 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-06-21 20:09:37 INFO Queuing job for member 8...
2026-06-21 20:09:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:09:37 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-06-21 20:09:37 INFO Found: ['4938026']
2026-06-21 20:09:42 INFO [TGCC-IRENE] Submitted job with ID:['4938026']
2026-06-21 20:09:42 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:09:42 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9
2026-06-21 20:09:42 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/end.2020021508_2_ENS9.nc
2026-06-21 20:09:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-06-21 20:09:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:09:43 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-06-21 20:09:43 INFO Queuing job for member 9...
2026-06-21 20:09:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:09:43 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-06-21 20:09:44 INFO Found: ['4938027']
2026-06-21 20:09:49 INFO [TGCC-IRENE] Submitted job with ID:['4938027']
2026-06-21 20:09:49 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:09:49 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10
2026-06-21 20:09:49 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/end.2020021508_2_ENS10.nc
2026-06-21 20:09:49 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-06-21 20:09:49 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:09:49 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-06-21 20:09:49 INFO Queuing job for member 10...
2026-06-21 20:09:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:09:49 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-06-21 20:09:49 INFO Found: ['4938028']
2026-06-21 20:09:54 INFO [TGCC-IRENE] Submitted job with ID:['4938028']
2026-06-21 20:09:54 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:09:54 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11
2026-06-21 20:09:54 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/end.2020021508_2_ENS11.nc
2026-06-21 20:09:54 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-06-21 20:09:54 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:09:54 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-06-21 20:09:54 INFO Queuing job for member 11...
2026-06-21 20:09:54 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:09:54 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-06-21 20:09:55 INFO Found: ['4938029']
2026-06-21 20:10:00 INFO [TGCC-IRENE] Submitted job with ID:['4938029']
2026-06-21 20:10:00 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:10:00 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12
2026-06-21 20:10:00 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/end.2020021508_2_ENS12.nc
2026-06-21 20:10:00 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-06-21 20:10:00 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:10:00 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-06-21 20:10:00 INFO Queuing job for member 12...
2026-06-21 20:10:00 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:10:00 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-06-21 20:10:01 INFO Found: ['4938031']
2026-06-21 20:10:06 INFO [TGCC-IRENE] Submitted job with ID:['4938031']
2026-06-21 20:10:06 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:10:06 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13
2026-06-21 20:10:06 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/end.2020021508_2_ENS13.nc
2026-06-21 20:10:06 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-06-21 20:10:06 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:10:06 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-06-21 20:10:06 INFO Queuing job for member 13...
2026-06-21 20:10:06 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:10:06 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-06-21 20:10:07 INFO Found: ['4938034']
2026-06-21 20:10:12 INFO [TGCC-IRENE] Submitted job with ID:['4938034']
2026-06-21 20:10:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:10:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14
2026-06-21 20:10:12 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/end.2020021508_2_ENS14.nc
2026-06-21 20:10:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-06-21 20:10:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:10:12 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-06-21 20:10:12 INFO Queuing job for member 14...
2026-06-21 20:10:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:10:12 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-06-21 20:10:13 INFO Found: ['4938036']
2026-06-21 20:10:18 INFO [TGCC-IRENE] Submitted job with ID:['4938036']
2026-06-21 20:10:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:10:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15
2026-06-21 20:10:18 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/end.2020021508_2_ENS15.nc
2026-06-21 20:10:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-06-21 20:10:18 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:10:18 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-06-21 20:10:18 INFO Queuing job for member 15...
2026-06-21 20:10:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:10:18 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-06-21 20:10:19 INFO Found: ['4938037']
2026-06-21 20:10:24 INFO [TGCC-IRENE] Submitted job with ID:['4938037']
2026-06-21 20:10:24 INFO Checking job status ...
2026-06-21 20:10:24 INFO None 4938018: status RUNNING/PENDING
2026-06-21 20:10:24 INFO None 4938020: status RUNNING/PENDING
2026-06-21 20:10:24 INFO None 4938021: status RUNNING/PENDING
2026-06-21 20:10:24 INFO None 4938022: status RUNNING/PENDING
2026-06-21 20:10:24 INFO None 4938023: status RUNNING/PENDING
2026-06-21 20:10:24 INFO None 4938024: status RUNNING/PENDING
2026-06-21 20:10:24 INFO None 4938025: status RUNNING/PENDING
2026-06-21 20:10:24 INFO None 4938026: status RUNNING/PENDING
2026-06-21 20:10:24 INFO None 4938027: status RUNNING/PENDING
2026-06-21 20:10:24 INFO None 4938028: status RUNNING/PENDING
2026-06-21 20:10:24 INFO None 4938029: status RUNNING/PENDING
2026-06-21 20:10:24 INFO None 4938031: status RUNNING/PENDING
2026-06-21 20:10:24 INFO None 4938034: status RUNNING/PENDING
2026-06-21 20:10:24 INFO None 4938036: status RUNNING/PENDING
2026-06-21 20:10:24 INFO None 4938037: status RUNNING/PENDING
2026-06-21 20:10:24 INFO Jobs still running: ['4938018', '4938020', '4938021', '4938022', '4938023', '4938024', '4938025', '4938026', '4938027', '4938028', '4938029', '4938031', '4938034', '4938036', '4938037']. Waiting...
2026-06-21 20:10:39 INFO None 4938018: status RUNNING/PENDING
2026-06-21 20:10:39 INFO None 4938020: status RUNNING/PENDING
2026-06-21 20:10:39 INFO None 4938021: status RUNNING/PENDING
2026-06-21 20:10:39 INFO None 4938022: status RUNNING/PENDING
2026-06-21 20:10:39 INFO None 4938023: status RUNNING/PENDING
2026-06-21 20:10:39 INFO None 4938024: status RUNNING/PENDING
2026-06-21 20:10:39 INFO None 4938025: status RUNNING/PENDING
2026-06-21 20:10:39 INFO None 4938026: status RUNNING/PENDING
2026-06-21 20:10:39 INFO None 4938027: status RUNNING/PENDING
2026-06-21 20:10:39 INFO None 4938028: status RUNNING/PENDING
2026-06-21 20:10:39 INFO None 4938029: status RUNNING/PENDING
2026-06-21 20:10:39 INFO None 4938031: status RUNNING/PENDING
2026-06-21 20:10:39 INFO None 4938034: status RUNNING/PENDING
2026-06-21 20:10:39 INFO None 4938036: status RUNNING/PENDING
2026-06-21 20:10:39 INFO None 4938037: status RUNNING/PENDING
2026-06-21 20:10:39 INFO Jobs still running: ['4938018', '4938020', '4938021', '4938022', '4938023', '4938024', '4938025', '4938026', '4938027', '4938028', '4938029', '4938031', '4938034', '4938036', '4938037']. Waiting...
2026-06-21 20:10:54 INFO None 4938018: status RUNNING/PENDING
2026-06-21 20:10:54 INFO None 4938020: status RUNNING/PENDING
2026-06-21 20:10:54 INFO None 4938021: status RUNNING/PENDING
2026-06-21 20:10:54 INFO None 4938022: status RUNNING/PENDING
2026-06-21 20:10:54 INFO None 4938023: status RUNNING/PENDING
2026-06-21 20:10:54 INFO None 4938024: status RUNNING/PENDING
2026-06-21 20:10:54 INFO None 4938025: status RUNNING/PENDING
2026-06-21 20:10:54 INFO None 4938026: status RUNNING/PENDING
2026-06-21 20:10:54 INFO None 4938027: status RUNNING/PENDING
2026-06-21 20:10:55 INFO None 4938028: status RUNNING/PENDING
2026-06-21 20:10:55 INFO None 4938029: status RUNNING/PENDING
2026-06-21 20:10:55 INFO None 4938031: status RUNNING/PENDING
2026-06-21 20:10:55 INFO None 4938034: status RUNNING/PENDING
2026-06-21 20:10:55 INFO None 4938036: status RUNNING/PENDING
2026-06-21 20:10:55 INFO None 4938037: status RUNNING/PENDING
2026-06-21 20:10:55 INFO Jobs still running: ['4938018', '4938020', '4938021', '4938022', '4938023', '4938024', '4938025', '4938026', '4938027', '4938028', '4938029', '4938031', '4938034', '4938036', '4938037']. Waiting...
2026-06-21 20:11:10 INFO None 4938018: status RUNNING/PENDING
2026-06-21 20:11:10 INFO None 4938020: status RUNNING/PENDING
2026-06-21 20:11:10 INFO None 4938021: status FINISHED
2026-06-21 20:11:10 INFO None 4938022: status RUNNING/PENDING
2026-06-21 20:11:10 INFO None 4938023: status RUNNING/PENDING
2026-06-21 20:11:10 INFO None 4938024: status RUNNING/PENDING
2026-06-21 20:11:10 INFO None 4938025: status RUNNING/PENDING
2026-06-21 20:11:10 INFO None 4938026: status RUNNING/PENDING
2026-06-21 20:11:10 INFO None 4938027: status RUNNING/PENDING
2026-06-21 20:11:10 INFO None 4938028: status RUNNING/PENDING
2026-06-21 20:11:10 INFO None 4938029: status RUNNING/PENDING
2026-06-21 20:11:10 INFO None 4938031: status RUNNING/PENDING
2026-06-21 20:11:10 INFO None 4938034: status RUNNING/PENDING
2026-06-21 20:11:10 INFO None 4938036: status RUNNING/PENDING
2026-06-21 20:11:10 INFO None 4938037: status RUNNING/PENDING
2026-06-21 20:11:10 INFO Jobs still running: ['4938018', '4938020', '4938022', '4938023', '4938024', '4938025', '4938026', '4938027', '4938028', '4938029', '4938031', '4938034', '4938036', '4938037']. Waiting...
2026-06-21 20:11:25 INFO None 4938018: status RUNNING/PENDING
2026-06-21 20:11:25 INFO None 4938020: status RUNNING/PENDING
2026-06-21 20:11:25 INFO None 4938021: status FINISHED
2026-06-21 20:11:25 INFO None 4938022: status FINISHED
2026-06-21 20:11:25 INFO None 4938023: status RUNNING/PENDING
2026-06-21 20:11:25 INFO None 4938024: status RUNNING/PENDING
2026-06-21 20:11:25 INFO None 4938025: status RUNNING/PENDING
2026-06-21 20:11:25 INFO None 4938026: status RUNNING/PENDING
2026-06-21 20:11:25 INFO None 4938027: status RUNNING/PENDING
2026-06-21 20:11:25 INFO None 4938028: status RUNNING/PENDING
2026-06-21 20:11:25 INFO None 4938029: status RUNNING/PENDING
2026-06-21 20:11:25 INFO None 4938031: status RUNNING/PENDING
2026-06-21 20:11:25 INFO None 4938034: status RUNNING/PENDING
2026-06-21 20:11:25 INFO None 4938036: status RUNNING/PENDING
2026-06-21 20:11:25 INFO None 4938037: status RUNNING/PENDING
2026-06-21 20:11:25 INFO Jobs still running: ['4938018', '4938020', '4938023', '4938024', '4938025', '4938026', '4938027', '4938028', '4938029', '4938031', '4938034', '4938036', '4938037']. Waiting...
2026-06-21 20:11:40 INFO None 4938018: status RUNNING/PENDING
2026-06-21 20:11:40 INFO None 4938020: status RUNNING/PENDING
2026-06-21 20:11:40 INFO None 4938021: status FINISHED
2026-06-21 20:11:40 INFO None 4938022: status FINISHED
2026-06-21 20:11:40 INFO None 4938023: status RUNNING/PENDING
2026-06-21 20:11:40 INFO None 4938024: status FINISHED
2026-06-21 20:11:40 INFO None 4938025: status RUNNING/PENDING
2026-06-21 20:11:40 INFO None 4938026: status RUNNING/PENDING
2026-06-21 20:11:40 INFO None 4938027: status RUNNING/PENDING
2026-06-21 20:11:40 INFO None 4938028: status RUNNING/PENDING
2026-06-21 20:11:40 INFO None 4938029: status RUNNING/PENDING
2026-06-21 20:11:40 INFO None 4938031: status RUNNING/PENDING
2026-06-21 20:11:40 INFO None 4938034: status RUNNING/PENDING
2026-06-21 20:11:40 INFO None 4938036: status RUNNING/PENDING
2026-06-21 20:11:41 INFO None 4938037: status RUNNING/PENDING
2026-06-21 20:11:41 INFO Jobs still running: ['4938018', '4938020', '4938023', '4938025', '4938026', '4938027', '4938028', '4938029', '4938031', '4938034', '4938036', '4938037']. Waiting...
2026-06-21 20:11:56 INFO None 4938018: status RUNNING/PENDING
2026-06-21 20:11:56 INFO None 4938020: status FINISHED
2026-06-21 20:11:56 INFO None 4938021: status FINISHED
2026-06-21 20:11:56 INFO None 4938022: status FINISHED
2026-06-21 20:11:56 INFO None 4938023: status RUNNING/PENDING
2026-06-21 20:11:56 INFO None 4938024: status FINISHED
2026-06-21 20:11:56 INFO None 4938025: status RUNNING/PENDING
2026-06-21 20:11:56 INFO None 4938026: status RUNNING/PENDING
2026-06-21 20:11:56 INFO None 4938027: status RUNNING/PENDING
2026-06-21 20:11:56 INFO None 4938028: status RUNNING/PENDING
2026-06-21 20:11:56 INFO None 4938029: status RUNNING/PENDING
2026-06-21 20:11:56 INFO None 4938031: status RUNNING/PENDING
2026-06-21 20:11:56 INFO None 4938034: status RUNNING/PENDING
2026-06-21 20:11:56 INFO None 4938036: status RUNNING/PENDING
2026-06-21 20:11:56 INFO None 4938037: status RUNNING/PENDING
2026-06-21 20:11:56 INFO Jobs still running: ['4938018', '4938023', '4938025', '4938026', '4938027', '4938028', '4938029', '4938031', '4938034', '4938036', '4938037']. Waiting...
2026-06-21 20:12:11 INFO None 4938018: status RUNNING/PENDING
2026-06-21 20:12:11 INFO None 4938020: status FINISHED
2026-06-21 20:12:11 INFO None 4938021: status FINISHED
2026-06-21 20:12:11 INFO None 4938022: status FINISHED
2026-06-21 20:12:11 INFO None 4938023: status FINISHED
2026-06-21 20:12:11 INFO None 4938024: status FINISHED
2026-06-21 20:12:11 INFO None 4938025: status RUNNING/PENDING
2026-06-21 20:12:11 INFO None 4938026: status RUNNING/PENDING
2026-06-21 20:12:11 INFO None 4938027: status RUNNING/PENDING
2026-06-21 20:12:11 INFO None 4938028: status RUNNING/PENDING
2026-06-21 20:12:11 INFO None 4938029: status RUNNING/PENDING
2026-06-21 20:12:11 INFO None 4938031: status RUNNING/PENDING
2026-06-21 20:12:11 INFO None 4938034: status RUNNING/PENDING
2026-06-21 20:12:11 INFO None 4938036: status FINISHED
2026-06-21 20:12:11 INFO None 4938037: status RUNNING/PENDING
2026-06-21 20:12:11 INFO Jobs still running: ['4938018', '4938025', '4938026', '4938027', '4938028', '4938029', '4938031', '4938034', '4938037']. Waiting...
2026-06-21 20:12:26 INFO None 4938018: status FINISHED
2026-06-21 20:12:26 INFO None 4938020: status FINISHED
2026-06-21 20:12:26 INFO None 4938021: status FINISHED
2026-06-21 20:12:26 INFO None 4938022: status FINISHED
2026-06-21 20:12:26 INFO None 4938023: status FINISHED
2026-06-21 20:12:26 INFO None 4938024: status FINISHED
2026-06-21 20:12:26 INFO None 4938025: status RUNNING/PENDING
2026-06-21 20:12:26 INFO None 4938026: status RUNNING/PENDING
2026-06-21 20:12:26 INFO None 4938027: status RUNNING/PENDING
2026-06-21 20:12:26 INFO None 4938028: status RUNNING/PENDING
2026-06-21 20:12:26 INFO None 4938029: status RUNNING/PENDING
2026-06-21 20:12:26 INFO None 4938031: status RUNNING/PENDING
2026-06-21 20:12:26 INFO None 4938034: status FINISHED
2026-06-21 20:12:26 INFO None 4938036: status FINISHED
2026-06-21 20:12:26 INFO None 4938037: status RUNNING/PENDING
2026-06-21 20:12:26 INFO Jobs still running: ['4938025', '4938026', '4938027', '4938028', '4938029', '4938031', '4938037']. Waiting...
2026-06-21 20:12:41 INFO None 4938018: status FINISHED
2026-06-21 20:12:41 INFO None 4938020: status FINISHED
2026-06-21 20:12:41 INFO None 4938021: status FINISHED
2026-06-21 20:12:41 INFO None 4938022: status FINISHED
2026-06-21 20:12:42 INFO None 4938023: status FINISHED
2026-06-21 20:12:42 INFO None 4938024: status FINISHED
2026-06-21 20:12:42 INFO None 4938025: status RUNNING/PENDING
2026-06-21 20:12:42 INFO None 4938026: status RUNNING/PENDING
2026-06-21 20:12:42 INFO None 4938027: status RUNNING/PENDING
2026-06-21 20:12:42 INFO None 4938028: status RUNNING/PENDING
2026-06-21 20:12:42 INFO None 4938029: status RUNNING/PENDING
2026-06-21 20:12:42 INFO None 4938031: status RUNNING/PENDING
2026-06-21 20:12:42 INFO None 4938034: status FINISHED
2026-06-21 20:12:42 INFO None 4938036: status FINISHED
2026-06-21 20:12:42 INFO None 4938037: status RUNNING/PENDING
2026-06-21 20:12:42 INFO Jobs still running: ['4938025', '4938026', '4938027', '4938028', '4938029', '4938031', '4938037']. Waiting...
2026-06-21 20:12:57 INFO None 4938018: status FINISHED
2026-06-21 20:12:57 INFO None 4938020: status FINISHED
2026-06-21 20:12:57 INFO None 4938021: status FINISHED
2026-06-21 20:12:57 INFO None 4938022: status FINISHED
2026-06-21 20:12:57 INFO None 4938023: status FINISHED
2026-06-21 20:12:57 INFO None 4938024: status FINISHED
2026-06-21 20:12:57 INFO None 4938025: status RUNNING/PENDING
2026-06-21 20:12:57 INFO None 4938026: status FINISHED
2026-06-21 20:12:57 INFO None 4938027: status RUNNING/PENDING
2026-06-21 20:12:57 INFO None 4938028: status RUNNING/PENDING
2026-06-21 20:12:57 INFO None 4938029: status RUNNING/PENDING
2026-06-21 20:12:57 INFO None 4938031: status FINISHED
2026-06-21 20:12:57 INFO None 4938034: status FINISHED
2026-06-21 20:12:57 INFO None 4938036: status FINISHED
2026-06-21 20:12:57 INFO None 4938037: status RUNNING/PENDING
2026-06-21 20:12:57 INFO Jobs still running: ['4938025', '4938027', '4938028', '4938029', '4938037']. Waiting...
2026-06-21 20:13:12 INFO None 4938018: status FINISHED
2026-06-21 20:13:12 INFO None 4938020: status FINISHED
2026-06-21 20:13:12 INFO None 4938021: status FINISHED
2026-06-21 20:13:12 INFO None 4938022: status FINISHED
2026-06-21 20:13:12 INFO None 4938023: status FINISHED
2026-06-21 20:13:12 INFO None 4938024: status FINISHED
2026-06-21 20:13:12 INFO None 4938025: status FINISHED
2026-06-21 20:13:12 INFO None 4938026: status FINISHED
2026-06-21 20:13:12 INFO None 4938027: status FINISHED
2026-06-21 20:13:12 INFO None 4938028: status FINISHED
2026-06-21 20:13:12 INFO None 4938029: status FINISHED
2026-06-21 20:13:12 INFO None 4938031: status FINISHED
2026-06-21 20:13:12 INFO None 4938034: status FINISHED
2026-06-21 20:13:12 INFO None 4938036: status FINISHED
2026-06-21 20:13:12 INFO None 4938037: status RUNNING/PENDING
2026-06-21 20:13:12 INFO Jobs still running: ['4938037']. Waiting...
2026-06-21 20:13:27 INFO None 4938018: status FINISHED
2026-06-21 20:13:27 INFO None 4938020: status FINISHED
2026-06-21 20:13:27 INFO None 4938021: status FINISHED
2026-06-21 20:13:27 INFO None 4938022: status FINISHED
2026-06-21 20:13:27 INFO None 4938023: status FINISHED
2026-06-21 20:13:27 INFO None 4938024: status FINISHED
2026-06-21 20:13:27 INFO None 4938025: status FINISHED
2026-06-21 20:13:27 INFO None 4938026: status FINISHED
2026-06-21 20:13:27 INFO None 4938027: status FINISHED
2026-06-21 20:13:27 INFO None 4938028: status FINISHED
2026-06-21 20:13:28 INFO None 4938029: status FINISHED
2026-06-21 20:13:28 INFO None 4938031: status FINISHED
2026-06-21 20:13:28 INFO None 4938034: status FINISHED
2026-06-21 20:13:28 INFO None 4938036: status FINISHED
2026-06-21 20:13:28 INFO None 4938037: status FINISHED
2026-06-21 20:13:28 INFO Jobs ['4938018', '4938020', '4938021', '4938022', '4938023', '4938024', '4938025', '4938026', '4938027', '4938028', '4938029', '4938031', '4938034', '4938036', '4938037'] have finished
2026-06-21 20:13:28 INFO Checking restart files were created ...
2026-06-21 20:13:28 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/end.2020021510_1_ENS1.nc(668832435 bytes)
2026-06-21 20:13:28 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/end.2020021510_1_ENS2.nc(668832435 bytes)
2026-06-21 20:13:28 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/end.2020021510_1_ENS3.nc(668832435 bytes)
2026-06-21 20:13:28 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/end.2020021510_1_ENS4.nc(668832435 bytes)
2026-06-21 20:13:28 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/end.2020021510_1_ENS5.nc(668832435 bytes)
2026-06-21 20:13:28 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/end.2020021510_1_ENS6.nc(668832435 bytes)
2026-06-21 20:13:28 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/end.2020021510_1_ENS7.nc(668832435 bytes)
2026-06-21 20:13:28 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/end.2020021510_1_ENS8.nc(668832435 bytes)
2026-06-21 20:13:28 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/end.2020021510_1_ENS9.nc(668832435 bytes)
2026-06-21 20:13:28 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/end.2020021510_1_ENS10.nc(668832435 bytes)
2026-06-21 20:13:28 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/end.2020021510_1_ENS11.nc(668832435 bytes)
2026-06-21 20:13:28 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/end.2020021510_1_ENS12.nc(668832435 bytes)
2026-06-21 20:13:28 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/end.2020021510_1_ENS13.nc(668832435 bytes)
2026-06-21 20:13:28 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/end.2020021510_1_ENS14.nc(668832435 bytes)
2026-06-21 20:13:28 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/end.2020021510_1_ENS15.nc(668832435 bytes)
2026-06-21 20:13:28 INFO  Run_model() completed successfully.
2026-06-21 20:13:28 INFO [TIME] after_model_set_simulated_time current_time=2020-02-15 10:00:00 simulated_time=2020-02-15 11:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 20:13:28 INFO [TIME] gregorian_conversion simulated_time=2020-02-15 11:00:00 days=153081 seconds=39600
2026-06-21 20:13:28 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-06-21 20:13:28 INFO [TIME] increment current_time 2020-02-15 10:00:00 -> 2020-02-15 11:00:00
2026-06-21 20:13:28 INFO [TIME] after_increment_before_assimilation current_time=2020-02-15 11:00:00 simulated_time=2020-02-15 11:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 20:13:28 INFO ---------->>> Running process_satellite_data()
2026-06-21 20:13:28 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12131.nc
2026-06-21 20:13:28 INFO ---------->>> Running run_obs_converter()
2026-06-21 20:13:28 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_40949_153081.out
2026-06-21 20:13:28 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_40949_153081.out
2026-06-21 20:13:28 INFO ---------->>> Running DART
2026-06-21 20:13:28 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-06-21 20:13:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/chim_ENS1_2020021510_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/chim_ENS1_2020021511_1_out_toDART.nc
2026-06-21 20:13:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/chim_ENS2_2020021510_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/chim_ENS2_2020021511_1_out_toDART.nc
2026-06-21 20:13:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/chim_ENS3_2020021510_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/chim_ENS3_2020021511_1_out_toDART.nc
2026-06-21 20:13:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/chim_ENS4_2020021510_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/chim_ENS4_2020021511_1_out_toDART.nc
2026-06-21 20:13:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/chim_ENS5_2020021510_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/chim_ENS5_2020021511_1_out_toDART.nc
2026-06-21 20:13:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/chim_ENS6_2020021510_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/chim_ENS6_2020021511_1_out_toDART.nc
2026-06-21 20:13:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/chim_ENS7_2020021510_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/chim_ENS7_2020021511_1_out_toDART.nc
2026-06-21 20:13:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/chim_ENS8_2020021510_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/chim_ENS8_2020021511_1_out_toDART.nc
2026-06-21 20:13:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/chim_ENS9_2020021510_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/chim_ENS9_2020021511_1_out_toDART.nc
2026-06-21 20:13:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/chim_ENS10_2020021510_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/chim_ENS10_2020021511_1_out_toDART.nc
2026-06-21 20:13:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/chim_ENS11_2020021510_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/chim_ENS11_2020021511_1_out_toDART.nc
2026-06-21 20:13:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/chim_ENS12_2020021510_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/chim_ENS12_2020021511_1_out_toDART.nc
2026-06-21 20:13:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/chim_ENS13_2020021510_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/chim_ENS13_2020021511_1_out_toDART.nc
2026-06-21 20:13:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/chim_ENS14_2020021510_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/chim_ENS14_2020021511_1_out_toDART.nc
2026-06-21 20:13:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/chim_ENS15_2020021510_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/chim_ENS15_2020021511_1_out_toDART.nc
2026-06-21 20:13:33 INFO Replacement input_template.nml → input.nml completed successfully.
2026-06-21 20:13:33 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-06-21 20:13:33 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-06-21 20:13:33 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-06-21 20:13:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-06-21 20:13:33 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-06-21 20:13:45 INFO Found: []
2026-06-21 20:13:45 INFO No job id returned by command ./run_filter.bsh
2026-06-21 20:13:45 INFO No monitoring will be performed
2026-06-21 20:13:45 INFO Moving DART output files to analysis and preassim directories for date 2020021511 if present ...
2026-06-21 20:13:45 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021511'
2026-06-21 20:13:45 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021511'
2026-06-21 20:13:45 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021511'
2026-06-21 20:13:45 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021511'
2026-06-21 20:13:45 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021511'
2026-06-21 20:13:45 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021511'
2026-06-21 20:13:45 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021511'
2026-06-21 20:13:45 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021511'
2026-06-21 20:13:45 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021511'
2026-06-21 20:13:45 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021511'
2026-06-21 20:13:45 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021511'
2026-06-21 20:13:45 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021511'
2026-06-21 20:13:45 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021511'
2026-06-21 20:13:45 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021511'
2026-06-21 20:13:45 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021511'
2026-06-21 20:13:45 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021511'
2026-06-21 20:13:45 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021511'
2026-06-21 20:13:45 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021511'
2026-06-21 20:13:45 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021511'
2026-06-21 20:13:45 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021511'
2026-06-21 20:13:45 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021511'
2026-06-21 20:13:45 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021511'
2026-06-21 20:13:45 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021511'
2026-06-21 20:13:45 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021511'
2026-06-21 20:13:45 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021511'
2026-06-21 20:13:45 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021511'
2026-06-21 20:13:45 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021511'
2026-06-21 20:13:45 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021511'
2026-06-21 20:13:45 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021511'
2026-06-21 20:13:45 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021511'
2026-06-21 20:13:45 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021511'
2026-06-21 20:13:45 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021511'
2026-06-21 20:13:45 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021511'
2026-06-21 20:13:45 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021511'
2026-06-21 20:13:45 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-06-21 20:13:45 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-06-21 20:13:45 INFO run_dart() is DONE.
2026-06-21 20:13:45 INFO ---------->>> Running update_pollutant_in_end()
2026-06-21 20:13:46 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/end.2020021510_1_ENS1.nc
2026-06-21 20:13:49 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:13:49 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:13:50 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:13:50 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:13:51 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:13:51 INFO [Posterior Diff ENS1] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021511/diff_posterior_ENS1_2020021511.nc
2026-06-21 20:13:51 INFO [Posterior Diff ENS1] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021511/diff_posterior_ENS1_2020021511.relative.nc
2026-06-21 20:13:51 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/end.2020021510_1_ENS2.nc
2026-06-21 20:13:55 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:13:55 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:13:56 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:13:56 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:13:56 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:13:56 INFO [Posterior Diff ENS2] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021511/diff_posterior_ENS2_2020021511.nc
2026-06-21 20:13:56 INFO [Posterior Diff ENS2] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021511/diff_posterior_ENS2_2020021511.relative.nc
2026-06-21 20:13:57 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/end.2020021510_1_ENS3.nc
2026-06-21 20:14:00 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:14:00 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:14:01 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:14:01 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:14:02 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:14:02 INFO [Posterior Diff ENS3] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021511/diff_posterior_ENS3_2020021511.nc
2026-06-21 20:14:02 INFO [Posterior Diff ENS3] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021511/diff_posterior_ENS3_2020021511.relative.nc
2026-06-21 20:14:02 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/end.2020021510_1_ENS4.nc
2026-06-21 20:14:05 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:14:05 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:14:06 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:14:06 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:14:07 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:14:07 INFO [Posterior Diff ENS4] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021511/diff_posterior_ENS4_2020021511.nc
2026-06-21 20:14:07 INFO [Posterior Diff ENS4] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021511/diff_posterior_ENS4_2020021511.relative.nc
2026-06-21 20:14:08 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/end.2020021510_1_ENS5.nc
2026-06-21 20:14:11 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:14:11 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:14:12 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:14:12 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:14:13 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:14:13 INFO [Posterior Diff ENS5] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021511/diff_posterior_ENS5_2020021511.nc
2026-06-21 20:14:13 INFO [Posterior Diff ENS5] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021511/diff_posterior_ENS5_2020021511.relative.nc
2026-06-21 20:14:13 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/end.2020021510_1_ENS6.nc
2026-06-21 20:14:17 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:14:17 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:14:18 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:14:18 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:14:18 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:14:18 INFO [Posterior Diff ENS6] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021511/diff_posterior_ENS6_2020021511.nc
2026-06-21 20:14:18 INFO [Posterior Diff ENS6] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021511/diff_posterior_ENS6_2020021511.relative.nc
2026-06-21 20:14:19 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/end.2020021510_1_ENS7.nc
2026-06-21 20:14:22 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:14:22 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:14:23 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:14:23 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:14:24 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:14:24 INFO [Posterior Diff ENS7] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021511/diff_posterior_ENS7_2020021511.nc
2026-06-21 20:14:24 INFO [Posterior Diff ENS7] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021511/diff_posterior_ENS7_2020021511.relative.nc
2026-06-21 20:14:24 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/end.2020021510_1_ENS8.nc
2026-06-21 20:14:28 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:14:28 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:14:29 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:14:29 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:14:29 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:14:29 INFO [Posterior Diff ENS8] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021511/diff_posterior_ENS8_2020021511.nc
2026-06-21 20:14:29 INFO [Posterior Diff ENS8] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021511/diff_posterior_ENS8_2020021511.relative.nc
2026-06-21 20:14:30 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/end.2020021510_1_ENS9.nc
2026-06-21 20:14:33 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:14:33 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:14:34 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:14:34 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:14:35 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:14:35 INFO [Posterior Diff ENS9] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021511/diff_posterior_ENS9_2020021511.nc
2026-06-21 20:14:35 INFO [Posterior Diff ENS9] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021511/diff_posterior_ENS9_2020021511.relative.nc
2026-06-21 20:14:35 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/end.2020021510_1_ENS10.nc
2026-06-21 20:14:39 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:14:39 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:14:40 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:14:40 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:14:41 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:14:41 INFO [Posterior Diff ENS10] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021511/diff_posterior_ENS10_2020021511.nc
2026-06-21 20:14:41 INFO [Posterior Diff ENS10] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021511/diff_posterior_ENS10_2020021511.relative.nc
2026-06-21 20:14:41 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/end.2020021510_1_ENS11.nc
2026-06-21 20:14:45 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:14:45 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:14:46 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:14:46 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:14:46 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:14:47 INFO [Posterior Diff ENS11] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021511/diff_posterior_ENS11_2020021511.nc
2026-06-21 20:14:47 INFO [Posterior Diff ENS11] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021511/diff_posterior_ENS11_2020021511.relative.nc
2026-06-21 20:14:47 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/end.2020021510_1_ENS12.nc
2026-06-21 20:14:51 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:14:51 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:14:52 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:14:52 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:14:52 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:14:52 INFO [Posterior Diff ENS12] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021511/diff_posterior_ENS12_2020021511.nc
2026-06-21 20:14:52 INFO [Posterior Diff ENS12] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021511/diff_posterior_ENS12_2020021511.relative.nc
2026-06-21 20:14:53 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/end.2020021510_1_ENS13.nc
2026-06-21 20:14:56 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:14:56 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:14:57 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:14:57 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:14:58 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:14:58 INFO [Posterior Diff ENS13] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021511/diff_posterior_ENS13_2020021511.nc
2026-06-21 20:14:58 INFO [Posterior Diff ENS13] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021511/diff_posterior_ENS13_2020021511.relative.nc
2026-06-21 20:14:58 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/end.2020021510_1_ENS14.nc
2026-06-21 20:15:02 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:15:02 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:15:03 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:15:03 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:15:03 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:15:04 INFO [Posterior Diff ENS14] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021511/diff_posterior_ENS14_2020021511.nc
2026-06-21 20:15:04 INFO [Posterior Diff ENS14] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021511/diff_posterior_ENS14_2020021511.relative.nc
2026-06-21 20:15:04 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/end.2020021510_1_ENS15.nc
2026-06-21 20:15:08 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:15:08 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:15:09 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:15:09 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:15:09 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:15:09 INFO [Posterior Diff ENS15] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021511/diff_posterior_ENS15_2020021511.nc
2026-06-21 20:15:09 INFO [Posterior Diff ENS15] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021511/diff_posterior_ENS15_2020021511.relative.nc
2026-06-21 20:15:09 INFO Next run starts from 2020-02-15 11:00:00
2026-06-21 20:15:09 INFO Cycle is DONE; starting a new loop!
2026-06-21 20:15:09 INFO [TIME] step_end current_time=2020-02-15 11:00:00 simulated_time=2020-02-15 11:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 20:15:09 INFO [TIME] step_start current_time=2020-02-15 11:00:00 simulated_time=2020-02-15 11:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 20:15:10 INFO [TIME] window start=2020-02-15 11:00:00 end=2020-02-15 13:00:00 run_hours=2 has_assimilation=True
2026-06-21 20:15:10 INFO Creating directories and links for ENS1 to run chimere's parallel part
2026-06-21 20:15:10 INFO Linking EMIS ...
2026-06-21 20:15:10 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens01.nc
2026-06-21 20:15:10 INFO >> Checking links...
2026-06-21 20:15:13 INFO >> All links are good for ENS1  ...
2026-06-21 20:15:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:15:14 INFO Hourly dataset computed and listing created
2026-06-21 20:15:19 INFO Hourly dataset computed
2026-06-21 20:15:19 INFO Creating directories and links for ENS2 to run chimere's parallel part
2026-06-21 20:15:19 INFO Linking EMIS ...
2026-06-21 20:15:20 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens02.nc
2026-06-21 20:15:20 INFO >> Checking links...
2026-06-21 20:15:22 INFO >> All links are good for ENS2  ...
2026-06-21 20:15:22 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:15:23 INFO Hourly dataset computed and listing created
2026-06-21 20:15:24 INFO Hourly dataset computed
2026-06-21 20:15:24 INFO Creating directories and links for ENS3 to run chimere's parallel part
2026-06-21 20:15:24 INFO Linking EMIS ...
2026-06-21 20:15:25 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens03.nc
2026-06-21 20:15:25 INFO >> Checking links...
2026-06-21 20:15:27 INFO >> All links are good for ENS3  ...
2026-06-21 20:15:27 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:15:28 INFO Hourly dataset computed and listing created
2026-06-21 20:15:29 INFO Hourly dataset computed
2026-06-21 20:15:29 INFO Creating directories and links for ENS4 to run chimere's parallel part
2026-06-21 20:15:29 INFO Linking EMIS ...
2026-06-21 20:15:29 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens04.nc
2026-06-21 20:15:29 INFO >> Checking links...
2026-06-21 20:15:32 INFO >> All links are good for ENS4  ...
2026-06-21 20:15:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:15:33 INFO Hourly dataset computed and listing created
2026-06-21 20:15:33 INFO Hourly dataset computed
2026-06-21 20:15:33 INFO Creating directories and links for ENS5 to run chimere's parallel part
2026-06-21 20:15:33 INFO Linking EMIS ...
2026-06-21 20:15:34 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens05.nc
2026-06-21 20:15:34 INFO >> Checking links...
2026-06-21 20:15:37 INFO >> All links are good for ENS5  ...
2026-06-21 20:15:37 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:15:38 INFO Hourly dataset computed and listing created
2026-06-21 20:15:38 INFO Hourly dataset computed
2026-06-21 20:15:38 INFO Creating directories and links for ENS6 to run chimere's parallel part
2026-06-21 20:15:38 INFO Linking EMIS ...
2026-06-21 20:15:39 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens06.nc
2026-06-21 20:15:39 INFO >> Checking links...
2026-06-21 20:15:41 INFO >> All links are good for ENS6  ...
2026-06-21 20:15:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:15:42 INFO Hourly dataset computed and listing created
2026-06-21 20:15:43 INFO Hourly dataset computed
2026-06-21 20:15:43 INFO Creating directories and links for ENS7 to run chimere's parallel part
2026-06-21 20:15:43 INFO Linking EMIS ...
2026-06-21 20:15:44 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens07.nc
2026-06-21 20:15:44 INFO >> Checking links...
2026-06-21 20:15:46 INFO >> All links are good for ENS7  ...
2026-06-21 20:15:46 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:15:47 INFO Hourly dataset computed and listing created
2026-06-21 20:15:48 INFO Hourly dataset computed
2026-06-21 20:15:48 INFO Creating directories and links for ENS8 to run chimere's parallel part
2026-06-21 20:15:48 INFO Linking EMIS ...
2026-06-21 20:15:48 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens08.nc
2026-06-21 20:15:48 INFO >> Checking links...
2026-06-21 20:15:51 INFO >> All links are good for ENS8  ...
2026-06-21 20:15:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:15:52 INFO Hourly dataset computed and listing created
2026-06-21 20:15:53 INFO Hourly dataset computed
2026-06-21 20:15:53 INFO Creating directories and links for ENS9 to run chimere's parallel part
2026-06-21 20:15:53 INFO Linking EMIS ...
2026-06-21 20:15:53 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens09.nc
2026-06-21 20:15:53 INFO >> Checking links...
2026-06-21 20:15:56 INFO >> All links are good for ENS9  ...
2026-06-21 20:15:56 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:15:57 INFO Hourly dataset computed and listing created
2026-06-21 20:15:58 INFO Hourly dataset computed
2026-06-21 20:15:58 INFO Creating directories and links for ENS10 to run chimere's parallel part
2026-06-21 20:15:58 INFO Linking EMIS ...
2026-06-21 20:15:58 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens10.nc
2026-06-21 20:15:58 INFO >> Checking links...
2026-06-21 20:16:00 INFO >> All links are good for ENS10  ...
2026-06-21 20:16:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:16:01 INFO Hourly dataset computed and listing created
2026-06-21 20:16:02 INFO Hourly dataset computed
2026-06-21 20:16:02 INFO Creating directories and links for ENS11 to run chimere's parallel part
2026-06-21 20:16:02 INFO Linking EMIS ...
2026-06-21 20:16:03 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens11.nc
2026-06-21 20:16:03 INFO >> Checking links...
2026-06-21 20:16:05 INFO >> All links are good for ENS11  ...
2026-06-21 20:16:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:16:06 INFO Hourly dataset computed and listing created
2026-06-21 20:16:07 INFO Hourly dataset computed
2026-06-21 20:16:07 INFO Creating directories and links for ENS12 to run chimere's parallel part
2026-06-21 20:16:07 INFO Linking EMIS ...
2026-06-21 20:16:07 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens12.nc
2026-06-21 20:16:07 INFO >> Checking links...
2026-06-21 20:16:10 INFO >> All links are good for ENS12  ...
2026-06-21 20:16:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:16:11 INFO Hourly dataset computed and listing created
2026-06-21 20:16:17 INFO Hourly dataset computed
2026-06-21 20:16:17 INFO Creating directories and links for ENS13 to run chimere's parallel part
2026-06-21 20:16:17 INFO Linking EMIS ...
2026-06-21 20:16:17 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens13.nc
2026-06-21 20:16:17 INFO >> Checking links...
2026-06-21 20:16:20 INFO >> All links are good for ENS13  ...
2026-06-21 20:16:20 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:16:21 INFO Hourly dataset computed and listing created
2026-06-21 20:16:22 INFO Hourly dataset computed
2026-06-21 20:16:22 INFO Creating directories and links for ENS14 to run chimere's parallel part
2026-06-21 20:16:22 INFO Linking EMIS ...
2026-06-21 20:16:22 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens14.nc
2026-06-21 20:16:22 INFO >> Checking links...
2026-06-21 20:16:25 INFO >> All links are good for ENS14  ...
2026-06-21 20:16:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:16:26 INFO Hourly dataset computed and listing created
2026-06-21 20:16:30 INFO Hourly dataset computed
2026-06-21 20:16:30 INFO Creating directories and links for ENS15 to run chimere's parallel part
2026-06-21 20:16:30 INFO Linking EMIS ...
2026-06-21 20:16:30 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens15.nc
2026-06-21 20:16:30 INFO >> Checking links...
2026-06-21 20:16:33 INFO >> All links are good for ENS15  ...
2026-06-21 20:16:33 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:16:34 INFO Hourly dataset computed and listing created
2026-06-21 20:16:36 INFO Hourly dataset computed
2026-06-21 20:16:36 INFO ---------->>> Running CHIMERE model from 2020-02-15 11:00:00 to 2020-02-15 13:00:00
2026-06-21 20:16:36 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:16:36 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1
2026-06-21 20:16:36 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/end.2020021510_1_ENS1.nc
2026-06-21 20:16:36 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-06-21 20:16:36 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:16:36 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-06-21 20:16:37 INFO Queuing job for member 1...
2026-06-21 20:16:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:16:37 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-06-21 20:16:37 INFO Found: ['4938086']
2026-06-21 20:16:42 INFO [TGCC-IRENE] Submitted job with ID:['4938086']
2026-06-21 20:16:42 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:16:42 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2
2026-06-21 20:16:42 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/end.2020021510_1_ENS2.nc
2026-06-21 20:16:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-06-21 20:16:42 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:16:42 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-06-21 20:16:42 INFO Queuing job for member 2...
2026-06-21 20:16:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:16:42 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-06-21 20:16:43 INFO Found: ['4938088']
2026-06-21 20:16:48 INFO [TGCC-IRENE] Submitted job with ID:['4938088']
2026-06-21 20:16:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:16:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3
2026-06-21 20:16:48 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/end.2020021510_1_ENS3.nc
2026-06-21 20:16:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-06-21 20:16:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:16:48 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-06-21 20:16:48 INFO Queuing job for member 3...
2026-06-21 20:16:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:16:48 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-06-21 20:16:49 INFO Found: ['4938089']
2026-06-21 20:16:54 INFO [TGCC-IRENE] Submitted job with ID:['4938089']
2026-06-21 20:16:54 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:16:54 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4
2026-06-21 20:16:54 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/end.2020021510_1_ENS4.nc
2026-06-21 20:16:54 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-06-21 20:16:54 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:16:54 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-06-21 20:16:54 INFO Queuing job for member 4...
2026-06-21 20:16:54 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:16:54 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-06-21 20:16:55 INFO Found: ['4938091']
2026-06-21 20:17:00 INFO [TGCC-IRENE] Submitted job with ID:['4938091']
2026-06-21 20:17:00 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:17:00 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5
2026-06-21 20:17:00 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/end.2020021510_1_ENS5.nc
2026-06-21 20:17:00 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-06-21 20:17:00 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:17:00 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-06-21 20:17:00 INFO Queuing job for member 5...
2026-06-21 20:17:00 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:17:00 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-06-21 20:17:00 INFO Found: ['4938092']
2026-06-21 20:17:05 INFO [TGCC-IRENE] Submitted job with ID:['4938092']
2026-06-21 20:17:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:17:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6
2026-06-21 20:17:05 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/end.2020021510_1_ENS6.nc
2026-06-21 20:17:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-06-21 20:17:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:17:05 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-06-21 20:17:05 INFO Queuing job for member 6...
2026-06-21 20:17:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:17:05 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-06-21 20:17:06 INFO Found: ['4938095']
2026-06-21 20:17:11 INFO [TGCC-IRENE] Submitted job with ID:['4938095']
2026-06-21 20:17:11 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:17:11 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7
2026-06-21 20:17:11 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/end.2020021510_1_ENS7.nc
2026-06-21 20:17:11 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-06-21 20:17:11 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:17:11 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-06-21 20:17:11 INFO Queuing job for member 7...
2026-06-21 20:17:11 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:17:11 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-06-21 20:17:12 INFO Found: ['4938096']
2026-06-21 20:17:17 INFO [TGCC-IRENE] Submitted job with ID:['4938096']
2026-06-21 20:17:17 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:17:17 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8
2026-06-21 20:17:17 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/end.2020021510_1_ENS8.nc
2026-06-21 20:17:17 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-06-21 20:17:17 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:17:17 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-06-21 20:17:17 INFO Queuing job for member 8...
2026-06-21 20:17:17 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:17:17 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-06-21 20:17:18 INFO Found: ['4938097']
2026-06-21 20:17:23 INFO [TGCC-IRENE] Submitted job with ID:['4938097']
2026-06-21 20:17:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:17:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9
2026-06-21 20:17:23 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/end.2020021510_1_ENS9.nc
2026-06-21 20:17:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-06-21 20:17:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:17:23 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-06-21 20:17:23 INFO Queuing job for member 9...
2026-06-21 20:17:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:17:23 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-06-21 20:17:24 INFO Found: ['4938099']
2026-06-21 20:17:29 INFO [TGCC-IRENE] Submitted job with ID:['4938099']
2026-06-21 20:17:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:17:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10
2026-06-21 20:17:29 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/end.2020021510_1_ENS10.nc
2026-06-21 20:17:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-06-21 20:17:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:17:29 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-06-21 20:17:29 INFO Queuing job for member 10...
2026-06-21 20:17:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:17:29 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-06-21 20:17:29 INFO Found: ['4938100']
2026-06-21 20:17:34 INFO [TGCC-IRENE] Submitted job with ID:['4938100']
2026-06-21 20:17:34 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:17:34 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11
2026-06-21 20:17:34 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/end.2020021510_1_ENS11.nc
2026-06-21 20:17:34 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-06-21 20:17:34 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:17:34 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-06-21 20:17:34 INFO Queuing job for member 11...
2026-06-21 20:17:34 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:17:34 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-06-21 20:17:35 INFO Found: ['4938101']
2026-06-21 20:17:40 INFO [TGCC-IRENE] Submitted job with ID:['4938101']
2026-06-21 20:17:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:17:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12
2026-06-21 20:17:40 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/end.2020021510_1_ENS12.nc
2026-06-21 20:17:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-06-21 20:17:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:17:40 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-06-21 20:17:40 INFO Queuing job for member 12...
2026-06-21 20:17:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:17:40 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-06-21 20:17:41 INFO Found: ['4938103']
2026-06-21 20:17:46 INFO [TGCC-IRENE] Submitted job with ID:['4938103']
2026-06-21 20:17:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:17:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13
2026-06-21 20:17:46 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/end.2020021510_1_ENS13.nc
2026-06-21 20:17:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-06-21 20:17:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:17:46 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-06-21 20:17:46 INFO Queuing job for member 13...
2026-06-21 20:17:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:17:46 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-06-21 20:17:47 INFO Found: ['4938105']
2026-06-21 20:17:52 INFO [TGCC-IRENE] Submitted job with ID:['4938105']
2026-06-21 20:17:52 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:17:52 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14
2026-06-21 20:17:52 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/end.2020021510_1_ENS14.nc
2026-06-21 20:17:52 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-06-21 20:17:52 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:17:52 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-06-21 20:17:52 INFO Queuing job for member 14...
2026-06-21 20:17:52 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:17:52 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-06-21 20:17:53 INFO Found: ['4938108']
2026-06-21 20:17:58 INFO [TGCC-IRENE] Submitted job with ID:['4938108']
2026-06-21 20:17:58 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:17:58 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15
2026-06-21 20:17:58 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/end.2020021510_1_ENS15.nc
2026-06-21 20:17:58 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-06-21 20:17:58 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:17:58 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-06-21 20:17:58 INFO Queuing job for member 15...
2026-06-21 20:17:58 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:17:58 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-06-21 20:17:59 INFO Found: ['4938110']
2026-06-21 20:18:04 INFO [TGCC-IRENE] Submitted job with ID:['4938110']
2026-06-21 20:18:04 INFO Checking job status ...
2026-06-21 20:18:04 INFO None 4938086: status RUNNING/PENDING
2026-06-21 20:18:04 INFO None 4938088: status RUNNING/PENDING
2026-06-21 20:18:04 INFO None 4938089: status RUNNING/PENDING
2026-06-21 20:18:04 INFO None 4938091: status RUNNING/PENDING
2026-06-21 20:18:04 INFO None 4938092: status RUNNING/PENDING
2026-06-21 20:18:04 INFO None 4938095: status RUNNING/PENDING
2026-06-21 20:18:04 INFO None 4938096: status RUNNING/PENDING
2026-06-21 20:18:04 INFO None 4938097: status RUNNING/PENDING
2026-06-21 20:18:04 INFO None 4938099: status RUNNING/PENDING
2026-06-21 20:18:04 INFO None 4938100: status RUNNING/PENDING
2026-06-21 20:18:04 INFO None 4938101: status RUNNING/PENDING
2026-06-21 20:18:04 INFO None 4938103: status RUNNING/PENDING
2026-06-21 20:18:04 INFO None 4938105: status RUNNING/PENDING
2026-06-21 20:18:04 INFO None 4938108: status RUNNING/PENDING
2026-06-21 20:18:04 INFO None 4938110: status RUNNING/PENDING
2026-06-21 20:18:04 INFO Jobs still running: ['4938086', '4938088', '4938089', '4938091', '4938092', '4938095', '4938096', '4938097', '4938099', '4938100', '4938101', '4938103', '4938105', '4938108', '4938110']. Waiting...
2026-06-21 20:18:19 INFO None 4938086: status RUNNING/PENDING
2026-06-21 20:18:19 INFO None 4938088: status RUNNING/PENDING
2026-06-21 20:18:19 INFO None 4938089: status RUNNING/PENDING
2026-06-21 20:18:19 INFO None 4938091: status RUNNING/PENDING
2026-06-21 20:18:19 INFO None 4938092: status RUNNING/PENDING
2026-06-21 20:18:19 INFO None 4938095: status RUNNING/PENDING
2026-06-21 20:18:19 INFO None 4938096: status RUNNING/PENDING
2026-06-21 20:18:19 INFO None 4938097: status RUNNING/PENDING
2026-06-21 20:18:19 INFO None 4938099: status RUNNING/PENDING
2026-06-21 20:18:19 INFO None 4938100: status RUNNING/PENDING
2026-06-21 20:18:19 INFO None 4938101: status RUNNING/PENDING
2026-06-21 20:18:19 INFO None 4938103: status RUNNING/PENDING
2026-06-21 20:18:19 INFO None 4938105: status RUNNING/PENDING
2026-06-21 20:18:19 INFO None 4938108: status RUNNING/PENDING
2026-06-21 20:18:19 INFO None 4938110: status RUNNING/PENDING
2026-06-21 20:18:19 INFO Jobs still running: ['4938086', '4938088', '4938089', '4938091', '4938092', '4938095', '4938096', '4938097', '4938099', '4938100', '4938101', '4938103', '4938105', '4938108', '4938110']. Waiting...
2026-06-21 20:18:34 INFO None 4938086: status RUNNING/PENDING
2026-06-21 20:18:34 INFO None 4938088: status RUNNING/PENDING
2026-06-21 20:18:34 INFO None 4938089: status RUNNING/PENDING
2026-06-21 20:18:34 INFO None 4938091: status RUNNING/PENDING
2026-06-21 20:18:34 INFO None 4938092: status RUNNING/PENDING
2026-06-21 20:18:34 INFO None 4938095: status RUNNING/PENDING
2026-06-21 20:18:34 INFO None 4938096: status RUNNING/PENDING
2026-06-21 20:18:34 INFO None 4938097: status RUNNING/PENDING
2026-06-21 20:18:34 INFO None 4938099: status RUNNING/PENDING
2026-06-21 20:18:35 INFO None 4938100: status RUNNING/PENDING
2026-06-21 20:18:35 INFO None 4938101: status RUNNING/PENDING
2026-06-21 20:18:35 INFO None 4938103: status RUNNING/PENDING
2026-06-21 20:18:35 INFO None 4938105: status RUNNING/PENDING
2026-06-21 20:18:35 INFO None 4938108: status RUNNING/PENDING
2026-06-21 20:18:35 INFO None 4938110: status RUNNING/PENDING
2026-06-21 20:18:35 INFO Jobs still running: ['4938086', '4938088', '4938089', '4938091', '4938092', '4938095', '4938096', '4938097', '4938099', '4938100', '4938101', '4938103', '4938105', '4938108', '4938110']. Waiting...
2026-06-21 20:18:50 INFO None 4938086: status RUNNING/PENDING
2026-06-21 20:18:50 INFO None 4938088: status RUNNING/PENDING
2026-06-21 20:18:50 INFO None 4938089: status RUNNING/PENDING
2026-06-21 20:18:50 INFO None 4938091: status RUNNING/PENDING
2026-06-21 20:18:50 INFO None 4938092: status RUNNING/PENDING
2026-06-21 20:18:50 INFO None 4938095: status RUNNING/PENDING
2026-06-21 20:18:50 INFO None 4938096: status RUNNING/PENDING
2026-06-21 20:18:50 INFO None 4938097: status RUNNING/PENDING
2026-06-21 20:18:50 INFO None 4938099: status RUNNING/PENDING
2026-06-21 20:18:50 INFO None 4938100: status RUNNING/PENDING
2026-06-21 20:18:50 INFO None 4938101: status RUNNING/PENDING
2026-06-21 20:18:50 INFO None 4938103: status RUNNING/PENDING
2026-06-21 20:18:50 INFO None 4938105: status RUNNING/PENDING
2026-06-21 20:18:50 INFO None 4938108: status RUNNING/PENDING
2026-06-21 20:18:50 INFO None 4938110: status RUNNING/PENDING
2026-06-21 20:18:50 INFO Jobs still running: ['4938086', '4938088', '4938089', '4938091', '4938092', '4938095', '4938096', '4938097', '4938099', '4938100', '4938101', '4938103', '4938105', '4938108', '4938110']. Waiting...
2026-06-21 20:19:05 INFO None 4938086: status RUNNING/PENDING
2026-06-21 20:19:05 INFO None 4938088: status RUNNING/PENDING
2026-06-21 20:19:05 INFO None 4938089: status RUNNING/PENDING
2026-06-21 20:19:05 INFO None 4938091: status RUNNING/PENDING
2026-06-21 20:19:05 INFO None 4938092: status RUNNING/PENDING
2026-06-21 20:19:05 INFO None 4938095: status RUNNING/PENDING
2026-06-21 20:19:05 INFO None 4938096: status RUNNING/PENDING
2026-06-21 20:19:05 INFO None 4938097: status RUNNING/PENDING
2026-06-21 20:19:05 INFO None 4938099: status RUNNING/PENDING
2026-06-21 20:19:05 INFO None 4938100: status RUNNING/PENDING
2026-06-21 20:19:05 INFO None 4938101: status RUNNING/PENDING
2026-06-21 20:19:05 INFO None 4938103: status RUNNING/PENDING
2026-06-21 20:19:05 INFO None 4938105: status RUNNING/PENDING
2026-06-21 20:19:05 INFO None 4938108: status RUNNING/PENDING
2026-06-21 20:19:05 INFO None 4938110: status RUNNING/PENDING
2026-06-21 20:19:05 INFO Jobs still running: ['4938086', '4938088', '4938089', '4938091', '4938092', '4938095', '4938096', '4938097', '4938099', '4938100', '4938101', '4938103', '4938105', '4938108', '4938110']. Waiting...
2026-06-21 20:19:20 INFO None 4938086: status RUNNING/PENDING
2026-06-21 20:19:20 INFO None 4938088: status RUNNING/PENDING
2026-06-21 20:19:20 INFO None 4938089: status RUNNING/PENDING
2026-06-21 20:19:20 INFO None 4938091: status RUNNING/PENDING
2026-06-21 20:19:20 INFO None 4938092: status RUNNING/PENDING
2026-06-21 20:19:20 INFO None 4938095: status RUNNING/PENDING
2026-06-21 20:19:20 INFO None 4938096: status RUNNING/PENDING
2026-06-21 20:19:20 INFO None 4938097: status RUNNING/PENDING
2026-06-21 20:19:20 INFO None 4938099: status RUNNING/PENDING
2026-06-21 20:19:20 INFO None 4938100: status RUNNING/PENDING
2026-06-21 20:19:20 INFO None 4938101: status RUNNING/PENDING
2026-06-21 20:19:20 INFO None 4938103: status RUNNING/PENDING
2026-06-21 20:19:21 INFO None 4938105: status RUNNING/PENDING
2026-06-21 20:19:21 INFO None 4938108: status RUNNING/PENDING
2026-06-21 20:19:21 INFO None 4938110: status RUNNING/PENDING
2026-06-21 20:19:21 INFO Jobs still running: ['4938086', '4938088', '4938089', '4938091', '4938092', '4938095', '4938096', '4938097', '4938099', '4938100', '4938101', '4938103', '4938105', '4938108', '4938110']. Waiting...
2026-06-21 20:19:36 INFO None 4938086: status RUNNING/PENDING
2026-06-21 20:19:36 INFO None 4938088: status RUNNING/PENDING
2026-06-21 20:19:36 INFO None 4938089: status RUNNING/PENDING
2026-06-21 20:19:36 INFO None 4938091: status RUNNING/PENDING
2026-06-21 20:19:36 INFO None 4938092: status RUNNING/PENDING
2026-06-21 20:19:36 INFO None 4938095: status RUNNING/PENDING
2026-06-21 20:19:36 INFO None 4938096: status RUNNING/PENDING
2026-06-21 20:19:36 INFO None 4938097: status RUNNING/PENDING
2026-06-21 20:19:36 INFO None 4938099: status RUNNING/PENDING
2026-06-21 20:19:36 INFO None 4938100: status RUNNING/PENDING
2026-06-21 20:19:36 INFO None 4938101: status RUNNING/PENDING
2026-06-21 20:19:36 INFO None 4938103: status RUNNING/PENDING
2026-06-21 20:19:36 INFO None 4938105: status RUNNING/PENDING
2026-06-21 20:19:36 INFO None 4938108: status RUNNING/PENDING
2026-06-21 20:19:36 INFO None 4938110: status RUNNING/PENDING
2026-06-21 20:19:36 INFO Jobs still running: ['4938086', '4938088', '4938089', '4938091', '4938092', '4938095', '4938096', '4938097', '4938099', '4938100', '4938101', '4938103', '4938105', '4938108', '4938110']. Waiting...
2026-06-21 20:19:51 INFO None 4938086: status FINISHED
2026-06-21 20:19:51 INFO None 4938088: status RUNNING/PENDING
2026-06-21 20:19:51 INFO None 4938089: status RUNNING/PENDING
2026-06-21 20:19:51 INFO None 4938091: status RUNNING/PENDING
2026-06-21 20:19:51 INFO None 4938092: status RUNNING/PENDING
2026-06-21 20:19:51 INFO None 4938095: status RUNNING/PENDING
2026-06-21 20:19:51 INFO None 4938096: status RUNNING/PENDING
2026-06-21 20:19:51 INFO None 4938097: status RUNNING/PENDING
2026-06-21 20:19:51 INFO None 4938099: status RUNNING/PENDING
2026-06-21 20:19:51 INFO None 4938100: status RUNNING/PENDING
2026-06-21 20:19:51 INFO None 4938101: status RUNNING/PENDING
2026-06-21 20:19:51 INFO None 4938103: status RUNNING/PENDING
2026-06-21 20:19:51 INFO None 4938105: status RUNNING/PENDING
2026-06-21 20:19:51 INFO None 4938108: status RUNNING/PENDING
2026-06-21 20:19:51 INFO None 4938110: status RUNNING/PENDING
2026-06-21 20:19:51 INFO Jobs still running: ['4938088', '4938089', '4938091', '4938092', '4938095', '4938096', '4938097', '4938099', '4938100', '4938101', '4938103', '4938105', '4938108', '4938110']. Waiting...
2026-06-21 20:20:06 INFO None 4938086: status FINISHED
2026-06-21 20:20:06 INFO None 4938088: status FINISHED
2026-06-21 20:20:06 INFO None 4938089: status FINISHED
2026-06-21 20:20:06 INFO None 4938091: status FINISHED
2026-06-21 20:20:06 INFO None 4938092: status FINISHED
2026-06-21 20:20:06 INFO None 4938095: status RUNNING/PENDING
2026-06-21 20:20:06 INFO None 4938096: status RUNNING/PENDING
2026-06-21 20:20:06 INFO None 4938097: status RUNNING/PENDING
2026-06-21 20:20:06 INFO None 4938099: status FINISHED
2026-06-21 20:20:06 INFO None 4938100: status RUNNING/PENDING
2026-06-21 20:20:06 INFO None 4938101: status RUNNING/PENDING
2026-06-21 20:20:06 INFO None 4938103: status RUNNING/PENDING
2026-06-21 20:20:06 INFO None 4938105: status RUNNING/PENDING
2026-06-21 20:20:06 INFO None 4938108: status RUNNING/PENDING
2026-06-21 20:20:07 INFO None 4938110: status RUNNING/PENDING
2026-06-21 20:20:07 INFO Jobs still running: ['4938095', '4938096', '4938097', '4938100', '4938101', '4938103', '4938105', '4938108', '4938110']. Waiting...
2026-06-21 20:20:22 INFO None 4938086: status FINISHED
2026-06-21 20:20:22 INFO None 4938088: status FINISHED
2026-06-21 20:20:22 INFO None 4938089: status FINISHED
2026-06-21 20:20:22 INFO None 4938091: status FINISHED
2026-06-21 20:20:22 INFO None 4938092: status FINISHED
2026-06-21 20:20:22 INFO None 4938095: status FINISHED
2026-06-21 20:20:22 INFO None 4938096: status RUNNING/PENDING
2026-06-21 20:20:22 INFO None 4938097: status RUNNING/PENDING
2026-06-21 20:20:22 INFO None 4938099: status FINISHED
2026-06-21 20:20:22 INFO None 4938100: status RUNNING/PENDING
2026-06-21 20:20:22 INFO None 4938101: status RUNNING/PENDING
2026-06-21 20:20:22 INFO None 4938103: status RUNNING/PENDING
2026-06-21 20:20:22 INFO None 4938105: status RUNNING/PENDING
2026-06-21 20:20:22 INFO None 4938108: status RUNNING/PENDING
2026-06-21 20:20:22 INFO None 4938110: status RUNNING/PENDING
2026-06-21 20:20:22 INFO Jobs still running: ['4938096', '4938097', '4938100', '4938101', '4938103', '4938105', '4938108', '4938110']. Waiting...
2026-06-21 20:20:37 INFO None 4938086: status FINISHED
2026-06-21 20:20:37 INFO None 4938088: status FINISHED
2026-06-21 20:20:37 INFO None 4938089: status FINISHED
2026-06-21 20:20:37 INFO None 4938091: status FINISHED
2026-06-21 20:20:37 INFO None 4938092: status FINISHED
2026-06-21 20:20:37 INFO None 4938095: status FINISHED
2026-06-21 20:20:37 INFO None 4938096: status RUNNING/PENDING
2026-06-21 20:20:37 INFO None 4938097: status RUNNING/PENDING
2026-06-21 20:20:37 INFO None 4938099: status FINISHED
2026-06-21 20:20:37 INFO None 4938100: status RUNNING/PENDING
2026-06-21 20:20:37 INFO None 4938101: status RUNNING/PENDING
2026-06-21 20:20:37 INFO None 4938103: status RUNNING/PENDING
2026-06-21 20:20:37 INFO None 4938105: status RUNNING/PENDING
2026-06-21 20:20:37 INFO None 4938108: status RUNNING/PENDING
2026-06-21 20:20:37 INFO None 4938110: status RUNNING/PENDING
2026-06-21 20:20:37 INFO Jobs still running: ['4938096', '4938097', '4938100', '4938101', '4938103', '4938105', '4938108', '4938110']. Waiting...
2026-06-21 20:20:52 INFO None 4938086: status FINISHED
2026-06-21 20:20:52 INFO None 4938088: status FINISHED
2026-06-21 20:20:52 INFO None 4938089: status FINISHED
2026-06-21 20:20:52 INFO None 4938091: status FINISHED
2026-06-21 20:20:52 INFO None 4938092: status FINISHED
2026-06-21 20:20:52 INFO None 4938095: status FINISHED
2026-06-21 20:20:52 INFO None 4938096: status RUNNING/PENDING
2026-06-21 20:20:52 INFO None 4938097: status FINISHED
2026-06-21 20:20:52 INFO None 4938099: status FINISHED
2026-06-21 20:20:52 INFO None 4938100: status RUNNING/PENDING
2026-06-21 20:20:52 INFO None 4938101: status RUNNING/PENDING
2026-06-21 20:20:52 INFO None 4938103: status RUNNING/PENDING
2026-06-21 20:20:52 INFO None 4938105: status RUNNING/PENDING
2026-06-21 20:20:52 INFO None 4938108: status RUNNING/PENDING
2026-06-21 20:20:52 INFO None 4938110: status RUNNING/PENDING
2026-06-21 20:20:52 INFO Jobs still running: ['4938096', '4938100', '4938101', '4938103', '4938105', '4938108', '4938110']. Waiting...
2026-06-21 20:21:07 INFO None 4938086: status FINISHED
2026-06-21 20:21:08 INFO None 4938088: status FINISHED
2026-06-21 20:21:08 INFO None 4938089: status FINISHED
2026-06-21 20:21:08 INFO None 4938091: status FINISHED
2026-06-21 20:21:08 INFO None 4938092: status FINISHED
2026-06-21 20:21:08 INFO None 4938095: status FINISHED
2026-06-21 20:21:08 INFO None 4938096: status RUNNING/PENDING
2026-06-21 20:21:08 INFO None 4938097: status FINISHED
2026-06-21 20:21:08 INFO None 4938099: status FINISHED
2026-06-21 20:21:08 INFO None 4938100: status RUNNING/PENDING
2026-06-21 20:21:08 INFO None 4938101: status RUNNING/PENDING
2026-06-21 20:21:08 INFO None 4938103: status RUNNING/PENDING
2026-06-21 20:21:08 INFO None 4938105: status RUNNING/PENDING
2026-06-21 20:21:08 INFO None 4938108: status RUNNING/PENDING
2026-06-21 20:21:08 INFO None 4938110: status RUNNING/PENDING
2026-06-21 20:21:08 INFO Jobs still running: ['4938096', '4938100', '4938101', '4938103', '4938105', '4938108', '4938110']. Waiting...
2026-06-21 20:21:23 INFO None 4938086: status FINISHED
2026-06-21 20:21:23 INFO None 4938088: status FINISHED
2026-06-21 20:21:23 INFO None 4938089: status FINISHED
2026-06-21 20:21:23 INFO None 4938091: status FINISHED
2026-06-21 20:21:23 INFO None 4938092: status FINISHED
2026-06-21 20:21:23 INFO None 4938095: status FINISHED
2026-06-21 20:21:23 INFO None 4938096: status FINISHED
2026-06-21 20:21:23 INFO None 4938097: status FINISHED
2026-06-21 20:21:23 INFO None 4938099: status FINISHED
2026-06-21 20:21:23 INFO None 4938100: status RUNNING/PENDING
2026-06-21 20:21:23 INFO None 4938101: status RUNNING/PENDING
2026-06-21 20:21:23 INFO None 4938103: status FINISHED
2026-06-21 20:21:23 INFO None 4938105: status FINISHED
2026-06-21 20:21:23 INFO None 4938108: status FINISHED
2026-06-21 20:21:23 INFO None 4938110: status RUNNING/PENDING
2026-06-21 20:21:23 INFO Jobs still running: ['4938100', '4938101', '4938110']. Waiting...
2026-06-21 20:21:38 INFO None 4938086: status FINISHED
2026-06-21 20:21:38 INFO None 4938088: status FINISHED
2026-06-21 20:21:38 INFO None 4938089: status FINISHED
2026-06-21 20:21:38 INFO None 4938091: status FINISHED
2026-06-21 20:21:38 INFO None 4938092: status FINISHED
2026-06-21 20:21:38 INFO None 4938095: status FINISHED
2026-06-21 20:21:38 INFO None 4938096: status FINISHED
2026-06-21 20:21:38 INFO None 4938097: status FINISHED
2026-06-21 20:21:38 INFO None 4938099: status FINISHED
2026-06-21 20:21:38 INFO None 4938100: status RUNNING/PENDING
2026-06-21 20:21:38 INFO None 4938101: status RUNNING/PENDING
2026-06-21 20:21:38 INFO None 4938103: status FINISHED
2026-06-21 20:21:38 INFO None 4938105: status FINISHED
2026-06-21 20:21:38 INFO None 4938108: status FINISHED
2026-06-21 20:21:38 INFO None 4938110: status RUNNING/PENDING
2026-06-21 20:21:38 INFO Jobs still running: ['4938100', '4938101', '4938110']. Waiting...
2026-06-21 20:21:53 INFO None 4938086: status FINISHED
2026-06-21 20:21:53 INFO None 4938088: status FINISHED
2026-06-21 20:21:53 INFO None 4938089: status FINISHED
2026-06-21 20:21:53 INFO None 4938091: status FINISHED
2026-06-21 20:21:53 INFO None 4938092: status FINISHED
2026-06-21 20:21:53 INFO None 4938095: status FINISHED
2026-06-21 20:21:53 INFO None 4938096: status FINISHED
2026-06-21 20:21:53 INFO None 4938097: status FINISHED
2026-06-21 20:21:53 INFO None 4938099: status FINISHED
2026-06-21 20:21:54 INFO None 4938100: status FINISHED
2026-06-21 20:21:54 INFO None 4938101: status FINISHED
2026-06-21 20:21:54 INFO None 4938103: status FINISHED
2026-06-21 20:21:54 INFO None 4938105: status FINISHED
2026-06-21 20:21:54 INFO None 4938108: status FINISHED
2026-06-21 20:21:54 INFO None 4938110: status RUNNING/PENDING
2026-06-21 20:21:54 INFO Jobs still running: ['4938110']. Waiting...
2026-06-21 20:22:09 INFO None 4938086: status FINISHED
2026-06-21 20:22:09 INFO None 4938088: status FINISHED
2026-06-21 20:22:09 INFO None 4938089: status FINISHED
2026-06-21 20:22:09 INFO None 4938091: status FINISHED
2026-06-21 20:22:09 INFO None 4938092: status FINISHED
2026-06-21 20:22:09 INFO None 4938095: status FINISHED
2026-06-21 20:22:09 INFO None 4938096: status FINISHED
2026-06-21 20:22:09 INFO None 4938097: status FINISHED
2026-06-21 20:22:09 INFO None 4938099: status FINISHED
2026-06-21 20:22:09 INFO None 4938100: status FINISHED
2026-06-21 20:22:09 INFO None 4938101: status FINISHED
2026-06-21 20:22:09 INFO None 4938103: status FINISHED
2026-06-21 20:22:09 INFO None 4938105: status FINISHED
2026-06-21 20:22:09 INFO None 4938108: status FINISHED
2026-06-21 20:22:09 INFO None 4938110: status RUNNING/PENDING
2026-06-21 20:22:09 INFO Jobs still running: ['4938110']. Waiting...
2026-06-21 20:22:24 INFO None 4938086: status FINISHED
2026-06-21 20:22:24 INFO None 4938088: status FINISHED
2026-06-21 20:22:24 INFO None 4938089: status FINISHED
2026-06-21 20:22:24 INFO None 4938091: status FINISHED
2026-06-21 20:22:24 INFO None 4938092: status FINISHED
2026-06-21 20:22:24 INFO None 4938095: status FINISHED
2026-06-21 20:22:24 INFO None 4938096: status FINISHED
2026-06-21 20:22:24 INFO None 4938097: status FINISHED
2026-06-21 20:22:24 INFO None 4938099: status FINISHED
2026-06-21 20:22:24 INFO None 4938100: status FINISHED
2026-06-21 20:22:24 INFO None 4938101: status FINISHED
2026-06-21 20:22:24 INFO None 4938103: status FINISHED
2026-06-21 20:22:24 INFO None 4938105: status FINISHED
2026-06-21 20:22:24 INFO None 4938108: status FINISHED
2026-06-21 20:22:24 INFO None 4938110: status RUNNING/PENDING
2026-06-21 20:22:24 INFO Jobs still running: ['4938110']. Waiting...
2026-06-21 20:22:39 INFO None 4938086: status FINISHED
2026-06-21 20:22:39 INFO None 4938088: status FINISHED
2026-06-21 20:22:39 INFO None 4938089: status FINISHED
2026-06-21 20:22:39 INFO None 4938091: status FINISHED
2026-06-21 20:22:39 INFO None 4938092: status FINISHED
2026-06-21 20:22:39 INFO None 4938095: status FINISHED
2026-06-21 20:22:39 INFO None 4938096: status FINISHED
2026-06-21 20:22:39 INFO None 4938097: status FINISHED
2026-06-21 20:22:39 INFO None 4938099: status FINISHED
2026-06-21 20:22:39 INFO None 4938100: status FINISHED
2026-06-21 20:22:39 INFO None 4938101: status FINISHED
2026-06-21 20:22:39 INFO None 4938103: status FINISHED
2026-06-21 20:22:39 INFO None 4938105: status FINISHED
2026-06-21 20:22:39 INFO None 4938108: status FINISHED
2026-06-21 20:22:39 INFO None 4938110: status FINISHED
2026-06-21 20:22:39 INFO Jobs ['4938086', '4938088', '4938089', '4938091', '4938092', '4938095', '4938096', '4938097', '4938099', '4938100', '4938101', '4938103', '4938105', '4938108', '4938110'] have finished
2026-06-21 20:22:39 INFO Checking restart files were created ...
2026-06-21 20:22:39 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/end.2020021511_2_ENS1.nc(1002685915 bytes)
2026-06-21 20:22:39 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/end.2020021511_2_ENS2.nc(1002685915 bytes)
2026-06-21 20:22:39 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/end.2020021511_2_ENS3.nc(1002685915 bytes)
2026-06-21 20:22:39 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/end.2020021511_2_ENS4.nc(1002685915 bytes)
2026-06-21 20:22:39 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/end.2020021511_2_ENS5.nc(1002685915 bytes)
2026-06-21 20:22:39 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/end.2020021511_2_ENS6.nc(1002685915 bytes)
2026-06-21 20:22:39 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/end.2020021511_2_ENS7.nc(1002685915 bytes)
2026-06-21 20:22:39 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/end.2020021511_2_ENS8.nc(1002685915 bytes)
2026-06-21 20:22:39 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/end.2020021511_2_ENS9.nc(1002685915 bytes)
2026-06-21 20:22:39 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/end.2020021511_2_ENS10.nc(1002685915 bytes)
2026-06-21 20:22:39 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/end.2020021511_2_ENS11.nc(1002685915 bytes)
2026-06-21 20:22:39 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/end.2020021511_2_ENS12.nc(1002685915 bytes)
2026-06-21 20:22:39 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/end.2020021511_2_ENS13.nc(1002685915 bytes)
2026-06-21 20:22:39 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/end.2020021511_2_ENS14.nc(1002685915 bytes)
2026-06-21 20:22:39 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/end.2020021511_2_ENS15.nc(1002685915 bytes)
2026-06-21 20:22:39 INFO  Run_model() completed successfully.
2026-06-21 20:22:39 INFO [TIME] after_model_set_simulated_time current_time=2020-02-15 11:00:00 simulated_time=2020-02-15 13:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 20:22:39 INFO [TIME] gregorian_conversion simulated_time=2020-02-15 13:00:00 days=153081 seconds=46800
2026-06-21 20:22:39 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-06-21 20:22:39 INFO [TIME] increment current_time 2020-02-15 11:00:00 -> 2020-02-15 13:00:00
2026-06-21 20:22:39 INFO [TIME] after_increment_before_assimilation current_time=2020-02-15 13:00:00 simulated_time=2020-02-15 13:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 20:22:39 INFO ---------->>> Running process_satellite_data()
2026-06-21 20:22:40 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12132.nc
2026-06-21 20:22:40 INFO ---------->>> Running run_obs_converter()
2026-06-21 20:22:40 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_47038_153081.out
2026-06-21 20:22:40 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_47038_153081.out
2026-06-21 20:22:40 INFO ---------->>> Running DART
2026-06-21 20:22:40 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-06-21 20:22:40 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/chim_ENS1_2020021511_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/chim_ENS1_2020021513_1_out_toDART.nc
2026-06-21 20:22:40 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/chim_ENS2_2020021511_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/chim_ENS2_2020021513_1_out_toDART.nc
2026-06-21 20:22:40 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/chim_ENS3_2020021511_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/chim_ENS3_2020021513_1_out_toDART.nc
2026-06-21 20:22:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/chim_ENS4_2020021511_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/chim_ENS4_2020021513_1_out_toDART.nc
2026-06-21 20:22:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/chim_ENS5_2020021511_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/chim_ENS5_2020021513_1_out_toDART.nc
2026-06-21 20:22:42 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/chim_ENS6_2020021511_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/chim_ENS6_2020021513_1_out_toDART.nc
2026-06-21 20:22:42 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/chim_ENS7_2020021511_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/chim_ENS7_2020021513_1_out_toDART.nc
2026-06-21 20:22:42 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/chim_ENS8_2020021511_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/chim_ENS8_2020021513_1_out_toDART.nc
2026-06-21 20:22:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/chim_ENS9_2020021511_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/chim_ENS9_2020021513_1_out_toDART.nc
2026-06-21 20:22:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/chim_ENS10_2020021511_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/chim_ENS10_2020021513_1_out_toDART.nc
2026-06-21 20:22:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/chim_ENS11_2020021511_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/chim_ENS11_2020021513_1_out_toDART.nc
2026-06-21 20:22:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/chim_ENS12_2020021511_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/chim_ENS12_2020021513_1_out_toDART.nc
2026-06-21 20:22:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/chim_ENS13_2020021511_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/chim_ENS13_2020021513_1_out_toDART.nc
2026-06-21 20:22:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/chim_ENS14_2020021511_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/chim_ENS14_2020021513_1_out_toDART.nc
2026-06-21 20:22:45 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/chim_ENS15_2020021511_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/chim_ENS15_2020021513_1_out_toDART.nc
2026-06-21 20:22:45 INFO Replacement input_template.nml → input.nml completed successfully.
2026-06-21 20:22:45 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-06-21 20:22:45 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-06-21 20:22:45 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-06-21 20:22:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-06-21 20:22:45 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-06-21 20:22:56 INFO Found: []
2026-06-21 20:22:56 INFO No job id returned by command ./run_filter.bsh
2026-06-21 20:22:56 INFO No monitoring will be performed
2026-06-21 20:22:56 INFO Moving DART output files to analysis and preassim directories for date 2020021513 if present ...
2026-06-21 20:22:56 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021513'
2026-06-21 20:22:56 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021513'
2026-06-21 20:22:56 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021513'
2026-06-21 20:22:56 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021513'
2026-06-21 20:22:56 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021513'
2026-06-21 20:22:56 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021513'
2026-06-21 20:22:56 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021513'
2026-06-21 20:22:56 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021513'
2026-06-21 20:22:56 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021513'
2026-06-21 20:22:56 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021513'
2026-06-21 20:22:56 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021513'
2026-06-21 20:22:56 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021513'
2026-06-21 20:22:56 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021513'
2026-06-21 20:22:56 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021513'
2026-06-21 20:22:56 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021513'
2026-06-21 20:22:56 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021513'
2026-06-21 20:22:56 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021513'
2026-06-21 20:22:56 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021513'
2026-06-21 20:22:56 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021513'
2026-06-21 20:22:56 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021513'
2026-06-21 20:22:56 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021513'
2026-06-21 20:22:56 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021513'
2026-06-21 20:22:56 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021513'
2026-06-21 20:22:56 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021513'
2026-06-21 20:22:56 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021513'
2026-06-21 20:22:56 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021513'
2026-06-21 20:22:56 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021513'
2026-06-21 20:22:56 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021513'
2026-06-21 20:22:56 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021513'
2026-06-21 20:22:56 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021513'
2026-06-21 20:22:56 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021513'
2026-06-21 20:22:56 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021513'
2026-06-21 20:22:56 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021513'
2026-06-21 20:22:56 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021513'
2026-06-21 20:22:56 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-06-21 20:22:56 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-06-21 20:22:56 INFO run_dart() is DONE.
2026-06-21 20:22:56 INFO ---------->>> Running update_pollutant_in_end()
2026-06-21 20:22:56 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/end.2020021511_2_ENS1.nc
2026-06-21 20:23:02 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:23:02 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:23:03 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:23:03 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:23:04 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:23:04 INFO [Posterior Diff ENS1] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021513/diff_posterior_ENS1_2020021513.nc
2026-06-21 20:23:04 INFO [Posterior Diff ENS1] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021513/diff_posterior_ENS1_2020021513.relative.nc
2026-06-21 20:23:04 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/end.2020021511_2_ENS2.nc
2026-06-21 20:23:10 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:23:10 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:23:11 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:23:11 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:23:11 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:23:12 INFO [Posterior Diff ENS2] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021513/diff_posterior_ENS2_2020021513.nc
2026-06-21 20:23:12 INFO [Posterior Diff ENS2] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021513/diff_posterior_ENS2_2020021513.relative.nc
2026-06-21 20:23:12 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/end.2020021511_2_ENS3.nc
2026-06-21 20:23:17 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:23:17 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:23:18 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:23:18 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:23:19 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:23:19 INFO [Posterior Diff ENS3] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021513/diff_posterior_ENS3_2020021513.nc
2026-06-21 20:23:19 INFO [Posterior Diff ENS3] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021513/diff_posterior_ENS3_2020021513.relative.nc
2026-06-21 20:23:19 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/end.2020021511_2_ENS4.nc
2026-06-21 20:23:24 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:23:25 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:23:25 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:23:25 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:23:26 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:23:26 INFO [Posterior Diff ENS4] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021513/diff_posterior_ENS4_2020021513.nc
2026-06-21 20:23:26 INFO [Posterior Diff ENS4] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021513/diff_posterior_ENS4_2020021513.relative.nc
2026-06-21 20:23:27 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/end.2020021511_2_ENS5.nc
2026-06-21 20:23:32 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:23:32 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:23:33 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:23:33 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:23:33 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:23:34 INFO [Posterior Diff ENS5] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021513/diff_posterior_ENS5_2020021513.nc
2026-06-21 20:23:34 INFO [Posterior Diff ENS5] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021513/diff_posterior_ENS5_2020021513.relative.nc
2026-06-21 20:23:34 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/end.2020021511_2_ENS6.nc
2026-06-21 20:23:39 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:23:39 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:23:40 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:23:40 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:23:41 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:23:41 INFO [Posterior Diff ENS6] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021513/diff_posterior_ENS6_2020021513.nc
2026-06-21 20:23:41 INFO [Posterior Diff ENS6] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021513/diff_posterior_ENS6_2020021513.relative.nc
2026-06-21 20:23:41 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/end.2020021511_2_ENS7.nc
2026-06-21 20:23:47 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:23:47 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:23:48 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:23:48 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:23:48 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:23:49 INFO [Posterior Diff ENS7] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021513/diff_posterior_ENS7_2020021513.nc
2026-06-21 20:23:49 INFO [Posterior Diff ENS7] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021513/diff_posterior_ENS7_2020021513.relative.nc
2026-06-21 20:23:49 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/end.2020021511_2_ENS8.nc
2026-06-21 20:23:54 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:23:54 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:23:55 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:23:55 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:23:56 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:23:56 INFO [Posterior Diff ENS8] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021513/diff_posterior_ENS8_2020021513.nc
2026-06-21 20:23:56 INFO [Posterior Diff ENS8] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021513/diff_posterior_ENS8_2020021513.relative.nc
2026-06-21 20:23:56 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/end.2020021511_2_ENS9.nc
2026-06-21 20:24:02 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:24:02 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:24:03 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:24:03 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:24:03 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:24:03 INFO [Posterior Diff ENS9] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021513/diff_posterior_ENS9_2020021513.nc
2026-06-21 20:24:03 INFO [Posterior Diff ENS9] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021513/diff_posterior_ENS9_2020021513.relative.nc
2026-06-21 20:24:04 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/end.2020021511_2_ENS10.nc
2026-06-21 20:24:09 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:24:09 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:24:10 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:24:10 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:24:10 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:24:10 INFO [Posterior Diff ENS10] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021513/diff_posterior_ENS10_2020021513.nc
2026-06-21 20:24:10 INFO [Posterior Diff ENS10] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021513/diff_posterior_ENS10_2020021513.relative.nc
2026-06-21 20:24:11 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/end.2020021511_2_ENS11.nc
2026-06-21 20:24:16 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:24:16 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:24:17 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:24:17 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:24:18 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:24:18 INFO [Posterior Diff ENS11] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021513/diff_posterior_ENS11_2020021513.nc
2026-06-21 20:24:18 INFO [Posterior Diff ENS11] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021513/diff_posterior_ENS11_2020021513.relative.nc
2026-06-21 20:24:18 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/end.2020021511_2_ENS12.nc
2026-06-21 20:24:23 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:24:23 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:24:24 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:24:24 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:24:25 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:24:25 INFO [Posterior Diff ENS12] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021513/diff_posterior_ENS12_2020021513.nc
2026-06-21 20:24:25 INFO [Posterior Diff ENS12] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021513/diff_posterior_ENS12_2020021513.relative.nc
2026-06-21 20:24:25 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/end.2020021511_2_ENS13.nc
2026-06-21 20:24:31 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:24:31 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:24:32 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:24:32 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:24:32 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:24:33 INFO [Posterior Diff ENS13] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021513/diff_posterior_ENS13_2020021513.nc
2026-06-21 20:24:33 INFO [Posterior Diff ENS13] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021513/diff_posterior_ENS13_2020021513.relative.nc
2026-06-21 20:24:33 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/end.2020021511_2_ENS14.nc
2026-06-21 20:24:38 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:24:38 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:24:39 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:24:39 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:24:40 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:24:40 INFO [Posterior Diff ENS14] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021513/diff_posterior_ENS14_2020021513.nc
2026-06-21 20:24:40 INFO [Posterior Diff ENS14] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021513/diff_posterior_ENS14_2020021513.relative.nc
2026-06-21 20:24:41 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/end.2020021511_2_ENS15.nc
2026-06-21 20:24:46 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:24:46 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:24:47 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:24:47 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:24:48 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:24:48 INFO [Posterior Diff ENS15] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021513/diff_posterior_ENS15_2020021513.nc
2026-06-21 20:24:48 INFO [Posterior Diff ENS15] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021513/diff_posterior_ENS15_2020021513.relative.nc
2026-06-21 20:24:48 INFO Next run starts from 2020-02-15 13:00:00
2026-06-21 20:24:48 INFO Cycle is DONE; starting a new loop!
2026-06-21 20:24:48 INFO [TIME] step_end current_time=2020-02-15 13:00:00 simulated_time=2020-02-15 13:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 20:24:48 INFO [TIME] step_start current_time=2020-02-15 13:00:00 simulated_time=2020-02-15 13:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 20:24:48 INFO [TIME] window start=2020-02-15 13:00:00 end=2020-02-15 15:00:00 run_hours=2 has_assimilation=True
2026-06-21 20:24:48 INFO Creating directories and links for ENS1 to run chimere's parallel part
2026-06-21 20:24:48 INFO Linking EMIS ...
2026-06-21 20:24:48 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens01.nc
2026-06-21 20:24:48 INFO >> Checking links...
2026-06-21 20:24:51 INFO >> All links are good for ENS1  ...
2026-06-21 20:24:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:24:52 INFO Hourly dataset computed and listing created
2026-06-21 20:24:56 INFO Hourly dataset computed
2026-06-21 20:24:56 INFO Creating directories and links for ENS2 to run chimere's parallel part
2026-06-21 20:24:56 INFO Linking EMIS ...
2026-06-21 20:24:56 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens02.nc
2026-06-21 20:24:56 INFO >> Checking links...
2026-06-21 20:24:59 INFO >> All links are good for ENS2  ...
2026-06-21 20:24:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:25:00 INFO Hourly dataset computed and listing created
2026-06-21 20:25:01 INFO Hourly dataset computed
2026-06-21 20:25:01 INFO Creating directories and links for ENS3 to run chimere's parallel part
2026-06-21 20:25:01 INFO Linking EMIS ...
2026-06-21 20:25:01 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens03.nc
2026-06-21 20:25:01 INFO >> Checking links...
2026-06-21 20:25:03 INFO >> All links are good for ENS3  ...
2026-06-21 20:25:03 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:25:04 INFO Hourly dataset computed and listing created
2026-06-21 20:25:10 INFO Hourly dataset computed
2026-06-21 20:25:10 INFO Creating directories and links for ENS4 to run chimere's parallel part
2026-06-21 20:25:10 INFO Linking EMIS ...
2026-06-21 20:25:11 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens04.nc
2026-06-21 20:25:11 INFO >> Checking links...
2026-06-21 20:25:14 INFO >> All links are good for ENS4  ...
2026-06-21 20:25:14 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:25:15 INFO Hourly dataset computed and listing created
2026-06-21 20:25:21 INFO Hourly dataset computed
2026-06-21 20:25:22 INFO Creating directories and links for ENS5 to run chimere's parallel part
2026-06-21 20:25:22 INFO Linking EMIS ...
2026-06-21 20:25:22 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens05.nc
2026-06-21 20:25:22 INFO >> Checking links...
2026-06-21 20:25:24 INFO >> All links are good for ENS5  ...
2026-06-21 20:25:24 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:25:25 INFO Hourly dataset computed and listing created
2026-06-21 20:25:32 INFO Hourly dataset computed
2026-06-21 20:25:32 INFO Creating directories and links for ENS6 to run chimere's parallel part
2026-06-21 20:25:32 INFO Linking EMIS ...
2026-06-21 20:25:33 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens06.nc
2026-06-21 20:25:33 INFO >> Checking links...
2026-06-21 20:25:35 INFO >> All links are good for ENS6  ...
2026-06-21 20:25:35 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:25:36 INFO Hourly dataset computed and listing created
2026-06-21 20:25:44 INFO Hourly dataset computed
2026-06-21 20:25:44 INFO Creating directories and links for ENS7 to run chimere's parallel part
2026-06-21 20:25:44 INFO Linking EMIS ...
2026-06-21 20:25:44 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens07.nc
2026-06-21 20:25:44 INFO >> Checking links...
2026-06-21 20:25:47 INFO >> All links are good for ENS7  ...
2026-06-21 20:25:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:25:48 INFO Hourly dataset computed and listing created
2026-06-21 20:25:53 INFO Hourly dataset computed
2026-06-21 20:25:53 INFO Creating directories and links for ENS8 to run chimere's parallel part
2026-06-21 20:25:53 INFO Linking EMIS ...
2026-06-21 20:25:54 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens08.nc
2026-06-21 20:25:54 INFO >> Checking links...
2026-06-21 20:25:56 INFO >> All links are good for ENS8  ...
2026-06-21 20:25:56 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:25:57 INFO Hourly dataset computed and listing created
2026-06-21 20:25:58 INFO Hourly dataset computed
2026-06-21 20:25:58 INFO Creating directories and links for ENS9 to run chimere's parallel part
2026-06-21 20:25:58 INFO Linking EMIS ...
2026-06-21 20:25:58 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens09.nc
2026-06-21 20:25:58 INFO >> Checking links...
2026-06-21 20:26:01 INFO >> All links are good for ENS9  ...
2026-06-21 20:26:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:26:02 INFO Hourly dataset computed and listing created
2026-06-21 20:26:03 INFO Hourly dataset computed
2026-06-21 20:26:03 INFO Creating directories and links for ENS10 to run chimere's parallel part
2026-06-21 20:26:03 INFO Linking EMIS ...
2026-06-21 20:26:03 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens10.nc
2026-06-21 20:26:03 INFO >> Checking links...
2026-06-21 20:26:05 INFO >> All links are good for ENS10  ...
2026-06-21 20:26:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:26:06 INFO Hourly dataset computed and listing created
2026-06-21 20:26:07 INFO Hourly dataset computed
2026-06-21 20:26:07 INFO Creating directories and links for ENS11 to run chimere's parallel part
2026-06-21 20:26:07 INFO Linking EMIS ...
2026-06-21 20:26:08 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens11.nc
2026-06-21 20:26:08 INFO >> Checking links...
2026-06-21 20:26:10 INFO >> All links are good for ENS11  ...
2026-06-21 20:26:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:26:11 INFO Hourly dataset computed and listing created
2026-06-21 20:26:12 INFO Hourly dataset computed
2026-06-21 20:26:12 INFO Creating directories and links for ENS12 to run chimere's parallel part
2026-06-21 20:26:12 INFO Linking EMIS ...
2026-06-21 20:26:13 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens12.nc
2026-06-21 20:26:13 INFO >> Checking links...
2026-06-21 20:26:15 INFO >> All links are good for ENS12  ...
2026-06-21 20:26:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:26:16 INFO Hourly dataset computed and listing created
2026-06-21 20:26:30 INFO Hourly dataset computed
2026-06-21 20:26:30 INFO Creating directories and links for ENS13 to run chimere's parallel part
2026-06-21 20:26:30 INFO Linking EMIS ...
2026-06-21 20:26:30 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens13.nc
2026-06-21 20:26:30 INFO >> Checking links...
2026-06-21 20:26:33 INFO >> All links are good for ENS13  ...
2026-06-21 20:26:33 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:26:34 INFO Hourly dataset computed and listing created
2026-06-21 20:26:35 INFO Hourly dataset computed
2026-06-21 20:26:35 INFO Creating directories and links for ENS14 to run chimere's parallel part
2026-06-21 20:26:35 INFO Linking EMIS ...
2026-06-21 20:26:35 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens14.nc
2026-06-21 20:26:35 INFO >> Checking links...
2026-06-21 20:26:38 INFO >> All links are good for ENS14  ...
2026-06-21 20:26:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:26:39 INFO Hourly dataset computed and listing created
2026-06-21 20:26:45 INFO Hourly dataset computed
2026-06-21 20:26:45 INFO Creating directories and links for ENS15 to run chimere's parallel part
2026-06-21 20:26:45 INFO Linking EMIS ...
2026-06-21 20:26:45 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens15.nc
2026-06-21 20:26:45 INFO >> Checking links...
2026-06-21 20:26:48 INFO >> All links are good for ENS15  ...
2026-06-21 20:26:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:26:49 INFO Hourly dataset computed and listing created
2026-06-21 20:26:49 INFO Hourly dataset computed
2026-06-21 20:26:49 INFO ---------->>> Running CHIMERE model from 2020-02-15 13:00:00 to 2020-02-15 15:00:00
2026-06-21 20:26:49 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:26:49 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1
2026-06-21 20:26:49 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/end.2020021511_2_ENS1.nc
2026-06-21 20:26:49 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-06-21 20:26:49 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:26:49 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-06-21 20:26:49 INFO Queuing job for member 1...
2026-06-21 20:26:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:26:49 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-06-21 20:26:50 INFO Found: ['4938155']
2026-06-21 20:26:55 INFO [TGCC-IRENE] Submitted job with ID:['4938155']
2026-06-21 20:26:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:26:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2
2026-06-21 20:26:55 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/end.2020021511_2_ENS2.nc
2026-06-21 20:26:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-06-21 20:26:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:26:55 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-06-21 20:26:55 INFO Queuing job for member 2...
2026-06-21 20:26:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:26:55 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-06-21 20:26:56 INFO Found: ['4938156']
2026-06-21 20:27:01 INFO [TGCC-IRENE] Submitted job with ID:['4938156']
2026-06-21 20:27:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:27:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3
2026-06-21 20:27:01 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/end.2020021511_2_ENS3.nc
2026-06-21 20:27:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-06-21 20:27:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:27:01 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-06-21 20:27:01 INFO Queuing job for member 3...
2026-06-21 20:27:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:27:01 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-06-21 20:27:02 INFO Found: ['4938159']
2026-06-21 20:27:07 INFO [TGCC-IRENE] Submitted job with ID:['4938159']
2026-06-21 20:27:07 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:27:07 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4
2026-06-21 20:27:07 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/end.2020021511_2_ENS4.nc
2026-06-21 20:27:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-06-21 20:27:07 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:27:07 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-06-21 20:27:07 INFO Queuing job for member 4...
2026-06-21 20:27:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:27:07 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-06-21 20:27:08 INFO Found: ['4938160']
2026-06-21 20:27:13 INFO [TGCC-IRENE] Submitted job with ID:['4938160']
2026-06-21 20:27:13 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:27:13 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5
2026-06-21 20:27:13 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/end.2020021511_2_ENS5.nc
2026-06-21 20:27:13 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-06-21 20:27:13 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:27:13 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-06-21 20:27:13 INFO Queuing job for member 5...
2026-06-21 20:27:13 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:27:13 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-06-21 20:27:13 INFO Found: ['4938162']
2026-06-21 20:27:18 INFO [TGCC-IRENE] Submitted job with ID:['4938162']
2026-06-21 20:27:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:27:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6
2026-06-21 20:27:18 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/end.2020021511_2_ENS6.nc
2026-06-21 20:27:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-06-21 20:27:18 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:27:18 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-06-21 20:27:18 INFO Queuing job for member 6...
2026-06-21 20:27:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:27:18 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-06-21 20:27:19 INFO Found: ['4938163']
2026-06-21 20:27:24 INFO [TGCC-IRENE] Submitted job with ID:['4938163']
2026-06-21 20:27:24 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:27:24 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7
2026-06-21 20:27:24 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/end.2020021511_2_ENS7.nc
2026-06-21 20:27:24 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-06-21 20:27:24 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:27:24 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-06-21 20:27:24 INFO Queuing job for member 7...
2026-06-21 20:27:24 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:27:24 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-06-21 20:27:25 INFO Found: ['4938165']
2026-06-21 20:27:30 INFO [TGCC-IRENE] Submitted job with ID:['4938165']
2026-06-21 20:27:30 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:27:30 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8
2026-06-21 20:27:30 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/end.2020021511_2_ENS8.nc
2026-06-21 20:27:30 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-06-21 20:27:30 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:27:30 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-06-21 20:27:30 INFO Queuing job for member 8...
2026-06-21 20:27:30 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:27:30 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-06-21 20:27:31 INFO Found: ['4938166']
2026-06-21 20:27:36 INFO [TGCC-IRENE] Submitted job with ID:['4938166']
2026-06-21 20:27:36 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:27:36 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9
2026-06-21 20:27:36 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/end.2020021511_2_ENS9.nc
2026-06-21 20:27:36 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-06-21 20:27:36 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:27:36 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-06-21 20:27:36 INFO Queuing job for member 9...
2026-06-21 20:27:36 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:27:36 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-06-21 20:27:36 INFO Found: ['4938167']
2026-06-21 20:27:41 INFO [TGCC-IRENE] Submitted job with ID:['4938167']
2026-06-21 20:27:41 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:27:41 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10
2026-06-21 20:27:41 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/end.2020021511_2_ENS10.nc
2026-06-21 20:27:41 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-06-21 20:27:41 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:27:41 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-06-21 20:27:41 INFO Queuing job for member 10...
2026-06-21 20:27:41 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:27:41 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-06-21 20:27:42 INFO Found: ['4938169']
2026-06-21 20:27:47 INFO [TGCC-IRENE] Submitted job with ID:['4938169']
2026-06-21 20:27:47 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:27:47 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11
2026-06-21 20:27:47 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/end.2020021511_2_ENS11.nc
2026-06-21 20:27:47 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-06-21 20:27:47 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:27:47 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-06-21 20:27:47 INFO Queuing job for member 11...
2026-06-21 20:27:47 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:27:47 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-06-21 20:27:48 INFO Found: ['4938176']
2026-06-21 20:27:53 INFO [TGCC-IRENE] Submitted job with ID:['4938176']
2026-06-21 20:27:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:27:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12
2026-06-21 20:27:53 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/end.2020021511_2_ENS12.nc
2026-06-21 20:27:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-06-21 20:27:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:27:53 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-06-21 20:27:53 INFO Queuing job for member 12...
2026-06-21 20:27:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:27:53 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-06-21 20:27:54 INFO Found: ['4938179']
2026-06-21 20:27:59 INFO [TGCC-IRENE] Submitted job with ID:['4938179']
2026-06-21 20:27:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:27:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13
2026-06-21 20:27:59 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/end.2020021511_2_ENS13.nc
2026-06-21 20:27:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-06-21 20:27:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:27:59 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-06-21 20:27:59 INFO Queuing job for member 13...
2026-06-21 20:27:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:27:59 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-06-21 20:28:00 INFO Found: ['4938180']
2026-06-21 20:28:05 INFO [TGCC-IRENE] Submitted job with ID:['4938180']
2026-06-21 20:28:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:28:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14
2026-06-21 20:28:05 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/end.2020021511_2_ENS14.nc
2026-06-21 20:28:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-06-21 20:28:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:28:05 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-06-21 20:28:05 INFO Queuing job for member 14...
2026-06-21 20:28:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:28:05 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-06-21 20:28:05 INFO Found: ['4938182']
2026-06-21 20:28:10 INFO [TGCC-IRENE] Submitted job with ID:['4938182']
2026-06-21 20:28:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:28:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15
2026-06-21 20:28:10 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/end.2020021511_2_ENS15.nc
2026-06-21 20:28:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-06-21 20:28:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:28:10 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-06-21 20:28:10 INFO Queuing job for member 15...
2026-06-21 20:28:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:28:10 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-06-21 20:28:11 INFO Found: ['4938183']
2026-06-21 20:28:16 INFO [TGCC-IRENE] Submitted job with ID:['4938183']
2026-06-21 20:28:16 INFO Checking job status ...
2026-06-21 20:28:16 INFO None 4938155: status RUNNING/PENDING
2026-06-21 20:28:16 INFO None 4938156: status RUNNING/PENDING
2026-06-21 20:28:16 INFO None 4938159: status RUNNING/PENDING
2026-06-21 20:28:16 INFO None 4938160: status RUNNING/PENDING
2026-06-21 20:28:16 INFO None 4938162: status RUNNING/PENDING
2026-06-21 20:28:16 INFO None 4938163: status RUNNING/PENDING
2026-06-21 20:28:16 INFO None 4938165: status RUNNING/PENDING
2026-06-21 20:28:16 INFO None 4938166: status RUNNING/PENDING
2026-06-21 20:28:16 INFO None 4938167: status RUNNING/PENDING
2026-06-21 20:28:16 INFO None 4938169: status RUNNING/PENDING
2026-06-21 20:28:16 INFO None 4938176: status RUNNING/PENDING
2026-06-21 20:28:16 INFO None 4938179: status RUNNING/PENDING
2026-06-21 20:28:16 INFO None 4938180: status RUNNING/PENDING
2026-06-21 20:28:16 INFO None 4938182: status RUNNING/PENDING
2026-06-21 20:28:16 INFO None 4938183: status RUNNING/PENDING
2026-06-21 20:28:16 INFO Jobs still running: ['4938155', '4938156', '4938159', '4938160', '4938162', '4938163', '4938165', '4938166', '4938167', '4938169', '4938176', '4938179', '4938180', '4938182', '4938183']. Waiting...
2026-06-21 20:28:31 INFO None 4938155: status RUNNING/PENDING
2026-06-21 20:28:31 INFO None 4938156: status RUNNING/PENDING
2026-06-21 20:28:32 INFO None 4938159: status RUNNING/PENDING
2026-06-21 20:28:32 INFO None 4938160: status RUNNING/PENDING
2026-06-21 20:28:32 INFO None 4938162: status RUNNING/PENDING
2026-06-21 20:28:32 INFO None 4938163: status RUNNING/PENDING
2026-06-21 20:28:32 INFO None 4938165: status RUNNING/PENDING
2026-06-21 20:28:32 INFO None 4938166: status RUNNING/PENDING
2026-06-21 20:28:32 INFO None 4938167: status RUNNING/PENDING
2026-06-21 20:28:32 INFO None 4938169: status RUNNING/PENDING
2026-06-21 20:28:32 INFO None 4938176: status RUNNING/PENDING
2026-06-21 20:28:32 INFO None 4938179: status RUNNING/PENDING
2026-06-21 20:28:32 INFO None 4938180: status RUNNING/PENDING
2026-06-21 20:28:32 INFO None 4938182: status RUNNING/PENDING
2026-06-21 20:28:32 INFO None 4938183: status RUNNING/PENDING
2026-06-21 20:28:32 INFO Jobs still running: ['4938155', '4938156', '4938159', '4938160', '4938162', '4938163', '4938165', '4938166', '4938167', '4938169', '4938176', '4938179', '4938180', '4938182', '4938183']. Waiting...
2026-06-21 20:28:47 INFO None 4938155: status RUNNING/PENDING
2026-06-21 20:28:47 INFO None 4938156: status RUNNING/PENDING
2026-06-21 20:28:47 INFO None 4938159: status RUNNING/PENDING
2026-06-21 20:28:47 INFO None 4938160: status RUNNING/PENDING
2026-06-21 20:28:47 INFO None 4938162: status RUNNING/PENDING
2026-06-21 20:28:47 INFO None 4938163: status RUNNING/PENDING
2026-06-21 20:28:47 INFO None 4938165: status RUNNING/PENDING
2026-06-21 20:28:47 INFO None 4938166: status RUNNING/PENDING
2026-06-21 20:28:47 INFO None 4938167: status RUNNING/PENDING
2026-06-21 20:28:47 INFO None 4938169: status RUNNING/PENDING
2026-06-21 20:28:47 INFO None 4938176: status RUNNING/PENDING
2026-06-21 20:28:47 INFO None 4938179: status RUNNING/PENDING
2026-06-21 20:28:47 INFO None 4938180: status RUNNING/PENDING
2026-06-21 20:28:47 INFO None 4938182: status RUNNING/PENDING
2026-06-21 20:28:47 INFO None 4938183: status RUNNING/PENDING
2026-06-21 20:28:47 INFO Jobs still running: ['4938155', '4938156', '4938159', '4938160', '4938162', '4938163', '4938165', '4938166', '4938167', '4938169', '4938176', '4938179', '4938180', '4938182', '4938183']. Waiting...
2026-06-21 20:29:02 INFO None 4938155: status RUNNING/PENDING
2026-06-21 20:29:02 INFO None 4938156: status RUNNING/PENDING
2026-06-21 20:29:02 INFO None 4938159: status RUNNING/PENDING
2026-06-21 20:29:02 INFO None 4938160: status RUNNING/PENDING
2026-06-21 20:29:02 INFO None 4938162: status RUNNING/PENDING
2026-06-21 20:29:02 INFO None 4938163: status RUNNING/PENDING
2026-06-21 20:29:02 INFO None 4938165: status RUNNING/PENDING
2026-06-21 20:29:02 INFO None 4938166: status RUNNING/PENDING
2026-06-21 20:29:02 INFO None 4938167: status RUNNING/PENDING
2026-06-21 20:29:02 INFO None 4938169: status RUNNING/PENDING
2026-06-21 20:29:02 INFO None 4938176: status RUNNING/PENDING
2026-06-21 20:29:02 INFO None 4938179: status RUNNING/PENDING
2026-06-21 20:29:02 INFO None 4938180: status RUNNING/PENDING
2026-06-21 20:29:02 INFO None 4938182: status RUNNING/PENDING
2026-06-21 20:29:02 INFO None 4938183: status RUNNING/PENDING
2026-06-21 20:29:02 INFO Jobs still running: ['4938155', '4938156', '4938159', '4938160', '4938162', '4938163', '4938165', '4938166', '4938167', '4938169', '4938176', '4938179', '4938180', '4938182', '4938183']. Waiting...
2026-06-21 20:29:17 INFO None 4938155: status RUNNING/PENDING
2026-06-21 20:29:17 INFO None 4938156: status RUNNING/PENDING
2026-06-21 20:29:17 INFO None 4938159: status RUNNING/PENDING
2026-06-21 20:29:17 INFO None 4938160: status RUNNING/PENDING
2026-06-21 20:29:17 INFO None 4938162: status RUNNING/PENDING
2026-06-21 20:29:17 INFO None 4938163: status RUNNING/PENDING
2026-06-21 20:29:17 INFO None 4938165: status RUNNING/PENDING
2026-06-21 20:29:17 INFO None 4938166: status RUNNING/PENDING
2026-06-21 20:29:17 INFO None 4938167: status RUNNING/PENDING
2026-06-21 20:29:18 INFO None 4938169: status RUNNING/PENDING
2026-06-21 20:29:18 INFO None 4938176: status RUNNING/PENDING
2026-06-21 20:29:18 INFO None 4938179: status RUNNING/PENDING
2026-06-21 20:29:18 INFO None 4938180: status RUNNING/PENDING
2026-06-21 20:29:18 INFO None 4938182: status RUNNING/PENDING
2026-06-21 20:29:18 INFO None 4938183: status RUNNING/PENDING
2026-06-21 20:29:18 INFO Jobs still running: ['4938155', '4938156', '4938159', '4938160', '4938162', '4938163', '4938165', '4938166', '4938167', '4938169', '4938176', '4938179', '4938180', '4938182', '4938183']. Waiting...
2026-06-21 20:29:33 INFO None 4938155: status RUNNING/PENDING
2026-06-21 20:29:33 INFO None 4938156: status FINISHED
2026-06-21 20:29:33 INFO None 4938159: status RUNNING/PENDING
2026-06-21 20:29:33 INFO None 4938160: status RUNNING/PENDING
2026-06-21 20:29:33 INFO None 4938162: status RUNNING/PENDING
2026-06-21 20:29:33 INFO None 4938163: status RUNNING/PENDING
2026-06-21 20:29:33 INFO None 4938165: status RUNNING/PENDING
2026-06-21 20:29:33 INFO None 4938166: status RUNNING/PENDING
2026-06-21 20:29:33 INFO None 4938167: status RUNNING/PENDING
2026-06-21 20:29:33 INFO None 4938169: status RUNNING/PENDING
2026-06-21 20:29:33 INFO None 4938176: status RUNNING/PENDING
2026-06-21 20:29:33 INFO None 4938179: status RUNNING/PENDING
2026-06-21 20:29:33 INFO None 4938180: status RUNNING/PENDING
2026-06-21 20:29:33 INFO None 4938182: status RUNNING/PENDING
2026-06-21 20:29:33 INFO None 4938183: status RUNNING/PENDING
2026-06-21 20:29:33 INFO Jobs still running: ['4938155', '4938159', '4938160', '4938162', '4938163', '4938165', '4938166', '4938167', '4938169', '4938176', '4938179', '4938180', '4938182', '4938183']. Waiting...
2026-06-21 20:29:48 INFO None 4938155: status RUNNING/PENDING
2026-06-21 20:29:48 INFO None 4938156: status FINISHED
2026-06-21 20:29:48 INFO None 4938159: status FINISHED
2026-06-21 20:29:48 INFO None 4938160: status RUNNING/PENDING
2026-06-21 20:29:48 INFO None 4938162: status RUNNING/PENDING
2026-06-21 20:29:48 INFO None 4938163: status RUNNING/PENDING
2026-06-21 20:29:48 INFO None 4938165: status RUNNING/PENDING
2026-06-21 20:29:48 INFO None 4938166: status RUNNING/PENDING
2026-06-21 20:29:48 INFO None 4938167: status RUNNING/PENDING
2026-06-21 20:29:48 INFO None 4938169: status RUNNING/PENDING
2026-06-21 20:29:48 INFO None 4938176: status RUNNING/PENDING
2026-06-21 20:29:48 INFO None 4938179: status RUNNING/PENDING
2026-06-21 20:29:48 INFO None 4938180: status RUNNING/PENDING
2026-06-21 20:29:48 INFO None 4938182: status RUNNING/PENDING
2026-06-21 20:29:48 INFO None 4938183: status RUNNING/PENDING
2026-06-21 20:29:48 INFO Jobs still running: ['4938155', '4938160', '4938162', '4938163', '4938165', '4938166', '4938167', '4938169', '4938176', '4938179', '4938180', '4938182', '4938183']. Waiting...
2026-06-21 20:30:03 INFO None 4938155: status RUNNING/PENDING
2026-06-21 20:30:03 INFO None 4938156: status FINISHED
2026-06-21 20:30:03 INFO None 4938159: status FINISHED
2026-06-21 20:30:03 INFO None 4938160: status FINISHED
2026-06-21 20:30:03 INFO None 4938162: status RUNNING/PENDING
2026-06-21 20:30:03 INFO None 4938163: status FINISHED
2026-06-21 20:30:03 INFO None 4938165: status FINISHED
2026-06-21 20:30:03 INFO None 4938166: status RUNNING/PENDING
2026-06-21 20:30:04 INFO None 4938167: status RUNNING/PENDING
2026-06-21 20:30:04 INFO None 4938169: status RUNNING/PENDING
2026-06-21 20:30:04 INFO None 4938176: status RUNNING/PENDING
2026-06-21 20:30:04 INFO None 4938179: status RUNNING/PENDING
2026-06-21 20:30:04 INFO None 4938180: status RUNNING/PENDING
2026-06-21 20:30:04 INFO None 4938182: status RUNNING/PENDING
2026-06-21 20:30:04 INFO None 4938183: status RUNNING/PENDING
2026-06-21 20:30:04 INFO Jobs still running: ['4938155', '4938162', '4938166', '4938167', '4938169', '4938176', '4938179', '4938180', '4938182', '4938183']. Waiting...
2026-06-21 20:30:19 INFO None 4938155: status RUNNING/PENDING
2026-06-21 20:30:19 INFO None 4938156: status FINISHED
2026-06-21 20:30:19 INFO None 4938159: status FINISHED
2026-06-21 20:30:19 INFO None 4938160: status FINISHED
2026-06-21 20:30:19 INFO None 4938162: status RUNNING/PENDING
2026-06-21 20:30:19 INFO None 4938163: status FINISHED
2026-06-21 20:30:19 INFO None 4938165: status FINISHED
2026-06-21 20:30:19 INFO None 4938166: status FINISHED
2026-06-21 20:30:19 INFO None 4938167: status FINISHED
2026-06-21 20:30:19 INFO None 4938169: status FINISHED
2026-06-21 20:30:19 INFO None 4938176: status RUNNING/PENDING
2026-06-21 20:30:19 INFO None 4938179: status RUNNING/PENDING
2026-06-21 20:30:19 INFO None 4938180: status RUNNING/PENDING
2026-06-21 20:30:19 INFO None 4938182: status RUNNING/PENDING
2026-06-21 20:30:19 INFO None 4938183: status RUNNING/PENDING
2026-06-21 20:30:19 INFO Jobs still running: ['4938155', '4938162', '4938176', '4938179', '4938180', '4938182', '4938183']. Waiting...
2026-06-21 20:30:34 INFO None 4938155: status RUNNING/PENDING
2026-06-21 20:30:34 INFO None 4938156: status FINISHED
2026-06-21 20:30:34 INFO None 4938159: status FINISHED
2026-06-21 20:30:34 INFO None 4938160: status FINISHED
2026-06-21 20:30:34 INFO None 4938162: status FINISHED
2026-06-21 20:30:34 INFO None 4938163: status FINISHED
2026-06-21 20:30:34 INFO None 4938165: status FINISHED
2026-06-21 20:30:34 INFO None 4938166: status FINISHED
2026-06-21 20:30:34 INFO None 4938167: status FINISHED
2026-06-21 20:30:34 INFO None 4938169: status FINISHED
2026-06-21 20:30:34 INFO None 4938176: status FINISHED
2026-06-21 20:30:34 INFO None 4938179: status RUNNING/PENDING
2026-06-21 20:30:34 INFO None 4938180: status RUNNING/PENDING
2026-06-21 20:30:34 INFO None 4938182: status RUNNING/PENDING
2026-06-21 20:30:34 INFO None 4938183: status RUNNING/PENDING
2026-06-21 20:30:34 INFO Jobs still running: ['4938155', '4938179', '4938180', '4938182', '4938183']. Waiting...
2026-06-21 20:30:49 INFO None 4938155: status RUNNING/PENDING
2026-06-21 20:30:49 INFO None 4938156: status FINISHED
2026-06-21 20:30:49 INFO None 4938159: status FINISHED
2026-06-21 20:30:49 INFO None 4938160: status FINISHED
2026-06-21 20:30:49 INFO None 4938162: status FINISHED
2026-06-21 20:30:49 INFO None 4938163: status FINISHED
2026-06-21 20:30:49 INFO None 4938165: status FINISHED
2026-06-21 20:30:49 INFO None 4938166: status FINISHED
2026-06-21 20:30:49 INFO None 4938167: status FINISHED
2026-06-21 20:30:49 INFO None 4938169: status FINISHED
2026-06-21 20:30:49 INFO None 4938176: status FINISHED
2026-06-21 20:30:49 INFO None 4938179: status RUNNING/PENDING
2026-06-21 20:30:49 INFO None 4938180: status RUNNING/PENDING
2026-06-21 20:30:49 INFO None 4938182: status RUNNING/PENDING
2026-06-21 20:30:49 INFO None 4938183: status RUNNING/PENDING
2026-06-21 20:30:49 INFO Jobs still running: ['4938155', '4938179', '4938180', '4938182', '4938183']. Waiting...
2026-06-21 20:31:04 INFO None 4938155: status RUNNING/PENDING
2026-06-21 20:31:05 INFO None 4938156: status FINISHED
2026-06-21 20:31:05 INFO None 4938159: status FINISHED
2026-06-21 20:31:05 INFO None 4938160: status FINISHED
2026-06-21 20:31:05 INFO None 4938162: status FINISHED
2026-06-21 20:31:05 INFO None 4938163: status FINISHED
2026-06-21 20:31:05 INFO None 4938165: status FINISHED
2026-06-21 20:31:05 INFO None 4938166: status FINISHED
2026-06-21 20:31:05 INFO None 4938167: status FINISHED
2026-06-21 20:31:05 INFO None 4938169: status FINISHED
2026-06-21 20:31:05 INFO None 4938176: status FINISHED
2026-06-21 20:31:05 INFO None 4938179: status FINISHED
2026-06-21 20:31:05 INFO None 4938180: status FINISHED
2026-06-21 20:31:05 INFO None 4938182: status FINISHED
2026-06-21 20:31:05 INFO None 4938183: status FINISHED
2026-06-21 20:31:05 INFO Jobs still running: ['4938155']. Waiting...
2026-06-21 20:31:20 INFO None 4938155: status FINISHED
2026-06-21 20:31:20 INFO None 4938156: status FINISHED
2026-06-21 20:31:20 INFO None 4938159: status FINISHED
2026-06-21 20:31:20 INFO None 4938160: status FINISHED
2026-06-21 20:31:20 INFO None 4938162: status FINISHED
2026-06-21 20:31:20 INFO None 4938163: status FINISHED
2026-06-21 20:31:20 INFO None 4938165: status FINISHED
2026-06-21 20:31:20 INFO None 4938166: status FINISHED
2026-06-21 20:31:20 INFO None 4938167: status FINISHED
2026-06-21 20:31:20 INFO None 4938169: status FINISHED
2026-06-21 20:31:20 INFO None 4938176: status FINISHED
2026-06-21 20:31:20 INFO None 4938179: status FINISHED
2026-06-21 20:31:20 INFO None 4938180: status FINISHED
2026-06-21 20:31:20 INFO None 4938182: status FINISHED
2026-06-21 20:31:20 INFO None 4938183: status FINISHED
2026-06-21 20:31:20 INFO Jobs ['4938155', '4938156', '4938159', '4938160', '4938162', '4938163', '4938165', '4938166', '4938167', '4938169', '4938176', '4938179', '4938180', '4938182', '4938183'] have finished
2026-06-21 20:31:20 INFO Checking restart files were created ...
2026-06-21 20:31:20 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/end.2020021513_2_ENS1.nc(1002685915 bytes)
2026-06-21 20:31:20 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/end.2020021513_2_ENS2.nc(1002685915 bytes)
2026-06-21 20:31:20 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/end.2020021513_2_ENS3.nc(1002685915 bytes)
2026-06-21 20:31:20 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/end.2020021513_2_ENS4.nc(1002685915 bytes)
2026-06-21 20:31:20 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/end.2020021513_2_ENS5.nc(1002685915 bytes)
2026-06-21 20:31:20 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/end.2020021513_2_ENS6.nc(1002685915 bytes)
2026-06-21 20:31:20 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/end.2020021513_2_ENS7.nc(1002685915 bytes)
2026-06-21 20:31:20 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/end.2020021513_2_ENS8.nc(1002685915 bytes)
2026-06-21 20:31:20 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/end.2020021513_2_ENS9.nc(1002685915 bytes)
2026-06-21 20:31:20 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/end.2020021513_2_ENS10.nc(1002685915 bytes)
2026-06-21 20:31:20 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/end.2020021513_2_ENS11.nc(1002685915 bytes)
2026-06-21 20:31:20 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/end.2020021513_2_ENS12.nc(1002685915 bytes)
2026-06-21 20:31:20 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/end.2020021513_2_ENS13.nc(1002685915 bytes)
2026-06-21 20:31:20 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/end.2020021513_2_ENS14.nc(1002685915 bytes)
2026-06-21 20:31:20 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/end.2020021513_2_ENS15.nc(1002685915 bytes)
2026-06-21 20:31:20 INFO  Run_model() completed successfully.
2026-06-21 20:31:20 INFO [TIME] after_model_set_simulated_time current_time=2020-02-15 13:00:00 simulated_time=2020-02-15 15:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 20:31:20 INFO [TIME] gregorian_conversion simulated_time=2020-02-15 15:00:00 days=153081 seconds=54000
2026-06-21 20:31:20 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-06-21 20:31:20 INFO [TIME] increment current_time 2020-02-15 13:00:00 -> 2020-02-15 15:00:00
2026-06-21 20:31:20 INFO [TIME] after_increment_before_assimilation current_time=2020-02-15 15:00:00 simulated_time=2020-02-15 15:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 20:31:20 INFO ---------->>> Running process_satellite_data()
2026-06-21 20:31:20 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12133.nc
2026-06-21 20:31:20 INFO ---------->>> Running run_obs_converter()
2026-06-21 20:31:20 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_53128_153081.out
2026-06-21 20:31:20 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_53128_153081.out
2026-06-21 20:31:20 INFO ---------->>> Running DART
2026-06-21 20:31:20 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-06-21 20:31:20 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/chim_ENS1_2020021513_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/chim_ENS1_2020021515_1_out_toDART.nc
2026-06-21 20:31:20 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/chim_ENS2_2020021513_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/chim_ENS2_2020021515_1_out_toDART.nc
2026-06-21 20:31:21 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/chim_ENS3_2020021513_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/chim_ENS3_2020021515_1_out_toDART.nc
2026-06-21 20:31:21 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/chim_ENS4_2020021513_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/chim_ENS4_2020021515_1_out_toDART.nc
2026-06-21 20:31:21 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/chim_ENS5_2020021513_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/chim_ENS5_2020021515_1_out_toDART.nc
2026-06-21 20:31:22 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/chim_ENS6_2020021513_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/chim_ENS6_2020021515_1_out_toDART.nc
2026-06-21 20:31:22 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/chim_ENS7_2020021513_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/chim_ENS7_2020021515_1_out_toDART.nc
2026-06-21 20:31:23 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/chim_ENS8_2020021513_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/chim_ENS8_2020021515_1_out_toDART.nc
2026-06-21 20:31:23 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/chim_ENS9_2020021513_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/chim_ENS9_2020021515_1_out_toDART.nc
2026-06-21 20:31:23 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/chim_ENS10_2020021513_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/chim_ENS10_2020021515_1_out_toDART.nc
2026-06-21 20:31:24 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/chim_ENS11_2020021513_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/chim_ENS11_2020021515_1_out_toDART.nc
2026-06-21 20:31:24 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/chim_ENS12_2020021513_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/chim_ENS12_2020021515_1_out_toDART.nc
2026-06-21 20:31:24 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/chim_ENS13_2020021513_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/chim_ENS13_2020021515_1_out_toDART.nc
2026-06-21 20:31:25 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/chim_ENS14_2020021513_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/chim_ENS14_2020021515_1_out_toDART.nc
2026-06-21 20:31:25 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/chim_ENS15_2020021513_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/chim_ENS15_2020021515_1_out_toDART.nc
2026-06-21 20:31:25 INFO Replacement input_template.nml → input.nml completed successfully.
2026-06-21 20:31:25 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-06-21 20:31:25 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-06-21 20:31:25 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-06-21 20:31:25 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-06-21 20:31:25 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-06-21 20:31:32 INFO Found: []
2026-06-21 20:31:32 INFO No job id returned by command ./run_filter.bsh
2026-06-21 20:31:32 INFO No monitoring will be performed
2026-06-21 20:31:32 INFO Moving DART output files to analysis and preassim directories for date 2020021515 if present ...
2026-06-21 20:31:32 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021515'
2026-06-21 20:31:32 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021515'
2026-06-21 20:31:32 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021515'
2026-06-21 20:31:32 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021515'
2026-06-21 20:31:32 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021515'
2026-06-21 20:31:32 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021515'
2026-06-21 20:31:32 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021515'
2026-06-21 20:31:32 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021515'
2026-06-21 20:31:32 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021515'
2026-06-21 20:31:32 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021515'
2026-06-21 20:31:32 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021515'
2026-06-21 20:31:32 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021515'
2026-06-21 20:31:32 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021515'
2026-06-21 20:31:32 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021515'
2026-06-21 20:31:32 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021515'
2026-06-21 20:31:32 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021515'
2026-06-21 20:31:32 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021515'
2026-06-21 20:31:32 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021515'
2026-06-21 20:31:32 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021515'
2026-06-21 20:31:32 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021515'
2026-06-21 20:31:33 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021515'
2026-06-21 20:31:33 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021515'
2026-06-21 20:31:33 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021515'
2026-06-21 20:31:33 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021515'
2026-06-21 20:31:33 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021515'
2026-06-21 20:31:33 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021515'
2026-06-21 20:31:33 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021515'
2026-06-21 20:31:33 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021515'
2026-06-21 20:31:33 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021515'
2026-06-21 20:31:33 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021515'
2026-06-21 20:31:33 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021515'
2026-06-21 20:31:33 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020021515'
2026-06-21 20:31:33 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021515'
2026-06-21 20:31:33 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020021515'
2026-06-21 20:31:33 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-06-21 20:31:33 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-06-21 20:31:33 INFO run_dart() is DONE.
2026-06-21 20:31:33 INFO ---------->>> Running update_pollutant_in_end()
2026-06-21 20:31:33 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/end.2020021513_2_ENS1.nc
2026-06-21 20:31:38 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:31:38 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:31:39 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:31:39 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:31:40 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:31:40 INFO [Posterior Diff ENS1] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021515/diff_posterior_ENS1_2020021515.nc
2026-06-21 20:31:40 INFO [Posterior Diff ENS1] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021515/diff_posterior_ENS1_2020021515.relative.nc
2026-06-21 20:31:41 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/end.2020021513_2_ENS2.nc
2026-06-21 20:31:46 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:31:46 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:31:47 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:31:47 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:31:47 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:31:48 INFO [Posterior Diff ENS2] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021515/diff_posterior_ENS2_2020021515.nc
2026-06-21 20:31:48 INFO [Posterior Diff ENS2] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021515/diff_posterior_ENS2_2020021515.relative.nc
2026-06-21 20:31:48 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/end.2020021513_2_ENS3.nc
2026-06-21 20:31:53 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:31:53 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:31:54 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:31:54 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:31:55 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:31:55 INFO [Posterior Diff ENS3] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021515/diff_posterior_ENS3_2020021515.nc
2026-06-21 20:31:55 INFO [Posterior Diff ENS3] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021515/diff_posterior_ENS3_2020021515.relative.nc
2026-06-21 20:31:56 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/end.2020021513_2_ENS4.nc
2026-06-21 20:32:01 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:32:01 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:32:02 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:32:02 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:32:02 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:32:02 INFO [Posterior Diff ENS4] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021515/diff_posterior_ENS4_2020021515.nc
2026-06-21 20:32:02 INFO [Posterior Diff ENS4] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021515/diff_posterior_ENS4_2020021515.relative.nc
2026-06-21 20:32:03 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/end.2020021513_2_ENS5.nc
2026-06-21 20:32:08 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:32:08 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:32:09 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:32:09 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:32:10 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:32:10 INFO [Posterior Diff ENS5] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021515/diff_posterior_ENS5_2020021515.nc
2026-06-21 20:32:10 INFO [Posterior Diff ENS5] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021515/diff_posterior_ENS5_2020021515.relative.nc
2026-06-21 20:32:10 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/end.2020021513_2_ENS6.nc
2026-06-21 20:32:16 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:32:16 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:32:17 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:32:17 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:32:17 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:32:17 INFO [Posterior Diff ENS6] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021515/diff_posterior_ENS6_2020021515.nc
2026-06-21 20:32:17 INFO [Posterior Diff ENS6] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021515/diff_posterior_ENS6_2020021515.relative.nc
2026-06-21 20:32:18 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/end.2020021513_2_ENS7.nc
2026-06-21 20:32:23 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:32:23 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:32:24 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:32:24 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:32:24 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:32:25 INFO [Posterior Diff ENS7] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021515/diff_posterior_ENS7_2020021515.nc
2026-06-21 20:32:25 INFO [Posterior Diff ENS7] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021515/diff_posterior_ENS7_2020021515.relative.nc
2026-06-21 20:32:25 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/end.2020021513_2_ENS8.nc
2026-06-21 20:32:31 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:32:31 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:32:32 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:32:32 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:32:32 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:32:33 INFO [Posterior Diff ENS8] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021515/diff_posterior_ENS8_2020021515.nc
2026-06-21 20:32:33 INFO [Posterior Diff ENS8] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021515/diff_posterior_ENS8_2020021515.relative.nc
2026-06-21 20:32:33 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/end.2020021513_2_ENS9.nc
2026-06-21 20:32:38 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:32:38 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:32:39 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:32:39 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:32:40 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:32:40 INFO [Posterior Diff ENS9] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021515/diff_posterior_ENS9_2020021515.nc
2026-06-21 20:32:40 INFO [Posterior Diff ENS9] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021515/diff_posterior_ENS9_2020021515.relative.nc
2026-06-21 20:32:40 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/end.2020021513_2_ENS10.nc
2026-06-21 20:32:46 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:32:46 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:32:47 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:32:47 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:32:47 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:32:48 INFO [Posterior Diff ENS10] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021515/diff_posterior_ENS10_2020021515.nc
2026-06-21 20:32:48 INFO [Posterior Diff ENS10] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021515/diff_posterior_ENS10_2020021515.relative.nc
2026-06-21 20:32:48 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/end.2020021513_2_ENS11.nc
2026-06-21 20:32:53 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:32:53 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:32:54 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:32:54 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:32:55 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:32:55 INFO [Posterior Diff ENS11] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021515/diff_posterior_ENS11_2020021515.nc
2026-06-21 20:32:55 INFO [Posterior Diff ENS11] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021515/diff_posterior_ENS11_2020021515.relative.nc
2026-06-21 20:32:55 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/end.2020021513_2_ENS12.nc
2026-06-21 20:33:01 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:33:01 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:33:02 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:33:02 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:33:02 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:33:02 INFO [Posterior Diff ENS12] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021515/diff_posterior_ENS12_2020021515.nc
2026-06-21 20:33:02 INFO [Posterior Diff ENS12] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021515/diff_posterior_ENS12_2020021515.relative.nc
2026-06-21 20:33:03 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/end.2020021513_2_ENS13.nc
2026-06-21 20:33:08 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:33:08 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:33:09 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:33:09 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:33:10 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:33:10 INFO [Posterior Diff ENS13] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021515/diff_posterior_ENS13_2020021515.nc
2026-06-21 20:33:10 INFO [Posterior Diff ENS13] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021515/diff_posterior_ENS13_2020021515.relative.nc
2026-06-21 20:33:10 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/end.2020021513_2_ENS14.nc
2026-06-21 20:33:15 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:33:15 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:33:16 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:33:16 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:33:17 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:33:17 INFO [Posterior Diff ENS14] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021515/diff_posterior_ENS14_2020021515.nc
2026-06-21 20:33:17 INFO [Posterior Diff ENS14] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021515/diff_posterior_ENS14_2020021515.relative.nc
2026-06-21 20:33:17 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/end.2020021513_2_ENS15.nc
2026-06-21 20:33:23 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:33:23 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:33:24 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-21 20:33:24 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc using posterior/prior ratio.
2026-06-21 20:33:24 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-21 20:33:25 INFO [Posterior Diff ENS15] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021515/diff_posterior_ENS15_2020021515.nc
2026-06-21 20:33:25 INFO [Posterior Diff ENS15] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020021515/diff_posterior_ENS15_2020021515.relative.nc
2026-06-21 20:33:25 INFO Next run starts from 2020-02-15 15:00:00
2026-06-21 20:33:25 INFO Cycle is DONE; starting a new loop!
2026-06-21 20:33:25 INFO [TIME] step_end current_time=2020-02-15 15:00:00 simulated_time=2020-02-15 15:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 20:33:25 INFO [TIME] step_start current_time=2020-02-15 15:00:00 simulated_time=2020-02-15 15:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 20:33:25 INFO [TIME] window start=2020-02-15 15:00:00 end=2020-02-16 00:00:00 run_hours=9 has_assimilation=False
2026-06-21 20:33:25 INFO Creating directories and links for ENS1 to run chimere's parallel part
2026-06-21 20:33:25 INFO Linking EMIS ...
2026-06-21 20:33:25 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens01.nc
2026-06-21 20:33:25 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Sunday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Sunday.s.ens01.nc
2026-06-21 20:33:25 INFO >> Checking links...
2026-06-21 20:33:28 INFO >> All links are good for ENS1  ...
2026-06-21 20:33:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:33:30 INFO Hourly dataset computed and listing created
2026-06-21 20:33:46 INFO Hourly dataset computed
2026-06-21 20:33:46 INFO Creating directories and links for ENS2 to run chimere's parallel part
2026-06-21 20:33:46 INFO Linking EMIS ...
2026-06-21 20:33:46 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens02.nc
2026-06-21 20:33:47 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Sunday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Sunday.s.ens02.nc
2026-06-21 20:33:47 INFO >> Checking links...
2026-06-21 20:33:49 INFO >> All links are good for ENS2  ...
2026-06-21 20:33:49 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:33:51 INFO Hourly dataset computed and listing created
2026-06-21 20:34:00 INFO Hourly dataset computed
2026-06-21 20:34:00 INFO Creating directories and links for ENS3 to run chimere's parallel part
2026-06-21 20:34:00 INFO Linking EMIS ...
2026-06-21 20:34:00 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens03.nc
2026-06-21 20:34:01 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Sunday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Sunday.s.ens03.nc
2026-06-21 20:34:01 INFO >> Checking links...
2026-06-21 20:34:03 INFO >> All links are good for ENS3  ...
2026-06-21 20:34:03 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:34:05 INFO Hourly dataset computed and listing created
2026-06-21 20:34:14 INFO Hourly dataset computed
2026-06-21 20:34:14 INFO Creating directories and links for ENS4 to run chimere's parallel part
2026-06-21 20:34:14 INFO Linking EMIS ...
2026-06-21 20:34:14 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens04.nc
2026-06-21 20:34:14 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Sunday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Sunday.s.ens04.nc
2026-06-21 20:34:14 INFO >> Checking links...
2026-06-21 20:34:17 INFO >> All links are good for ENS4  ...
2026-06-21 20:34:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:34:18 INFO Hourly dataset computed and listing created
2026-06-21 20:34:21 INFO Hourly dataset computed
2026-06-21 20:34:21 INFO Creating directories and links for ENS5 to run chimere's parallel part
2026-06-21 20:34:21 INFO Linking EMIS ...
2026-06-21 20:34:21 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens05.nc
2026-06-21 20:34:22 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Sunday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Sunday.s.ens05.nc
2026-06-21 20:34:22 INFO >> Checking links...
2026-06-21 20:34:24 INFO >> All links are good for ENS5  ...
2026-06-21 20:34:24 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:34:25 INFO Hourly dataset computed and listing created
2026-06-21 20:34:27 INFO Hourly dataset computed
2026-06-21 20:34:28 INFO Creating directories and links for ENS6 to run chimere's parallel part
2026-06-21 20:34:28 INFO Linking EMIS ...
2026-06-21 20:34:28 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens06.nc
2026-06-21 20:34:28 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Sunday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Sunday.s.ens06.nc
2026-06-21 20:34:28 INFO >> Checking links...
2026-06-21 20:34:31 INFO >> All links are good for ENS6  ...
2026-06-21 20:34:31 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:34:32 INFO Hourly dataset computed and listing created
2026-06-21 20:34:35 INFO Hourly dataset computed
2026-06-21 20:34:35 INFO Creating directories and links for ENS7 to run chimere's parallel part
2026-06-21 20:34:35 INFO Linking EMIS ...
2026-06-21 20:34:35 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens07.nc
2026-06-21 20:34:35 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Sunday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Sunday.s.ens07.nc
2026-06-21 20:34:35 INFO >> Checking links...
2026-06-21 20:34:38 INFO >> All links are good for ENS7  ...
2026-06-21 20:34:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:34:39 INFO Hourly dataset computed and listing created
2026-06-21 20:34:42 INFO Hourly dataset computed
2026-06-21 20:34:42 INFO Creating directories and links for ENS8 to run chimere's parallel part
2026-06-21 20:34:42 INFO Linking EMIS ...
2026-06-21 20:34:42 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens08.nc
2026-06-21 20:34:42 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Sunday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Sunday.s.ens08.nc
2026-06-21 20:34:42 INFO >> Checking links...
2026-06-21 20:34:45 INFO >> All links are good for ENS8  ...
2026-06-21 20:34:45 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:34:46 INFO Hourly dataset computed and listing created
2026-06-21 20:34:49 INFO Hourly dataset computed
2026-06-21 20:34:49 INFO Creating directories and links for ENS9 to run chimere's parallel part
2026-06-21 20:34:49 INFO Linking EMIS ...
2026-06-21 20:34:49 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens09.nc
2026-06-21 20:34:50 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Sunday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Sunday.s.ens09.nc
2026-06-21 20:34:50 INFO >> Checking links...
2026-06-21 20:34:52 INFO >> All links are good for ENS9  ...
2026-06-21 20:34:52 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:34:54 INFO Hourly dataset computed and listing created
2026-06-21 20:34:56 INFO Hourly dataset computed
2026-06-21 20:34:56 INFO Creating directories and links for ENS10 to run chimere's parallel part
2026-06-21 20:34:56 INFO Linking EMIS ...
2026-06-21 20:34:56 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens10.nc
2026-06-21 20:34:57 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Sunday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Sunday.s.ens10.nc
2026-06-21 20:34:57 INFO >> Checking links...
2026-06-21 20:34:59 INFO >> All links are good for ENS10  ...
2026-06-21 20:34:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:35:00 INFO Hourly dataset computed and listing created
2026-06-21 20:35:03 INFO Hourly dataset computed
2026-06-21 20:35:03 INFO Creating directories and links for ENS11 to run chimere's parallel part
2026-06-21 20:35:03 INFO Linking EMIS ...
2026-06-21 20:35:03 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens11.nc
2026-06-21 20:35:04 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Sunday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Sunday.s.ens11.nc
2026-06-21 20:35:04 INFO >> Checking links...
2026-06-21 20:35:06 INFO >> All links are good for ENS11  ...
2026-06-21 20:35:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:35:07 INFO Hourly dataset computed and listing created
2026-06-21 20:35:09 INFO Hourly dataset computed
2026-06-21 20:35:09 INFO Creating directories and links for ENS12 to run chimere's parallel part
2026-06-21 20:35:09 INFO Linking EMIS ...
2026-06-21 20:35:10 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens12.nc
2026-06-21 20:35:10 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Sunday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Sunday.s.ens12.nc
2026-06-21 20:35:10 INFO >> Checking links...
2026-06-21 20:35:13 INFO >> All links are good for ENS12  ...
2026-06-21 20:35:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:35:14 INFO Hourly dataset computed and listing created
2026-06-21 20:35:26 INFO Hourly dataset computed
2026-06-21 20:35:26 INFO Creating directories and links for ENS13 to run chimere's parallel part
2026-06-21 20:35:26 INFO Linking EMIS ...
2026-06-21 20:35:26 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens13.nc
2026-06-21 20:35:27 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Sunday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Sunday.s.ens13.nc
2026-06-21 20:35:27 INFO >> Checking links...
2026-06-21 20:35:29 INFO >> All links are good for ENS13  ...
2026-06-21 20:35:29 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:35:30 INFO Hourly dataset computed and listing created
2026-06-21 20:35:45 INFO Hourly dataset computed
2026-06-21 20:35:45 INFO Creating directories and links for ENS14 to run chimere's parallel part
2026-06-21 20:35:45 INFO Linking EMIS ...
2026-06-21 20:35:46 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens14.nc
2026-06-21 20:35:46 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Sunday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Sunday.s.ens14.nc
2026-06-21 20:35:46 INFO >> Checking links...
2026-06-21 20:35:49 INFO >> All links are good for ENS14  ...
2026-06-21 20:35:49 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:35:50 INFO Hourly dataset computed and listing created
2026-06-21 20:36:01 INFO Hourly dataset computed
2026-06-21 20:36:01 INFO Creating directories and links for ENS15 to run chimere's parallel part
2026-06-21 20:36:01 INFO Linking EMIS ...
2026-06-21 20:36:02 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens15.nc
2026-06-21 20:36:02 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Sunday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Sunday.s.ens15.nc
2026-06-21 20:36:02 INFO >> Checking links...
2026-06-21 20:36:05 INFO >> All links are good for ENS15  ...
2026-06-21 20:36:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 20:36:06 INFO Hourly dataset computed and listing created
2026-06-21 20:36:17 INFO Hourly dataset computed
2026-06-21 20:36:17 INFO ---------->>> Running CHIMERE model from 2020-02-15 15:00:00 to 2020-02-16 00:00:00
2026-06-21 20:36:17 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:36:17 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1
2026-06-21 20:36:17 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/end.2020021513_2_ENS1.nc
2026-06-21 20:36:17 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-06-21 20:36:17 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:36:17 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-06-21 20:36:17 INFO Queuing job for member 1...
2026-06-21 20:36:17 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:36:17 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-06-21 20:36:18 INFO Found: ['4938219']
2026-06-21 20:36:23 INFO [TGCC-IRENE] Submitted job with ID:['4938219']
2026-06-21 20:36:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:36:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2
2026-06-21 20:36:23 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/end.2020021513_2_ENS2.nc
2026-06-21 20:36:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-06-21 20:36:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:36:23 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-06-21 20:36:23 INFO Queuing job for member 2...
2026-06-21 20:36:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:36:23 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-06-21 20:36:24 INFO Found: ['4938220']
2026-06-21 20:36:29 INFO [TGCC-IRENE] Submitted job with ID:['4938220']
2026-06-21 20:36:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:36:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3
2026-06-21 20:36:29 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/end.2020021513_2_ENS3.nc
2026-06-21 20:36:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-06-21 20:36:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:36:29 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-06-21 20:36:29 INFO Queuing job for member 3...
2026-06-21 20:36:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:36:29 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-06-21 20:36:30 INFO Found: ['4938222']
2026-06-21 20:36:35 INFO [TGCC-IRENE] Submitted job with ID:['4938222']
2026-06-21 20:36:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:36:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4
2026-06-21 20:36:35 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/end.2020021513_2_ENS4.nc
2026-06-21 20:36:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-06-21 20:36:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:36:35 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-06-21 20:36:35 INFO Queuing job for member 4...
2026-06-21 20:36:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:36:35 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-06-21 20:36:35 INFO Found: ['4938223']
2026-06-21 20:36:40 INFO [TGCC-IRENE] Submitted job with ID:['4938223']
2026-06-21 20:36:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:36:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5
2026-06-21 20:36:40 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/end.2020021513_2_ENS5.nc
2026-06-21 20:36:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-06-21 20:36:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:36:40 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-06-21 20:36:40 INFO Queuing job for member 5...
2026-06-21 20:36:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:36:40 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-06-21 20:36:41 INFO Found: ['4938224']
2026-06-21 20:36:46 INFO [TGCC-IRENE] Submitted job with ID:['4938224']
2026-06-21 20:36:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:36:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6
2026-06-21 20:36:46 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/end.2020021513_2_ENS6.nc
2026-06-21 20:36:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-06-21 20:36:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:36:46 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-06-21 20:36:46 INFO Queuing job for member 6...
2026-06-21 20:36:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:36:46 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-06-21 20:36:47 INFO Found: ['4938225']
2026-06-21 20:36:52 INFO [TGCC-IRENE] Submitted job with ID:['4938225']
2026-06-21 20:36:52 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:36:52 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7
2026-06-21 20:36:52 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/end.2020021513_2_ENS7.nc
2026-06-21 20:36:52 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-06-21 20:36:52 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:36:52 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-06-21 20:36:52 INFO Queuing job for member 7...
2026-06-21 20:36:52 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:36:52 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-06-21 20:36:53 INFO Found: ['4938226']
2026-06-21 20:36:58 INFO [TGCC-IRENE] Submitted job with ID:['4938226']
2026-06-21 20:36:58 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:36:58 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8
2026-06-21 20:36:58 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/end.2020021513_2_ENS8.nc
2026-06-21 20:36:58 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-06-21 20:36:58 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:36:58 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-06-21 20:36:58 INFO Queuing job for member 8...
2026-06-21 20:36:58 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:36:58 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-06-21 20:36:59 INFO Found: ['4938227']
2026-06-21 20:37:04 INFO [TGCC-IRENE] Submitted job with ID:['4938227']
2026-06-21 20:37:04 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:37:04 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9
2026-06-21 20:37:04 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/end.2020021513_2_ENS9.nc
2026-06-21 20:37:04 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-06-21 20:37:04 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:37:04 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-06-21 20:37:04 INFO Queuing job for member 9...
2026-06-21 20:37:04 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:37:04 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-06-21 20:37:05 INFO Found: ['4938229']
2026-06-21 20:37:10 INFO [TGCC-IRENE] Submitted job with ID:['4938229']
2026-06-21 20:37:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:37:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10
2026-06-21 20:37:10 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/end.2020021513_2_ENS10.nc
2026-06-21 20:37:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-06-21 20:37:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:37:10 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-06-21 20:37:10 INFO Queuing job for member 10...
2026-06-21 20:37:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:37:10 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-06-21 20:37:10 INFO Found: ['4938230']
2026-06-21 20:37:15 INFO [TGCC-IRENE] Submitted job with ID:['4938230']
2026-06-21 20:37:15 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:37:15 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11
2026-06-21 20:37:15 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/end.2020021513_2_ENS11.nc
2026-06-21 20:37:15 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-06-21 20:37:15 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:37:15 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-06-21 20:37:15 INFO Queuing job for member 11...
2026-06-21 20:37:15 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:37:15 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-06-21 20:37:16 INFO Found: ['4938231']
2026-06-21 20:37:21 INFO [TGCC-IRENE] Submitted job with ID:['4938231']
2026-06-21 20:37:21 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:37:21 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12
2026-06-21 20:37:21 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/end.2020021513_2_ENS12.nc
2026-06-21 20:37:21 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-06-21 20:37:21 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:37:21 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-06-21 20:37:21 INFO Queuing job for member 12...
2026-06-21 20:37:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:37:21 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-06-21 20:37:22 INFO Found: ['4938232']
2026-06-21 20:37:27 INFO [TGCC-IRENE] Submitted job with ID:['4938232']
2026-06-21 20:37:27 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:37:27 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13
2026-06-21 20:37:27 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/end.2020021513_2_ENS13.nc
2026-06-21 20:37:27 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-06-21 20:37:27 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:37:27 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-06-21 20:37:27 INFO Queuing job for member 13...
2026-06-21 20:37:27 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:37:27 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-06-21 20:37:28 INFO Found: ['4938233']
2026-06-21 20:37:33 INFO [TGCC-IRENE] Submitted job with ID:['4938233']
2026-06-21 20:37:33 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:37:33 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14
2026-06-21 20:37:33 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/end.2020021513_2_ENS14.nc
2026-06-21 20:37:33 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-06-21 20:37:33 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:37:33 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-06-21 20:37:33 INFO Queuing job for member 14...
2026-06-21 20:37:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:37:33 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-06-21 20:37:34 INFO Found: ['4938236']
2026-06-21 20:37:39 INFO [TGCC-IRENE] Submitted job with ID:['4938236']
2026-06-21 20:37:39 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 20:37:39 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15
2026-06-21 20:37:39 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/end.2020021513_2_ENS15.nc
2026-06-21 20:37:39 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-06-21 20:37:39 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 20:37:39 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-06-21 20:37:39 INFO Queuing job for member 15...
2026-06-21 20:37:39 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 20:37:39 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-06-21 20:37:39 INFO Found: ['4938238']
2026-06-21 20:37:44 INFO [TGCC-IRENE] Submitted job with ID:['4938238']
2026-06-21 20:37:44 INFO Checking job status ...
2026-06-21 20:37:44 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:37:44 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:37:45 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:37:45 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:37:45 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:37:45 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:37:45 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:37:45 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:37:45 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:37:45 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:37:45 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:37:45 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:37:45 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:37:45 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:37:45 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:37:45 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:38:00 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:38:00 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:38:00 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:38:00 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:38:00 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:38:00 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:38:00 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:38:00 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:38:00 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:38:00 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:38:00 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:38:00 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:38:00 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:38:00 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:38:00 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:38:00 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:38:15 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:38:15 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:38:15 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:38:15 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:38:15 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:38:15 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:38:15 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:38:15 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:38:15 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:38:15 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:38:15 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:38:15 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:38:15 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:38:15 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:38:15 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:38:15 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:38:30 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:38:30 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:38:31 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:38:31 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:38:31 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:38:31 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:38:31 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:38:31 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:38:31 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:38:31 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:38:31 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:38:31 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:38:31 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:38:31 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:38:31 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:38:31 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:38:46 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:38:46 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:38:46 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:38:46 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:38:46 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:38:46 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:38:46 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:38:46 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:38:46 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:38:46 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:38:46 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:38:46 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:38:46 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:38:46 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:38:46 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:38:46 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:39:01 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:39:01 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:39:01 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:39:01 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:39:01 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:39:01 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:39:01 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:39:01 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:39:01 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:39:01 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:39:01 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:39:01 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:39:01 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:39:01 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:39:01 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:39:01 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:39:16 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:39:16 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:39:16 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:39:16 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:39:16 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:39:16 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:39:17 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:39:17 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:39:17 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:39:17 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:39:17 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:39:17 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:39:17 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:39:17 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:39:17 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:39:17 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:39:32 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:39:32 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:39:32 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:39:32 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:39:32 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:39:32 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:39:32 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:39:32 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:39:32 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:39:32 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:39:32 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:39:32 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:39:32 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:39:32 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:39:32 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:39:32 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:39:47 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:39:47 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:39:47 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:39:47 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:39:47 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:39:47 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:39:47 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:39:47 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:39:47 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:39:47 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:39:47 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:39:47 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:39:47 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:39:47 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:39:47 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:39:47 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:40:02 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:40:02 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:40:02 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:40:02 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:40:02 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:40:02 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:40:02 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:40:02 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:40:03 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:40:03 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:40:03 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:40:03 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:40:03 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:40:03 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:40:03 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:40:03 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:40:18 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:40:18 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:40:18 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:40:18 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:40:18 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:40:18 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:40:18 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:40:18 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:40:18 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:40:18 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:40:18 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:40:18 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:40:18 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:40:18 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:40:18 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:40:18 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:40:33 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:40:33 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:40:33 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:40:33 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:40:33 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:40:33 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:40:33 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:40:33 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:40:33 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:40:33 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:40:33 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:40:33 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:40:33 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:40:33 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:40:33 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:40:33 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:40:48 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:40:48 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:40:48 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:40:48 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:40:48 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:40:48 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:40:48 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:40:49 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:40:49 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:40:49 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:40:49 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:40:49 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:40:49 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:40:49 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:40:49 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:40:49 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:41:04 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:41:04 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:41:04 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:41:04 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:41:04 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:41:04 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:41:04 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:41:04 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:41:04 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:41:04 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:41:04 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:41:04 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:41:04 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:41:04 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:41:04 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:41:04 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:41:19 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:41:19 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:41:19 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:41:19 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:41:19 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:41:19 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:41:19 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:41:19 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:41:19 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:41:19 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:41:19 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:41:19 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:41:19 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:41:19 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:41:19 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:41:19 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:41:34 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:41:34 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:41:34 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:41:34 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:41:34 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:41:35 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:41:35 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:41:35 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:41:35 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:41:35 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:41:35 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:41:35 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:41:35 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:41:35 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:41:35 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:41:35 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:41:50 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:41:50 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:41:50 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:41:50 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:41:50 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:41:50 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:41:50 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:41:50 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:41:50 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:41:50 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:41:50 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:41:50 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:41:50 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:41:50 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:41:50 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:41:50 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:42:05 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:42:05 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:42:05 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:42:05 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:42:05 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:42:05 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:42:05 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:42:05 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:42:06 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:42:06 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:42:06 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:42:06 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:42:06 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:42:06 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:42:06 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:42:06 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:42:21 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:42:21 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:42:21 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:42:21 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:42:21 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:42:21 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:42:21 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:42:21 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:42:21 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:42:21 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:42:21 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:42:21 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:42:21 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:42:21 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:42:21 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:42:21 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:42:36 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:42:36 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:42:36 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:42:36 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:42:36 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:42:36 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:42:36 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:42:36 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:42:36 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:42:36 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:42:36 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:42:36 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:42:36 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:42:36 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:42:36 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:42:36 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:42:51 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:42:51 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:42:51 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:42:51 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:42:51 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:42:51 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:42:51 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:42:51 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:42:51 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:42:51 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:42:51 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:42:51 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:42:52 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:42:52 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:42:52 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:42:52 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:43:07 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:43:07 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:43:07 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:43:07 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:43:07 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:43:07 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:43:07 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:43:07 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:43:07 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:43:07 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:43:07 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:43:07 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:43:07 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:43:07 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:43:07 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:43:07 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:43:22 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:43:22 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:43:22 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:43:22 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:43:22 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:43:22 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:43:22 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:43:22 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:43:22 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:43:22 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:43:22 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:43:22 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:43:22 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:43:22 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:43:22 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:43:22 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:43:37 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:43:37 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:43:37 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:43:37 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:43:37 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:43:37 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:43:37 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:43:37 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:43:37 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:43:37 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:43:37 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:43:37 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:43:37 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:43:37 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:43:37 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:43:37 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:43:53 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:43:53 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:43:53 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:43:53 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:43:53 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:43:53 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:43:53 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:43:53 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:43:53 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:43:53 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:43:53 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:43:53 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:43:53 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:43:53 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:43:53 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:43:53 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:44:08 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:44:08 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:44:08 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:44:08 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:44:08 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:44:08 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:44:08 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:44:08 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:44:08 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:44:08 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:44:08 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:44:08 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:44:08 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:44:08 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:44:08 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:44:08 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:44:23 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:44:23 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:44:23 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:44:23 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:44:23 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:44:23 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:44:23 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:44:23 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:44:23 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:44:23 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:44:23 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:44:23 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:44:23 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:44:23 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:44:23 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:44:23 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:44:39 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:44:39 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:44:39 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:44:39 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:44:39 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:44:39 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:44:39 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:44:39 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:44:39 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:44:39 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:44:39 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:44:39 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:44:39 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:44:39 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:44:39 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:44:39 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:44:54 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:44:54 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:44:54 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:44:54 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:44:54 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:44:54 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:44:54 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:44:54 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:44:54 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:44:54 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:44:54 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:44:54 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:44:54 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:44:54 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:44:54 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:44:54 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:45:09 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:45:09 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:45:09 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:45:09 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:45:09 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:45:09 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:45:09 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:45:09 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:45:09 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:45:09 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:45:09 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:45:09 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:45:09 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:45:09 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:45:09 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:45:09 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:45:24 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:45:25 INFO None 4938220: status RUNNING/PENDING
2026-06-21 20:45:25 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:45:25 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:45:25 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:45:25 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:45:25 INFO None 4938226: status RUNNING/PENDING
2026-06-21 20:45:25 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:45:25 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:45:25 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:45:25 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:45:25 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:45:25 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:45:25 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:45:25 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:45:25 INFO Jobs still running: ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:45:40 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:45:40 INFO None 4938220: status FINISHED
2026-06-21 20:45:40 INFO None 4938222: status RUNNING/PENDING
2026-06-21 20:45:40 INFO None 4938223: status RUNNING/PENDING
2026-06-21 20:45:40 INFO None 4938224: status RUNNING/PENDING
2026-06-21 20:45:40 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:45:40 INFO None 4938226: status FINISHED
2026-06-21 20:45:40 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:45:40 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:45:40 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:45:40 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:45:40 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:45:40 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:45:40 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:45:40 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:45:40 INFO Jobs still running: ['4938219', '4938222', '4938223', '4938224', '4938225', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:45:55 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:45:55 INFO None 4938220: status FINISHED
2026-06-21 20:45:55 INFO None 4938222: status FINISHED
2026-06-21 20:45:55 INFO None 4938223: status FINISHED
2026-06-21 20:45:55 INFO None 4938224: status FINISHED
2026-06-21 20:45:55 INFO None 4938225: status RUNNING/PENDING
2026-06-21 20:45:55 INFO None 4938226: status FINISHED
2026-06-21 20:45:55 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:45:55 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:45:55 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:45:55 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:45:55 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:45:55 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:45:55 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:45:55 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:45:55 INFO Jobs still running: ['4938219', '4938225', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:46:10 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:46:10 INFO None 4938220: status FINISHED
2026-06-21 20:46:10 INFO None 4938222: status FINISHED
2026-06-21 20:46:10 INFO None 4938223: status FINISHED
2026-06-21 20:46:10 INFO None 4938224: status FINISHED
2026-06-21 20:46:10 INFO None 4938225: status FINISHED
2026-06-21 20:46:10 INFO None 4938226: status FINISHED
2026-06-21 20:46:11 INFO None 4938227: status RUNNING/PENDING
2026-06-21 20:46:11 INFO None 4938229: status RUNNING/PENDING
2026-06-21 20:46:11 INFO None 4938230: status RUNNING/PENDING
2026-06-21 20:46:11 INFO None 4938231: status RUNNING/PENDING
2026-06-21 20:46:11 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:46:12 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:46:12 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:46:12 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:46:12 INFO Jobs still running: ['4938219', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:46:27 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:46:27 INFO None 4938220: status FINISHED
2026-06-21 20:46:27 INFO None 4938222: status FINISHED
2026-06-21 20:46:27 INFO None 4938223: status FINISHED
2026-06-21 20:46:27 INFO None 4938224: status FINISHED
2026-06-21 20:46:27 INFO None 4938225: status FINISHED
2026-06-21 20:46:27 INFO None 4938226: status FINISHED
2026-06-21 20:46:27 INFO None 4938227: status FINISHED
2026-06-21 20:46:27 INFO None 4938229: status FINISHED
2026-06-21 20:46:27 INFO None 4938230: status FINISHED
2026-06-21 20:46:27 INFO None 4938231: status FINISHED
2026-06-21 20:46:27 INFO None 4938232: status RUNNING/PENDING
2026-06-21 20:46:27 INFO None 4938233: status RUNNING/PENDING
2026-06-21 20:46:27 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:46:27 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:46:27 INFO Jobs still running: ['4938219', '4938232', '4938233', '4938236', '4938238']. Waiting...
2026-06-21 20:46:42 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:46:42 INFO None 4938220: status FINISHED
2026-06-21 20:46:42 INFO None 4938222: status FINISHED
2026-06-21 20:46:42 INFO None 4938223: status FINISHED
2026-06-21 20:46:42 INFO None 4938224: status FINISHED
2026-06-21 20:46:42 INFO None 4938225: status FINISHED
2026-06-21 20:46:42 INFO None 4938226: status FINISHED
2026-06-21 20:46:42 INFO None 4938227: status FINISHED
2026-06-21 20:46:42 INFO None 4938229: status FINISHED
2026-06-21 20:46:42 INFO None 4938230: status FINISHED
2026-06-21 20:46:42 INFO None 4938231: status FINISHED
2026-06-21 20:46:42 INFO None 4938232: status FINISHED
2026-06-21 20:46:42 INFO None 4938233: status FINISHED
2026-06-21 20:46:42 INFO None 4938236: status RUNNING/PENDING
2026-06-21 20:46:42 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:46:42 INFO Jobs still running: ['4938219', '4938236', '4938238']. Waiting...
2026-06-21 20:46:57 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:46:57 INFO None 4938220: status FINISHED
2026-06-21 20:46:57 INFO None 4938222: status FINISHED
2026-06-21 20:46:57 INFO None 4938223: status FINISHED
2026-06-21 20:46:57 INFO None 4938224: status FINISHED
2026-06-21 20:46:57 INFO None 4938225: status FINISHED
2026-06-21 20:46:57 INFO None 4938226: status FINISHED
2026-06-21 20:46:57 INFO None 4938227: status FINISHED
2026-06-21 20:46:57 INFO None 4938229: status FINISHED
2026-06-21 20:46:57 INFO None 4938230: status FINISHED
2026-06-21 20:46:57 INFO None 4938231: status FINISHED
2026-06-21 20:46:58 INFO None 4938232: status FINISHED
2026-06-21 20:46:58 INFO None 4938233: status FINISHED
2026-06-21 20:46:58 INFO None 4938236: status FINISHED
2026-06-21 20:46:58 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:46:58 INFO Jobs still running: ['4938219', '4938238']. Waiting...
2026-06-21 20:47:13 INFO None 4938219: status RUNNING/PENDING
2026-06-21 20:47:13 INFO None 4938220: status FINISHED
2026-06-21 20:47:13 INFO None 4938222: status FINISHED
2026-06-21 20:47:13 INFO None 4938223: status FINISHED
2026-06-21 20:47:13 INFO None 4938224: status FINISHED
2026-06-21 20:47:13 INFO None 4938225: status FINISHED
2026-06-21 20:47:13 INFO None 4938226: status FINISHED
2026-06-21 20:47:13 INFO None 4938227: status FINISHED
2026-06-21 20:47:13 INFO None 4938229: status FINISHED
2026-06-21 20:47:13 INFO None 4938230: status FINISHED
2026-06-21 20:47:13 INFO None 4938231: status FINISHED
2026-06-21 20:47:13 INFO None 4938232: status FINISHED
2026-06-21 20:47:13 INFO None 4938233: status FINISHED
2026-06-21 20:47:13 INFO None 4938236: status FINISHED
2026-06-21 20:47:13 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:47:13 INFO Jobs still running: ['4938219', '4938238']. Waiting...
2026-06-21 20:47:28 INFO None 4938219: status FINISHED
2026-06-21 20:47:28 INFO None 4938220: status FINISHED
2026-06-21 20:47:28 INFO None 4938222: status FINISHED
2026-06-21 20:47:28 INFO None 4938223: status FINISHED
2026-06-21 20:47:28 INFO None 4938224: status FINISHED
2026-06-21 20:47:28 INFO None 4938225: status FINISHED
2026-06-21 20:47:28 INFO None 4938226: status FINISHED
2026-06-21 20:47:28 INFO None 4938227: status FINISHED
2026-06-21 20:47:28 INFO None 4938229: status FINISHED
2026-06-21 20:47:28 INFO None 4938230: status FINISHED
2026-06-21 20:47:28 INFO None 4938231: status FINISHED
2026-06-21 20:47:28 INFO None 4938232: status FINISHED
2026-06-21 20:47:28 INFO None 4938233: status FINISHED
2026-06-21 20:47:28 INFO None 4938236: status FINISHED
2026-06-21 20:47:28 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:47:28 INFO Jobs still running: ['4938238']. Waiting...
2026-06-21 20:47:43 INFO None 4938219: status FINISHED
2026-06-21 20:47:43 INFO None 4938220: status FINISHED
2026-06-21 20:47:43 INFO None 4938222: status FINISHED
2026-06-21 20:47:43 INFO None 4938223: status FINISHED
2026-06-21 20:47:43 INFO None 4938224: status FINISHED
2026-06-21 20:47:43 INFO None 4938225: status FINISHED
2026-06-21 20:47:43 INFO None 4938226: status FINISHED
2026-06-21 20:47:43 INFO None 4938227: status FINISHED
2026-06-21 20:47:43 INFO None 4938229: status FINISHED
2026-06-21 20:47:43 INFO None 4938230: status FINISHED
2026-06-21 20:47:43 INFO None 4938231: status FINISHED
2026-06-21 20:47:43 INFO None 4938232: status FINISHED
2026-06-21 20:47:43 INFO None 4938233: status FINISHED
2026-06-21 20:47:43 INFO None 4938236: status FINISHED
2026-06-21 20:47:43 INFO None 4938238: status RUNNING/PENDING
2026-06-21 20:47:43 INFO Jobs still running: ['4938238']. Waiting...
2026-06-21 20:47:58 INFO None 4938219: status FINISHED
2026-06-21 20:47:58 INFO None 4938220: status FINISHED
2026-06-21 20:47:58 INFO None 4938222: status FINISHED
2026-06-21 20:47:59 INFO None 4938223: status FINISHED
2026-06-21 20:47:59 INFO None 4938224: status FINISHED
2026-06-21 20:47:59 INFO None 4938225: status FINISHED
2026-06-21 20:47:59 INFO None 4938226: status FINISHED
2026-06-21 20:47:59 INFO None 4938227: status FINISHED
2026-06-21 20:47:59 INFO None 4938229: status FINISHED
2026-06-21 20:47:59 INFO None 4938230: status FINISHED
2026-06-21 20:47:59 INFO None 4938231: status FINISHED
2026-06-21 20:47:59 INFO None 4938232: status FINISHED
2026-06-21 20:47:59 INFO None 4938233: status FINISHED
2026-06-21 20:47:59 INFO None 4938236: status FINISHED
2026-06-21 20:47:59 INFO None 4938238: status FINISHED
2026-06-21 20:47:59 INFO Jobs ['4938219', '4938220', '4938222', '4938223', '4938224', '4938225', '4938226', '4938227', '4938229', '4938230', '4938231', '4938232', '4938233', '4938236', '4938238'] have finished
2026-06-21 20:47:59 INFO Checking restart files were created ...
2026-06-21 20:47:59 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/end.2020021515_9_ENS1.nc(3339660275 bytes)
2026-06-21 20:47:59 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/end.2020021515_9_ENS2.nc(3339660275 bytes)
2026-06-21 20:47:59 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/end.2020021515_9_ENS3.nc(3339660275 bytes)
2026-06-21 20:47:59 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/end.2020021515_9_ENS4.nc(3339660275 bytes)
2026-06-21 20:47:59 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/end.2020021515_9_ENS5.nc(3339660275 bytes)
2026-06-21 20:47:59 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/end.2020021515_9_ENS6.nc(3339660275 bytes)
2026-06-21 20:47:59 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/end.2020021515_9_ENS7.nc(3339660275 bytes)
2026-06-21 20:47:59 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/end.2020021515_9_ENS8.nc(3339660275 bytes)
2026-06-21 20:47:59 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/end.2020021515_9_ENS9.nc(3339660275 bytes)
2026-06-21 20:47:59 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/end.2020021515_9_ENS10.nc(3339660275 bytes)
2026-06-21 20:47:59 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/end.2020021515_9_ENS11.nc(3339660275 bytes)
2026-06-21 20:47:59 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/end.2020021515_9_ENS12.nc(3339660275 bytes)
2026-06-21 20:47:59 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/end.2020021515_9_ENS13.nc(3339660275 bytes)
2026-06-21 20:47:59 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/end.2020021515_9_ENS14.nc(3339660275 bytes)
2026-06-21 20:47:59 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/end.2020021515_9_ENS15.nc(3339660275 bytes)
2026-06-21 20:47:59 INFO  Run_model() completed successfully.
2026-06-21 20:47:59 INFO [TIME] after_model_set_simulated_time current_time=2020-02-15 15:00:00 simulated_time=2020-02-16 00:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 20:47:59 INFO [TIME] gregorian_conversion simulated_time=2020-02-16 00:00:00 days=153082 seconds=0
2026-06-21 20:47:59 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-06-21 20:47:59 INFO [TIME] increment current_time 2020-02-15 15:00:00 -> 2020-02-16 00:00:00
2026-06-21 20:47:59 INFO [TIME] after_increment_before_assimilation current_time=2020-02-16 00:00:00 simulated_time=2020-02-16 00:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 20:47:59 INFO ---------->>> Running process_satellite_data()
2026-06-21 20:47:59 INFO [DART] No satellite data found, skipping assimilation
2026-06-21 20:47:59 INFO after_assimilation() skipped
2026-06-21 20:47:59 INFO Next run starts from 2020-02-16 00:00:00
2026-06-21 20:47:59 INFO Cycle is DONE; starting a new loop!
2026-06-21 20:47:59 INFO [TIME] step_end current_time=2020-02-16 00:00:00 simulated_time=2020-02-16 00:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 20:47:59 INFO [PIPELINE] ---- TIME LOOP END ----
+ exit 0
