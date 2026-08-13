+ SCRIPT_PID=2757416
+ /bin/bash -x /tmp/tmp.gJForznf79
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
2026-05-27 17:39:38 INFO 
███    ███ ██ ███    ███ ███████ ███████ ██
████  ████ ██ ████  ████ ██      ██      ██
██ ████ ██ ██ ██ ████ ██ █████   ███████ ██
██  ██  ██ ██ ██  ██  ██ ██           ██ ██
██      ██ ██ ██      ██ ███████ ███████ ██



2026-05-27 17:39:38 INFO [PIPELINE] =======================================
2026-05-27 17:39:38 INFO [PIPELINE] Starting chimere–DART orchestrator
2026-05-27 17:39:38 INFO [PIPELINE] Config file: config/config_irene_IM.yaml
2026-05-27 17:39:38 INFO [PIPELINE] Run dir: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart
2026-05-27 17:39:38 INFO [PIPELINE] Log file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/mimesi_orchestrator_logs/chimere_DART_20260527_173938.log
2026-05-27 17:39:38 INFO [PIPELINE] =======================================
2026-05-27 17:39:38 INFO Running assimilation with model_type=ModelType.CHIMERE
2026-05-27 17:39:38 INFO Using scheduler=Scheduler.SLURM, queue=rome
2026-05-27 17:39:38 INFO [STEP] ---- TIME LOOP START ----
2026-05-27 17:39:38 INFO [TIME] step_start current_time=2020-02-06 14:00:00 simulated_time=None dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-05-27 17:39:38 INFO [TIME] window start=2020-02-06 14:00:00 end=2020-02-07 00:00:00 run_hours=10 has_assimilation=False
2026-05-27 17:39:38 INFO Creating directories and links for ENS1 to run chimere's parallel part
2026-05-27 17:39:38 INFO Linking EMIS ...
2026-05-27 17:39:38 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/nocy_0615_15m_high_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc
2026-05-27 17:39:38 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/nocy_0615_15m_high_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc -> /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_tripled2/EMIS.EUROCOMEX3.02.Friday.s.ens01.nc
2026-05-27 17:39:38 INFO >> Checking links...
2026-05-27 17:39:39 INFO >> All links are good for ENS1  ...
2026-05-27 17:39:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-05-27 17:39:43 INFO Hourly dataset computed and listing created
2026-05-27 17:39:59 INFO Hourly dataset computed
2026-05-27 17:39:59 INFO Creating directories and links for ENS2 to run chimere's parallel part
2026-05-27 17:39:59 INFO Linking EMIS ...
2026-05-27 17:39:59 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/nocy_0615_15m_high_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc
2026-05-27 17:39:59 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/nocy_0615_15m_high_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc -> /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_tripled2/EMIS.EUROCOMEX3.02.Friday.s.ens02.nc
2026-05-27 17:39:59 INFO >> Checking links...
2026-05-27 17:39:59 INFO >> All links are good for ENS2  ...
2026-05-27 17:39:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-05-27 17:40:01 INFO Hourly dataset computed and listing created
2026-05-27 17:40:15 INFO Hourly d