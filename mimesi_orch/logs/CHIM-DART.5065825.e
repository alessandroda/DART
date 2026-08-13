+ SCRIPT_PID=404285
+ /bin/bash -x /tmp/tmp.SNjC4E7429
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
2026-07-02 12:53:21 INFO 
███    ███ ██ ███    ███ ███████ ███████ ██
████  ████ ██ ████  ████ ██      ██      ██
██ ████ ██ ██ ██ ████ ██ █████   ███████ ██
██  ██  ██ ██ ██  ██  ██ ██           ██ ██
██      ██ ██ ██      ██ ███████ ███████ ██



2026-07-02 12:53:21 INFO [PIPELINE] =======================================
2026-07-02 12:53:21 INFO [PIPELINE] Starting chimere–DART orchestrator
2026-07-02 12:53:21 INFO [PIPELINE] Config file: config/config_irene_IM.yaml
2026-07-02 12:53:21 INFO [PIPELINE] Run dir: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart
2026-07-02 12:53:21 INFO [PIPELINE] Log file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/mimesi_orchestrator_logs/chimere_DART_20260702_125321.log
2026-07-02 12:53:21 INFO [PIPELINE] =======================================
2026-07-02 12:53:21 INFO Running assimilation with model_type=ModelType.CHIMERE
2026-07-02 12:53:21 INFO Using scheduler=Scheduler.SLURM, queue=rome
2026-07-02 12:53:21 INFO [STEP] ---- TIME LOOP START ----
2026-07-02 12:53:21 INFO [TIME] step_start current_time=2020-02-06 00:00:00 simulated_time=None dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 12:53:21 INFO [TIME] window start=2020-02-06 00:00:00 end=2020-02-06 01:00:00 run_hours=1 has_assimilation=False
2026-07-02 12:53:21 INFO Creating directories and links for ENS1 to run chimere's parallel part
2026-07-02 12:53:21 INFO Copying EMIS ...
2026-07-02 12:53:21 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens01.nc
2026-07-02 12:53:21 INFO Linking END ...
2026-07-02 12:53:21 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020500_24_ENS1.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-02 12:53:21 INFO >> Checking links...
2026-07-02 12:53:21 INFO >> All links are good for ENS1  ...
2026-07-02 12:53:21 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 12:53:29 INFO Hourly dataset computed and listing created
2026-07-02 12:53:32 INFO Hourly dataset computed
2026-07-02 12:53:32 INFO Creating directories and links for ENS2 to run chimere's parallel part
2026-07-02 12:53:32 INFO Copying EMIS ...
2026-07-02 12:53:32 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens02.nc
2026-07-02 12:53:32 INFO Linking END ...
2026-07-02 12:53:32 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020500_24_ENS2.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-02 12:53:32 INFO >> Checking links...
2026-07-02 12:53:32 INFO >> All links are good for ENS2  ...
2026-07-02 12:53:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 12:53:33 INFO Hourly dataset computed and listing created
2026-07-02 12:53:34 INFO Hourly dataset computed
2026-07-02 12:53:34 INFO Creating directories and links for ENS3 to run chimere's parallel part
2026-07-02 12:53:34 INFO Copying EMIS ...
2026-07-02 12:53:34 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens03.nc
2026-07-02 12:53:34 INFO Linking END ...
2026-07-02 12:53:34 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020500_24_ENS3.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-02 12:53:34 INFO >> Checking links...
2026-07-02 12:53:34 INFO >> All links are good for ENS3  ...
2026-07-02 12:53:34 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 12:53:35 INFO Hourly dataset computed and listing created
2026-07-02 12:53:35 INFO Hourly dataset computed
2026-07-02 12:53:36 INFO Creating directories and links for ENS4 to run chimere's parallel part
2026-07-02 12:53:36 INFO Copying EMIS ...
2026-07-02 12:53:36 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens04.nc
2026-07-02 12:53:36 INFO Linking END ...
2026-07-02 12:53:36 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020500_24_ENS4.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-02 12:53:36 INFO >> Checking links...
2026-07-02 12:53:36 INFO >> All links are good for ENS4  ...
2026-07-02 12:53:36 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 12:53:37 INFO Hourly dataset computed and listing created
2026-07-02 12:53:37 INFO Hourly dataset computed
2026-07-02 12:53:37 INFO Creating directories and links for ENS5 to run chimere's parallel part
2026-07-02 12:53:37 INFO Copying EMIS ...
2026-07-02 12:53:38 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens05.nc
2026-07-02 12:53:38 INFO Linking END ...
2026-07-02 12:53:38 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020500_24_ENS5.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-02 12:53:38 INFO >> Checking links...
2026-07-02 12:53:38 INFO >> All links are good for ENS5  ...
2026-07-02 12:53:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 12:53:39 INFO Hourly dataset computed and listing created
2026-07-02 12:53:39 INFO Hourly dataset computed
2026-07-02 12:53:39 INFO Creating directories and links for ENS6 to run chimere's parallel part
2026-07-02 12:53:39 INFO Copying EMIS ...
2026-07-02 12:53:40 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens06.nc
2026-07-02 12:53:40 INFO Linking END ...
2026-07-02 12:53:40 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020500_24_ENS6.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-02 12:53:40 INFO >> Checking links...
2026-07-02 12:53:40 INFO >> All links are good for ENS6  ...
2026-07-02 12:53:40 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 12:53:41 INFO Hourly dataset computed and listing created
2026-07-02 12:53:41 INFO Hourly dataset computed
2026-07-02 12:53:41 INFO Creating directories and links for ENS7 to run chimere's parallel part
2026-07-02 12:53:41 INFO Copying EMIS ...
2026-07-02 12:53:42 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens07.nc
2026-07-02 12:53:42 INFO Linking END ...
2026-07-02 12:53:42 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020500_24_ENS7.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-02 12:53:42 INFO >> Checking links...
2026-07-02 12:53:42 INFO >> All links are good for ENS7  ...
2026-07-02 12:53:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 12:53:43 INFO Hourly dataset computed and listing created
2026-07-02 12:53:43 INFO Hourly dataset computed
2026-07-02 12:53:43 INFO Creating directories and links for ENS8 to run chimere's parallel part
2026-07-02 12:53:43 INFO Copying EMIS ...
2026-07-02 12:53:44 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens08.nc
2026-07-02 12:53:44 INFO Linking END ...
2026-07-02 12:53:44 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020500_24_ENS8.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-02 12:53:44 INFO >> Checking links...
2026-07-02 12:53:44 INFO >> All links are good for ENS8  ...
2026-07-02 12:53:44 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 12:53:45 INFO Hourly dataset computed and listing created
2026-07-02 12:53:47 INFO Hourly dataset computed
2026-07-02 12:53:47 INFO Creating directories and links for ENS9 to run chimere's parallel part
2026-07-02 12:53:47 INFO Copying EMIS ...
2026-07-02 12:53:48 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens09.nc
2026-07-02 12:53:48 INFO Linking END ...
2026-07-02 12:53:48 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020500_24_ENS9.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-02 12:53:48 INFO >> Checking links...
2026-07-02 12:53:48 INFO >> All links are good for ENS9  ...
2026-07-02 12:53:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 12:53:49 INFO Hourly dataset computed and listing created
2026-07-02 12:53:51 INFO Hourly dataset computed
2026-07-02 12:53:51 INFO Creating directories and links for ENS10 to run chimere's parallel part
2026-07-02 12:53:51 INFO Copying EMIS ...
2026-07-02 12:53:51 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens10.nc
2026-07-02 12:53:51 INFO Linking END ...
2026-07-02 12:53:51 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020500_24_ENS10.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-02 12:53:51 INFO >> Checking links...
2026-07-02 12:53:51 INFO >> All links are good for ENS10  ...
2026-07-02 12:53:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 12:53:52 INFO Hourly dataset computed and listing created
2026-07-02 12:53:52 INFO Hourly dataset computed
2026-07-02 12:53:53 INFO Creating directories and links for ENS11 to run chimere's parallel part
2026-07-02 12:53:53 INFO Copying EMIS ...
2026-07-02 12:53:53 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens11.nc
2026-07-02 12:53:53 INFO Linking END ...
2026-07-02 12:53:53 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020500_24_ENS11.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-02 12:53:53 INFO >> Checking links...
2026-07-02 12:53:53 INFO >> All links are good for ENS11  ...
2026-07-02 12:53:53 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 12:53:54 INFO Hourly dataset computed and listing created
2026-07-02 12:53:54 INFO Hourly dataset computed
2026-07-02 12:53:54 INFO Creating directories and links for ENS12 to run chimere's parallel part
2026-07-02 12:53:54 INFO Copying EMIS ...
2026-07-02 12:53:55 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens12.nc
2026-07-02 12:53:55 INFO Linking END ...
2026-07-02 12:53:55 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020500_24_ENS12.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-02 12:53:55 INFO >> Checking links...
2026-07-02 12:53:55 INFO >> All links are good for ENS12  ...
2026-07-02 12:53:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 12:53:56 INFO Hourly dataset computed and listing created
2026-07-02 12:53:56 INFO Hourly dataset computed
2026-07-02 12:53:56 INFO Creating directories and links for ENS13 to run chimere's parallel part
2026-07-02 12:53:56 INFO Copying EMIS ...
2026-07-02 12:53:57 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens13.nc
2026-07-02 12:53:57 INFO Linking END ...
2026-07-02 12:53:57 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020500_24_ENS13.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-02 12:53:57 INFO >> Checking links...
2026-07-02 12:53:57 INFO >> All links are good for ENS13  ...
2026-07-02 12:53:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 12:53:58 INFO Hourly dataset computed and listing created
2026-07-02 12:53:58 INFO Hourly dataset computed
2026-07-02 12:53:58 INFO Creating directories and links for ENS14 to run chimere's parallel part
2026-07-02 12:53:58 INFO Copying EMIS ...
2026-07-02 12:53:59 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens14.nc
2026-07-02 12:53:59 INFO Linking END ...
2026-07-02 12:53:59 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020500_24_ENS14.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-02 12:53:59 INFO >> Checking links...
2026-07-02 12:53:59 INFO >> All links are good for ENS14  ...
2026-07-02 12:53:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 12:54:00 INFO Hourly dataset computed and listing created
2026-07-02 12:54:00 INFO Hourly dataset computed
2026-07-02 12:54:00 INFO Creating directories and links for ENS15 to run chimere's parallel part
2026-07-02 12:54:00 INFO Copying EMIS ...
2026-07-02 12:54:01 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens15.nc
2026-07-02 12:54:01 INFO Linking END ...
2026-07-02 12:54:01 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020500_24_ENS15.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-02 12:54:01 INFO >> Checking links...
2026-07-02 12:54:01 INFO >> All links are good for ENS15  ...
2026-07-02 12:54:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 12:54:01 INFO Hourly dataset computed and listing created
2026-07-02 12:54:02 INFO Hourly dataset computed
2026-07-02 12:54:02 INFO ---------->>> Running CHIMERE model from 2020-02-06 00:00:00 to 2020-02-06 01:00:00
2026-07-02 12:54:02 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 12:54:02 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1
2026-07-02 12:54:02 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020500_24_ENS1.nc
2026-07-02 12:54:02 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-02 12:54:02 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 12:54:02 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-02 12:54:02 INFO Queuing job for member 1...
2026-07-02 12:54:02 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 12:54:02 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-02 12:54:03 INFO Found: ['5065830']
2026-07-02 12:54:08 INFO [TGCC-IRENE] Submitted job with ID:['5065830']
2026-07-02 12:54:08 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 12:54:08 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2
2026-07-02 12:54:08 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020500_24_ENS2.nc
2026-07-02 12:54:08 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-02 12:54:08 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 12:54:08 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-02 12:54:08 INFO Queuing job for member 2...
2026-07-02 12:54:08 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 12:54:08 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-02 12:54:09 INFO Found: ['5065832']
2026-07-02 12:54:14 INFO [TGCC-IRENE] Submitted job with ID:['5065832']
2026-07-02 12:54:14 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 12:54:14 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3
2026-07-02 12:54:14 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020500_24_ENS3.nc
2026-07-02 12:54:14 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-02 12:54:14 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 12:54:14 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-02 12:54:14 INFO Queuing job for member 3...
2026-07-02 12:54:14 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 12:54:14 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-02 12:54:15 INFO Found: ['5065833']
2026-07-02 12:54:20 INFO [TGCC-IRENE] Submitted job with ID:['5065833']
2026-07-02 12:54:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 12:54:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4
2026-07-02 12:54:20 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020500_24_ENS4.nc
2026-07-02 12:54:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-02 12:54:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 12:54:20 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-02 12:54:20 INFO Queuing job for member 4...
2026-07-02 12:54:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 12:54:20 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-02 12:54:21 INFO Found: ['5065835']
2026-07-02 12:54:26 INFO [TGCC-IRENE] Submitted job with ID:['5065835']
2026-07-02 12:54:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 12:54:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5
2026-07-02 12:54:26 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020500_24_ENS5.nc
2026-07-02 12:54:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-02 12:54:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 12:54:26 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-02 12:54:26 INFO Queuing job for member 5...
2026-07-02 12:54:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 12:54:26 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-02 12:54:26 INFO Found: ['5065836']
2026-07-02 12:54:31 INFO [TGCC-IRENE] Submitted job with ID:['5065836']
2026-07-02 12:54:31 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 12:54:31 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6
2026-07-02 12:54:31 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020500_24_ENS6.nc
2026-07-02 12:54:31 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-02 12:54:31 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 12:54:31 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-02 12:54:31 INFO Queuing job for member 6...
2026-07-02 12:54:31 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 12:54:31 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-02 12:54:32 INFO Found: ['5065838']
2026-07-02 12:54:37 INFO [TGCC-IRENE] Submitted job with ID:['5065838']
2026-07-02 12:54:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 12:54:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7
2026-07-02 12:54:37 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020500_24_ENS7.nc
2026-07-02 12:54:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-02 12:54:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 12:54:37 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-02 12:54:37 INFO Queuing job for member 7...
2026-07-02 12:54:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 12:54:37 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-02 12:54:38 INFO Found: ['5065839']
2026-07-02 12:54:43 INFO [TGCC-IRENE] Submitted job with ID:['5065839']
2026-07-02 12:54:43 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 12:54:43 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8
2026-07-02 12:54:43 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020500_24_ENS8.nc
2026-07-02 12:54:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-02 12:54:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 12:54:43 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-02 12:54:43 INFO Queuing job for member 8...
2026-07-02 12:54:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 12:54:43 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-02 12:54:44 INFO Found: ['5065840']
2026-07-02 12:54:49 INFO [TGCC-IRENE] Submitted job with ID:['5065840']
2026-07-02 12:54:49 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 12:54:49 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9
2026-07-02 12:54:49 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020500_24_ENS9.nc
2026-07-02 12:54:49 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-02 12:54:49 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 12:54:49 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-02 12:54:49 INFO Queuing job for member 9...
2026-07-02 12:54:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 12:54:49 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-02 12:54:51 INFO Found: ['5065841']
2026-07-02 12:54:56 INFO [TGCC-IRENE] Submitted job with ID:['5065841']
2026-07-02 12:54:56 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 12:54:56 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10
2026-07-02 12:54:56 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020500_24_ENS10.nc
2026-07-02 12:54:56 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-02 12:54:56 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 12:54:56 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-02 12:54:56 INFO Queuing job for member 10...
2026-07-02 12:54:56 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 12:54:56 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-02 12:54:56 INFO Found: ['5065842']
2026-07-02 12:55:01 INFO [TGCC-IRENE] Submitted job with ID:['5065842']
2026-07-02 12:55:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 12:55:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11
2026-07-02 12:55:01 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020500_24_ENS11.nc
2026-07-02 12:55:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-02 12:55:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 12:55:01 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-02 12:55:02 INFO Queuing job for member 11...
2026-07-02 12:55:02 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 12:55:02 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-02 12:55:02 INFO Found: ['5065844']
2026-07-02 12:55:07 INFO [TGCC-IRENE] Submitted job with ID:['5065844']
2026-07-02 12:55:07 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 12:55:07 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12
2026-07-02 12:55:07 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020500_24_ENS12.nc
2026-07-02 12:55:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-02 12:55:07 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 12:55:07 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-02 12:55:07 INFO Queuing job for member 12...
2026-07-02 12:55:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 12:55:07 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-02 12:55:08 INFO Found: ['5065845']
2026-07-02 12:55:13 INFO [TGCC-IRENE] Submitted job with ID:['5065845']
2026-07-02 12:55:13 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 12:55:13 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13
2026-07-02 12:55:13 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020500_24_ENS13.nc
2026-07-02 12:55:13 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-02 12:55:13 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 12:55:13 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-02 12:55:13 INFO Queuing job for member 13...
2026-07-02 12:55:13 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 12:55:13 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-02 12:55:14 INFO Found: ['5065846']
2026-07-02 12:55:19 INFO [TGCC-IRENE] Submitted job with ID:['5065846']
2026-07-02 12:55:19 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 12:55:19 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14
2026-07-02 12:55:19 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020500_24_ENS14.nc
2026-07-02 12:55:19 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-02 12:55:19 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 12:55:19 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-02 12:55:19 INFO Queuing job for member 14...
2026-07-02 12:55:19 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 12:55:19 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-02 12:55:20 INFO Found: ['5065847']
2026-07-02 12:55:25 INFO [TGCC-IRENE] Submitted job with ID:['5065847']
2026-07-02 12:55:25 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 12:55:25 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15
2026-07-02 12:55:25 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020500_24_ENS15.nc
2026-07-02 12:55:25 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-02 12:55:25 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 12:55:25 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-02 12:55:25 INFO Queuing job for member 15...
2026-07-02 12:55:25 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 12:55:25 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-02 12:55:25 INFO Found: ['5065849']
2026-07-02 12:55:30 INFO [TGCC-IRENE] Submitted job with ID:['5065849']
2026-07-02 12:55:30 INFO Checking job status ...
2026-07-02 12:55:31 INFO None 5065830: status RUNNING/PENDING
2026-07-02 12:55:31 INFO None 5065832: status RUNNING/PENDING
2026-07-02 12:55:31 INFO None 5065833: status RUNNING/PENDING
2026-07-02 12:55:31 INFO None 5065835: status RUNNING/PENDING
2026-07-02 12:55:31 INFO None 5065836: status RUNNING/PENDING
2026-07-02 12:55:31 INFO None 5065838: status RUNNING/PENDING
2026-07-02 12:55:31 INFO None 5065839: status RUNNING/PENDING
2026-07-02 12:55:31 INFO None 5065840: status RUNNING/PENDING
2026-07-02 12:55:31 INFO None 5065841: status RUNNING/PENDING
2026-07-02 12:55:31 INFO None 5065842: status RUNNING/PENDING
2026-07-02 12:55:31 INFO None 5065844: status RUNNING/PENDING
2026-07-02 12:55:31 INFO None 5065845: status RUNNING/PENDING
2026-07-02 12:55:31 INFO None 5065846: status RUNNING/PENDING
2026-07-02 12:55:31 INFO None 5065847: status RUNNING/PENDING
2026-07-02 12:55:31 INFO None 5065849: status RUNNING/PENDING
2026-07-02 12:55:31 INFO Jobs still running: ['5065830', '5065832', '5065833', '5065835', '5065836', '5065838', '5065839', '5065840', '5065841', '5065842', '5065844', '5065845', '5065846', '5065847', '5065849']. Waiting...
2026-07-02 12:55:47 INFO None 5065830: status RUNNING/PENDING
2026-07-02 12:55:47 INFO None 5065832: status RUNNING/PENDING
2026-07-02 12:55:47 INFO None 5065833: status RUNNING/PENDING
2026-07-02 12:55:47 INFO None 5065835: status RUNNING/PENDING
2026-07-02 12:55:47 INFO None 5065836: status RUNNING/PENDING
2026-07-02 12:55:47 INFO None 5065838: status RUNNING/PENDING
2026-07-02 12:55:47 INFO None 5065839: status RUNNING/PENDING
2026-07-02 12:55:47 INFO None 5065840: status RUNNING/PENDING
2026-07-02 12:55:47 INFO None 5065841: status RUNNING/PENDING
2026-07-02 12:55:47 INFO None 5065842: status RUNNING/PENDING
2026-07-02 12:55:47 INFO None 5065844: status RUNNING/PENDING
2026-07-02 12:55:47 INFO None 5065845: status RUNNING/PENDING
2026-07-02 12:55:47 INFO None 5065846: status RUNNING/PENDING
2026-07-02 12:55:47 INFO None 5065847: status RUNNING/PENDING
2026-07-02 12:55:47 INFO None 5065849: status RUNNING/PENDING
2026-07-02 12:55:47 INFO Jobs still running: ['5065830', '5065832', '5065833', '5065835', '5065836', '5065838', '5065839', '5065840', '5065841', '5065842', '5065844', '5065845', '5065846', '5065847', '5065849']. Waiting...
2026-07-02 12:56:02 INFO None 5065830: status RUNNING/PENDING
2026-07-02 12:56:02 INFO None 5065832: status RUNNING/PENDING
2026-07-02 12:56:02 INFO None 5065833: status RUNNING/PENDING
2026-07-02 12:56:02 INFO None 5065835: status RUNNING/PENDING
2026-07-02 12:56:02 INFO None 5065836: status RUNNING/PENDING
2026-07-02 12:56:02 INFO None 5065838: status RUNNING/PENDING
2026-07-02 12:56:02 INFO None 5065839: status RUNNING/PENDING
2026-07-02 12:56:02 INFO None 5065840: status RUNNING/PENDING
2026-07-02 12:56:02 INFO None 5065841: status RUNNING/PENDING
2026-07-02 12:56:02 INFO None 5065842: status RUNNING/PENDING
2026-07-02 12:56:02 INFO None 5065844: status RUNNING/PENDING
2026-07-02 12:56:02 INFO None 5065845: status RUNNING/PENDING
2026-07-02 12:56:03 INFO None 5065846: status RUNNING/PENDING
2026-07-02 12:56:03 INFO None 5065847: status RUNNING/PENDING
2026-07-02 12:56:03 INFO None 5065849: status RUNNING/PENDING
2026-07-02 12:56:03 INFO Jobs still running: ['5065830', '5065832', '5065833', '5065835', '5065836', '5065838', '5065839', '5065840', '5065841', '5065842', '5065844', '5065845', '5065846', '5065847', '5065849']. Waiting...
2026-07-02 12:56:18 INFO None 5065830: status RUNNING/PENDING
2026-07-02 12:56:18 INFO None 5065832: status RUNNING/PENDING
2026-07-02 12:56:18 INFO None 5065833: status RUNNING/PENDING
2026-07-02 12:56:18 INFO None 5065835: status RUNNING/PENDING
2026-07-02 12:56:18 INFO None 5065836: status RUNNING/PENDING
2026-07-02 12:56:18 INFO None 5065838: status RUNNING/PENDING
2026-07-02 12:56:18 INFO None 5065839: status RUNNING/PENDING
2026-07-02 12:56:18 INFO None 5065840: status RUNNING/PENDING
2026-07-02 12:56:18 INFO None 5065841: status RUNNING/PENDING
2026-07-02 12:56:18 INFO None 5065842: status RUNNING/PENDING
2026-07-02 12:56:18 INFO None 5065844: status RUNNING/PENDING
2026-07-02 12:56:18 INFO None 5065845: status RUNNING/PENDING
2026-07-02 12:56:18 INFO None 5065846: status RUNNING/PENDING
2026-07-02 12:56:18 INFO None 5065847: status RUNNING/PENDING
2026-07-02 12:56:18 INFO None 5065849: status RUNNING/PENDING
2026-07-02 12:56:18 INFO Jobs still running: ['5065830', '5065832', '5065833', '5065835', '5065836', '5065838', '5065839', '5065840', '5065841', '5065842', '5065844', '5065845', '5065846', '5065847', '5065849']. Waiting...
2026-07-02 12:56:34 INFO None 5065830: status RUNNING/PENDING
2026-07-02 12:56:34 INFO None 5065832: status RUNNING/PENDING
2026-07-02 12:56:34 INFO None 5065833: status RUNNING/PENDING
2026-07-02 12:56:34 INFO None 5065835: status RUNNING/PENDING
2026-07-02 12:56:34 INFO None 5065836: status RUNNING/PENDING
2026-07-02 12:56:34 INFO None 5065838: status RUNNING/PENDING
2026-07-02 12:56:34 INFO None 5065839: status RUNNING/PENDING
2026-07-02 12:56:34 INFO None 5065840: status RUNNING/PENDING
2026-07-02 12:56:34 INFO None 5065841: status RUNNING/PENDING
2026-07-02 12:56:34 INFO None 5065842: status RUNNING/PENDING
2026-07-02 12:56:34 INFO None 5065844: status RUNNING/PENDING
2026-07-02 12:56:34 INFO None 5065845: status RUNNING/PENDING
2026-07-02 12:56:34 INFO None 5065846: status RUNNING/PENDING
2026-07-02 12:56:34 INFO None 5065847: status RUNNING/PENDING
2026-07-02 12:56:34 INFO None 5065849: status RUNNING/PENDING
2026-07-02 12:56:34 INFO Jobs still running: ['5065830', '5065832', '5065833', '5065835', '5065836', '5065838', '5065839', '5065840', '5065841', '5065842', '5065844', '5065845', '5065846', '5065847', '5065849']. Waiting...
2026-07-02 12:56:49 INFO None 5065830: status RUNNING/PENDING
2026-07-02 12:56:49 INFO None 5065832: status RUNNING/PENDING
2026-07-02 12:56:49 INFO None 5065833: status RUNNING/PENDING
2026-07-02 12:56:49 INFO None 5065835: status RUNNING/PENDING
2026-07-02 12:56:49 INFO None 5065836: status RUNNING/PENDING
2026-07-02 12:56:49 INFO None 5065838: status RUNNING/PENDING
2026-07-02 12:56:49 INFO None 5065839: status RUNNING/PENDING
2026-07-02 12:56:49 INFO None 5065840: status RUNNING/PENDING
2026-07-02 12:56:49 INFO None 5065841: status RUNNING/PENDING
2026-07-02 12:56:49 INFO None 5065842: status RUNNING/PENDING
2026-07-02 12:56:50 INFO None 5065844: status RUNNING/PENDING
2026-07-02 12:56:50 INFO None 5065845: status RUNNING/PENDING
2026-07-02 12:56:50 INFO None 5065846: status RUNNING/PENDING
2026-07-02 12:56:50 INFO None 5065847: status RUNNING/PENDING
2026-07-02 12:56:50 INFO None 5065849: status RUNNING/PENDING
2026-07-02 12:56:50 INFO Jobs still running: ['5065830', '5065832', '5065833', '5065835', '5065836', '5065838', '5065839', '5065840', '5065841', '5065842', '5065844', '5065845', '5065846', '5065847', '5065849']. Waiting...
2026-07-02 12:57:05 INFO None 5065830: status RUNNING/PENDING
2026-07-02 12:57:05 INFO None 5065832: status FINISHED
2026-07-02 12:57:05 INFO None 5065833: status RUNNING/PENDING
2026-07-02 12:57:05 INFO None 5065835: status RUNNING/PENDING
2026-07-02 12:57:05 INFO None 5065836: status FINISHED
2026-07-02 12:57:05 INFO None 5065838: status FINISHED
2026-07-02 12:57:05 INFO None 5065839: status RUNNING/PENDING
2026-07-02 12:57:05 INFO None 5065840: status RUNNING/PENDING
2026-07-02 12:57:05 INFO None 5065841: status RUNNING/PENDING
2026-07-02 12:57:05 INFO None 5065842: status RUNNING/PENDING
2026-07-02 12:57:05 INFO None 5065844: status RUNNING/PENDING
2026-07-02 12:57:05 INFO None 5065845: status RUNNING/PENDING
2026-07-02 12:57:05 INFO None 5065846: status RUNNING/PENDING
2026-07-02 12:57:05 INFO None 5065847: status RUNNING/PENDING
2026-07-02 12:57:05 INFO None 5065849: status RUNNING/PENDING
2026-07-02 12:57:05 INFO Jobs still running: ['5065830', '5065833', '5065835', '5065839', '5065840', '5065841', '5065842', '5065844', '5065845', '5065846', '5065847', '5065849']. Waiting...
2026-07-02 12:57:20 INFO None 5065830: status RUNNING/PENDING
2026-07-02 12:57:20 INFO None 5065832: status FINISHED
2026-07-02 12:57:20 INFO None 5065833: status RUNNING/PENDING
2026-07-02 12:57:20 INFO None 5065835: status RUNNING/PENDING
2026-07-02 12:57:20 INFO None 5065836: status FINISHED
2026-07-02 12:57:20 INFO None 5065838: status FINISHED
2026-07-02 12:57:20 INFO None 5065839: status RUNNING/PENDING
2026-07-02 12:57:20 INFO None 5065840: status FINISHED
2026-07-02 12:57:20 INFO None 5065841: status RUNNING/PENDING
2026-07-02 12:57:20 INFO None 5065842: status RUNNING/PENDING
2026-07-02 12:57:20 INFO None 5065844: status RUNNING/PENDING
2026-07-02 12:57:20 INFO None 5065845: status RUNNING/PENDING
2026-07-02 12:57:20 INFO None 5065846: status RUNNING/PENDING
2026-07-02 12:57:20 INFO None 5065847: status RUNNING/PENDING
2026-07-02 12:57:20 INFO None 5065849: status RUNNING/PENDING
2026-07-02 12:57:20 INFO Jobs still running: ['5065830', '5065833', '5065835', '5065839', '5065841', '5065842', '5065844', '5065845', '5065846', '5065847', '5065849']. Waiting...
2026-07-02 12:57:35 INFO None 5065830: status RUNNING/PENDING
2026-07-02 12:57:35 INFO None 5065832: status FINISHED
2026-07-02 12:57:35 INFO None 5065833: status FINISHED
2026-07-02 12:57:35 INFO None 5065835: status FINISHED
2026-07-02 12:57:35 INFO None 5065836: status FINISHED
2026-07-02 12:57:35 INFO None 5065838: status FINISHED
2026-07-02 12:57:35 INFO None 5065839: status FINISHED
2026-07-02 12:57:35 INFO None 5065840: status FINISHED
2026-07-02 12:57:35 INFO None 5065841: status RUNNING/PENDING
2026-07-02 12:57:35 INFO None 5065842: status RUNNING/PENDING
2026-07-02 12:57:35 INFO None 5065844: status RUNNING/PENDING
2026-07-02 12:57:35 INFO None 5065845: status RUNNING/PENDING
2026-07-02 12:57:35 INFO None 5065846: status RUNNING/PENDING
2026-07-02 12:57:36 INFO None 5065847: status RUNNING/PENDING
2026-07-02 12:57:36 INFO None 5065849: status RUNNING/PENDING
2026-07-02 12:57:36 INFO Jobs still running: ['5065830', '5065841', '5065842', '5065844', '5065845', '5065846', '5065847', '5065849']. Waiting...
2026-07-02 12:57:51 INFO None 5065830: status RUNNING/PENDING
2026-07-02 12:57:51 INFO None 5065832: status FINISHED
2026-07-02 12:57:51 INFO None 5065833: status FINISHED
2026-07-02 12:57:51 INFO None 5065835: status FINISHED
2026-07-02 12:57:51 INFO None 5065836: status FINISHED
2026-07-02 12:57:51 INFO None 5065838: status FINISHED
2026-07-02 12:57:51 INFO None 5065839: status FINISHED
2026-07-02 12:57:51 INFO None 5065840: status FINISHED
2026-07-02 12:57:51 INFO None 5065841: status FINISHED
2026-07-02 12:57:51 INFO None 5065842: status RUNNING/PENDING
2026-07-02 12:57:51 INFO None 5065844: status RUNNING/PENDING
2026-07-02 12:57:51 INFO None 5065845: status RUNNING/PENDING
2026-07-02 12:57:51 INFO None 5065846: status RUNNING/PENDING
2026-07-02 12:57:51 INFO None 5065847: status RUNNING/PENDING
2026-07-02 12:57:51 INFO None 5065849: status RUNNING/PENDING
2026-07-02 12:57:51 INFO Jobs still running: ['5065830', '5065842', '5065844', '5065845', '5065846', '5065847', '5065849']. Waiting...
2026-07-02 12:58:06 INFO None 5065830: status RUNNING/PENDING
2026-07-02 12:58:06 INFO None 5065832: status FINISHED
2026-07-02 12:58:06 INFO None 5065833: status FINISHED
2026-07-02 12:58:06 INFO None 5065835: status FINISHED
2026-07-02 12:58:06 INFO None 5065836: status FINISHED
2026-07-02 12:58:06 INFO None 5065838: status FINISHED
2026-07-02 12:58:06 INFO None 5065839: status FINISHED
2026-07-02 12:58:06 INFO None 5065840: status FINISHED
2026-07-02 12:58:06 INFO None 5065841: status FINISHED
2026-07-02 12:58:06 INFO None 5065842: status FINISHED
2026-07-02 12:58:06 INFO None 5065844: status RUNNING/PENDING
2026-07-02 12:58:06 INFO None 5065845: status RUNNING/PENDING
2026-07-02 12:58:06 INFO None 5065846: status FINISHED
2026-07-02 12:58:06 INFO None 5065847: status FINISHED
2026-07-02 12:58:06 INFO None 5065849: status RUNNING/PENDING
2026-07-02 12:58:06 INFO Jobs still running: ['5065830', '5065844', '5065845', '5065849']. Waiting...
2026-07-02 12:58:21 INFO None 5065830: status RUNNING/PENDING
2026-07-02 12:58:22 INFO None 5065832: status FINISHED
2026-07-02 12:58:22 INFO None 5065833: status FINISHED
2026-07-02 12:58:22 INFO None 5065835: status FINISHED
2026-07-02 12:58:22 INFO None 5065836: status FINISHED
2026-07-02 12:58:22 INFO None 5065838: status FINISHED
2026-07-02 12:58:22 INFO None 5065839: status FINISHED
2026-07-02 12:58:22 INFO None 5065840: status FINISHED
2026-07-02 12:58:22 INFO None 5065841: status FINISHED
2026-07-02 12:58:22 INFO None 5065842: status FINISHED
2026-07-02 12:58:22 INFO None 5065844: status FINISHED
2026-07-02 12:58:22 INFO None 5065845: status RUNNING/PENDING
2026-07-02 12:58:22 INFO None 5065846: status FINISHED
2026-07-02 12:58:22 INFO None 5065847: status FINISHED
2026-07-02 12:58:22 INFO None 5065849: status RUNNING/PENDING
2026-07-02 12:58:22 INFO Jobs still running: ['5065830', '5065845', '5065849']. Waiting...
2026-07-02 12:58:37 INFO None 5065830: status FINISHED
2026-07-02 12:58:37 INFO None 5065832: status FINISHED
2026-07-02 12:58:37 INFO None 5065833: status FINISHED
2026-07-02 12:58:37 INFO None 5065835: status FINISHED
2026-07-02 12:58:37 INFO None 5065836: status FINISHED
2026-07-02 12:58:37 INFO None 5065838: status FINISHED
2026-07-02 12:58:37 INFO None 5065839: status FINISHED
2026-07-02 12:58:37 INFO None 5065840: status FINISHED
2026-07-02 12:58:37 INFO None 5065841: status FINISHED
2026-07-02 12:58:37 INFO None 5065842: status FINISHED
2026-07-02 12:58:37 INFO None 5065844: status FINISHED
2026-07-02 12:58:37 INFO None 5065845: status FINISHED
2026-07-02 12:58:37 INFO None 5065846: status FINISHED
2026-07-02 12:58:37 INFO None 5065847: status FINISHED
2026-07-02 12:58:37 INFO None 5065849: status RUNNING/PENDING
2026-07-02 12:58:37 INFO Jobs still running: ['5065849']. Waiting...
2026-07-02 12:58:52 INFO None 5065830: status FINISHED
2026-07-02 12:58:52 INFO None 5065832: status FINISHED
2026-07-02 12:58:52 INFO None 5065833: status FINISHED
2026-07-02 12:58:52 INFO None 5065835: status FINISHED
2026-07-02 12:58:52 INFO None 5065836: status FINISHED
2026-07-02 12:58:52 INFO None 5065838: status FINISHED
2026-07-02 12:58:52 INFO None 5065839: status FINISHED
2026-07-02 12:58:52 INFO None 5065840: status FINISHED
2026-07-02 12:58:52 INFO None 5065841: status FINISHED
2026-07-02 12:58:52 INFO None 5065842: status FINISHED
2026-07-02 12:58:52 INFO None 5065844: status FINISHED
2026-07-02 12:58:52 INFO None 5065845: status FINISHED
2026-07-02 12:58:52 INFO None 5065846: status FINISHED
2026-07-02 12:58:52 INFO None 5065847: status FINISHED
2026-07-02 12:58:52 INFO None 5065849: status RUNNING/PENDING
2026-07-02 12:58:52 INFO Jobs still running: ['5065849']. Waiting...
2026-07-02 12:59:07 INFO None 5065830: status FINISHED
2026-07-02 12:59:07 INFO None 5065832: status FINISHED
2026-07-02 12:59:07 INFO None 5065833: status FINISHED
2026-07-02 12:59:07 INFO None 5065835: status FINISHED
2026-07-02 12:59:07 INFO None 5065836: status FINISHED
2026-07-02 12:59:07 INFO None 5065838: status FINISHED
2026-07-02 12:59:07 INFO None 5065839: status FINISHED
2026-07-02 12:59:07 INFO None 5065840: status FINISHED
2026-07-02 12:59:07 INFO None 5065841: status FINISHED
2026-07-02 12:59:07 INFO None 5065842: status FINISHED
2026-07-02 12:59:07 INFO None 5065844: status FINISHED
2026-07-02 12:59:08 INFO None 5065845: status FINISHED
2026-07-02 12:59:08 INFO None 5065846: status FINISHED
2026-07-02 12:59:08 INFO None 5065847: status FINISHED
2026-07-02 12:59:08 INFO None 5065849: status FINISHED
2026-07-02 12:59:08 INFO Jobs ['5065830', '5065832', '5065833', '5065835', '5065836', '5065838', '5065839', '5065840', '5065841', '5065842', '5065844', '5065845', '5065846', '5065847', '5065849'] have finished
2026-07-02 12:59:08 INFO Checking restart files were created ...
2026-07-02 12:59:08 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020600_1_ENS1.nc(668832435 bytes)
2026-07-02 12:59:08 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020600_1_ENS2.nc(668832435 bytes)
2026-07-02 12:59:08 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020600_1_ENS3.nc(668832435 bytes)
2026-07-02 12:59:08 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020600_1_ENS4.nc(668832435 bytes)
2026-07-02 12:59:08 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020600_1_ENS5.nc(668832435 bytes)
2026-07-02 12:59:08 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020600_1_ENS6.nc(668832435 bytes)
2026-07-02 12:59:08 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020600_1_ENS7.nc(668832435 bytes)
2026-07-02 12:59:08 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020600_1_ENS8.nc(668832435 bytes)
2026-07-02 12:59:08 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020600_1_ENS9.nc(668832435 bytes)
2026-07-02 12:59:08 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020600_1_ENS10.nc(668832435 bytes)
2026-07-02 12:59:08 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020600_1_ENS11.nc(668832435 bytes)
2026-07-02 12:59:08 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020600_1_ENS12.nc(668832435 bytes)
2026-07-02 12:59:08 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020600_1_ENS13.nc(668832435 bytes)
2026-07-02 12:59:08 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020600_1_ENS14.nc(668832435 bytes)
2026-07-02 12:59:08 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020600_1_ENS15.nc(668832435 bytes)
2026-07-02 12:59:08 INFO  Run_model() completed successfully.
2026-07-02 12:59:08 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 00:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 12:59:08 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 01:00:00 days=153072 seconds=3600
2026-07-02 12:59:08 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-02 12:59:08 INFO [TIME] increment current_time 2020-02-06 00:00:00 -> 2020-02-06 01:00:00
2026-07-02 12:59:08 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 12:59:08 INFO ---------->>> Running process_satellite_data()
2026-07-02 12:59:08 INFO [DART] No satellite data found, skipping assimilation
2026-07-02 12:59:08 INFO after_assimilation() skipped
2026-07-02 12:59:08 INFO Next run starts from 2020-02-06 01:00:00
2026-07-02 12:59:08 INFO Cycle is DONE; starting a new loop!
2026-07-02 12:59:08 INFO [TIME] step_end current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 12:59:08 INFO [TIME] step_start current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 12:59:08 INFO [TIME] window start=2020-02-06 01:00:00 end=2020-02-06 09:00:00 run_hours=8 has_assimilation=True
2026-07-02 12:59:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 12:59:09 INFO Hourly dataset computed and listing created
2026-07-02 12:59:19 INFO Hourly dataset computed
2026-07-02 12:59:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 12:59:20 INFO Hourly dataset computed and listing created
2026-07-02 12:59:22 INFO Hourly dataset computed
2026-07-02 12:59:22 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 12:59:23 INFO Hourly dataset computed and listing created
2026-07-02 12:59:25 INFO Hourly dataset computed
2026-07-02 12:59:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 12:59:26 INFO Hourly dataset computed and listing created
2026-07-02 12:59:28 INFO Hourly dataset computed
2026-07-02 12:59:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 12:59:29 INFO Hourly dataset computed and listing created
2026-07-02 12:59:31 INFO Hourly dataset computed
2026-07-02 12:59:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 12:59:33 INFO Hourly dataset computed and listing created
2026-07-02 12:59:35 INFO Hourly dataset computed
2026-07-02 12:59:35 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 12:59:36 INFO Hourly dataset computed and listing created
2026-07-02 12:59:38 INFO Hourly dataset computed
2026-07-02 12:59:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 12:59:39 INFO Hourly dataset computed and listing created
2026-07-02 12:59:41 INFO Hourly dataset computed
2026-07-02 12:59:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 12:59:42 INFO Hourly dataset computed and listing created
2026-07-02 12:59:44 INFO Hourly dataset computed
2026-07-02 12:59:44 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 12:59:45 INFO Hourly dataset computed and listing created
2026-07-02 12:59:47 INFO Hourly dataset computed
2026-07-02 12:59:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 12:59:49 INFO Hourly dataset computed and listing created
2026-07-02 12:59:51 INFO Hourly dataset computed
2026-07-02 12:59:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 12:59:52 INFO Hourly dataset computed and listing created
2026-07-02 12:59:54 INFO Hourly dataset computed
2026-07-02 12:59:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 12:59:55 INFO Hourly dataset computed and listing created
2026-07-02 12:59:57 INFO Hourly dataset computed
2026-07-02 12:59:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 12:59:58 INFO Hourly dataset computed and listing created
2026-07-02 13:00:00 INFO Hourly dataset computed
2026-07-02 13:00:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:00:01 INFO Hourly dataset computed and listing created
2026-07-02 13:00:03 INFO Hourly dataset computed
2026-07-02 13:00:03 INFO ---------->>> Running CHIMERE model from 2020-02-06 01:00:00 to 2020-02-06 09:00:00
2026-07-02 13:00:03 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:00:03 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1
2026-07-02 13:00:03 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020600_1_ENS1.nc
2026-07-02 13:00:03 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-02 13:00:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:00:03 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-02 13:00:03 INFO Queuing job for member 1...
2026-07-02 13:00:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:00:03 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-02 13:00:04 INFO Found: ['5065863']
2026-07-02 13:00:09 INFO [TGCC-IRENE] Submitted job with ID:['5065863']
2026-07-02 13:00:09 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:00:09 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2
2026-07-02 13:00:09 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020600_1_ENS2.nc
2026-07-02 13:00:09 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-02 13:00:09 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:00:09 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-02 13:00:09 INFO Queuing job for member 2...
2026-07-02 13:00:09 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:00:09 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-02 13:00:10 INFO Found: ['5065864']
2026-07-02 13:00:15 INFO [TGCC-IRENE] Submitted job with ID:['5065864']
2026-07-02 13:00:15 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:00:15 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3
2026-07-02 13:00:15 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020600_1_ENS3.nc
2026-07-02 13:00:15 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-02 13:00:15 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:00:15 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-02 13:00:15 INFO Queuing job for member 3...
2026-07-02 13:00:15 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:00:15 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-02 13:00:16 INFO Found: ['5065865']
2026-07-02 13:00:21 INFO [TGCC-IRENE] Submitted job with ID:['5065865']
2026-07-02 13:00:21 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:00:21 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4
2026-07-02 13:00:21 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020600_1_ENS4.nc
2026-07-02 13:00:21 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-02 13:00:21 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:00:21 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-02 13:00:21 INFO Queuing job for member 4...
2026-07-02 13:00:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:00:21 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-02 13:00:22 INFO Found: ['5065866']
2026-07-02 13:00:27 INFO [TGCC-IRENE] Submitted job with ID:['5065866']
2026-07-02 13:00:27 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:00:27 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5
2026-07-02 13:00:27 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020600_1_ENS5.nc
2026-07-02 13:00:27 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-02 13:00:27 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:00:27 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-02 13:00:27 INFO Queuing job for member 5...
2026-07-02 13:00:27 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:00:27 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-02 13:00:27 INFO Found: ['5065868']
2026-07-02 13:00:32 INFO [TGCC-IRENE] Submitted job with ID:['5065868']
2026-07-02 13:00:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:00:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6
2026-07-02 13:00:32 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020600_1_ENS6.nc
2026-07-02 13:00:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-02 13:00:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:00:32 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-02 13:00:32 INFO Queuing job for member 6...
2026-07-02 13:00:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:00:32 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-02 13:00:33 INFO Found: ['5065869']
2026-07-02 13:00:38 INFO [TGCC-IRENE] Submitted job with ID:['5065869']
2026-07-02 13:00:38 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:00:38 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7
2026-07-02 13:00:38 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020600_1_ENS7.nc
2026-07-02 13:00:38 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-02 13:00:38 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:00:38 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-02 13:00:38 INFO Queuing job for member 7...
2026-07-02 13:00:38 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:00:38 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-02 13:00:39 INFO Found: ['5065870']
2026-07-02 13:00:44 INFO [TGCC-IRENE] Submitted job with ID:['5065870']
2026-07-02 13:00:44 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:00:44 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8
2026-07-02 13:00:44 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020600_1_ENS8.nc
2026-07-02 13:00:44 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-02 13:00:44 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:00:44 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-02 13:00:44 INFO Queuing job for member 8...
2026-07-02 13:00:44 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:00:44 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-02 13:00:45 INFO Found: ['5065871']
2026-07-02 13:00:50 INFO [TGCC-IRENE] Submitted job with ID:['5065871']
2026-07-02 13:00:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:00:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9
2026-07-02 13:00:50 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020600_1_ENS9.nc
2026-07-02 13:00:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-02 13:00:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:00:50 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-02 13:00:50 INFO Queuing job for member 9...
2026-07-02 13:00:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:00:50 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-02 13:00:50 INFO Found: ['5065872']
2026-07-02 13:00:55 INFO [TGCC-IRENE] Submitted job with ID:['5065872']
2026-07-02 13:00:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:00:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10
2026-07-02 13:00:55 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020600_1_ENS10.nc
2026-07-02 13:00:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-02 13:00:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:00:55 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-02 13:00:55 INFO Queuing job for member 10...
2026-07-02 13:00:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:00:55 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-02 13:00:56 INFO Found: ['5065873']
2026-07-02 13:01:01 INFO [TGCC-IRENE] Submitted job with ID:['5065873']
2026-07-02 13:01:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:01:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11
2026-07-02 13:01:01 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020600_1_ENS11.nc
2026-07-02 13:01:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-02 13:01:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:01:01 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-02 13:01:01 INFO Queuing job for member 11...
2026-07-02 13:01:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:01:01 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-02 13:01:03 INFO Found: ['5065874']
2026-07-02 13:01:08 INFO [TGCC-IRENE] Submitted job with ID:['5065874']
2026-07-02 13:01:08 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:01:08 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12
2026-07-02 13:01:08 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020600_1_ENS12.nc
2026-07-02 13:01:08 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-02 13:01:08 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:01:08 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-02 13:01:08 INFO Queuing job for member 12...
2026-07-02 13:01:08 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:01:08 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-02 13:01:09 INFO Found: ['5065876']
2026-07-02 13:01:14 INFO [TGCC-IRENE] Submitted job with ID:['5065876']
2026-07-02 13:01:14 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:01:14 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13
2026-07-02 13:01:14 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020600_1_ENS13.nc
2026-07-02 13:01:14 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-02 13:01:14 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:01:14 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-02 13:01:14 INFO Queuing job for member 13...
2026-07-02 13:01:14 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:01:14 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-02 13:01:16 INFO Found: ['5065878']
2026-07-02 13:01:21 INFO [TGCC-IRENE] Submitted job with ID:['5065878']
2026-07-02 13:01:21 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:01:21 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14
2026-07-02 13:01:21 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020600_1_ENS14.nc
2026-07-02 13:01:21 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-02 13:01:21 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:01:21 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-02 13:01:21 INFO Queuing job for member 14...
2026-07-02 13:01:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:01:21 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-02 13:01:22 INFO Found: ['5065882']
2026-07-02 13:01:27 INFO [TGCC-IRENE] Submitted job with ID:['5065882']
2026-07-02 13:01:27 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:01:27 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15
2026-07-02 13:01:27 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020600_1_ENS15.nc
2026-07-02 13:01:27 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-02 13:01:27 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:01:27 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-02 13:01:27 INFO Queuing job for member 15...
2026-07-02 13:01:27 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:01:27 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-02 13:01:28 INFO Found: ['5065884']
2026-07-02 13:01:33 INFO [TGCC-IRENE] Submitted job with ID:['5065884']
2026-07-02 13:01:33 INFO Checking job status ...
2026-07-02 13:01:33 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:01:33 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:01:33 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:01:33 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:01:33 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:01:33 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:01:33 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:01:33 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:01:33 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:01:33 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:01:33 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:01:33 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:01:33 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:01:33 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:01:33 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:01:33 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:01:48 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:01:48 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:01:48 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:01:48 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:01:48 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:01:48 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:01:48 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:01:48 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:01:48 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:01:48 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:01:48 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:01:48 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:01:48 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:01:48 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:01:48 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:01:48 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:02:04 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:02:04 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:02:04 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:02:04 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:02:04 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:02:04 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:02:04 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:02:04 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:02:04 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:02:04 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:02:04 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:02:04 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:02:04 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:02:04 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:02:04 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:02:04 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:02:19 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:02:19 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:02:19 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:02:19 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:02:19 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:02:19 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:02:19 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:02:19 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:02:19 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:02:19 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:02:19 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:02:19 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:02:19 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:02:19 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:02:19 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:02:19 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:02:35 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:02:35 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:02:35 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:02:35 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:02:35 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:02:35 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:02:35 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:02:35 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:02:35 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:02:35 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:02:35 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:02:35 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:02:35 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:02:35 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:02:35 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:02:35 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:02:50 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:02:50 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:02:50 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:02:50 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:02:50 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:02:50 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:02:50 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:02:50 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:02:50 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:02:50 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:02:50 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:02:50 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:02:50 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:02:50 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:02:50 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:02:50 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:03:07 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:03:07 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:03:07 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:03:07 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:03:07 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:03:07 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:03:07 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:03:07 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:03:07 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:03:07 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:03:07 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:03:07 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:03:07 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:03:07 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:03:07 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:03:07 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:03:22 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:03:22 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:03:22 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:03:22 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:03:22 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:03:22 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:03:22 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:03:22 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:03:22 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:03:22 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:03:22 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:03:22 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:03:22 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:03:22 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:03:22 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:03:22 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:03:37 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:03:38 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:03:38 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:03:38 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:03:38 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:03:38 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:03:38 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:03:38 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:03:38 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:03:38 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:03:38 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:03:38 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:03:38 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:03:38 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:03:38 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:03:38 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:03:53 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:03:53 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:03:53 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:03:53 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:03:53 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:03:53 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:03:53 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:03:53 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:03:53 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:03:53 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:03:53 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:03:53 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:03:53 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:03:53 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:03:53 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:03:53 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:04:08 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:04:08 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:04:08 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:04:08 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:04:08 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:04:08 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:04:08 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:04:08 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:04:08 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:04:08 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:04:08 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:04:08 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:04:08 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:04:08 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:04:08 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:04:08 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:04:23 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:04:23 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:04:23 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:04:23 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:04:23 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:04:23 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:04:23 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:04:24 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:04:24 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:04:24 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:04:24 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:04:24 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:04:24 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:04:24 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:04:24 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:04:24 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:04:39 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:04:39 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:04:39 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:04:39 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:04:39 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:04:39 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:04:39 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:04:39 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:04:39 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:04:39 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:04:39 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:04:39 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:04:39 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:04:39 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:04:39 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:04:39 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:04:55 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:04:55 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:04:55 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:04:55 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:04:55 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:04:55 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:04:55 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:04:55 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:04:55 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:04:55 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:04:55 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:04:55 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:04:55 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:04:55 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:04:55 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:04:55 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:05:10 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:05:10 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:05:10 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:05:10 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:05:10 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:05:10 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:05:10 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:05:10 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:05:10 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:05:10 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:05:10 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:05:10 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:05:10 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:05:10 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:05:10 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:05:10 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:05:25 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:05:25 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:05:25 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:05:25 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:05:25 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:05:25 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:05:25 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:05:25 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:05:26 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:05:26 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:05:26 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:05:26 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:05:26 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:05:26 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:05:26 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:05:26 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:05:41 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:05:41 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:05:41 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:05:41 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:05:41 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:05:41 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:05:41 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:05:41 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:05:41 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:05:41 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:05:41 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:05:41 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:05:41 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:05:41 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:05:41 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:05:41 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:05:56 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:05:56 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:05:56 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:05:56 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:05:56 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:05:56 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:05:56 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:05:56 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:05:56 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:05:56 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:05:56 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:05:56 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:05:56 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:05:56 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:05:56 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:05:56 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:06:11 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:06:11 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:06:11 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:06:11 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:06:11 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:06:11 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:06:11 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:06:11 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:06:11 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:06:11 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:06:11 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:06:11 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:06:11 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:06:11 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:06:11 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:06:11 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:06:27 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:06:27 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:06:27 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:06:27 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:06:27 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:06:27 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:06:27 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:06:27 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:06:27 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:06:27 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:06:27 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:06:27 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:06:27 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:06:27 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:06:27 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:06:27 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:06:43 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:06:43 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:06:43 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:06:43 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:06:43 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:06:43 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:06:43 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:06:43 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:06:43 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:06:43 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:06:43 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:06:43 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:06:43 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:06:43 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:06:43 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:06:43 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:06:58 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:06:58 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:06:58 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:06:58 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:06:58 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:06:58 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:06:58 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:06:58 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:06:58 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:06:58 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:06:58 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:06:58 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:06:58 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:06:58 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:06:58 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:06:58 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:07:13 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:07:13 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:07:13 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:07:13 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:07:14 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:07:14 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:07:14 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:07:14 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:07:14 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:07:14 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:07:14 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:07:14 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:07:14 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:07:14 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:07:14 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:07:14 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:07:29 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:07:29 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:07:29 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:07:29 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:07:29 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:07:29 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:07:29 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:07:29 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:07:29 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:07:29 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:07:29 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:07:29 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:07:29 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:07:29 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:07:29 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:07:29 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:07:44 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:07:44 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:07:44 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:07:44 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:07:44 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:07:44 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:07:44 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:07:44 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:07:44 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:07:44 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:07:44 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:07:44 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:07:44 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:07:44 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:07:44 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:07:44 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:07:59 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:07:59 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:07:59 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:07:59 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:07:59 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:07:59 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:07:59 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:07:59 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:07:59 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:07:59 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:08:00 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:08:00 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:08:00 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:08:00 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:08:00 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:08:00 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:08:15 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:08:15 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:08:15 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:08:15 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:08:15 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:08:15 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:08:15 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:08:15 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:08:15 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:08:15 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:08:15 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:08:15 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:08:15 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:08:15 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:08:15 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:08:15 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:08:31 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:08:31 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:08:31 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:08:31 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:08:31 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:08:31 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:08:31 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:08:31 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:08:31 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:08:31 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:08:31 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:08:31 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:08:31 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:08:31 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:08:31 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:08:31 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:08:46 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:08:46 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:08:46 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:08:46 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:08:46 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:08:46 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:08:46 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:08:46 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:08:46 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:08:46 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:08:46 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:08:46 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:08:46 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:08:46 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:08:46 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:08:46 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:09:01 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:09:01 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:09:01 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:09:01 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:09:02 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:09:02 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:09:02 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:09:02 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:09:02 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:09:02 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:09:02 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:09:02 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:09:02 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:09:02 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:09:02 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:09:02 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:09:17 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:09:17 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:09:17 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:09:17 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:09:17 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:09:17 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:09:17 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:09:17 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:09:17 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:09:17 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:09:17 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:09:17 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:09:17 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:09:17 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:09:17 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:09:17 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:09:32 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:09:32 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:09:32 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:09:32 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:09:32 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:09:32 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:09:32 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:09:32 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:09:32 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:09:32 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:09:32 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:09:32 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:09:32 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:09:32 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:09:32 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:09:32 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:09:47 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:09:47 INFO None 5065864: status RUNNING/PENDING
2026-07-02 13:09:47 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:09:47 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:09:47 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:09:47 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:09:47 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:09:47 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:09:47 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:09:47 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:09:47 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:09:48 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:09:48 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:09:48 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:09:48 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:09:48 INFO Jobs still running: ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:10:03 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:10:03 INFO None 5065864: status FINISHED
2026-07-02 13:10:03 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:10:03 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:10:03 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:10:03 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:10:03 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:10:03 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:10:03 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:10:03 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:10:03 INFO None 5065874: status RUNNING/PENDING
2026-07-02 13:10:03 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:10:03 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:10:03 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:10:03 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:10:03 INFO Jobs still running: ['5065863', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:10:19 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:10:19 INFO None 5065864: status FINISHED
2026-07-02 13:10:19 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:10:19 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:10:19 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:10:19 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:10:19 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:10:19 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:10:19 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:10:19 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:10:19 INFO None 5065874: status FINISHED
2026-07-02 13:10:19 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:10:20 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:10:20 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:10:20 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:10:20 INFO Jobs still running: ['5065863', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:10:35 INFO None 5065863: status RUNNING/PENDING
2026-07-02 13:10:35 INFO None 5065864: status FINISHED
2026-07-02 13:10:35 INFO None 5065865: status RUNNING/PENDING
2026-07-02 13:10:35 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:10:35 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:10:35 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:10:35 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:10:35 INFO None 5065871: status RUNNING/PENDING
2026-07-02 13:10:35 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:10:35 INFO None 5065873: status RUNNING/PENDING
2026-07-02 13:10:35 INFO None 5065874: status FINISHED
2026-07-02 13:10:35 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:10:35 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:10:35 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:10:35 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:10:35 INFO Jobs still running: ['5065863', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:10:50 INFO None 5065863: status FINISHED
2026-07-02 13:10:50 INFO None 5065864: status FINISHED
2026-07-02 13:10:50 INFO None 5065865: status FINISHED
2026-07-02 13:10:50 INFO None 5065866: status RUNNING/PENDING
2026-07-02 13:10:50 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:10:50 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:10:50 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:10:50 INFO None 5065871: status FINISHED
2026-07-02 13:10:50 INFO None 5065872: status RUNNING/PENDING
2026-07-02 13:10:50 INFO None 5065873: status FINISHED
2026-07-02 13:10:50 INFO None 5065874: status FINISHED
2026-07-02 13:10:50 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:10:50 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:10:50 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:10:50 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:10:50 INFO Jobs still running: ['5065866', '5065868', '5065869', '5065870', '5065872', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:11:05 INFO None 5065863: status FINISHED
2026-07-02 13:11:05 INFO None 5065864: status FINISHED
2026-07-02 13:11:05 INFO None 5065865: status FINISHED
2026-07-02 13:11:05 INFO None 5065866: status FINISHED
2026-07-02 13:11:05 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:11:05 INFO None 5065869: status RUNNING/PENDING
2026-07-02 13:11:05 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:11:05 INFO None 5065871: status FINISHED
2026-07-02 13:11:05 INFO None 5065872: status FINISHED
2026-07-02 13:11:05 INFO None 5065873: status FINISHED
2026-07-02 13:11:05 INFO None 5065874: status FINISHED
2026-07-02 13:11:05 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:11:05 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:11:05 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:11:05 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:11:05 INFO Jobs still running: ['5065868', '5065869', '5065870', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:11:20 INFO None 5065863: status FINISHED
2026-07-02 13:11:20 INFO None 5065864: status FINISHED
2026-07-02 13:11:20 INFO None 5065865: status FINISHED
2026-07-02 13:11:20 INFO None 5065866: status FINISHED
2026-07-02 13:11:21 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:11:21 INFO None 5065869: status FINISHED
2026-07-02 13:11:21 INFO None 5065870: status RUNNING/PENDING
2026-07-02 13:11:21 INFO None 5065871: status FINISHED
2026-07-02 13:11:21 INFO None 5065872: status FINISHED
2026-07-02 13:11:21 INFO None 5065873: status FINISHED
2026-07-02 13:11:21 INFO None 5065874: status FINISHED
2026-07-02 13:11:21 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:11:21 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:11:21 INFO None 5065882: status RUNNING/PENDING
2026-07-02 13:11:21 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:11:21 INFO Jobs still running: ['5065868', '5065870', '5065876', '5065878', '5065882', '5065884']. Waiting...
2026-07-02 13:11:36 INFO None 5065863: status FINISHED
2026-07-02 13:11:36 INFO None 5065864: status FINISHED
2026-07-02 13:11:36 INFO None 5065865: status FINISHED
2026-07-02 13:11:36 INFO None 5065866: status FINISHED
2026-07-02 13:11:36 INFO None 5065868: status RUNNING/PENDING
2026-07-02 13:11:36 INFO None 5065869: status FINISHED
2026-07-02 13:11:36 INFO None 5065870: status FINISHED
2026-07-02 13:11:36 INFO None 5065871: status FINISHED
2026-07-02 13:11:36 INFO None 5065872: status FINISHED
2026-07-02 13:11:36 INFO None 5065873: status FINISHED
2026-07-02 13:11:36 INFO None 5065874: status FINISHED
2026-07-02 13:11:36 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:11:36 INFO None 5065878: status RUNNING/PENDING
2026-07-02 13:11:36 INFO None 5065882: status FINISHED
2026-07-02 13:11:36 INFO None 5065884: status RUNNING/PENDING
2026-07-02 13:11:36 INFO Jobs still running: ['5065868', '5065876', '5065878', '5065884']. Waiting...
2026-07-02 13:11:51 INFO None 5065863: status FINISHED
2026-07-02 13:11:51 INFO None 5065864: status FINISHED
2026-07-02 13:11:51 INFO None 5065865: status FINISHED
2026-07-02 13:11:51 INFO None 5065866: status FINISHED
2026-07-02 13:11:51 INFO None 5065868: status FINISHED
2026-07-02 13:11:51 INFO None 5065869: status FINISHED
2026-07-02 13:11:51 INFO None 5065870: status FINISHED
2026-07-02 13:11:51 INFO None 5065871: status FINISHED
2026-07-02 13:11:51 INFO None 5065872: status FINISHED
2026-07-02 13:11:51 INFO None 5065873: status FINISHED
2026-07-02 13:11:51 INFO None 5065874: status FINISHED
2026-07-02 13:11:51 INFO None 5065876: status RUNNING/PENDING
2026-07-02 13:11:51 INFO None 5065878: status FINISHED
2026-07-02 13:11:51 INFO None 5065882: status FINISHED
2026-07-02 13:11:51 INFO None 5065884: status FINISHED
2026-07-02 13:11:51 INFO Jobs still running: ['5065876']. Waiting...
2026-07-02 13:12:06 INFO None 5065863: status FINISHED
2026-07-02 13:12:06 INFO None 5065864: status FINISHED
2026-07-02 13:12:06 INFO None 5065865: status FINISHED
2026-07-02 13:12:06 INFO None 5065866: status FINISHED
2026-07-02 13:12:06 INFO None 5065868: status FINISHED
2026-07-02 13:12:06 INFO None 5065869: status FINISHED
2026-07-02 13:12:06 INFO None 5065870: status FINISHED
2026-07-02 13:12:06 INFO None 5065871: status FINISHED
2026-07-02 13:12:06 INFO None 5065872: status FINISHED
2026-07-02 13:12:06 INFO None 5065873: status FINISHED
2026-07-02 13:12:06 INFO None 5065874: status FINISHED
2026-07-02 13:12:06 INFO None 5065876: status FINISHED
2026-07-02 13:12:06 INFO None 5065878: status FINISHED
2026-07-02 13:12:06 INFO None 5065882: status FINISHED
2026-07-02 13:12:07 INFO None 5065884: status FINISHED
2026-07-02 13:12:07 INFO Jobs ['5065863', '5065864', '5065865', '5065866', '5065868', '5065869', '5065870', '5065871', '5065872', '5065873', '5065874', '5065876', '5065878', '5065882', '5065884'] have finished
2026-07-02 13:12:07 INFO Checking restart files were created ...
2026-07-02 13:12:07 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020601_8_ENS1.nc(3005806795 bytes)
2026-07-02 13:12:07 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020601_8_ENS2.nc(3005806795 bytes)
2026-07-02 13:12:07 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020601_8_ENS3.nc(3005806795 bytes)
2026-07-02 13:12:07 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020601_8_ENS4.nc(3005806795 bytes)
2026-07-02 13:12:07 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020601_8_ENS5.nc(3005806795 bytes)
2026-07-02 13:12:07 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020601_8_ENS6.nc(3005806795 bytes)
2026-07-02 13:12:07 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020601_8_ENS7.nc(3005806795 bytes)
2026-07-02 13:12:07 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020601_8_ENS8.nc(3005806795 bytes)
2026-07-02 13:12:07 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020601_8_ENS9.nc(3005806795 bytes)
2026-07-02 13:12:07 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020601_8_ENS10.nc(3005806795 bytes)
2026-07-02 13:12:07 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020601_8_ENS11.nc(3005806795 bytes)
2026-07-02 13:12:07 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020601_8_ENS12.nc(3005806795 bytes)
2026-07-02 13:12:07 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020601_8_ENS13.nc(3005806795 bytes)
2026-07-02 13:12:07 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020601_8_ENS14.nc(3005806795 bytes)
2026-07-02 13:12:07 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020601_8_ENS15.nc(3005806795 bytes)
2026-07-02 13:12:07 INFO  Run_model() completed successfully.
2026-07-02 13:12:07 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 13:12:07 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 09:00:00 days=153072 seconds=32400
2026-07-02 13:12:07 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-02 13:12:07 INFO [TIME] increment current_time 2020-02-06 01:00:00 -> 2020-02-06 09:00:00
2026-07-02 13:12:07 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 13:12:07 INFO ---------->>> Running process_satellite_data()
2026-07-02 13:12:07 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12002.nc
2026-07-02 13:12:07 INFO ---------->>> Running run_obs_converter()
2026-07-02 13:12:07 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_33067_153072.out
2026-07-02 13:12:07 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_33067_153072.out
2026-07-02 13:12:07 INFO ---------->>> Running DART
2026-07-02 13:12:07 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-02 13:12:07 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/chim_ENS1_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/chim_ENS1_2020020609_1_out_toDART.nc
2026-07-02 13:12:07 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/chim_ENS2_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/chim_ENS2_2020020609_1_out_toDART.nc
2026-07-02 13:12:07 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/chim_ENS3_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/chim_ENS3_2020020609_1_out_toDART.nc
2026-07-02 13:12:08 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/chim_ENS4_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/chim_ENS4_2020020609_1_out_toDART.nc
2026-07-02 13:12:08 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/chim_ENS5_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/chim_ENS5_2020020609_1_out_toDART.nc
2026-07-02 13:12:09 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/chim_ENS6_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/chim_ENS6_2020020609_1_out_toDART.nc
2026-07-02 13:12:09 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/chim_ENS7_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/chim_ENS7_2020020609_1_out_toDART.nc
2026-07-02 13:12:09 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/chim_ENS8_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/chim_ENS8_2020020609_1_out_toDART.nc
2026-07-02 13:12:10 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/chim_ENS9_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/chim_ENS9_2020020609_1_out_toDART.nc
2026-07-02 13:12:10 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/chim_ENS10_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/chim_ENS10_2020020609_1_out_toDART.nc
2026-07-02 13:12:10 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/chim_ENS11_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/chim_ENS11_2020020609_1_out_toDART.nc
2026-07-02 13:12:11 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/chim_ENS12_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/chim_ENS12_2020020609_1_out_toDART.nc
2026-07-02 13:12:11 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/chim_ENS13_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/chim_ENS13_2020020609_1_out_toDART.nc
2026-07-02 13:12:12 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/chim_ENS14_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/chim_ENS14_2020020609_1_out_toDART.nc
2026-07-02 13:12:12 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/chim_ENS15_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/chim_ENS15_2020020609_1_out_toDART.nc
2026-07-02 13:12:12 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-02 13:12:12 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-02 13:12:12 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-02 13:12:12 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-02 13:12:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-02 13:12:12 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-02 13:12:20 INFO Found: []
2026-07-02 13:12:20 INFO No job id returned by command ./run_filter.bsh
2026-07-02 13:12:20 INFO No monitoring will be performed
2026-07-02 13:12:20 INFO Moving DART output files to analysis and preassim directories for date 2020020609 if present ...
2026-07-02 13:12:20 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-07-02 13:12:20 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-07-02 13:12:20 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-07-02 13:12:20 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-07-02 13:12:20 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-07-02 13:12:20 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-07-02 13:12:20 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-07-02 13:12:20 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-07-02 13:12:20 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-07-02 13:12:20 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-07-02 13:12:20 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-07-02 13:12:20 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-07-02 13:12:20 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-07-02 13:12:20 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-07-02 13:12:20 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-07-02 13:12:20 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-07-02 13:12:20 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-07-02 13:12:20 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-07-02 13:12:20 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-07-02 13:12:20 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-07-02 13:12:20 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-07-02 13:12:20 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-07-02 13:12:20 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-07-02 13:12:20 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-07-02 13:12:20 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-07-02 13:12:20 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-07-02 13:12:20 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-07-02 13:12:20 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-07-02 13:12:20 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-07-02 13:12:20 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-07-02 13:12:20 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-07-02 13:12:20 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-07-02 13:12:20 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-07-02 13:12:20 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-07-02 13:12:20 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-02 13:12:23 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-02 13:12:23 INFO run_dart() is DONE.
2026-07-02 13:12:23 INFO ---------->>> Running update_pollutant_in_end()
2026-07-02 13:12:23 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:12:23 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:12:24 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:12:24 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:12:25 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:12:25 INFO [Posterior Diff ENS1] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS1_2020020609.nc
2026-07-02 13:12:25 INFO [Posterior Diff ENS1] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS1_2020020609.relative.nc
2026-07-02 13:12:25 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:12:25 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:12:26 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:12:26 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:12:27 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:12:27 INFO [Posterior Diff ENS2] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS2_2020020609.nc
2026-07-02 13:12:27 INFO [Posterior Diff ENS2] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS2_2020020609.relative.nc
2026-07-02 13:12:27 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:12:27 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:12:28 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:12:28 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:12:29 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:12:29 INFO [Posterior Diff ENS3] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS3_2020020609.nc
2026-07-02 13:12:29 INFO [Posterior Diff ENS3] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS3_2020020609.relative.nc
2026-07-02 13:12:29 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:12:29 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:12:30 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:12:30 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:12:31 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:12:31 INFO [Posterior Diff ENS4] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS4_2020020609.nc
2026-07-02 13:12:31 INFO [Posterior Diff ENS4] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS4_2020020609.relative.nc
2026-07-02 13:12:31 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:12:31 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:12:32 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:12:32 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:12:33 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:12:33 INFO [Posterior Diff ENS5] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS5_2020020609.nc
2026-07-02 13:12:33 INFO [Posterior Diff ENS5] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS5_2020020609.relative.nc
2026-07-02 13:12:33 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:12:33 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:12:34 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:12:34 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:12:35 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:12:35 INFO [Posterior Diff ENS6] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS6_2020020609.nc
2026-07-02 13:12:35 INFO [Posterior Diff ENS6] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS6_2020020609.relative.nc
2026-07-02 13:12:36 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:12:36 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:12:36 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:12:36 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:12:37 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:12:37 INFO [Posterior Diff ENS7] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS7_2020020609.nc
2026-07-02 13:12:37 INFO [Posterior Diff ENS7] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS7_2020020609.relative.nc
2026-07-02 13:12:38 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:12:38 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:12:39 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:12:39 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:12:39 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:12:39 INFO [Posterior Diff ENS8] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS8_2020020609.nc
2026-07-02 13:12:39 INFO [Posterior Diff ENS8] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS8_2020020609.relative.nc
2026-07-02 13:12:40 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:12:40 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:12:41 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:12:41 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:12:41 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:12:41 INFO [Posterior Diff ENS9] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS9_2020020609.nc
2026-07-02 13:12:41 INFO [Posterior Diff ENS9] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS9_2020020609.relative.nc
2026-07-02 13:12:42 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:12:42 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:12:43 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:12:43 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:12:43 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:12:43 INFO [Posterior Diff ENS10] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS10_2020020609.nc
2026-07-02 13:12:43 INFO [Posterior Diff ENS10] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS10_2020020609.relative.nc
2026-07-02 13:12:44 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:12:44 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:12:45 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:12:45 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:12:46 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:12:46 INFO [Posterior Diff ENS11] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS11_2020020609.nc
2026-07-02 13:12:46 INFO [Posterior Diff ENS11] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS11_2020020609.relative.nc
2026-07-02 13:12:46 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:12:46 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:12:47 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:12:47 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:12:48 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:12:48 INFO [Posterior Diff ENS12] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS12_2020020609.nc
2026-07-02 13:12:48 INFO [Posterior Diff ENS12] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS12_2020020609.relative.nc
2026-07-02 13:12:48 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:12:48 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:12:49 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:12:49 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:12:50 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:12:50 INFO [Posterior Diff ENS13] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS13_2020020609.nc
2026-07-02 13:12:50 INFO [Posterior Diff ENS13] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS13_2020020609.relative.nc
2026-07-02 13:12:50 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:12:51 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:12:51 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:12:51 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:12:52 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:12:52 INFO [Posterior Diff ENS14] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS14_2020020609.nc
2026-07-02 13:12:52 INFO [Posterior Diff ENS14] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS14_2020020609.relative.nc
2026-07-02 13:12:53 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:12:53 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:12:54 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:12:54 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:12:54 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:12:55 INFO [Posterior Diff ENS15] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS15_2020020609.nc
2026-07-02 13:12:55 INFO [Posterior Diff ENS15] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS15_2020020609.relative.nc
2026-07-02 13:12:55 INFO Next run starts from 2020-02-06 09:00:00
2026-07-02 13:12:55 INFO Cycle is DONE; starting a new loop!
2026-07-02 13:12:55 INFO [TIME] step_end current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 13:12:55 INFO [TIME] step_start current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 13:12:55 INFO [TIME] window start=2020-02-06 09:00:00 end=2020-02-06 11:00:00 run_hours=2 has_assimilation=True
2026-07-02 13:12:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:12:56 INFO Hourly dataset computed and listing created
2026-07-02 13:12:58 INFO Hourly dataset computed
2026-07-02 13:12:58 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:12:59 INFO Hourly dataset computed and listing created
2026-07-02 13:13:00 INFO Hourly dataset computed
2026-07-02 13:13:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:13:01 INFO Hourly dataset computed and listing created
2026-07-02 13:13:02 INFO Hourly dataset computed
2026-07-02 13:13:02 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:13:03 INFO Hourly dataset computed and listing created
2026-07-02 13:13:03 INFO Hourly dataset computed
2026-07-02 13:13:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:13:04 INFO Hourly dataset computed and listing created
2026-07-02 13:13:05 INFO Hourly dataset computed
2026-07-02 13:13:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:13:06 INFO Hourly dataset computed and listing created
2026-07-02 13:13:07 INFO Hourly dataset computed
2026-07-02 13:13:07 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:13:08 INFO Hourly dataset computed and listing created
2026-07-02 13:13:09 INFO Hourly dataset computed
2026-07-02 13:13:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:13:10 INFO Hourly dataset computed and listing created
2026-07-02 13:13:10 INFO Hourly dataset computed
2026-07-02 13:13:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:13:11 INFO Hourly dataset computed and listing created
2026-07-02 13:13:12 INFO Hourly dataset computed
2026-07-02 13:13:12 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:13:13 INFO Hourly dataset computed and listing created
2026-07-02 13:13:14 INFO Hourly dataset computed
2026-07-02 13:13:14 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:13:15 INFO Hourly dataset computed and listing created
2026-07-02 13:13:16 INFO Hourly dataset computed
2026-07-02 13:13:16 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:13:16 INFO Hourly dataset computed and listing created
2026-07-02 13:13:17 INFO Hourly dataset computed
2026-07-02 13:13:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:13:18 INFO Hourly dataset computed and listing created
2026-07-02 13:13:19 INFO Hourly dataset computed
2026-07-02 13:13:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:13:20 INFO Hourly dataset computed and listing created
2026-07-02 13:13:21 INFO Hourly dataset computed
2026-07-02 13:13:21 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:13:22 INFO Hourly dataset computed and listing created
2026-07-02 13:13:23 INFO Hourly dataset computed
2026-07-02 13:13:23 INFO ---------->>> Running CHIMERE model from 2020-02-06 09:00:00 to 2020-02-06 11:00:00
2026-07-02 13:13:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:13:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1
2026-07-02 13:13:23 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020601_8_ENS1.nc
2026-07-02 13:13:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-02 13:13:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:13:23 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-02 13:13:23 INFO Queuing job for member 1...
2026-07-02 13:13:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:13:23 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-02 13:13:23 INFO Found: ['5065928']
2026-07-02 13:13:28 INFO [TGCC-IRENE] Submitted job with ID:['5065928']
2026-07-02 13:13:28 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:13:28 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2
2026-07-02 13:13:28 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020601_8_ENS2.nc
2026-07-02 13:13:28 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-02 13:13:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:13:28 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-02 13:13:28 INFO Queuing job for member 2...
2026-07-02 13:13:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:13:28 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-02 13:13:29 INFO Found: ['5065929']
2026-07-02 13:13:34 INFO [TGCC-IRENE] Submitted job with ID:['5065929']
2026-07-02 13:13:34 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:13:34 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3
2026-07-02 13:13:34 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020601_8_ENS3.nc
2026-07-02 13:13:34 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-02 13:13:34 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:13:34 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-02 13:13:34 INFO Queuing job for member 3...
2026-07-02 13:13:34 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:13:34 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-02 13:13:35 INFO Found: ['5065930']
2026-07-02 13:13:40 INFO [TGCC-IRENE] Submitted job with ID:['5065930']
2026-07-02 13:13:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:13:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4
2026-07-02 13:13:40 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020601_8_ENS4.nc
2026-07-02 13:13:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-02 13:13:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:13:40 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-02 13:13:40 INFO Queuing job for member 4...
2026-07-02 13:13:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:13:40 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-02 13:13:41 INFO Found: ['5065932']
2026-07-02 13:13:46 INFO [TGCC-IRENE] Submitted job with ID:['5065932']
2026-07-02 13:13:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:13:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5
2026-07-02 13:13:46 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020601_8_ENS5.nc
2026-07-02 13:13:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-02 13:13:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:13:46 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-02 13:13:46 INFO Queuing job for member 5...
2026-07-02 13:13:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:13:46 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-02 13:13:46 INFO Found: ['5065933']
2026-07-02 13:13:51 INFO [TGCC-IRENE] Submitted job with ID:['5065933']
2026-07-02 13:13:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:13:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6
2026-07-02 13:13:51 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020601_8_ENS6.nc
2026-07-02 13:13:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-02 13:13:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:13:51 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-02 13:13:51 INFO Queuing job for member 6...
2026-07-02 13:13:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:13:51 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-02 13:13:52 INFO Found: ['5065937']
2026-07-02 13:13:57 INFO [TGCC-IRENE] Submitted job with ID:['5065937']
2026-07-02 13:13:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:13:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7
2026-07-02 13:13:57 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020601_8_ENS7.nc
2026-07-02 13:13:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-02 13:13:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:13:57 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-02 13:13:57 INFO Queuing job for member 7...
2026-07-02 13:13:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:13:57 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-02 13:13:58 INFO Found: ['5065940']
2026-07-02 13:14:03 INFO [TGCC-IRENE] Submitted job with ID:['5065940']
2026-07-02 13:14:03 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:14:03 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8
2026-07-02 13:14:03 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020601_8_ENS8.nc
2026-07-02 13:14:03 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-02 13:14:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:14:03 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-02 13:14:03 INFO Queuing job for member 8...
2026-07-02 13:14:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:14:03 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-02 13:14:03 INFO Found: ['5065944']
2026-07-02 13:14:08 INFO [TGCC-IRENE] Submitted job with ID:['5065944']
2026-07-02 13:14:08 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:14:08 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9
2026-07-02 13:14:08 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020601_8_ENS9.nc
2026-07-02 13:14:08 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-02 13:14:08 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:14:08 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-02 13:14:08 INFO Queuing job for member 9...
2026-07-02 13:14:08 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:14:08 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-02 13:14:11 INFO Found: ['5065948']
2026-07-02 13:14:16 INFO [TGCC-IRENE] Submitted job with ID:['5065948']
2026-07-02 13:14:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:14:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10
2026-07-02 13:14:16 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020601_8_ENS10.nc
2026-07-02 13:14:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-02 13:14:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:14:16 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-02 13:14:16 INFO Queuing job for member 10...
2026-07-02 13:14:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:14:16 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-02 13:14:17 INFO Found: ['5065953']
2026-07-02 13:14:22 INFO [TGCC-IRENE] Submitted job with ID:['5065953']
2026-07-02 13:14:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:14:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11
2026-07-02 13:14:22 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020601_8_ENS11.nc
2026-07-02 13:14:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-02 13:14:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:14:22 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-02 13:14:22 INFO Queuing job for member 11...
2026-07-02 13:14:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:14:22 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-02 13:14:23 INFO Found: ['5065957']
2026-07-02 13:14:28 INFO [TGCC-IRENE] Submitted job with ID:['5065957']
2026-07-02 13:14:28 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:14:28 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12
2026-07-02 13:14:28 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020601_8_ENS12.nc
2026-07-02 13:14:28 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-02 13:14:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:14:28 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-02 13:14:28 INFO Queuing job for member 12...
2026-07-02 13:14:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:14:28 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-02 13:14:29 INFO Found: ['5065961']
2026-07-02 13:14:34 INFO [TGCC-IRENE] Submitted job with ID:['5065961']
2026-07-02 13:14:34 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:14:34 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13
2026-07-02 13:14:34 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020601_8_ENS13.nc
2026-07-02 13:14:34 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-02 13:14:34 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:14:34 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-02 13:14:34 INFO Queuing job for member 13...
2026-07-02 13:14:34 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:14:34 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-02 13:14:34 INFO Found: ['5065963']
2026-07-02 13:14:39 INFO [TGCC-IRENE] Submitted job with ID:['5065963']
2026-07-02 13:14:39 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:14:39 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14
2026-07-02 13:14:39 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020601_8_ENS14.nc
2026-07-02 13:14:39 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-02 13:14:39 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:14:39 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-02 13:14:39 INFO Queuing job for member 14...
2026-07-02 13:14:39 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:14:39 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-02 13:14:40 INFO Found: ['5065965']
2026-07-02 13:14:45 INFO [TGCC-IRENE] Submitted job with ID:['5065965']
2026-07-02 13:14:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:14:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15
2026-07-02 13:14:45 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020601_8_ENS15.nc
2026-07-02 13:14:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-02 13:14:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:14:45 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-02 13:14:45 INFO Queuing job for member 15...
2026-07-02 13:14:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:14:45 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-02 13:14:46 INFO Found: ['5065966']
2026-07-02 13:14:51 INFO [TGCC-IRENE] Submitted job with ID:['5065966']
2026-07-02 13:14:51 INFO Checking job status ...
2026-07-02 13:14:51 INFO None 5065928: status RUNNING/PENDING
2026-07-02 13:14:51 INFO None 5065929: status RUNNING/PENDING
2026-07-02 13:14:51 INFO None 5065930: status RUNNING/PENDING
2026-07-02 13:14:51 INFO None 5065932: status RUNNING/PENDING
2026-07-02 13:14:51 INFO None 5065933: status RUNNING/PENDING
2026-07-02 13:14:51 INFO None 5065937: status RUNNING/PENDING
2026-07-02 13:14:51 INFO None 5065940: status RUNNING/PENDING
2026-07-02 13:14:51 INFO None 5065944: status RUNNING/PENDING
2026-07-02 13:14:51 INFO None 5065948: status RUNNING/PENDING
2026-07-02 13:14:51 INFO None 5065953: status RUNNING/PENDING
2026-07-02 13:14:51 INFO None 5065957: status RUNNING/PENDING
2026-07-02 13:14:51 INFO None 5065961: status RUNNING/PENDING
2026-07-02 13:14:51 INFO None 5065963: status RUNNING/PENDING
2026-07-02 13:14:51 INFO None 5065965: status RUNNING/PENDING
2026-07-02 13:14:51 INFO None 5065966: status RUNNING/PENDING
2026-07-02 13:14:51 INFO Jobs still running: ['5065928', '5065929', '5065930', '5065932', '5065933', '5065937', '5065940', '5065944', '5065948', '5065953', '5065957', '5065961', '5065963', '5065965', '5065966']. Waiting...
2026-07-02 13:15:07 INFO None 5065928: status RUNNING/PENDING
2026-07-02 13:15:07 INFO None 5065929: status RUNNING/PENDING
2026-07-02 13:15:07 INFO None 5065930: status RUNNING/PENDING
2026-07-02 13:15:07 INFO None 5065932: status RUNNING/PENDING
2026-07-02 13:15:07 INFO None 5065933: status RUNNING/PENDING
2026-07-02 13:15:07 INFO None 5065937: status RUNNING/PENDING
2026-07-02 13:15:07 INFO None 5065940: status RUNNING/PENDING
2026-07-02 13:15:07 INFO None 5065944: status RUNNING/PENDING
2026-07-02 13:15:07 INFO None 5065948: status RUNNING/PENDING
2026-07-02 13:15:07 INFO None 5065953: status RUNNING/PENDING
2026-07-02 13:15:07 INFO None 5065957: status RUNNING/PENDING
2026-07-02 13:15:07 INFO None 5065961: status RUNNING/PENDING
2026-07-02 13:15:07 INFO None 5065963: status RUNNING/PENDING
2026-07-02 13:15:07 INFO None 5065965: status RUNNING/PENDING
2026-07-02 13:15:07 INFO None 5065966: status RUNNING/PENDING
2026-07-02 13:15:07 INFO Jobs still running: ['5065928', '5065929', '5065930', '5065932', '5065933', '5065937', '5065940', '5065944', '5065948', '5065953', '5065957', '5065961', '5065963', '5065965', '5065966']. Waiting...
2026-07-02 13:15:22 INFO None 5065928: status RUNNING/PENDING
2026-07-02 13:15:22 INFO None 5065929: status RUNNING/PENDING
2026-07-02 13:15:22 INFO None 5065930: status RUNNING/PENDING
2026-07-02 13:15:22 INFO None 5065932: status RUNNING/PENDING
2026-07-02 13:15:22 INFO None 5065933: status RUNNING/PENDING
2026-07-02 13:15:22 INFO None 5065937: status RUNNING/PENDING
2026-07-02 13:15:22 INFO None 5065940: status RUNNING/PENDING
2026-07-02 13:15:22 INFO None 5065944: status RUNNING/PENDING
2026-07-02 13:15:22 INFO None 5065948: status RUNNING/PENDING
2026-07-02 13:15:22 INFO None 5065953: status RUNNING/PENDING
2026-07-02 13:15:22 INFO None 5065957: status RUNNING/PENDING
2026-07-02 13:15:22 INFO None 5065961: status RUNNING/PENDING
2026-07-02 13:15:22 INFO None 5065963: status RUNNING/PENDING
2026-07-02 13:15:22 INFO None 5065965: status RUNNING/PENDING
2026-07-02 13:15:22 INFO None 5065966: status RUNNING/PENDING
2026-07-02 13:15:22 INFO Jobs still running: ['5065928', '5065929', '5065930', '5065932', '5065933', '5065937', '5065940', '5065944', '5065948', '5065953', '5065957', '5065961', '5065963', '5065965', '5065966']. Waiting...
2026-07-02 13:15:37 INFO None 5065928: status RUNNING/PENDING
2026-07-02 13:15:37 INFO None 5065929: status RUNNING/PENDING
2026-07-02 13:15:37 INFO None 5065930: status RUNNING/PENDING
2026-07-02 13:15:37 INFO None 5065932: status RUNNING/PENDING
2026-07-02 13:15:37 INFO None 5065933: status RUNNING/PENDING
2026-07-02 13:15:37 INFO None 5065937: status RUNNING/PENDING
2026-07-02 13:15:37 INFO None 5065940: status RUNNING/PENDING
2026-07-02 13:15:37 INFO None 5065944: status RUNNING/PENDING
2026-07-02 13:15:37 INFO None 5065948: status RUNNING/PENDING
2026-07-02 13:15:37 INFO None 5065953: status RUNNING/PENDING
2026-07-02 13:15:37 INFO None 5065957: status RUNNING/PENDING
2026-07-02 13:15:37 INFO None 5065961: status RUNNING/PENDING
2026-07-02 13:15:37 INFO None 5065963: status RUNNING/PENDING
2026-07-02 13:15:37 INFO None 5065965: status RUNNING/PENDING
2026-07-02 13:15:37 INFO None 5065966: status RUNNING/PENDING
2026-07-02 13:15:37 INFO Jobs still running: ['5065928', '5065929', '5065930', '5065932', '5065933', '5065937', '5065940', '5065944', '5065948', '5065953', '5065957', '5065961', '5065963', '5065965', '5065966']. Waiting...
2026-07-02 13:15:52 INFO None 5065928: status RUNNING/PENDING
2026-07-02 13:15:52 INFO None 5065929: status RUNNING/PENDING
2026-07-02 13:15:52 INFO None 5065930: status RUNNING/PENDING
2026-07-02 13:15:52 INFO None 5065932: status RUNNING/PENDING
2026-07-02 13:15:52 INFO None 5065933: status RUNNING/PENDING
2026-07-02 13:15:53 INFO None 5065937: status RUNNING/PENDING
2026-07-02 13:15:53 INFO None 5065940: status RUNNING/PENDING
2026-07-02 13:15:53 INFO None 5065944: status RUNNING/PENDING
2026-07-02 13:15:53 INFO None 5065948: status RUNNING/PENDING
2026-07-02 13:15:53 INFO None 5065953: status RUNNING/PENDING
2026-07-02 13:15:53 INFO None 5065957: status RUNNING/PENDING
2026-07-02 13:15:53 INFO None 5065961: status RUNNING/PENDING
2026-07-02 13:15:53 INFO None 5065963: status RUNNING/PENDING
2026-07-02 13:15:53 INFO None 5065965: status RUNNING/PENDING
2026-07-02 13:15:53 INFO None 5065966: status RUNNING/PENDING
2026-07-02 13:15:53 INFO Jobs still running: ['5065928', '5065929', '5065930', '5065932', '5065933', '5065937', '5065940', '5065944', '5065948', '5065953', '5065957', '5065961', '5065963', '5065965', '5065966']. Waiting...
2026-07-02 13:16:08 INFO None 5065928: status RUNNING/PENDING
2026-07-02 13:16:08 INFO None 5065929: status RUNNING/PENDING
2026-07-02 13:16:08 INFO None 5065930: status RUNNING/PENDING
2026-07-02 13:16:08 INFO None 5065932: status RUNNING/PENDING
2026-07-02 13:16:08 INFO None 5065933: status RUNNING/PENDING
2026-07-02 13:16:08 INFO None 5065937: status RUNNING/PENDING
2026-07-02 13:16:08 INFO None 5065940: status RUNNING/PENDING
2026-07-02 13:16:08 INFO None 5065944: status RUNNING/PENDING
2026-07-02 13:16:08 INFO None 5065948: status RUNNING/PENDING
2026-07-02 13:16:08 INFO None 5065953: status RUNNING/PENDING
2026-07-02 13:16:08 INFO None 5065957: status RUNNING/PENDING
2026-07-02 13:16:08 INFO None 5065961: status RUNNING/PENDING
2026-07-02 13:16:08 INFO None 5065963: status RUNNING/PENDING
2026-07-02 13:16:08 INFO None 5065965: status RUNNING/PENDING
2026-07-02 13:16:08 INFO None 5065966: status RUNNING/PENDING
2026-07-02 13:16:08 INFO Jobs still running: ['5065928', '5065929', '5065930', '5065932', '5065933', '5065937', '5065940', '5065944', '5065948', '5065953', '5065957', '5065961', '5065963', '5065965', '5065966']. Waiting...
2026-07-02 13:16:23 INFO None 5065928: status RUNNING/PENDING
2026-07-02 13:16:23 INFO None 5065929: status RUNNING/PENDING
2026-07-02 13:16:23 INFO None 5065930: status RUNNING/PENDING
2026-07-02 13:16:23 INFO None 5065932: status RUNNING/PENDING
2026-07-02 13:16:23 INFO None 5065933: status RUNNING/PENDING
2026-07-02 13:16:23 INFO None 5065937: status RUNNING/PENDING
2026-07-02 13:16:23 INFO None 5065940: status RUNNING/PENDING
2026-07-02 13:16:23 INFO None 5065944: status RUNNING/PENDING
2026-07-02 13:16:23 INFO None 5065948: status RUNNING/PENDING
2026-07-02 13:16:23 INFO None 5065953: status RUNNING/PENDING
2026-07-02 13:16:23 INFO None 5065957: status RUNNING/PENDING
2026-07-02 13:16:23 INFO None 5065961: status RUNNING/PENDING
2026-07-02 13:16:23 INFO None 5065963: status RUNNING/PENDING
2026-07-02 13:16:23 INFO None 5065965: status RUNNING/PENDING
2026-07-02 13:16:23 INFO None 5065966: status RUNNING/PENDING
2026-07-02 13:16:23 INFO Jobs still running: ['5065928', '5065929', '5065930', '5065932', '5065933', '5065937', '5065940', '5065944', '5065948', '5065953', '5065957', '5065961', '5065963', '5065965', '5065966']. Waiting...
2026-07-02 13:16:38 INFO None 5065928: status RUNNING/PENDING
2026-07-02 13:16:39 INFO None 5065929: status RUNNING/PENDING
2026-07-02 13:16:39 INFO None 5065930: status RUNNING/PENDING
2026-07-02 13:16:39 INFO None 5065932: status RUNNING/PENDING
2026-07-02 13:16:39 INFO None 5065933: status RUNNING/PENDING
2026-07-02 13:16:39 INFO None 5065937: status RUNNING/PENDING
2026-07-02 13:16:39 INFO None 5065940: status RUNNING/PENDING
2026-07-02 13:16:39 INFO None 5065944: status RUNNING/PENDING
2026-07-02 13:16:39 INFO None 5065948: status RUNNING/PENDING
2026-07-02 13:16:39 INFO None 5065953: status RUNNING/PENDING
2026-07-02 13:16:39 INFO None 5065957: status RUNNING/PENDING
2026-07-02 13:16:39 INFO None 5065961: status RUNNING/PENDING
2026-07-02 13:16:39 INFO None 5065963: status RUNNING/PENDING
2026-07-02 13:16:39 INFO None 5065965: status RUNNING/PENDING
2026-07-02 13:16:39 INFO None 5065966: status RUNNING/PENDING
2026-07-02 13:16:39 INFO Jobs still running: ['5065928', '5065929', '5065930', '5065932', '5065933', '5065937', '5065940', '5065944', '5065948', '5065953', '5065957', '5065961', '5065963', '5065965', '5065966']. Waiting...
2026-07-02 13:16:54 INFO None 5065928: status FINISHED
2026-07-02 13:16:54 INFO None 5065929: status RUNNING/PENDING
2026-07-02 13:16:54 INFO None 5065930: status RUNNING/PENDING
2026-07-02 13:16:54 INFO None 5065932: status FINISHED
2026-07-02 13:16:54 INFO None 5065933: status FINISHED
2026-07-02 13:16:54 INFO None 5065937: status RUNNING/PENDING
2026-07-02 13:16:54 INFO None 5065940: status RUNNING/PENDING
2026-07-02 13:16:54 INFO None 5065944: status RUNNING/PENDING
2026-07-02 13:16:54 INFO None 5065948: status RUNNING/PENDING
2026-07-02 13:16:54 INFO None 5065953: status RUNNING/PENDING
2026-07-02 13:16:54 INFO None 5065957: status RUNNING/PENDING
2026-07-02 13:16:54 INFO None 5065961: status RUNNING/PENDING
2026-07-02 13:16:54 INFO None 5065963: status RUNNING/PENDING
2026-07-02 13:16:54 INFO None 5065965: status RUNNING/PENDING
2026-07-02 13:16:54 INFO None 5065966: status RUNNING/PENDING
2026-07-02 13:16:54 INFO Jobs still running: ['5065929', '5065930', '5065937', '5065940', '5065944', '5065948', '5065953', '5065957', '5065961', '5065963', '5065965', '5065966']. Waiting...
2026-07-02 13:17:09 INFO None 5065928: status FINISHED
2026-07-02 13:17:09 INFO None 5065929: status FINISHED
2026-07-02 13:17:09 INFO None 5065930: status RUNNING/PENDING
2026-07-02 13:17:10 INFO None 5065932: status FINISHED
2026-07-02 13:17:10 INFO None 5065933: status FINISHED
2026-07-02 13:17:10 INFO None 5065937: status FINISHED
2026-07-02 13:17:10 INFO None 5065940: status RUNNING/PENDING
2026-07-02 13:17:10 INFO None 5065944: status RUNNING/PENDING
2026-07-02 13:17:10 INFO None 5065948: status RUNNING/PENDING
2026-07-02 13:17:10 INFO None 5065953: status RUNNING/PENDING
2026-07-02 13:17:10 INFO None 5065957: status RUNNING/PENDING
2026-07-02 13:17:10 INFO None 5065961: status RUNNING/PENDING
2026-07-02 13:17:10 INFO None 5065963: status RUNNING/PENDING
2026-07-02 13:17:10 INFO None 5065965: status RUNNING/PENDING
2026-07-02 13:17:10 INFO None 5065966: status RUNNING/PENDING
2026-07-02 13:17:10 INFO Jobs still running: ['5065930', '5065940', '5065944', '5065948', '5065953', '5065957', '5065961', '5065963', '5065965', '5065966']. Waiting...
2026-07-02 13:17:25 INFO None 5065928: status FINISHED
2026-07-02 13:17:25 INFO None 5065929: status FINISHED
2026-07-02 13:17:25 INFO None 5065930: status RUNNING/PENDING
2026-07-02 13:17:25 INFO None 5065932: status FINISHED
2026-07-02 13:17:25 INFO None 5065933: status FINISHED
2026-07-02 13:17:25 INFO None 5065937: status FINISHED
2026-07-02 13:17:25 INFO None 5065940: status FINISHED
2026-07-02 13:17:25 INFO None 5065944: status FINISHED
2026-07-02 13:17:25 INFO None 5065948: status RUNNING/PENDING
2026-07-02 13:17:25 INFO None 5065953: status RUNNING/PENDING
2026-07-02 13:17:25 INFO None 5065957: status RUNNING/PENDING
2026-07-02 13:17:25 INFO None 5065961: status RUNNING/PENDING
2026-07-02 13:17:25 INFO None 5065963: status RUNNING/PENDING
2026-07-02 13:17:25 INFO None 5065965: status RUNNING/PENDING
2026-07-02 13:17:25 INFO None 5065966: status RUNNING/PENDING
2026-07-02 13:17:25 INFO Jobs still running: ['5065930', '5065948', '5065953', '5065957', '5065961', '5065963', '5065965', '5065966']. Waiting...
2026-07-02 13:17:40 INFO None 5065928: status FINISHED
2026-07-02 13:17:40 INFO None 5065929: status FINISHED
2026-07-02 13:17:40 INFO None 5065930: status FINISHED
2026-07-02 13:17:40 INFO None 5065932: status FINISHED
2026-07-02 13:17:40 INFO None 5065933: status FINISHED
2026-07-02 13:17:40 INFO None 5065937: status FINISHED
2026-07-02 13:17:40 INFO None 5065940: status FINISHED
2026-07-02 13:17:40 INFO None 5065944: status FINISHED
2026-07-02 13:17:40 INFO None 5065948: status RUNNING/PENDING
2026-07-02 13:17:40 INFO None 5065953: status RUNNING/PENDING
2026-07-02 13:17:40 INFO None 5065957: status RUNNING/PENDING
2026-07-02 13:17:40 INFO None 5065961: status RUNNING/PENDING
2026-07-02 13:17:40 INFO None 5065963: status RUNNING/PENDING
2026-07-02 13:17:40 INFO None 5065965: status RUNNING/PENDING
2026-07-02 13:17:40 INFO None 5065966: status RUNNING/PENDING
2026-07-02 13:17:40 INFO Jobs still running: ['5065948', '5065953', '5065957', '5065961', '5065963', '5065965', '5065966']. Waiting...
2026-07-02 13:17:55 INFO None 5065928: status FINISHED
2026-07-02 13:17:55 INFO None 5065929: status FINISHED
2026-07-02 13:17:55 INFO None 5065930: status FINISHED
2026-07-02 13:17:55 INFO None 5065932: status FINISHED
2026-07-02 13:17:55 INFO None 5065933: status FINISHED
2026-07-02 13:17:55 INFO None 5065937: status FINISHED
2026-07-02 13:17:55 INFO None 5065940: status FINISHED
2026-07-02 13:17:55 INFO None 5065944: status FINISHED
2026-07-02 13:17:55 INFO None 5065948: status RUNNING/PENDING
2026-07-02 13:17:55 INFO None 5065953: status RUNNING/PENDING
2026-07-02 13:17:55 INFO None 5065957: status RUNNING/PENDING
2026-07-02 13:17:56 INFO None 5065961: status RUNNING/PENDING
2026-07-02 13:17:56 INFO None 5065963: status FINISHED
2026-07-02 13:17:56 INFO None 5065965: status RUNNING/PENDING
2026-07-02 13:17:56 INFO None 5065966: status FINISHED
2026-07-02 13:17:56 INFO Jobs still running: ['5065948', '5065953', '5065957', '5065961', '5065965']. Waiting...
2026-07-02 13:18:11 INFO None 5065928: status FINISHED
2026-07-02 13:18:11 INFO None 5065929: status FINISHED
2026-07-02 13:18:11 INFO None 5065930: status FINISHED
2026-07-02 13:18:11 INFO None 5065932: status FINISHED
2026-07-02 13:18:11 INFO None 5065933: status FINISHED
2026-07-02 13:18:11 INFO None 5065937: status FINISHED
2026-07-02 13:18:11 INFO None 5065940: status FINISHED
2026-07-02 13:18:11 INFO None 5065944: status FINISHED
2026-07-02 13:18:11 INFO None 5065948: status RUNNING/PENDING
2026-07-02 13:18:11 INFO None 5065953: status RUNNING/PENDING
2026-07-02 13:18:11 INFO None 5065957: status FINISHED
2026-07-02 13:18:11 INFO None 5065961: status FINISHED
2026-07-02 13:18:11 INFO None 5065963: status FINISHED
2026-07-02 13:18:11 INFO None 5065965: status FINISHED
2026-07-02 13:18:11 INFO None 5065966: status FINISHED
2026-07-02 13:18:11 INFO Jobs still running: ['5065948', '5065953']. Waiting...
2026-07-02 13:18:26 INFO None 5065928: status FINISHED
2026-07-02 13:18:26 INFO None 5065929: status FINISHED
2026-07-02 13:18:26 INFO None 5065930: status FINISHED
2026-07-02 13:18:26 INFO None 5065932: status FINISHED
2026-07-02 13:18:26 INFO None 5065933: status FINISHED
2026-07-02 13:18:26 INFO None 5065937: status FINISHED
2026-07-02 13:18:26 INFO None 5065940: status FINISHED
2026-07-02 13:18:26 INFO None 5065944: status FINISHED
2026-07-02 13:18:26 INFO None 5065948: status FINISHED
2026-07-02 13:18:26 INFO None 5065953: status FINISHED
2026-07-02 13:18:26 INFO None 5065957: status FINISHED
2026-07-02 13:18:26 INFO None 5065961: status FINISHED
2026-07-02 13:18:26 INFO None 5065963: status FINISHED
2026-07-02 13:18:26 INFO None 5065965: status FINISHED
2026-07-02 13:18:27 INFO None 5065966: status FINISHED
2026-07-02 13:18:27 INFO Jobs ['5065928', '5065929', '5065930', '5065932', '5065933', '5065937', '5065940', '5065944', '5065948', '5065953', '5065957', '5065961', '5065963', '5065965', '5065966'] have finished
2026-07-02 13:18:27 INFO Checking restart files were created ...
2026-07-02 13:18:27 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc(1002685915 bytes)
2026-07-02 13:18:27 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc(1002685915 bytes)
2026-07-02 13:18:27 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc(1002685915 bytes)
2026-07-02 13:18:27 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc(1002685915 bytes)
2026-07-02 13:18:27 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc(1002685915 bytes)
2026-07-02 13:18:27 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc(1002685915 bytes)
2026-07-02 13:18:27 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc(1002685915 bytes)
2026-07-02 13:18:27 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc(1002685915 bytes)
2026-07-02 13:18:27 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc(1002685915 bytes)
2026-07-02 13:18:27 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc(1002685915 bytes)
2026-07-02 13:18:27 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc(1002685915 bytes)
2026-07-02 13:18:27 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc(1002685915 bytes)
2026-07-02 13:18:27 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc(1002685915 bytes)
2026-07-02 13:18:27 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc(1002685915 bytes)
2026-07-02 13:18:27 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc(1002685915 bytes)
2026-07-02 13:18:27 INFO  Run_model() completed successfully.
2026-07-02 13:18:27 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 13:18:27 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 11:00:00 days=153072 seconds=39600
2026-07-02 13:18:27 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-02 13:18:27 INFO [TIME] increment current_time 2020-02-06 09:00:00 -> 2020-02-06 11:00:00
2026-07-02 13:18:27 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 13:18:27 INFO ---------->>> Running process_satellite_data()
2026-07-02 13:18:27 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12003.nc
2026-07-02 13:18:27 INFO ---------->>> Running run_obs_converter()
2026-07-02 13:18:27 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_39045_153072.out
2026-07-02 13:18:27 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_39045_153072.out
2026-07-02 13:18:27 INFO ---------->>> Running DART
2026-07-02 13:18:27 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-02 13:18:27 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/chim_ENS1_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/chim_ENS1_2020020611_1_out_toDART.nc
2026-07-02 13:18:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/chim_ENS2_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/chim_ENS2_2020020611_1_out_toDART.nc
2026-07-02 13:18:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/chim_ENS3_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/chim_ENS3_2020020611_1_out_toDART.nc
2026-07-02 13:18:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/chim_ENS4_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/chim_ENS4_2020020611_1_out_toDART.nc
2026-07-02 13:18:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/chim_ENS5_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/chim_ENS5_2020020611_1_out_toDART.nc
2026-07-02 13:18:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/chim_ENS6_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/chim_ENS6_2020020611_1_out_toDART.nc
2026-07-02 13:18:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/chim_ENS7_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/chim_ENS7_2020020611_1_out_toDART.nc
2026-07-02 13:18:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/chim_ENS8_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/chim_ENS8_2020020611_1_out_toDART.nc
2026-07-02 13:18:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/chim_ENS9_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/chim_ENS9_2020020611_1_out_toDART.nc
2026-07-02 13:18:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/chim_ENS10_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/chim_ENS10_2020020611_1_out_toDART.nc
2026-07-02 13:18:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/chim_ENS11_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/chim_ENS11_2020020611_1_out_toDART.nc
2026-07-02 13:18:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/chim_ENS12_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/chim_ENS12_2020020611_1_out_toDART.nc
2026-07-02 13:18:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/chim_ENS13_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/chim_ENS13_2020020611_1_out_toDART.nc
2026-07-02 13:18:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/chim_ENS14_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/chim_ENS14_2020020611_1_out_toDART.nc
2026-07-02 13:18:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/chim_ENS15_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/chim_ENS15_2020020611_1_out_toDART.nc
2026-07-02 13:18:32 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-02 13:18:32 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-02 13:18:32 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-02 13:18:32 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-02 13:18:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-02 13:18:32 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-02 13:18:44 INFO Found: []
2026-07-02 13:18:44 INFO No job id returned by command ./run_filter.bsh
2026-07-02 13:18:44 INFO No monitoring will be performed
2026-07-02 13:18:44 INFO Moving DART output files to analysis and preassim directories for date 2020020611 if present ...
2026-07-02 13:18:44 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-07-02 13:18:44 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-07-02 13:18:44 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-07-02 13:18:44 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-07-02 13:18:44 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-07-02 13:18:44 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-07-02 13:18:44 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-07-02 13:18:44 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-07-02 13:18:44 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-07-02 13:18:44 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-07-02 13:18:45 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-07-02 13:18:45 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-07-02 13:18:45 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-07-02 13:18:45 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-07-02 13:18:45 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-07-02 13:18:45 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-07-02 13:18:45 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-07-02 13:18:45 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-07-02 13:18:45 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-07-02 13:18:45 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-07-02 13:18:45 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-07-02 13:18:45 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-07-02 13:18:45 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-07-02 13:18:45 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-07-02 13:18:45 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-07-02 13:18:45 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-07-02 13:18:45 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-07-02 13:18:45 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-07-02 13:18:45 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-07-02 13:18:45 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-07-02 13:18:45 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-07-02 13:18:45 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-07-02 13:18:45 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-07-02 13:18:45 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-07-02 13:18:45 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-02 13:18:45 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-02 13:18:45 INFO run_dart() is DONE.
2026-07-02 13:18:45 INFO ---------->>> Running update_pollutant_in_end()
2026-07-02 13:18:45 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:18:45 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:18:46 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:18:46 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:18:46 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:18:47 INFO [Posterior Diff ENS1] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS1_2020020611.nc
2026-07-02 13:18:47 INFO [Posterior Diff ENS1] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS1_2020020611.relative.nc
2026-07-02 13:18:47 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:18:47 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:18:48 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:18:48 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:18:49 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:18:49 INFO [Posterior Diff ENS2] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS2_2020020611.nc
2026-07-02 13:18:49 INFO [Posterior Diff ENS2] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS2_2020020611.relative.nc
2026-07-02 13:18:49 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:18:49 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:18:50 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:18:50 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:18:51 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:18:51 INFO [Posterior Diff ENS3] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS3_2020020611.nc
2026-07-02 13:18:51 INFO [Posterior Diff ENS3] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS3_2020020611.relative.nc
2026-07-02 13:18:51 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:18:51 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:18:52 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:18:52 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:18:53 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:18:53 INFO [Posterior Diff ENS4] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS4_2020020611.nc
2026-07-02 13:18:53 INFO [Posterior Diff ENS4] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS4_2020020611.relative.nc
2026-07-02 13:18:53 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:18:53 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:18:54 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:18:54 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:18:55 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:18:55 INFO [Posterior Diff ENS5] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS5_2020020611.nc
2026-07-02 13:18:55 INFO [Posterior Diff ENS5] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS5_2020020611.relative.nc
2026-07-02 13:18:55 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:18:55 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:18:56 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:18:56 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:18:57 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:18:57 INFO [Posterior Diff ENS6] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS6_2020020611.nc
2026-07-02 13:18:57 INFO [Posterior Diff ENS6] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS6_2020020611.relative.nc
2026-07-02 13:18:57 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:18:57 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:18:58 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:18:58 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:18:59 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:18:59 INFO [Posterior Diff ENS7] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS7_2020020611.nc
2026-07-02 13:18:59 INFO [Posterior Diff ENS7] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS7_2020020611.relative.nc
2026-07-02 13:19:00 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:19:00 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:19:01 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:19:01 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:19:01 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:19:01 INFO [Posterior Diff ENS8] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS8_2020020611.nc
2026-07-02 13:19:01 INFO [Posterior Diff ENS8] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS8_2020020611.relative.nc
2026-07-02 13:19:02 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:19:02 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:19:03 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:19:03 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:19:03 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:19:03 INFO [Posterior Diff ENS9] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS9_2020020611.nc
2026-07-02 13:19:03 INFO [Posterior Diff ENS9] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS9_2020020611.relative.nc
2026-07-02 13:19:04 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:19:04 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:19:05 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:19:05 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:19:05 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:19:06 INFO [Posterior Diff ENS10] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS10_2020020611.nc
2026-07-02 13:19:06 INFO [Posterior Diff ENS10] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS10_2020020611.relative.nc
2026-07-02 13:19:06 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:19:06 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:19:07 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:19:07 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:19:07 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:19:08 INFO [Posterior Diff ENS11] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS11_2020020611.nc
2026-07-02 13:19:08 INFO [Posterior Diff ENS11] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS11_2020020611.relative.nc
2026-07-02 13:19:08 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:19:08 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:19:09 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:19:09 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:19:09 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:19:10 INFO [Posterior Diff ENS12] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS12_2020020611.nc
2026-07-02 13:19:10 INFO [Posterior Diff ENS12] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS12_2020020611.relative.nc
2026-07-02 13:19:10 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:19:10 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:19:11 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:19:11 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:19:11 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:19:12 INFO [Posterior Diff ENS13] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS13_2020020611.nc
2026-07-02 13:19:12 INFO [Posterior Diff ENS13] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS13_2020020611.relative.nc
2026-07-02 13:19:12 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:19:12 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:19:13 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:19:13 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:19:14 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:19:14 INFO [Posterior Diff ENS14] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS14_2020020611.nc
2026-07-02 13:19:14 INFO [Posterior Diff ENS14] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS14_2020020611.relative.nc
2026-07-02 13:19:14 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:19:14 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:19:15 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:19:15 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:19:16 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:19:16 INFO [Posterior Diff ENS15] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS15_2020020611.nc
2026-07-02 13:19:16 INFO [Posterior Diff ENS15] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS15_2020020611.relative.nc
2026-07-02 13:19:16 INFO Next run starts from 2020-02-06 11:00:00
2026-07-02 13:19:16 INFO Cycle is DONE; starting a new loop!
2026-07-02 13:19:16 INFO [TIME] step_end current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 13:19:16 INFO [TIME] step_start current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 13:19:16 INFO [TIME] window start=2020-02-06 11:00:00 end=2020-02-06 13:00:00 run_hours=2 has_assimilation=True
2026-07-02 13:19:16 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:19:17 INFO Hourly dataset computed and listing created
2026-07-02 13:19:20 INFO Hourly dataset computed
2026-07-02 13:19:20 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:19:21 INFO Hourly dataset computed and listing created
2026-07-02 13:19:22 INFO Hourly dataset computed
2026-07-02 13:19:22 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:19:23 INFO Hourly dataset computed and listing created
2026-07-02 13:19:24 INFO Hourly dataset computed
2026-07-02 13:19:24 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:19:24 INFO Hourly dataset computed and listing created
2026-07-02 13:19:25 INFO Hourly dataset computed
2026-07-02 13:19:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:19:26 INFO Hourly dataset computed and listing created
2026-07-02 13:19:27 INFO Hourly dataset computed
2026-07-02 13:19:27 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:19:28 INFO Hourly dataset computed and listing created
2026-07-02 13:19:29 INFO Hourly dataset computed
2026-07-02 13:19:29 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:19:30 INFO Hourly dataset computed and listing created
2026-07-02 13:19:31 INFO Hourly dataset computed
2026-07-02 13:19:31 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:19:32 INFO Hourly dataset computed and listing created
2026-07-02 13:19:32 INFO Hourly dataset computed
2026-07-02 13:19:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:19:33 INFO Hourly dataset computed and listing created
2026-07-02 13:19:34 INFO Hourly dataset computed
2026-07-02 13:19:34 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:19:35 INFO Hourly dataset computed and listing created
2026-07-02 13:19:36 INFO Hourly dataset computed
2026-07-02 13:19:36 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:19:37 INFO Hourly dataset computed and listing created
2026-07-02 13:19:38 INFO Hourly dataset computed
2026-07-02 13:19:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:19:38 INFO Hourly dataset computed and listing created
2026-07-02 13:19:39 INFO Hourly dataset computed
2026-07-02 13:19:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:19:40 INFO Hourly dataset computed and listing created
2026-07-02 13:19:41 INFO Hourly dataset computed
2026-07-02 13:19:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:19:42 INFO Hourly dataset computed and listing created
2026-07-02 13:19:43 INFO Hourly dataset computed
2026-07-02 13:19:43 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:19:43 INFO Hourly dataset computed and listing created
2026-07-02 13:19:44 INFO Hourly dataset computed
2026-07-02 13:19:44 INFO ---------->>> Running CHIMERE model from 2020-02-06 11:00:00 to 2020-02-06 13:00:00
2026-07-02 13:19:44 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:19:44 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1
2026-07-02 13:19:44 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc
2026-07-02 13:19:44 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-02 13:19:44 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:19:44 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-02 13:19:44 INFO Queuing job for member 1...
2026-07-02 13:19:44 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:19:44 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-02 13:19:45 INFO Found: ['5065981']
2026-07-02 13:19:50 INFO [TGCC-IRENE] Submitted job with ID:['5065981']
2026-07-02 13:19:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:19:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2
2026-07-02 13:19:50 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc
2026-07-02 13:19:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-02 13:19:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:19:50 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-02 13:19:50 INFO Queuing job for member 2...
2026-07-02 13:19:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:19:50 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-02 13:19:51 INFO Found: ['5065982']
2026-07-02 13:19:56 INFO [TGCC-IRENE] Submitted job with ID:['5065982']
2026-07-02 13:19:56 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:19:56 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3
2026-07-02 13:19:56 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc
2026-07-02 13:19:56 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-02 13:19:56 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:19:56 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-02 13:19:56 INFO Queuing job for member 3...
2026-07-02 13:19:56 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:19:56 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-02 13:19:57 INFO Found: ['5065983']
2026-07-02 13:20:02 INFO [TGCC-IRENE] Submitted job with ID:['5065983']
2026-07-02 13:20:02 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:20:02 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4
2026-07-02 13:20:02 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc
2026-07-02 13:20:02 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-02 13:20:02 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:20:02 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-02 13:20:02 INFO Queuing job for member 4...
2026-07-02 13:20:02 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:20:02 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-02 13:20:02 INFO Found: ['5065985']
2026-07-02 13:20:07 INFO [TGCC-IRENE] Submitted job with ID:['5065985']
2026-07-02 13:20:07 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:20:07 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5
2026-07-02 13:20:07 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc
2026-07-02 13:20:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-02 13:20:07 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:20:07 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-02 13:20:07 INFO Queuing job for member 5...
2026-07-02 13:20:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:20:07 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-02 13:20:08 INFO Found: ['5065986']
2026-07-02 13:20:13 INFO [TGCC-IRENE] Submitted job with ID:['5065986']
2026-07-02 13:20:13 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:20:13 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6
2026-07-02 13:20:13 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc
2026-07-02 13:20:13 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-02 13:20:13 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:20:13 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-02 13:20:13 INFO Queuing job for member 6...
2026-07-02 13:20:13 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:20:13 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-02 13:20:14 INFO Found: ['5065987']
2026-07-02 13:20:19 INFO [TGCC-IRENE] Submitted job with ID:['5065987']
2026-07-02 13:20:19 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:20:19 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7
2026-07-02 13:20:19 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc
2026-07-02 13:20:19 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-02 13:20:19 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:20:19 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-02 13:20:19 INFO Queuing job for member 7...
2026-07-02 13:20:19 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:20:19 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-02 13:20:22 INFO Found: ['5065988']
2026-07-02 13:20:27 INFO [TGCC-IRENE] Submitted job with ID:['5065988']
2026-07-02 13:20:27 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:20:27 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8
2026-07-02 13:20:27 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc
2026-07-02 13:20:27 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-02 13:20:27 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:20:27 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-02 13:20:27 INFO Queuing job for member 8...
2026-07-02 13:20:27 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:20:27 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-02 13:20:29 INFO Found: ['5065989']
2026-07-02 13:20:34 INFO [TGCC-IRENE] Submitted job with ID:['5065989']
2026-07-02 13:20:34 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:20:34 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9
2026-07-02 13:20:34 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc
2026-07-02 13:20:34 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-02 13:20:34 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:20:34 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-02 13:20:34 INFO Queuing job for member 9...
2026-07-02 13:20:34 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:20:34 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-02 13:20:35 INFO Found: ['5065994']
2026-07-02 13:20:40 INFO [TGCC-IRENE] Submitted job with ID:['5065994']
2026-07-02 13:20:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:20:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10
2026-07-02 13:20:40 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc
2026-07-02 13:20:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-02 13:20:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:20:40 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-02 13:20:40 INFO Queuing job for member 10...
2026-07-02 13:20:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:20:40 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-02 13:20:41 INFO Found: ['5065995']
2026-07-02 13:20:46 INFO [TGCC-IRENE] Submitted job with ID:['5065995']
2026-07-02 13:20:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:20:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11
2026-07-02 13:20:46 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc
2026-07-02 13:20:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-02 13:20:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:20:46 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-02 13:20:46 INFO Queuing job for member 11...
2026-07-02 13:20:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:20:46 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-02 13:20:47 INFO Found: ['5065997']
2026-07-02 13:20:52 INFO [TGCC-IRENE] Submitted job with ID:['5065997']
2026-07-02 13:20:52 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:20:52 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12
2026-07-02 13:20:52 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc
2026-07-02 13:20:52 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-02 13:20:52 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:20:52 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-02 13:20:52 INFO Queuing job for member 12...
2026-07-02 13:20:52 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:20:52 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-02 13:20:52 INFO Found: ['5065998']
2026-07-02 13:20:57 INFO [TGCC-IRENE] Submitted job with ID:['5065998']
2026-07-02 13:20:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:20:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13
2026-07-02 13:20:57 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc
2026-07-02 13:20:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-02 13:20:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:20:57 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-02 13:20:57 INFO Queuing job for member 13...
2026-07-02 13:20:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:20:57 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-02 13:20:58 INFO Found: ['5065999']
2026-07-02 13:21:03 INFO [TGCC-IRENE] Submitted job with ID:['5065999']
2026-07-02 13:21:03 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:21:03 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14
2026-07-02 13:21:03 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc
2026-07-02 13:21:03 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-02 13:21:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:21:03 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-02 13:21:03 INFO Queuing job for member 14...
2026-07-02 13:21:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:21:03 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-02 13:21:04 INFO Found: ['5066001']
2026-07-02 13:21:09 INFO [TGCC-IRENE] Submitted job with ID:['5066001']
2026-07-02 13:21:09 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:21:09 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15
2026-07-02 13:21:09 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc
2026-07-02 13:21:09 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-02 13:21:09 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:21:09 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-02 13:21:09 INFO Queuing job for member 15...
2026-07-02 13:21:09 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:21:09 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-02 13:21:10 INFO Found: ['5066002']
2026-07-02 13:21:15 INFO [TGCC-IRENE] Submitted job with ID:['5066002']
2026-07-02 13:21:15 INFO Checking job status ...
2026-07-02 13:21:15 INFO None 5065981: status RUNNING/PENDING
2026-07-02 13:21:15 INFO None 5065982: status RUNNING/PENDING
2026-07-02 13:21:15 INFO None 5065983: status RUNNING/PENDING
2026-07-02 13:21:15 INFO None 5065985: status RUNNING/PENDING
2026-07-02 13:21:15 INFO None 5065986: status RUNNING/PENDING
2026-07-02 13:21:15 INFO None 5065987: status RUNNING/PENDING
2026-07-02 13:21:15 INFO None 5065988: status RUNNING/PENDING
2026-07-02 13:21:15 INFO None 5065989: status RUNNING/PENDING
2026-07-02 13:21:15 INFO None 5065994: status RUNNING/PENDING
2026-07-02 13:21:15 INFO None 5065995: status RUNNING/PENDING
2026-07-02 13:21:15 INFO None 5065997: status RUNNING/PENDING
2026-07-02 13:21:15 INFO None 5065998: status RUNNING/PENDING
2026-07-02 13:21:15 INFO None 5065999: status RUNNING/PENDING
2026-07-02 13:21:15 INFO None 5066001: status RUNNING/PENDING
2026-07-02 13:21:15 INFO None 5066002: status RUNNING/PENDING
2026-07-02 13:21:15 INFO Jobs still running: ['5065981', '5065982', '5065983', '5065985', '5065986', '5065987', '5065988', '5065989', '5065994', '5065995', '5065997', '5065998', '5065999', '5066001', '5066002']. Waiting...
2026-07-02 13:21:30 INFO None 5065981: status RUNNING/PENDING
2026-07-02 13:21:30 INFO None 5065982: status RUNNING/PENDING
2026-07-02 13:21:30 INFO None 5065983: status RUNNING/PENDING
2026-07-02 13:21:30 INFO None 5065985: status RUNNING/PENDING
2026-07-02 13:21:30 INFO None 5065986: status RUNNING/PENDING
2026-07-02 13:21:30 INFO None 5065987: status RUNNING/PENDING
2026-07-02 13:21:30 INFO None 5065988: status RUNNING/PENDING
2026-07-02 13:21:31 INFO None 5065989: status RUNNING/PENDING
2026-07-02 13:21:31 INFO None 5065994: status RUNNING/PENDING
2026-07-02 13:21:31 INFO None 5065995: status RUNNING/PENDING
2026-07-02 13:21:31 INFO None 5065997: status RUNNING/PENDING
2026-07-02 13:21:31 INFO None 5065998: status RUNNING/PENDING
2026-07-02 13:21:31 INFO None 5065999: status RUNNING/PENDING
2026-07-02 13:21:31 INFO None 5066001: status RUNNING/PENDING
2026-07-02 13:21:31 INFO None 5066002: status RUNNING/PENDING
2026-07-02 13:21:31 INFO Jobs still running: ['5065981', '5065982', '5065983', '5065985', '5065986', '5065987', '5065988', '5065989', '5065994', '5065995', '5065997', '5065998', '5065999', '5066001', '5066002']. Waiting...
2026-07-02 13:21:46 INFO None 5065981: status RUNNING/PENDING
2026-07-02 13:21:46 INFO None 5065982: status RUNNING/PENDING
2026-07-02 13:21:46 INFO None 5065983: status RUNNING/PENDING
2026-07-02 13:21:46 INFO None 5065985: status RUNNING/PENDING
2026-07-02 13:21:46 INFO None 5065986: status RUNNING/PENDING
2026-07-02 13:21:46 INFO None 5065987: status RUNNING/PENDING
2026-07-02 13:21:46 INFO None 5065988: status RUNNING/PENDING
2026-07-02 13:21:46 INFO None 5065989: status RUNNING/PENDING
2026-07-02 13:21:46 INFO None 5065994: status RUNNING/PENDING
2026-07-02 13:21:46 INFO None 5065995: status RUNNING/PENDING
2026-07-02 13:21:46 INFO None 5065997: status RUNNING/PENDING
2026-07-02 13:21:46 INFO None 5065998: status RUNNING/PENDING
2026-07-02 13:21:46 INFO None 5065999: status RUNNING/PENDING
2026-07-02 13:21:46 INFO None 5066001: status RUNNING/PENDING
2026-07-02 13:21:46 INFO None 5066002: status RUNNING/PENDING
2026-07-02 13:21:46 INFO Jobs still running: ['5065981', '5065982', '5065983', '5065985', '5065986', '5065987', '5065988', '5065989', '5065994', '5065995', '5065997', '5065998', '5065999', '5066001', '5066002']. Waiting...
2026-07-02 13:22:01 INFO None 5065981: status RUNNING/PENDING
2026-07-02 13:22:01 INFO None 5065982: status RUNNING/PENDING
2026-07-02 13:22:01 INFO None 5065983: status RUNNING/PENDING
2026-07-02 13:22:01 INFO None 5065985: status RUNNING/PENDING
2026-07-02 13:22:01 INFO None 5065986: status RUNNING/PENDING
2026-07-02 13:22:01 INFO None 5065987: status RUNNING/PENDING
2026-07-02 13:22:02 INFO None 5065988: status RUNNING/PENDING
2026-07-02 13:22:02 INFO None 5065989: status RUNNING/PENDING
2026-07-02 13:22:02 INFO None 5065994: status RUNNING/PENDING
2026-07-02 13:22:02 INFO None 5065995: status RUNNING/PENDING
2026-07-02 13:22:02 INFO None 5065997: status RUNNING/PENDING
2026-07-02 13:22:02 INFO None 5065998: status RUNNING/PENDING
2026-07-02 13:22:02 INFO None 5065999: status RUNNING/PENDING
2026-07-02 13:22:02 INFO None 5066001: status RUNNING/PENDING
2026-07-02 13:22:02 INFO None 5066002: status RUNNING/PENDING
2026-07-02 13:22:02 INFO Jobs still running: ['5065981', '5065982', '5065983', '5065985', '5065986', '5065987', '5065988', '5065989', '5065994', '5065995', '5065997', '5065998', '5065999', '5066001', '5066002']. Waiting...
2026-07-02 13:22:18 INFO None 5065981: status RUNNING/PENDING
2026-07-02 13:22:18 INFO None 5065982: status RUNNING/PENDING
2026-07-02 13:22:18 INFO None 5065983: status RUNNING/PENDING
2026-07-02 13:22:18 INFO None 5065985: status RUNNING/PENDING
2026-07-02 13:22:18 INFO None 5065986: status RUNNING/PENDING
2026-07-02 13:22:18 INFO None 5065987: status RUNNING/PENDING
2026-07-02 13:22:18 INFO None 5065988: status RUNNING/PENDING
2026-07-02 13:22:18 INFO None 5065989: status RUNNING/PENDING
2026-07-02 13:22:18 INFO None 5065994: status RUNNING/PENDING
2026-07-02 13:22:18 INFO None 5065995: status RUNNING/PENDING
2026-07-02 13:22:18 INFO None 5065997: status RUNNING/PENDING
2026-07-02 13:22:18 INFO None 5065998: status RUNNING/PENDING
2026-07-02 13:22:18 INFO None 5065999: status RUNNING/PENDING
2026-07-02 13:22:18 INFO None 5066001: status RUNNING/PENDING
2026-07-02 13:22:18 INFO None 5066002: status RUNNING/PENDING
2026-07-02 13:22:18 INFO Jobs still running: ['5065981', '5065982', '5065983', '5065985', '5065986', '5065987', '5065988', '5065989', '5065994', '5065995', '5065997', '5065998', '5065999', '5066001', '5066002']. Waiting...
2026-07-02 13:22:33 INFO None 5065981: status RUNNING/PENDING
2026-07-02 13:22:33 INFO None 5065982: status RUNNING/PENDING
2026-07-02 13:22:33 INFO None 5065983: status RUNNING/PENDING
2026-07-02 13:22:33 INFO None 5065985: status RUNNING/PENDING
2026-07-02 13:22:33 INFO None 5065986: status RUNNING/PENDING
2026-07-02 13:22:33 INFO None 5065987: status RUNNING/PENDING
2026-07-02 13:22:33 INFO None 5065988: status RUNNING/PENDING
2026-07-02 13:22:33 INFO None 5065989: status RUNNING/PENDING
2026-07-02 13:22:33 INFO None 5065994: status RUNNING/PENDING
2026-07-02 13:22:34 INFO None 5065995: status RUNNING/PENDING
2026-07-02 13:22:34 INFO None 5065997: status RUNNING/PENDING
2026-07-02 13:22:34 INFO None 5065998: status RUNNING/PENDING
2026-07-02 13:22:34 INFO None 5065999: status RUNNING/PENDING
2026-07-02 13:22:34 INFO None 5066001: status RUNNING/PENDING
2026-07-02 13:22:34 INFO None 5066002: status RUNNING/PENDING
2026-07-02 13:22:34 INFO Jobs still running: ['5065981', '5065982', '5065983', '5065985', '5065986', '5065987', '5065988', '5065989', '5065994', '5065995', '5065997', '5065998', '5065999', '5066001', '5066002']. Waiting...
2026-07-02 13:22:49 INFO None 5065981: status RUNNING/PENDING
2026-07-02 13:22:49 INFO None 5065982: status RUNNING/PENDING
2026-07-02 13:22:49 INFO None 5065983: status RUNNING/PENDING
2026-07-02 13:22:49 INFO None 5065985: status RUNNING/PENDING
2026-07-02 13:22:49 INFO None 5065986: status RUNNING/PENDING
2026-07-02 13:22:49 INFO None 5065987: status RUNNING/PENDING
2026-07-02 13:22:49 INFO None 5065988: status RUNNING/PENDING
2026-07-02 13:22:49 INFO None 5065989: status RUNNING/PENDING
2026-07-02 13:22:49 INFO None 5065994: status RUNNING/PENDING
2026-07-02 13:22:49 INFO None 5065995: status RUNNING/PENDING
2026-07-02 13:22:49 INFO None 5065997: status RUNNING/PENDING
2026-07-02 13:22:49 INFO None 5065998: status RUNNING/PENDING
2026-07-02 13:22:49 INFO None 5065999: status RUNNING/PENDING
2026-07-02 13:22:49 INFO None 5066001: status RUNNING/PENDING
2026-07-02 13:22:49 INFO None 5066002: status RUNNING/PENDING
2026-07-02 13:22:49 INFO Jobs still running: ['5065981', '5065982', '5065983', '5065985', '5065986', '5065987', '5065988', '5065989', '5065994', '5065995', '5065997', '5065998', '5065999', '5066001', '5066002']. Waiting...
2026-07-02 13:23:04 INFO None 5065981: status RUNNING/PENDING
2026-07-02 13:23:04 INFO None 5065982: status FINISHED
2026-07-02 13:23:04 INFO None 5065983: status FINISHED
2026-07-02 13:23:04 INFO None 5065985: status RUNNING/PENDING
2026-07-02 13:23:04 INFO None 5065986: status RUNNING/PENDING
2026-07-02 13:23:04 INFO None 5065987: status RUNNING/PENDING
2026-07-02 13:23:04 INFO None 5065988: status RUNNING/PENDING
2026-07-02 13:23:04 INFO None 5065989: status RUNNING/PENDING
2026-07-02 13:23:04 INFO None 5065994: status RUNNING/PENDING
2026-07-02 13:23:04 INFO None 5065995: status RUNNING/PENDING
2026-07-02 13:23:04 INFO None 5065997: status RUNNING/PENDING
2026-07-02 13:23:04 INFO None 5065998: status RUNNING/PENDING
2026-07-02 13:23:04 INFO None 5065999: status RUNNING/PENDING
2026-07-02 13:23:04 INFO None 5066001: status RUNNING/PENDING
2026-07-02 13:23:04 INFO None 5066002: status RUNNING/PENDING
2026-07-02 13:23:04 INFO Jobs still running: ['5065981', '5065985', '5065986', '5065987', '5065988', '5065989', '5065994', '5065995', '5065997', '5065998', '5065999', '5066001', '5066002']. Waiting...
2026-07-02 13:23:21 INFO None 5065981: status FINISHED
2026-07-02 13:23:21 INFO None 5065982: status FINISHED
2026-07-02 13:23:21 INFO None 5065983: status FINISHED
2026-07-02 13:23:21 INFO None 5065985: status RUNNING/PENDING
2026-07-02 13:23:21 INFO None 5065986: status RUNNING/PENDING
2026-07-02 13:23:21 INFO None 5065987: status RUNNING/PENDING
2026-07-02 13:23:21 INFO None 5065988: status RUNNING/PENDING
2026-07-02 13:23:21 INFO None 5065989: status RUNNING/PENDING
2026-07-02 13:23:21 INFO None 5065994: status RUNNING/PENDING
2026-07-02 13:23:21 INFO None 5065995: status RUNNING/PENDING
2026-07-02 13:23:21 INFO None 5065997: status RUNNING/PENDING
2026-07-02 13:23:21 INFO None 5065998: status RUNNING/PENDING
2026-07-02 13:23:21 INFO None 5065999: status RUNNING/PENDING
2026-07-02 13:23:21 INFO None 5066001: status RUNNING/PENDING
2026-07-02 13:23:21 INFO None 5066002: status RUNNING/PENDING
2026-07-02 13:23:21 INFO Jobs still running: ['5065985', '5065986', '5065987', '5065988', '5065989', '5065994', '5065995', '5065997', '5065998', '5065999', '5066001', '5066002']. Waiting...
2026-07-02 13:23:37 INFO None 5065981: status FINISHED
2026-07-02 13:23:37 INFO None 5065982: status FINISHED
2026-07-02 13:23:37 INFO None 5065983: status FINISHED
2026-07-02 13:23:37 INFO None 5065985: status RUNNING/PENDING
2026-07-02 13:23:37 INFO None 5065986: status FINISHED
2026-07-02 13:23:37 INFO None 5065987: status RUNNING/PENDING
2026-07-02 13:23:37 INFO None 5065988: status RUNNING/PENDING
2026-07-02 13:23:37 INFO None 5065989: status RUNNING/PENDING
2026-07-02 13:23:37 INFO None 5065994: status RUNNING/PENDING
2026-07-02 13:23:37 INFO None 5065995: status RUNNING/PENDING
2026-07-02 13:23:37 INFO None 5065997: status RUNNING/PENDING
2026-07-02 13:23:37 INFO None 5065998: status RUNNING/PENDING
2026-07-02 13:23:37 INFO None 5065999: status RUNNING/PENDING
2026-07-02 13:23:37 INFO None 5066001: status RUNNING/PENDING
2026-07-02 13:23:37 INFO None 5066002: status RUNNING/PENDING
2026-07-02 13:23:37 INFO Jobs still running: ['5065985', '5065987', '5065988', '5065989', '5065994', '5065995', '5065997', '5065998', '5065999', '5066001', '5066002']. Waiting...
2026-07-02 13:23:52 INFO None 5065981: status FINISHED
2026-07-02 13:23:52 INFO None 5065982: status FINISHED
2026-07-02 13:23:52 INFO None 5065983: status FINISHED
2026-07-02 13:23:52 INFO None 5065985: status FINISHED
2026-07-02 13:23:52 INFO None 5065986: status FINISHED
2026-07-02 13:23:52 INFO None 5065987: status RUNNING/PENDING
2026-07-02 13:23:52 INFO None 5065988: status RUNNING/PENDING
2026-07-02 13:23:52 INFO None 5065989: status RUNNING/PENDING
2026-07-02 13:23:52 INFO None 5065994: status RUNNING/PENDING
2026-07-02 13:23:52 INFO None 5065995: status RUNNING/PENDING
2026-07-02 13:23:52 INFO None 5065997: status RUNNING/PENDING
2026-07-02 13:23:52 INFO None 5065998: status RUNNING/PENDING
2026-07-02 13:23:52 INFO None 5065999: status RUNNING/PENDING
2026-07-02 13:23:52 INFO None 5066001: status RUNNING/PENDING
2026-07-02 13:23:52 INFO None 5066002: status RUNNING/PENDING
2026-07-02 13:23:52 INFO Jobs still running: ['5065987', '5065988', '5065989', '5065994', '5065995', '5065997', '5065998', '5065999', '5066001', '5066002']. Waiting...
2026-07-02 13:24:07 INFO None 5065981: status FINISHED
2026-07-02 13:24:07 INFO None 5065982: status FINISHED
2026-07-02 13:24:07 INFO None 5065983: status FINISHED
2026-07-02 13:24:07 INFO None 5065985: status FINISHED
2026-07-02 13:24:07 INFO None 5065986: status FINISHED
2026-07-02 13:24:07 INFO None 5065987: status RUNNING/PENDING
2026-07-02 13:24:07 INFO None 5065988: status RUNNING/PENDING
2026-07-02 13:24:07 INFO None 5065989: status RUNNING/PENDING
2026-07-02 13:24:07 INFO None 5065994: status RUNNING/PENDING
2026-07-02 13:24:07 INFO None 5065995: status RUNNING/PENDING
2026-07-02 13:24:07 INFO None 5065997: status FINISHED
2026-07-02 13:24:07 INFO None 5065998: status RUNNING/PENDING
2026-07-02 13:24:07 INFO None 5065999: status RUNNING/PENDING
2026-07-02 13:24:07 INFO None 5066001: status RUNNING/PENDING
2026-07-02 13:24:07 INFO None 5066002: status FINISHED
2026-07-02 13:24:07 INFO Jobs still running: ['5065987', '5065988', '5065989', '5065994', '5065995', '5065998', '5065999', '5066001']. Waiting...
2026-07-02 13:24:22 INFO None 5065981: status FINISHED
2026-07-02 13:24:22 INFO None 5065982: status FINISHED
2026-07-02 13:24:23 INFO None 5065983: status FINISHED
2026-07-02 13:24:23 INFO None 5065985: status FINISHED
2026-07-02 13:24:23 INFO None 5065986: status FINISHED
2026-07-02 13:24:23 INFO None 5065987: status RUNNING/PENDING
2026-07-02 13:24:23 INFO None 5065988: status RUNNING/PENDING
2026-07-02 13:24:23 INFO None 5065989: status RUNNING/PENDING
2026-07-02 13:24:23 INFO None 5065994: status FINISHED
2026-07-02 13:24:23 INFO None 5065995: status RUNNING/PENDING
2026-07-02 13:24:23 INFO None 5065997: status FINISHED
2026-07-02 13:24:23 INFO None 5065998: status RUNNING/PENDING
2026-07-02 13:24:23 INFO None 5065999: status FINISHED
2026-07-02 13:24:23 INFO None 5066001: status FINISHED
2026-07-02 13:24:23 INFO None 5066002: status FINISHED
2026-07-02 13:24:23 INFO Jobs still running: ['5065987', '5065988', '5065989', '5065995', '5065998']. Waiting...
2026-07-02 13:24:38 INFO None 5065981: status FINISHED
2026-07-02 13:24:38 INFO None 5065982: status FINISHED
2026-07-02 13:24:38 INFO None 5065983: status FINISHED
2026-07-02 13:24:38 INFO None 5065985: status FINISHED
2026-07-02 13:24:38 INFO None 5065986: status FINISHED
2026-07-02 13:24:38 INFO None 5065987: status FINISHED
2026-07-02 13:24:38 INFO None 5065988: status FINISHED
2026-07-02 13:24:38 INFO None 5065989: status FINISHED
2026-07-02 13:24:38 INFO None 5065994: status FINISHED
2026-07-02 13:24:38 INFO None 5065995: status RUNNING/PENDING
2026-07-02 13:24:38 INFO None 5065997: status FINISHED
2026-07-02 13:24:38 INFO None 5065998: status FINISHED
2026-07-02 13:24:38 INFO None 5065999: status FINISHED
2026-07-02 13:24:38 INFO None 5066001: status FINISHED
2026-07-02 13:24:38 INFO None 5066002: status FINISHED
2026-07-02 13:24:38 INFO Jobs still running: ['5065995']. Waiting...
2026-07-02 13:24:53 INFO None 5065981: status FINISHED
2026-07-02 13:24:53 INFO None 5065982: status FINISHED
2026-07-02 13:24:53 INFO None 5065983: status FINISHED
2026-07-02 13:24:53 INFO None 5065985: status FINISHED
2026-07-02 13:24:53 INFO None 5065986: status FINISHED
2026-07-02 13:24:53 INFO None 5065987: status FINISHED
2026-07-02 13:24:53 INFO None 5065988: status FINISHED
2026-07-02 13:24:53 INFO None 5065989: status FINISHED
2026-07-02 13:24:53 INFO None 5065994: status FINISHED
2026-07-02 13:24:53 INFO None 5065995: status RUNNING/PENDING
2026-07-02 13:24:53 INFO None 5065997: status FINISHED
2026-07-02 13:24:53 INFO None 5065998: status FINISHED
2026-07-02 13:24:53 INFO None 5065999: status FINISHED
2026-07-02 13:24:53 INFO None 5066001: status FINISHED
2026-07-02 13:24:53 INFO None 5066002: status FINISHED
2026-07-02 13:24:53 INFO Jobs still running: ['5065995']. Waiting...
2026-07-02 13:25:08 INFO None 5065981: status FINISHED
2026-07-02 13:25:08 INFO None 5065982: status FINISHED
2026-07-02 13:25:08 INFO None 5065983: status FINISHED
2026-07-02 13:25:08 INFO None 5065985: status FINISHED
2026-07-02 13:25:08 INFO None 5065986: status FINISHED
2026-07-02 13:25:08 INFO None 5065987: status FINISHED
2026-07-02 13:25:08 INFO None 5065988: status FINISHED
2026-07-02 13:25:08 INFO None 5065989: status FINISHED
2026-07-02 13:25:09 INFO None 5065994: status FINISHED
2026-07-02 13:25:09 INFO None 5065995: status RUNNING/PENDING
2026-07-02 13:25:09 INFO None 5065997: status FINISHED
2026-07-02 13:25:09 INFO None 5065998: status FINISHED
2026-07-02 13:25:09 INFO None 5065999: status FINISHED
2026-07-02 13:25:09 INFO None 5066001: status FINISHED
2026-07-02 13:25:09 INFO None 5066002: status FINISHED
2026-07-02 13:25:09 INFO Jobs still running: ['5065995']. Waiting...
2026-07-02 13:25:24 INFO None 5065981: status FINISHED
2026-07-02 13:25:24 INFO None 5065982: status FINISHED
2026-07-02 13:25:24 INFO None 5065983: status FINISHED
2026-07-02 13:25:24 INFO None 5065985: status FINISHED
2026-07-02 13:25:24 INFO None 5065986: status FINISHED
2026-07-02 13:25:24 INFO None 5065987: status FINISHED
2026-07-02 13:25:24 INFO None 5065988: status FINISHED
2026-07-02 13:25:24 INFO None 5065989: status FINISHED
2026-07-02 13:25:24 INFO None 5065994: status FINISHED
2026-07-02 13:25:24 INFO None 5065995: status RUNNING/PENDING
2026-07-02 13:25:24 INFO None 5065997: status FINISHED
2026-07-02 13:25:24 INFO None 5065998: status FINISHED
2026-07-02 13:25:24 INFO None 5065999: status FINISHED
2026-07-02 13:25:24 INFO None 5066001: status FINISHED
2026-07-02 13:25:24 INFO None 5066002: status FINISHED
2026-07-02 13:25:24 INFO Jobs still running: ['5065995']. Waiting...
2026-07-02 13:25:39 INFO None 5065981: status FINISHED
2026-07-02 13:25:39 INFO None 5065982: status FINISHED
2026-07-02 13:25:39 INFO None 5065983: status FINISHED
2026-07-02 13:25:39 INFO None 5065985: status FINISHED
2026-07-02 13:25:39 INFO None 5065986: status FINISHED
2026-07-02 13:25:39 INFO None 5065987: status FINISHED
2026-07-02 13:25:39 INFO None 5065988: status FINISHED
2026-07-02 13:25:39 INFO None 5065989: status FINISHED
2026-07-02 13:25:39 INFO None 5065994: status FINISHED
2026-07-02 13:25:39 INFO None 5065995: status RUNNING/PENDING
2026-07-02 13:25:39 INFO None 5065997: status FINISHED
2026-07-02 13:25:39 INFO None 5065998: status FINISHED
2026-07-02 13:25:39 INFO None 5065999: status FINISHED
2026-07-02 13:25:39 INFO None 5066001: status FINISHED
2026-07-02 13:25:39 INFO None 5066002: status FINISHED
2026-07-02 13:25:39 INFO Jobs still running: ['5065995']. Waiting...
2026-07-02 13:25:54 INFO None 5065981: status FINISHED
2026-07-02 13:25:54 INFO None 5065982: status FINISHED
2026-07-02 13:25:54 INFO None 5065983: status FINISHED
2026-07-02 13:25:54 INFO None 5065985: status FINISHED
2026-07-02 13:25:54 INFO None 5065986: status FINISHED
2026-07-02 13:25:54 INFO None 5065987: status FINISHED
2026-07-02 13:25:54 INFO None 5065988: status FINISHED
2026-07-02 13:25:54 INFO None 5065989: status FINISHED
2026-07-02 13:25:54 INFO None 5065994: status FINISHED
2026-07-02 13:25:54 INFO None 5065995: status FINISHED
2026-07-02 13:25:54 INFO None 5065997: status FINISHED
2026-07-02 13:25:54 INFO None 5065998: status FINISHED
2026-07-02 13:25:54 INFO None 5065999: status FINISHED
2026-07-02 13:25:54 INFO None 5066001: status FINISHED
2026-07-02 13:25:54 INFO None 5066002: status FINISHED
2026-07-02 13:25:54 INFO Jobs ['5065981', '5065982', '5065983', '5065985', '5065986', '5065987', '5065988', '5065989', '5065994', '5065995', '5065997', '5065998', '5065999', '5066001', '5066002'] have finished
2026-07-02 13:25:54 INFO Checking restart files were created ...
2026-07-02 13:25:54 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020611_2_ENS1.nc(1002685915 bytes)
2026-07-02 13:25:54 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020611_2_ENS2.nc(1002685915 bytes)
2026-07-02 13:25:54 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020611_2_ENS3.nc(1002685915 bytes)
2026-07-02 13:25:54 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020611_2_ENS4.nc(1002685915 bytes)
2026-07-02 13:25:54 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020611_2_ENS5.nc(1002685915 bytes)
2026-07-02 13:25:54 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020611_2_ENS6.nc(1002685915 bytes)
2026-07-02 13:25:54 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020611_2_ENS7.nc(1002685915 bytes)
2026-07-02 13:25:54 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020611_2_ENS8.nc(1002685915 bytes)
2026-07-02 13:25:54 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020611_2_ENS9.nc(1002685915 bytes)
2026-07-02 13:25:54 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020611_2_ENS10.nc(1002685915 bytes)
2026-07-02 13:25:54 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020611_2_ENS11.nc(1002685915 bytes)
2026-07-02 13:25:54 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020611_2_ENS12.nc(1002685915 bytes)
2026-07-02 13:25:54 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020611_2_ENS13.nc(1002685915 bytes)
2026-07-02 13:25:54 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020611_2_ENS14.nc(1002685915 bytes)
2026-07-02 13:25:54 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020611_2_ENS15.nc(1002685915 bytes)
2026-07-02 13:25:54 INFO  Run_model() completed successfully.
2026-07-02 13:25:54 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 13:25:54 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 13:00:00 days=153072 seconds=46800
2026-07-02 13:25:54 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-02 13:25:54 INFO [TIME] increment current_time 2020-02-06 11:00:00 -> 2020-02-06 13:00:00
2026-07-02 13:25:54 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 13:25:54 INFO ---------->>> Running process_satellite_data()
2026-07-02 13:25:55 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12004.nc
2026-07-02 13:25:55 INFO ---------->>> Running run_obs_converter()
2026-07-02 13:25:55 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_45135_153072.out
2026-07-02 13:25:55 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_45135_153072.out
2026-07-02 13:25:55 INFO ---------->>> Running DART
2026-07-02 13:25:55 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-02 13:25:55 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/chim_ENS1_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/chim_ENS1_2020020613_1_out_toDART.nc
2026-07-02 13:25:55 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/chim_ENS2_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/chim_ENS2_2020020613_1_out_toDART.nc
2026-07-02 13:25:55 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/chim_ENS3_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/chim_ENS3_2020020613_1_out_toDART.nc
2026-07-02 13:25:56 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/chim_ENS4_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/chim_ENS4_2020020613_1_out_toDART.nc
2026-07-02 13:25:56 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/chim_ENS5_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/chim_ENS5_2020020613_1_out_toDART.nc
2026-07-02 13:25:56 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/chim_ENS6_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/chim_ENS6_2020020613_1_out_toDART.nc
2026-07-02 13:25:57 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/chim_ENS7_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/chim_ENS7_2020020613_1_out_toDART.nc
2026-07-02 13:25:57 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/chim_ENS8_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/chim_ENS8_2020020613_1_out_toDART.nc
2026-07-02 13:25:57 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/chim_ENS9_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/chim_ENS9_2020020613_1_out_toDART.nc
2026-07-02 13:25:57 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/chim_ENS10_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/chim_ENS10_2020020613_1_out_toDART.nc
2026-07-02 13:25:58 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/chim_ENS11_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/chim_ENS11_2020020613_1_out_toDART.nc
2026-07-02 13:25:58 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/chim_ENS12_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/chim_ENS12_2020020613_1_out_toDART.nc
2026-07-02 13:25:58 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/chim_ENS13_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/chim_ENS13_2020020613_1_out_toDART.nc
2026-07-02 13:25:59 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/chim_ENS14_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/chim_ENS14_2020020613_1_out_toDART.nc
2026-07-02 13:25:59 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/chim_ENS15_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/chim_ENS15_2020020613_1_out_toDART.nc
2026-07-02 13:25:59 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-02 13:25:59 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-02 13:25:59 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-02 13:25:59 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-02 13:25:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-02 13:25:59 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-02 13:26:13 INFO Found: []
2026-07-02 13:26:13 INFO No job id returned by command ./run_filter.bsh
2026-07-02 13:26:13 INFO No monitoring will be performed
2026-07-02 13:26:13 INFO Moving DART output files to analysis and preassim directories for date 2020020613 if present ...
2026-07-02 13:26:13 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-07-02 13:26:13 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-07-02 13:26:13 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-07-02 13:26:13 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-07-02 13:26:13 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-07-02 13:26:13 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-07-02 13:26:13 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-07-02 13:26:13 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-07-02 13:26:13 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-07-02 13:26:13 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-07-02 13:26:13 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-07-02 13:26:13 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-07-02 13:26:13 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-07-02 13:26:13 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-07-02 13:26:13 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-07-02 13:26:13 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-07-02 13:26:13 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-07-02 13:26:13 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-07-02 13:26:13 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-07-02 13:26:13 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-07-02 13:26:13 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-07-02 13:26:13 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-07-02 13:26:13 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-07-02 13:26:13 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-07-02 13:26:13 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-07-02 13:26:13 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-07-02 13:26:13 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-07-02 13:26:13 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-07-02 13:26:13 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-07-02 13:26:13 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-07-02 13:26:13 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-07-02 13:26:13 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-07-02 13:26:13 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-07-02 13:26:13 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-07-02 13:26:13 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-02 13:26:14 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-02 13:26:14 INFO run_dart() is DONE.
2026-07-02 13:26:14 INFO ---------->>> Running update_pollutant_in_end()
2026-07-02 13:26:14 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:26:14 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:26:15 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:26:15 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:26:15 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:26:16 INFO [Posterior Diff ENS1] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS1_2020020613.nc
2026-07-02 13:26:16 INFO [Posterior Diff ENS1] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS1_2020020613.relative.nc
2026-07-02 13:26:16 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:26:16 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:26:17 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:26:17 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:26:18 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:26:18 INFO [Posterior Diff ENS2] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS2_2020020613.nc
2026-07-02 13:26:18 INFO [Posterior Diff ENS2] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS2_2020020613.relative.nc
2026-07-02 13:26:18 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:26:18 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:26:19 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:26:19 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:26:20 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:26:20 INFO [Posterior Diff ENS3] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS3_2020020613.nc
2026-07-02 13:26:20 INFO [Posterior Diff ENS3] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS3_2020020613.relative.nc
2026-07-02 13:26:20 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:26:20 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:26:21 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:26:21 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:26:22 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:26:22 INFO [Posterior Diff ENS4] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS4_2020020613.nc
2026-07-02 13:26:22 INFO [Posterior Diff ENS4] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS4_2020020613.relative.nc
2026-07-02 13:26:22 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:26:22 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:26:23 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:26:23 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:26:24 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:26:24 INFO [Posterior Diff ENS5] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS5_2020020613.nc
2026-07-02 13:26:24 INFO [Posterior Diff ENS5] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS5_2020020613.relative.nc
2026-07-02 13:26:25 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:26:25 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:26:26 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:26:26 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:26:26 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:26:26 INFO [Posterior Diff ENS6] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS6_2020020613.nc
2026-07-02 13:26:26 INFO [Posterior Diff ENS6] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS6_2020020613.relative.nc
2026-07-02 13:26:27 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:26:27 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:26:28 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:26:28 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:26:28 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:26:29 INFO [Posterior Diff ENS7] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS7_2020020613.nc
2026-07-02 13:26:29 INFO [Posterior Diff ENS7] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS7_2020020613.relative.nc
2026-07-02 13:26:29 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:26:29 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:26:30 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:26:30 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:26:31 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:26:31 INFO [Posterior Diff ENS8] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS8_2020020613.nc
2026-07-02 13:26:31 INFO [Posterior Diff ENS8] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS8_2020020613.relative.nc
2026-07-02 13:26:31 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:26:31 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:26:32 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:26:32 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:26:33 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:26:33 INFO [Posterior Diff ENS9] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS9_2020020613.nc
2026-07-02 13:26:33 INFO [Posterior Diff ENS9] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS9_2020020613.relative.nc
2026-07-02 13:26:33 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:26:33 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:26:34 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:26:34 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:26:35 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:26:35 INFO [Posterior Diff ENS10] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS10_2020020613.nc
2026-07-02 13:26:35 INFO [Posterior Diff ENS10] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS10_2020020613.relative.nc
2026-07-02 13:26:36 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:26:36 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:26:36 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:26:36 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:26:37 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:26:37 INFO [Posterior Diff ENS11] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS11_2020020613.nc
2026-07-02 13:26:37 INFO [Posterior Diff ENS11] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS11_2020020613.relative.nc
2026-07-02 13:26:38 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:26:38 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:26:39 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:26:39 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:26:39 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:26:39 INFO [Posterior Diff ENS12] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS12_2020020613.nc
2026-07-02 13:26:39 INFO [Posterior Diff ENS12] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS12_2020020613.relative.nc
2026-07-02 13:26:40 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:26:40 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:26:41 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:26:41 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:26:41 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:26:42 INFO [Posterior Diff ENS13] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS13_2020020613.nc
2026-07-02 13:26:42 INFO [Posterior Diff ENS13] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS13_2020020613.relative.nc
2026-07-02 13:26:42 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:26:42 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:26:43 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:26:43 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:26:44 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:26:44 INFO [Posterior Diff ENS14] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS14_2020020613.nc
2026-07-02 13:26:44 INFO [Posterior Diff ENS14] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS14_2020020613.relative.nc
2026-07-02 13:26:44 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:26:44 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:26:45 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:26:45 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:26:46 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:26:46 INFO [Posterior Diff ENS15] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS15_2020020613.nc
2026-07-02 13:26:46 INFO [Posterior Diff ENS15] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS15_2020020613.relative.nc
2026-07-02 13:26:46 INFO Next run starts from 2020-02-06 13:00:00
2026-07-02 13:26:46 INFO Cycle is DONE; starting a new loop!
2026-07-02 13:26:46 INFO [TIME] step_end current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 13:26:46 INFO [TIME] step_start current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 13:26:46 INFO [TIME] window start=2020-02-06 13:00:00 end=2020-02-06 14:00:00 run_hours=1 has_assimilation=True
2026-07-02 13:26:46 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:26:47 INFO Hourly dataset computed and listing created
2026-07-02 13:26:49 INFO Hourly dataset computed
2026-07-02 13:26:49 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:26:50 INFO Hourly dataset computed and listing created
2026-07-02 13:26:50 INFO Hourly dataset computed
2026-07-02 13:26:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:26:51 INFO Hourly dataset computed and listing created
2026-07-02 13:26:52 INFO Hourly dataset computed
2026-07-02 13:26:52 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:26:53 INFO Hourly dataset computed and listing created
2026-07-02 13:26:53 INFO Hourly dataset computed
2026-07-02 13:26:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:26:54 INFO Hourly dataset computed and listing created
2026-07-02 13:26:55 INFO Hourly dataset computed
2026-07-02 13:26:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:26:56 INFO Hourly dataset computed and listing created
2026-07-02 13:26:56 INFO Hourly dataset computed
2026-07-02 13:26:56 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:26:57 INFO Hourly dataset computed and listing created
2026-07-02 13:26:58 INFO Hourly dataset computed
2026-07-02 13:26:58 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:26:59 INFO Hourly dataset computed and listing created
2026-07-02 13:26:59 INFO Hourly dataset computed
2026-07-02 13:26:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:27:00 INFO Hourly dataset computed and listing created
2026-07-02 13:27:01 INFO Hourly dataset computed
2026-07-02 13:27:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:27:02 INFO Hourly dataset computed and listing created
2026-07-02 13:27:02 INFO Hourly dataset computed
2026-07-02 13:27:02 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:27:03 INFO Hourly dataset computed and listing created
2026-07-02 13:27:04 INFO Hourly dataset computed
2026-07-02 13:27:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:27:05 INFO Hourly dataset computed and listing created
2026-07-02 13:27:05 INFO Hourly dataset computed
2026-07-02 13:27:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:27:06 INFO Hourly dataset computed and listing created
2026-07-02 13:27:07 INFO Hourly dataset computed
2026-07-02 13:27:07 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:27:08 INFO Hourly dataset computed and listing created
2026-07-02 13:27:08 INFO Hourly dataset computed
2026-07-02 13:27:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:27:09 INFO Hourly dataset computed and listing created
2026-07-02 13:27:10 INFO Hourly dataset computed
2026-07-02 13:27:10 INFO ---------->>> Running CHIMERE model from 2020-02-06 13:00:00 to 2020-02-06 14:00:00
2026-07-02 13:27:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:27:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1
2026-07-02 13:27:10 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020611_2_ENS1.nc
2026-07-02 13:27:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-02 13:27:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:27:10 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-02 13:27:10 INFO Queuing job for member 1...
2026-07-02 13:27:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:27:10 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-02 13:27:11 INFO Found: ['5066022']
2026-07-02 13:27:16 INFO [TGCC-IRENE] Submitted job with ID:['5066022']
2026-07-02 13:27:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:27:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2
2026-07-02 13:27:16 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020611_2_ENS2.nc
2026-07-02 13:27:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-02 13:27:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:27:16 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-02 13:27:16 INFO Queuing job for member 2...
2026-07-02 13:27:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:27:16 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-02 13:27:16 INFO Found: ['5066024']
2026-07-02 13:27:21 INFO [TGCC-IRENE] Submitted job with ID:['5066024']
2026-07-02 13:27:21 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:27:21 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3
2026-07-02 13:27:21 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020611_2_ENS3.nc
2026-07-02 13:27:21 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-02 13:27:21 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:27:21 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-02 13:27:21 INFO Queuing job for member 3...
2026-07-02 13:27:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:27:21 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-02 13:27:22 INFO Found: ['5066025']
2026-07-02 13:27:27 INFO [TGCC-IRENE] Submitted job with ID:['5066025']
2026-07-02 13:27:27 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:27:27 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4
2026-07-02 13:27:27 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020611_2_ENS4.nc
2026-07-02 13:27:27 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-02 13:27:27 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:27:27 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-02 13:27:27 INFO Queuing job for member 4...
2026-07-02 13:27:27 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:27:27 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-02 13:27:28 INFO Found: ['5066026']
2026-07-02 13:27:33 INFO [TGCC-IRENE] Submitted job with ID:['5066026']
2026-07-02 13:27:33 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:27:33 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5
2026-07-02 13:27:33 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020611_2_ENS5.nc
2026-07-02 13:27:33 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-02 13:27:33 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:27:33 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-02 13:27:33 INFO Queuing job for member 5...
2026-07-02 13:27:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:27:33 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-02 13:27:34 INFO Found: ['5066027']
2026-07-02 13:27:39 INFO [TGCC-IRENE] Submitted job with ID:['5066027']
2026-07-02 13:27:39 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:27:39 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6
2026-07-02 13:27:39 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020611_2_ENS6.nc
2026-07-02 13:27:39 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-02 13:27:39 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:27:39 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-02 13:27:39 INFO Queuing job for member 6...
2026-07-02 13:27:39 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:27:39 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-02 13:27:39 INFO Found: ['5066028']
2026-07-02 13:27:44 INFO [TGCC-IRENE] Submitted job with ID:['5066028']
2026-07-02 13:27:44 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:27:44 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7
2026-07-02 13:27:44 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020611_2_ENS7.nc
2026-07-02 13:27:44 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-02 13:27:44 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:27:44 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-02 13:27:44 INFO Queuing job for member 7...
2026-07-02 13:27:44 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:27:44 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-02 13:27:45 INFO Found: ['5066029']
2026-07-02 13:27:50 INFO [TGCC-IRENE] Submitted job with ID:['5066029']
2026-07-02 13:27:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:27:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8
2026-07-02 13:27:50 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020611_2_ENS8.nc
2026-07-02 13:27:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-02 13:27:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:27:50 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-02 13:27:50 INFO Queuing job for member 8...
2026-07-02 13:27:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:27:50 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-02 13:27:51 INFO Found: ['5066030']
2026-07-02 13:27:56 INFO [TGCC-IRENE] Submitted job with ID:['5066030']
2026-07-02 13:27:56 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:27:56 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9
2026-07-02 13:27:56 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020611_2_ENS9.nc
2026-07-02 13:27:56 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-02 13:27:56 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:27:56 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-02 13:27:56 INFO Queuing job for member 9...
2026-07-02 13:27:56 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:27:56 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-02 13:27:58 INFO Found: ['5066031']
2026-07-02 13:28:03 INFO [TGCC-IRENE] Submitted job with ID:['5066031']
2026-07-02 13:28:03 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:28:03 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10
2026-07-02 13:28:03 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020611_2_ENS10.nc
2026-07-02 13:28:03 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-02 13:28:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:28:03 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-02 13:28:03 INFO Queuing job for member 10...
2026-07-02 13:28:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:28:03 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-02 13:28:05 INFO Found: ['5066033']
2026-07-02 13:28:10 INFO [TGCC-IRENE] Submitted job with ID:['5066033']
2026-07-02 13:28:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:28:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11
2026-07-02 13:28:10 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020611_2_ENS11.nc
2026-07-02 13:28:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-02 13:28:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:28:10 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-02 13:28:10 INFO Queuing job for member 11...
2026-07-02 13:28:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:28:10 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-02 13:28:10 INFO Found: ['5066034']
2026-07-02 13:28:15 INFO [TGCC-IRENE] Submitted job with ID:['5066034']
2026-07-02 13:28:15 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:28:15 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12
2026-07-02 13:28:15 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020611_2_ENS12.nc
2026-07-02 13:28:15 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-02 13:28:15 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:28:15 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-02 13:28:15 INFO Queuing job for member 12...
2026-07-02 13:28:15 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:28:15 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-02 13:28:16 INFO Found: ['5066035']
2026-07-02 13:28:21 INFO [TGCC-IRENE] Submitted job with ID:['5066035']
2026-07-02 13:28:21 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:28:21 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13
2026-07-02 13:28:21 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020611_2_ENS13.nc
2026-07-02 13:28:21 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-02 13:28:21 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:28:21 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-02 13:28:21 INFO Queuing job for member 13...
2026-07-02 13:28:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:28:21 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-02 13:28:22 INFO Found: ['5066036']
2026-07-02 13:28:27 INFO [TGCC-IRENE] Submitted job with ID:['5066036']
2026-07-02 13:28:27 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:28:27 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14
2026-07-02 13:28:27 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020611_2_ENS14.nc
2026-07-02 13:28:27 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-02 13:28:27 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:28:27 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-02 13:28:27 INFO Queuing job for member 14...
2026-07-02 13:28:27 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:28:27 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-02 13:28:28 INFO Found: ['5066037']
2026-07-02 13:28:33 INFO [TGCC-IRENE] Submitted job with ID:['5066037']
2026-07-02 13:28:33 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:28:33 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15
2026-07-02 13:28:33 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020611_2_ENS15.nc
2026-07-02 13:28:33 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-02 13:28:33 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:28:33 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-02 13:28:33 INFO Queuing job for member 15...
2026-07-02 13:28:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:28:33 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-02 13:28:33 INFO Found: ['5066038']
2026-07-02 13:28:38 INFO [TGCC-IRENE] Submitted job with ID:['5066038']
2026-07-02 13:28:38 INFO Checking job status ...
2026-07-02 13:28:38 INFO None 5066022: status RUNNING/PENDING
2026-07-02 13:28:38 INFO None 5066024: status RUNNING/PENDING
2026-07-02 13:28:38 INFO None 5066025: status RUNNING/PENDING
2026-07-02 13:28:38 INFO None 5066026: status RUNNING/PENDING
2026-07-02 13:28:38 INFO None 5066027: status RUNNING/PENDING
2026-07-02 13:28:38 INFO None 5066028: status RUNNING/PENDING
2026-07-02 13:28:39 INFO None 5066029: status RUNNING/PENDING
2026-07-02 13:28:39 INFO None 5066030: status RUNNING/PENDING
2026-07-02 13:28:39 INFO None 5066031: status RUNNING/PENDING
2026-07-02 13:28:39 INFO None 5066033: status RUNNING/PENDING
2026-07-02 13:28:39 INFO None 5066034: status RUNNING/PENDING
2026-07-02 13:28:39 INFO None 5066035: status RUNNING/PENDING
2026-07-02 13:28:39 INFO None 5066036: status RUNNING/PENDING
2026-07-02 13:28:39 INFO None 5066037: status RUNNING/PENDING
2026-07-02 13:28:39 INFO None 5066038: status RUNNING/PENDING
2026-07-02 13:28:39 INFO Jobs still running: ['5066022', '5066024', '5066025', '5066026', '5066027', '5066028', '5066029', '5066030', '5066031', '5066033', '5066034', '5066035', '5066036', '5066037', '5066038']. Waiting...
2026-07-02 13:28:54 INFO None 5066022: status RUNNING/PENDING
2026-07-02 13:28:54 INFO None 5066024: status RUNNING/PENDING
2026-07-02 13:28:54 INFO None 5066025: status RUNNING/PENDING
2026-07-02 13:28:54 INFO None 5066026: status RUNNING/PENDING
2026-07-02 13:28:54 INFO None 5066027: status RUNNING/PENDING
2026-07-02 13:28:54 INFO None 5066028: status RUNNING/PENDING
2026-07-02 13:28:54 INFO None 5066029: status RUNNING/PENDING
2026-07-02 13:28:54 INFO None 5066030: status RUNNING/PENDING
2026-07-02 13:28:54 INFO None 5066031: status RUNNING/PENDING
2026-07-02 13:28:54 INFO None 5066033: status RUNNING/PENDING
2026-07-02 13:28:54 INFO None 5066034: status RUNNING/PENDING
2026-07-02 13:28:54 INFO None 5066035: status RUNNING/PENDING
2026-07-02 13:28:54 INFO None 5066036: status RUNNING/PENDING
2026-07-02 13:28:54 INFO None 5066037: status RUNNING/PENDING
2026-07-02 13:28:54 INFO None 5066038: status RUNNING/PENDING
2026-07-02 13:28:54 INFO Jobs still running: ['5066022', '5066024', '5066025', '5066026', '5066027', '5066028', '5066029', '5066030', '5066031', '5066033', '5066034', '5066035', '5066036', '5066037', '5066038']. Waiting...
2026-07-02 13:29:09 INFO None 5066022: status RUNNING/PENDING
2026-07-02 13:29:09 INFO None 5066024: status RUNNING/PENDING
2026-07-02 13:29:09 INFO None 5066025: status RUNNING/PENDING
2026-07-02 13:29:09 INFO None 5066026: status RUNNING/PENDING
2026-07-02 13:29:09 INFO None 5066027: status RUNNING/PENDING
2026-07-02 13:29:09 INFO None 5066028: status RUNNING/PENDING
2026-07-02 13:29:09 INFO None 5066029: status RUNNING/PENDING
2026-07-02 13:29:09 INFO None 5066030: status RUNNING/PENDING
2026-07-02 13:29:09 INFO None 5066031: status RUNNING/PENDING
2026-07-02 13:29:09 INFO None 5066033: status RUNNING/PENDING
2026-07-02 13:29:09 INFO None 5066034: status RUNNING/PENDING
2026-07-02 13:29:10 INFO None 5066035: status RUNNING/PENDING
2026-07-02 13:29:10 INFO None 5066036: status RUNNING/PENDING
2026-07-02 13:29:10 INFO None 5066037: status RUNNING/PENDING
2026-07-02 13:29:10 INFO None 5066038: status RUNNING/PENDING
2026-07-02 13:29:10 INFO Jobs still running: ['5066022', '5066024', '5066025', '5066026', '5066027', '5066028', '5066029', '5066030', '5066031', '5066033', '5066034', '5066035', '5066036', '5066037', '5066038']. Waiting...
2026-07-02 13:29:25 INFO None 5066022: status RUNNING/PENDING
2026-07-02 13:29:25 INFO None 5066024: status RUNNING/PENDING
2026-07-02 13:29:25 INFO None 5066025: status RUNNING/PENDING
2026-07-02 13:29:25 INFO None 5066026: status RUNNING/PENDING
2026-07-02 13:29:25 INFO None 5066027: status RUNNING/PENDING
2026-07-02 13:29:25 INFO None 5066028: status RUNNING/PENDING
2026-07-02 13:29:25 INFO None 5066029: status RUNNING/PENDING
2026-07-02 13:29:25 INFO None 5066030: status RUNNING/PENDING
2026-07-02 13:29:25 INFO None 5066031: status RUNNING/PENDING
2026-07-02 13:29:25 INFO None 5066033: status RUNNING/PENDING
2026-07-02 13:29:25 INFO None 5066034: status RUNNING/PENDING
2026-07-02 13:29:25 INFO None 5066035: status RUNNING/PENDING
2026-07-02 13:29:25 INFO None 5066036: status RUNNING/PENDING
2026-07-02 13:29:25 INFO None 5066037: status RUNNING/PENDING
2026-07-02 13:29:25 INFO None 5066038: status RUNNING/PENDING
2026-07-02 13:29:25 INFO Jobs still running: ['5066022', '5066024', '5066025', '5066026', '5066027', '5066028', '5066029', '5066030', '5066031', '5066033', '5066034', '5066035', '5066036', '5066037', '5066038']. Waiting...
2026-07-02 13:29:40 INFO None 5066022: status RUNNING/PENDING
2026-07-02 13:29:40 INFO None 5066024: status RUNNING/PENDING
2026-07-02 13:29:40 INFO None 5066025: status RUNNING/PENDING
2026-07-02 13:29:40 INFO None 5066026: status RUNNING/PENDING
2026-07-02 13:29:40 INFO None 5066027: status RUNNING/PENDING
2026-07-02 13:29:40 INFO None 5066028: status RUNNING/PENDING
2026-07-02 13:29:40 INFO None 5066029: status RUNNING/PENDING
2026-07-02 13:29:40 INFO None 5066030: status RUNNING/PENDING
2026-07-02 13:29:40 INFO None 5066031: status RUNNING/PENDING
2026-07-02 13:29:40 INFO None 5066033: status RUNNING/PENDING
2026-07-02 13:29:40 INFO None 5066034: status RUNNING/PENDING
2026-07-02 13:29:40 INFO None 5066035: status RUNNING/PENDING
2026-07-02 13:29:40 INFO None 5066036: status RUNNING/PENDING
2026-07-02 13:29:40 INFO None 5066037: status RUNNING/PENDING
2026-07-02 13:29:40 INFO None 5066038: status RUNNING/PENDING
2026-07-02 13:29:40 INFO Jobs still running: ['5066022', '5066024', '5066025', '5066026', '5066027', '5066028', '5066029', '5066030', '5066031', '5066033', '5066034', '5066035', '5066036', '5066037', '5066038']. Waiting...
2026-07-02 13:29:55 INFO None 5066022: status FINISHED
2026-07-02 13:29:55 INFO None 5066024: status RUNNING/PENDING
2026-07-02 13:29:55 INFO None 5066025: status RUNNING/PENDING
2026-07-02 13:29:55 INFO None 5066026: status RUNNING/PENDING
2026-07-02 13:29:55 INFO None 5066027: status RUNNING/PENDING
2026-07-02 13:29:55 INFO None 5066028: status RUNNING/PENDING
2026-07-02 13:29:55 INFO None 5066029: status RUNNING/PENDING
2026-07-02 13:29:55 INFO None 5066030: status RUNNING/PENDING
2026-07-02 13:29:55 INFO None 5066031: status RUNNING/PENDING
2026-07-02 13:29:55 INFO None 5066033: status RUNNING/PENDING
2026-07-02 13:29:55 INFO None 5066034: status RUNNING/PENDING
2026-07-02 13:29:55 INFO None 5066035: status RUNNING/PENDING
2026-07-02 13:29:56 INFO None 5066036: status RUNNING/PENDING
2026-07-02 13:29:56 INFO None 5066037: status RUNNING/PENDING
2026-07-02 13:29:56 INFO None 5066038: status RUNNING/PENDING
2026-07-02 13:29:56 INFO Jobs still running: ['5066024', '5066025', '5066026', '5066027', '5066028', '5066029', '5066030', '5066031', '5066033', '5066034', '5066035', '5066036', '5066037', '5066038']. Waiting...
2026-07-02 13:30:11 INFO None 5066022: status FINISHED
2026-07-02 13:30:11 INFO None 5066024: status FINISHED
2026-07-02 13:30:11 INFO None 5066025: status FINISHED
2026-07-02 13:30:11 INFO None 5066026: status RUNNING/PENDING
2026-07-02 13:30:11 INFO None 5066027: status RUNNING/PENDING
2026-07-02 13:30:11 INFO None 5066028: status RUNNING/PENDING
2026-07-02 13:30:11 INFO None 5066029: status FINISHED
2026-07-02 13:30:11 INFO None 5066030: status FINISHED
2026-07-02 13:30:11 INFO None 5066031: status RUNNING/PENDING
2026-07-02 13:30:11 INFO None 5066033: status RUNNING/PENDING
2026-07-02 13:30:11 INFO None 5066034: status RUNNING/PENDING
2026-07-02 13:30:11 INFO None 5066035: status RUNNING/PENDING
2026-07-02 13:30:11 INFO None 5066036: status RUNNING/PENDING
2026-07-02 13:30:11 INFO None 5066037: status RUNNING/PENDING
2026-07-02 13:30:11 INFO None 5066038: status RUNNING/PENDING
2026-07-02 13:30:11 INFO Jobs still running: ['5066026', '5066027', '5066028', '5066031', '5066033', '5066034', '5066035', '5066036', '5066037', '5066038']. Waiting...
2026-07-02 13:30:26 INFO None 5066022: status FINISHED
2026-07-02 13:30:26 INFO None 5066024: status FINISHED
2026-07-02 13:30:26 INFO None 5066025: status FINISHED
2026-07-02 13:30:26 INFO None 5066026: status FINISHED
2026-07-02 13:30:26 INFO None 5066027: status FINISHED
2026-07-02 13:30:26 INFO None 5066028: status FINISHED
2026-07-02 13:30:27 INFO None 5066029: status FINISHED
2026-07-02 13:30:27 INFO None 5066030: status FINISHED
2026-07-02 13:30:27 INFO None 5066031: status RUNNING/PENDING
2026-07-02 13:30:27 INFO None 5066033: status RUNNING/PENDING
2026-07-02 13:30:27 INFO None 5066034: status RUNNING/PENDING
2026-07-02 13:30:27 INFO None 5066035: status RUNNING/PENDING
2026-07-02 13:30:27 INFO None 5066036: status RUNNING/PENDING
2026-07-02 13:30:27 INFO None 5066037: status RUNNING/PENDING
2026-07-02 13:30:27 INFO None 5066038: status RUNNING/PENDING
2026-07-02 13:30:27 INFO Jobs still running: ['5066031', '5066033', '5066034', '5066035', '5066036', '5066037', '5066038']. Waiting...
2026-07-02 13:30:43 INFO None 5066022: status FINISHED
2026-07-02 13:30:43 INFO None 5066024: status FINISHED
2026-07-02 13:30:43 INFO None 5066025: status FINISHED
2026-07-02 13:30:43 INFO None 5066026: status FINISHED
2026-07-02 13:30:43 INFO None 5066027: status FINISHED
2026-07-02 13:30:43 INFO None 5066028: status FINISHED
2026-07-02 13:30:43 INFO None 5066029: status FINISHED
2026-07-02 13:30:43 INFO None 5066030: status FINISHED
2026-07-02 13:30:43 INFO None 5066031: status RUNNING/PENDING
2026-07-02 13:30:43 INFO None 5066033: status RUNNING/PENDING
2026-07-02 13:30:43 INFO None 5066034: status RUNNING/PENDING
2026-07-02 13:30:43 INFO None 5066035: status RUNNING/PENDING
2026-07-02 13:30:43 INFO None 5066036: status RUNNING/PENDING
2026-07-02 13:30:43 INFO None 5066037: status FINISHED
2026-07-02 13:30:43 INFO None 5066038: status RUNNING/PENDING
2026-07-02 13:30:43 INFO Jobs still running: ['5066031', '5066033', '5066034', '5066035', '5066036', '5066038']. Waiting...
2026-07-02 13:30:58 INFO None 5066022: status FINISHED
2026-07-02 13:30:58 INFO None 5066024: status FINISHED
2026-07-02 13:30:58 INFO None 5066025: status FINISHED
2026-07-02 13:30:58 INFO None 5066026: status FINISHED
2026-07-02 13:30:58 INFO None 5066027: status FINISHED
2026-07-02 13:30:58 INFO None 5066028: status FINISHED
2026-07-02 13:30:58 INFO None 5066029: status FINISHED
2026-07-02 13:30:58 INFO None 5066030: status FINISHED
2026-07-02 13:30:58 INFO None 5066031: status RUNNING/PENDING
2026-07-02 13:30:58 INFO None 5066033: status RUNNING/PENDING
2026-07-02 13:30:59 INFO None 5066034: status RUNNING/PENDING
2026-07-02 13:30:59 INFO None 5066035: status RUNNING/PENDING
2026-07-02 13:30:59 INFO None 5066036: status RUNNING/PENDING
2026-07-02 13:30:59 INFO None 5066037: status FINISHED
2026-07-02 13:30:59 INFO None 5066038: status FINISHED
2026-07-02 13:30:59 INFO Jobs still running: ['5066031', '5066033', '5066034', '5066035', '5066036']. Waiting...
2026-07-02 13:31:14 INFO None 5066022: status FINISHED
2026-07-02 13:31:14 INFO None 5066024: status FINISHED
2026-07-02 13:31:14 INFO None 5066025: status FINISHED
2026-07-02 13:31:14 INFO None 5066026: status FINISHED
2026-07-02 13:31:14 INFO None 5066027: status FINISHED
2026-07-02 13:31:14 INFO None 5066028: status FINISHED
2026-07-02 13:31:14 INFO None 5066029: status FINISHED
2026-07-02 13:31:14 INFO None 5066030: status FINISHED
2026-07-02 13:31:14 INFO None 5066031: status FINISHED
2026-07-02 13:31:14 INFO None 5066033: status FINISHED
2026-07-02 13:31:14 INFO None 5066034: status FINISHED
2026-07-02 13:31:14 INFO None 5066035: status RUNNING/PENDING
2026-07-02 13:31:14 INFO None 5066036: status RUNNING/PENDING
2026-07-02 13:31:14 INFO None 5066037: status FINISHED
2026-07-02 13:31:14 INFO None 5066038: status FINISHED
2026-07-02 13:31:14 INFO Jobs still running: ['5066035', '5066036']. Waiting...
2026-07-02 13:31:29 INFO None 5066022: status FINISHED
2026-07-02 13:31:29 INFO None 5066024: status FINISHED
2026-07-02 13:31:29 INFO None 5066025: status FINISHED
2026-07-02 13:31:29 INFO None 5066026: status FINISHED
2026-07-02 13:31:29 INFO None 5066027: status FINISHED
2026-07-02 13:31:29 INFO None 5066028: status FINISHED
2026-07-02 13:31:29 INFO None 5066029: status FINISHED
2026-07-02 13:31:29 INFO None 5066030: status FINISHED
2026-07-02 13:31:29 INFO None 5066031: status FINISHED
2026-07-02 13:31:29 INFO None 5066033: status FINISHED
2026-07-02 13:31:29 INFO None 5066034: status FINISHED
2026-07-02 13:31:29 INFO None 5066035: status RUNNING/PENDING
2026-07-02 13:31:29 INFO None 5066036: status RUNNING/PENDING
2026-07-02 13:31:29 INFO None 5066037: status FINISHED
2026-07-02 13:31:29 INFO None 5066038: status FINISHED
2026-07-02 13:31:29 INFO Jobs still running: ['5066035', '5066036']. Waiting...
2026-07-02 13:31:44 INFO None 5066022: status FINISHED
2026-07-02 13:31:44 INFO None 5066024: status FINISHED
2026-07-02 13:31:44 INFO None 5066025: status FINISHED
2026-07-02 13:31:44 INFO None 5066026: status FINISHED
2026-07-02 13:31:44 INFO None 5066027: status FINISHED
2026-07-02 13:31:44 INFO None 5066028: status FINISHED
2026-07-02 13:31:44 INFO None 5066029: status FINISHED
2026-07-02 13:31:44 INFO None 5066030: status FINISHED
2026-07-02 13:31:44 INFO None 5066031: status FINISHED
2026-07-02 13:31:44 INFO None 5066033: status FINISHED
2026-07-02 13:31:44 INFO None 5066034: status FINISHED
2026-07-02 13:31:44 INFO None 5066035: status RUNNING/PENDING
2026-07-02 13:31:44 INFO None 5066036: status RUNNING/PENDING
2026-07-02 13:31:44 INFO None 5066037: status FINISHED
2026-07-02 13:31:45 INFO None 5066038: status FINISHED
2026-07-02 13:31:45 INFO Jobs still running: ['5066035', '5066036']. Waiting...
2026-07-02 13:32:00 INFO None 5066022: status FINISHED
2026-07-02 13:32:00 INFO None 5066024: status FINISHED
2026-07-02 13:32:00 INFO None 5066025: status FINISHED
2026-07-02 13:32:00 INFO None 5066026: status FINISHED
2026-07-02 13:32:00 INFO None 5066027: status FINISHED
2026-07-02 13:32:00 INFO None 5066028: status FINISHED
2026-07-02 13:32:00 INFO None 5066029: status FINISHED
2026-07-02 13:32:00 INFO None 5066030: status FINISHED
2026-07-02 13:32:00 INFO None 5066031: status FINISHED
2026-07-02 13:32:00 INFO None 5066033: status FINISHED
2026-07-02 13:32:00 INFO None 5066034: status FINISHED
2026-07-02 13:32:00 INFO None 5066035: status RUNNING/PENDING
2026-07-02 13:32:00 INFO None 5066036: status RUNNING/PENDING
2026-07-02 13:32:00 INFO None 5066037: status FINISHED
2026-07-02 13:32:00 INFO None 5066038: status FINISHED
2026-07-02 13:32:00 INFO Jobs still running: ['5066035', '5066036']. Waiting...
2026-07-02 13:32:15 INFO None 5066022: status FINISHED
2026-07-02 13:32:15 INFO None 5066024: status FINISHED
2026-07-02 13:32:15 INFO None 5066025: status FINISHED
2026-07-02 13:32:15 INFO None 5066026: status FINISHED
2026-07-02 13:32:15 INFO None 5066027: status FINISHED
2026-07-02 13:32:15 INFO None 5066028: status FINISHED
2026-07-02 13:32:15 INFO None 5066029: status FINISHED
2026-07-02 13:32:15 INFO None 5066030: status FINISHED
2026-07-02 13:32:15 INFO None 5066031: status FINISHED
2026-07-02 13:32:15 INFO None 5066033: status FINISHED
2026-07-02 13:32:15 INFO None 5066034: status FINISHED
2026-07-02 13:32:15 INFO None 5066035: status FINISHED
2026-07-02 13:32:15 INFO None 5066036: status FINISHED
2026-07-02 13:32:15 INFO None 5066037: status FINISHED
2026-07-02 13:32:15 INFO None 5066038: status FINISHED
2026-07-02 13:32:15 INFO Jobs ['5066022', '5066024', '5066025', '5066026', '5066027', '5066028', '5066029', '5066030', '5066031', '5066033', '5066034', '5066035', '5066036', '5066037', '5066038'] have finished
2026-07-02 13:32:15 INFO Checking restart files were created ...
2026-07-02 13:32:15 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020613_1_ENS1.nc(668832435 bytes)
2026-07-02 13:32:15 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020613_1_ENS2.nc(668832435 bytes)
2026-07-02 13:32:15 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020613_1_ENS3.nc(668832435 bytes)
2026-07-02 13:32:15 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020613_1_ENS4.nc(668832435 bytes)
2026-07-02 13:32:15 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020613_1_ENS5.nc(668832435 bytes)
2026-07-02 13:32:15 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020613_1_ENS6.nc(668832435 bytes)
2026-07-02 13:32:15 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020613_1_ENS7.nc(668832435 bytes)
2026-07-02 13:32:15 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020613_1_ENS8.nc(668832435 bytes)
2026-07-02 13:32:15 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020613_1_ENS9.nc(668832435 bytes)
2026-07-02 13:32:15 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020613_1_ENS10.nc(668832435 bytes)
2026-07-02 13:32:15 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020613_1_ENS11.nc(668832435 bytes)
2026-07-02 13:32:15 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020613_1_ENS12.nc(668832435 bytes)
2026-07-02 13:32:15 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020613_1_ENS13.nc(668832435 bytes)
2026-07-02 13:32:15 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020613_1_ENS14.nc(668832435 bytes)
2026-07-02 13:32:15 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020613_1_ENS15.nc(668832435 bytes)
2026-07-02 13:32:15 INFO  Run_model() completed successfully.
2026-07-02 13:32:15 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 13:32:15 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 14:00:00 days=153072 seconds=50400
2026-07-02 13:32:15 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-02 13:32:15 INFO [TIME] increment current_time 2020-02-06 13:00:00 -> 2020-02-06 14:00:00
2026-07-02 13:32:15 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 13:32:15 INFO ---------->>> Running process_satellite_data()
2026-07-02 13:32:15 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12005.nc
2026-07-02 13:32:15 INFO ---------->>> Running run_obs_converter()
2026-07-02 13:32:15 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_51225_153072.out
2026-07-02 13:32:15 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_51225_153072.out
2026-07-02 13:32:15 INFO ---------->>> Running DART
2026-07-02 13:32:15 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-02 13:32:15 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/chim_ENS1_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/chim_ENS1_2020020614_1_out_toDART.nc
2026-07-02 13:32:16 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/chim_ENS2_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/chim_ENS2_2020020614_1_out_toDART.nc
2026-07-02 13:32:16 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/chim_ENS3_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/chim_ENS3_2020020614_1_out_toDART.nc
2026-07-02 13:32:16 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/chim_ENS4_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/chim_ENS4_2020020614_1_out_toDART.nc
2026-07-02 13:32:17 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/chim_ENS5_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/chim_ENS5_2020020614_1_out_toDART.nc
2026-07-02 13:32:17 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/chim_ENS6_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/chim_ENS6_2020020614_1_out_toDART.nc
2026-07-02 13:32:17 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/chim_ENS7_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/chim_ENS7_2020020614_1_out_toDART.nc
2026-07-02 13:32:18 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/chim_ENS8_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/chim_ENS8_2020020614_1_out_toDART.nc
2026-07-02 13:32:18 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/chim_ENS9_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/chim_ENS9_2020020614_1_out_toDART.nc
2026-07-02 13:32:18 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/chim_ENS10_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/chim_ENS10_2020020614_1_out_toDART.nc
2026-07-02 13:32:18 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/chim_ENS11_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/chim_ENS11_2020020614_1_out_toDART.nc
2026-07-02 13:32:19 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/chim_ENS12_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/chim_ENS12_2020020614_1_out_toDART.nc
2026-07-02 13:32:19 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/chim_ENS13_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/chim_ENS13_2020020614_1_out_toDART.nc
2026-07-02 13:32:19 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/chim_ENS14_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/chim_ENS14_2020020614_1_out_toDART.nc
2026-07-02 13:32:20 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/chim_ENS15_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/chim_ENS15_2020020614_1_out_toDART.nc
2026-07-02 13:32:20 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-02 13:32:20 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-02 13:32:20 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-02 13:32:20 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-02 13:32:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-02 13:32:20 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-02 13:32:28 INFO Found: []
2026-07-02 13:32:28 INFO No job id returned by command ./run_filter.bsh
2026-07-02 13:32:28 INFO No monitoring will be performed
2026-07-02 13:32:28 INFO Moving DART output files to analysis and preassim directories for date 2020020614 if present ...
2026-07-02 13:32:28 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-07-02 13:32:28 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-07-02 13:32:28 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-07-02 13:32:28 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-07-02 13:32:28 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-07-02 13:32:28 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-07-02 13:32:28 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-07-02 13:32:28 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-07-02 13:32:28 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-07-02 13:32:28 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-07-02 13:32:28 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-07-02 13:32:28 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-07-02 13:32:28 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-07-02 13:32:28 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-07-02 13:32:28 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-07-02 13:32:28 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-07-02 13:32:28 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-07-02 13:32:28 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-07-02 13:32:29 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-07-02 13:32:29 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-07-02 13:32:29 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-07-02 13:32:29 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-07-02 13:32:29 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-07-02 13:32:29 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-07-02 13:32:29 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-07-02 13:32:29 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-07-02 13:32:29 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-07-02 13:32:29 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-07-02 13:32:29 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-07-02 13:32:29 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-07-02 13:32:29 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-07-02 13:32:29 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-07-02 13:32:29 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-07-02 13:32:29 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-07-02 13:32:29 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-02 13:32:29 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-02 13:32:29 INFO run_dart() is DONE.
2026-07-02 13:32:29 INFO ---------->>> Running update_pollutant_in_end()
2026-07-02 13:32:29 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:32:29 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:32:30 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:32:30 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:32:31 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:32:31 INFO [Posterior Diff ENS1] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS1_2020020614.nc
2026-07-02 13:32:31 INFO [Posterior Diff ENS1] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS1_2020020614.relative.nc
2026-07-02 13:32:32 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:32:32 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:32:32 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:32:32 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:32:33 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:32:33 INFO [Posterior Diff ENS2] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS2_2020020614.nc
2026-07-02 13:32:33 INFO [Posterior Diff ENS2] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS2_2020020614.relative.nc
2026-07-02 13:32:34 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:32:34 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:32:35 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:32:35 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:32:36 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:32:36 INFO [Posterior Diff ENS3] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS3_2020020614.nc
2026-07-02 13:32:36 INFO [Posterior Diff ENS3] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS3_2020020614.relative.nc
2026-07-02 13:32:37 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:32:37 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:32:38 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:32:38 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:32:38 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:32:39 INFO [Posterior Diff ENS4] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS4_2020020614.nc
2026-07-02 13:32:39 INFO [Posterior Diff ENS4] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS4_2020020614.relative.nc
2026-07-02 13:32:39 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:32:39 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:32:40 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:32:40 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:32:41 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:32:41 INFO [Posterior Diff ENS5] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS5_2020020614.nc
2026-07-02 13:32:41 INFO [Posterior Diff ENS5] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS5_2020020614.relative.nc
2026-07-02 13:32:42 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:32:42 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:32:43 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:32:43 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:32:43 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:32:43 INFO [Posterior Diff ENS6] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS6_2020020614.nc
2026-07-02 13:32:43 INFO [Posterior Diff ENS6] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS6_2020020614.relative.nc
2026-07-02 13:32:44 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:32:44 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:32:45 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:32:45 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:32:46 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:32:46 INFO [Posterior Diff ENS7] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS7_2020020614.nc
2026-07-02 13:32:46 INFO [Posterior Diff ENS7] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS7_2020020614.relative.nc
2026-07-02 13:32:46 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:32:46 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:32:47 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:32:47 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:32:48 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:32:48 INFO [Posterior Diff ENS8] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS8_2020020614.nc
2026-07-02 13:32:48 INFO [Posterior Diff ENS8] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS8_2020020614.relative.nc
2026-07-02 13:32:48 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:32:49 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:32:50 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:32:50 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:32:50 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:32:51 INFO [Posterior Diff ENS9] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS9_2020020614.nc
2026-07-02 13:32:51 INFO [Posterior Diff ENS9] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS9_2020020614.relative.nc
2026-07-02 13:32:51 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:32:51 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:32:52 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:32:52 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:32:53 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:32:53 INFO [Posterior Diff ENS10] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS10_2020020614.nc
2026-07-02 13:32:53 INFO [Posterior Diff ENS10] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS10_2020020614.relative.nc
2026-07-02 13:32:53 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:32:53 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:32:54 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:32:54 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:32:55 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:32:55 INFO [Posterior Diff ENS11] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS11_2020020614.nc
2026-07-02 13:32:55 INFO [Posterior Diff ENS11] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS11_2020020614.relative.nc
2026-07-02 13:32:55 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:32:55 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:32:56 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:32:56 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:32:57 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:32:57 INFO [Posterior Diff ENS12] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS12_2020020614.nc
2026-07-02 13:32:57 INFO [Posterior Diff ENS12] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS12_2020020614.relative.nc
2026-07-02 13:32:58 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:32:58 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:32:59 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:32:59 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:32:59 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:32:59 INFO [Posterior Diff ENS13] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS13_2020020614.nc
2026-07-02 13:32:59 INFO [Posterior Diff ENS13] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS13_2020020614.relative.nc
2026-07-02 13:33:00 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:33:00 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:33:01 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:33:01 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:33:02 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:33:02 INFO [Posterior Diff ENS14] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS14_2020020614.nc
2026-07-02 13:33:02 INFO [Posterior Diff ENS14] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS14_2020020614.relative.nc
2026-07-02 13:33:02 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:33:02 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:33:03 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 13:33:03 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-07-02 13:33:04 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 13:33:04 INFO [Posterior Diff ENS15] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS15_2020020614.nc
2026-07-02 13:33:04 INFO [Posterior Diff ENS15] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS15_2020020614.relative.nc
2026-07-02 13:33:04 INFO Next run starts from 2020-02-06 14:00:00
2026-07-02 13:33:04 INFO Cycle is DONE; starting a new loop!
2026-07-02 13:33:04 INFO [TIME] step_end current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 13:33:04 INFO [TIME] step_start current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 13:33:04 INFO [TIME] window start=2020-02-06 14:00:00 end=2020-02-07 00:00:00 run_hours=10 has_assimilation=False
2026-07-02 13:33:04 INFO Copying EMIS of next day ...
2026-07-02 13:33:05 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens01.nc
2026-07-02 13:33:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:33:06 INFO Hourly dataset computed and listing created
2026-07-02 13:33:17 INFO Hourly dataset computed
2026-07-02 13:33:17 INFO Copying EMIS of next day ...
2026-07-02 13:33:17 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens02.nc
2026-07-02 13:33:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:33:19 INFO Hourly dataset computed and listing created
2026-07-02 13:33:22 INFO Hourly dataset computed
2026-07-02 13:33:22 INFO Copying EMIS of next day ...
2026-07-02 13:33:22 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens03.nc
2026-07-02 13:33:22 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:33:23 INFO Hourly dataset computed and listing created
2026-07-02 13:33:27 INFO Hourly dataset computed
2026-07-02 13:33:27 INFO Copying EMIS of next day ...
2026-07-02 13:33:27 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens04.nc
2026-07-02 13:33:27 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:33:28 INFO Hourly dataset computed and listing created
2026-07-02 13:33:31 INFO Hourly dataset computed
2026-07-02 13:33:31 INFO Copying EMIS of next day ...
2026-07-02 13:33:32 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens05.nc
2026-07-02 13:33:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:33:33 INFO Hourly dataset computed and listing created
2026-07-02 13:33:37 INFO Hourly dataset computed
2026-07-02 13:33:38 INFO Copying EMIS of next day ...
2026-07-02 13:33:38 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens06.nc
2026-07-02 13:33:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:33:39 INFO Hourly dataset computed and listing created
2026-07-02 13:33:42 INFO Hourly dataset computed
2026-07-02 13:33:42 INFO Copying EMIS of next day ...
2026-07-02 13:33:42 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens07.nc
2026-07-02 13:33:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:33:43 INFO Hourly dataset computed and listing created
2026-07-02 13:33:46 INFO Hourly dataset computed
2026-07-02 13:33:46 INFO Copying EMIS of next day ...
2026-07-02 13:33:47 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens08.nc
2026-07-02 13:33:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:33:48 INFO Hourly dataset computed and listing created
2026-07-02 13:33:50 INFO Hourly dataset computed
2026-07-02 13:33:50 INFO Copying EMIS of next day ...
2026-07-02 13:33:51 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens09.nc
2026-07-02 13:33:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:33:52 INFO Hourly dataset computed and listing created
2026-07-02 13:33:55 INFO Hourly dataset computed
2026-07-02 13:33:55 INFO Copying EMIS of next day ...
2026-07-02 13:33:55 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens10.nc
2026-07-02 13:33:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:33:57 INFO Hourly dataset computed and listing created
2026-07-02 13:33:59 INFO Hourly dataset computed
2026-07-02 13:33:59 INFO Copying EMIS of next day ...
2026-07-02 13:34:00 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens11.nc
2026-07-02 13:34:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:34:01 INFO Hourly dataset computed and listing created
2026-07-02 13:34:03 INFO Hourly dataset computed
2026-07-02 13:34:03 INFO Copying EMIS of next day ...
2026-07-02 13:34:04 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens12.nc
2026-07-02 13:34:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:34:05 INFO Hourly dataset computed and listing created
2026-07-02 13:34:08 INFO Hourly dataset computed
2026-07-02 13:34:08 INFO Copying EMIS of next day ...
2026-07-02 13:34:08 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens13.nc
2026-07-02 13:34:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:34:09 INFO Hourly dataset computed and listing created
2026-07-02 13:34:12 INFO Hourly dataset computed
2026-07-02 13:34:12 INFO Copying EMIS of next day ...
2026-07-02 13:34:13 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens14.nc
2026-07-02 13:34:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:34:14 INFO Hourly dataset computed and listing created
2026-07-02 13:34:17 INFO Hourly dataset computed
2026-07-02 13:34:17 INFO Copying EMIS of next day ...
2026-07-02 13:34:17 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens15.nc
2026-07-02 13:34:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:34:18 INFO Hourly dataset computed and listing created
2026-07-02 13:34:21 INFO Hourly dataset computed
2026-07-02 13:34:21 INFO ---------->>> Running CHIMERE model from 2020-02-06 14:00:00 to 2020-02-07 00:00:00
2026-07-02 13:34:21 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:34:21 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1
2026-07-02 13:34:21 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020613_1_ENS1.nc
2026-07-02 13:34:21 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-02 13:34:21 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:34:21 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-02 13:34:21 INFO Queuing job for member 1...
2026-07-02 13:34:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:34:21 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-02 13:34:22 INFO Found: ['5066059']
2026-07-02 13:34:27 INFO [TGCC-IRENE] Submitted job with ID:['5066059']
2026-07-02 13:34:27 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:34:27 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2
2026-07-02 13:34:27 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020613_1_ENS2.nc
2026-07-02 13:34:27 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-02 13:34:27 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:34:27 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-02 13:34:27 INFO Queuing job for member 2...
2026-07-02 13:34:27 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:34:27 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-02 13:34:27 INFO Found: ['5066060']
2026-07-02 13:34:32 INFO [TGCC-IRENE] Submitted job with ID:['5066060']
2026-07-02 13:34:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:34:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3
2026-07-02 13:34:32 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020613_1_ENS3.nc
2026-07-02 13:34:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-02 13:34:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:34:32 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-02 13:34:32 INFO Queuing job for member 3...
2026-07-02 13:34:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:34:32 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-02 13:34:33 INFO Found: ['5066061']
2026-07-02 13:34:38 INFO [TGCC-IRENE] Submitted job with ID:['5066061']
2026-07-02 13:34:38 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:34:38 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4
2026-07-02 13:34:38 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020613_1_ENS4.nc
2026-07-02 13:34:38 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-02 13:34:38 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:34:38 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-02 13:34:38 INFO Queuing job for member 4...
2026-07-02 13:34:38 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:34:38 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-02 13:34:39 INFO Found: ['5066062']
2026-07-02 13:34:44 INFO [TGCC-IRENE] Submitted job with ID:['5066062']
2026-07-02 13:34:44 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:34:44 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5
2026-07-02 13:34:44 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020613_1_ENS5.nc
2026-07-02 13:34:44 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-02 13:34:44 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:34:44 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-02 13:34:44 INFO Queuing job for member 5...
2026-07-02 13:34:44 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:34:44 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-02 13:34:45 INFO Found: ['5066064']
2026-07-02 13:34:50 INFO [TGCC-IRENE] Submitted job with ID:['5066064']
2026-07-02 13:34:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:34:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6
2026-07-02 13:34:50 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020613_1_ENS6.nc
2026-07-02 13:34:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-02 13:34:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:34:50 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-02 13:34:50 INFO Queuing job for member 6...
2026-07-02 13:34:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:34:50 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-02 13:34:51 INFO Found: ['5066066']
2026-07-02 13:34:56 INFO [TGCC-IRENE] Submitted job with ID:['5066066']
2026-07-02 13:34:56 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:34:56 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7
2026-07-02 13:34:56 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020613_1_ENS7.nc
2026-07-02 13:34:56 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-02 13:34:56 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:34:56 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-02 13:34:56 INFO Queuing job for member 7...
2026-07-02 13:34:56 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:34:56 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-02 13:34:56 INFO Found: ['5066067']
2026-07-02 13:35:01 INFO [TGCC-IRENE] Submitted job with ID:['5066067']
2026-07-02 13:35:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:35:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8
2026-07-02 13:35:01 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020613_1_ENS8.nc
2026-07-02 13:35:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-02 13:35:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:35:01 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-02 13:35:01 INFO Queuing job for member 8...
2026-07-02 13:35:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:35:01 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-02 13:35:02 INFO Found: ['5066070']
2026-07-02 13:35:07 INFO [TGCC-IRENE] Submitted job with ID:['5066070']
2026-07-02 13:35:07 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:35:07 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9
2026-07-02 13:35:07 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020613_1_ENS9.nc
2026-07-02 13:35:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-02 13:35:07 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:35:07 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-02 13:35:07 INFO Queuing job for member 9...
2026-07-02 13:35:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:35:07 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-02 13:35:10 INFO Found: ['5066072']
2026-07-02 13:35:15 INFO [TGCC-IRENE] Submitted job with ID:['5066072']
2026-07-02 13:35:15 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:35:15 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10
2026-07-02 13:35:15 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020613_1_ENS10.nc
2026-07-02 13:35:15 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-02 13:35:15 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:35:15 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-02 13:35:15 INFO Queuing job for member 10...
2026-07-02 13:35:15 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:35:15 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-02 13:35:15 INFO Found: ['5066073']
2026-07-02 13:35:20 INFO [TGCC-IRENE] Submitted job with ID:['5066073']
2026-07-02 13:35:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:35:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11
2026-07-02 13:35:20 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020613_1_ENS11.nc
2026-07-02 13:35:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-02 13:35:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:35:20 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-02 13:35:20 INFO Queuing job for member 11...
2026-07-02 13:35:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:35:20 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-02 13:35:21 INFO Found: ['5066074']
2026-07-02 13:35:26 INFO [TGCC-IRENE] Submitted job with ID:['5066074']
2026-07-02 13:35:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:35:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12
2026-07-02 13:35:26 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020613_1_ENS12.nc
2026-07-02 13:35:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-02 13:35:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:35:26 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-02 13:35:26 INFO Queuing job for member 12...
2026-07-02 13:35:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:35:26 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-02 13:35:27 INFO Found: ['5066075']
2026-07-02 13:35:32 INFO [TGCC-IRENE] Submitted job with ID:['5066075']
2026-07-02 13:35:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:35:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13
2026-07-02 13:35:32 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020613_1_ENS13.nc
2026-07-02 13:35:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-02 13:35:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:35:32 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-02 13:35:32 INFO Queuing job for member 13...
2026-07-02 13:35:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:35:32 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-02 13:35:33 INFO Found: ['5066076']
2026-07-02 13:35:38 INFO [TGCC-IRENE] Submitted job with ID:['5066076']
2026-07-02 13:35:38 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:35:38 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14
2026-07-02 13:35:38 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020613_1_ENS14.nc
2026-07-02 13:35:38 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-02 13:35:38 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:35:38 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-02 13:35:38 INFO Queuing job for member 14...
2026-07-02 13:35:38 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:35:38 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-02 13:35:39 INFO Found: ['5066078']
2026-07-02 13:35:44 INFO [TGCC-IRENE] Submitted job with ID:['5066078']
2026-07-02 13:35:44 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:35:44 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15
2026-07-02 13:35:44 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020613_1_ENS15.nc
2026-07-02 13:35:44 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-02 13:35:44 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:35:44 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-02 13:35:44 INFO Queuing job for member 15...
2026-07-02 13:35:44 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:35:44 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-02 13:35:44 INFO Found: ['5066079']
2026-07-02 13:35:49 INFO [TGCC-IRENE] Submitted job with ID:['5066079']
2026-07-02 13:35:49 INFO Checking job status ...
2026-07-02 13:35:49 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:35:49 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:35:49 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:35:50 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:35:50 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:35:50 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:35:50 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:35:50 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:35:50 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:35:50 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:35:50 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:35:50 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:35:50 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:35:50 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:35:50 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:35:50 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:36:05 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:36:05 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:36:05 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:36:05 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:36:05 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:36:05 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:36:05 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:36:05 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:36:05 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:36:05 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:36:05 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:36:05 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:36:05 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:36:05 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:36:05 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:36:05 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:36:20 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:36:20 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:36:20 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:36:20 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:36:20 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:36:21 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:36:21 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:36:21 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:36:21 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:36:21 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:36:21 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:36:21 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:36:21 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:36:21 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:36:21 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:36:21 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:36:36 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:36:36 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:36:36 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:36:36 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:36:36 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:36:36 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:36:36 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:36:36 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:36:36 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:36:36 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:36:36 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:36:36 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:36:36 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:36:36 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:36:36 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:36:36 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:36:52 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:36:52 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:36:52 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:36:52 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:36:52 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:36:52 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:36:52 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:36:52 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:36:52 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:36:52 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:36:52 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:36:52 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:36:52 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:36:52 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:36:52 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:36:52 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:37:07 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:37:07 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:37:07 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:37:07 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:37:08 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:37:08 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:37:08 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:37:08 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:37:08 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:37:08 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:37:08 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:37:08 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:37:08 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:37:08 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:37:08 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:37:08 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:37:23 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:37:23 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:37:23 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:37:23 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:37:23 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:37:23 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:37:23 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:37:23 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:37:23 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:37:23 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:37:23 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:37:23 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:37:23 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:37:23 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:37:23 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:37:23 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:37:38 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:37:38 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:37:38 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:37:38 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:37:38 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:37:38 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:37:38 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:37:38 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:37:38 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:37:38 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:37:38 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:37:38 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:37:38 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:37:38 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:37:38 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:37:38 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:37:53 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:37:53 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:37:53 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:37:53 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:37:53 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:37:54 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:37:54 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:37:54 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:37:54 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:37:54 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:37:54 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:37:54 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:37:54 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:37:54 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:37:54 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:37:54 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:38:09 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:38:09 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:38:09 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:38:09 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:38:09 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:38:09 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:38:09 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:38:09 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:38:09 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:38:09 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:38:09 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:38:09 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:38:09 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:38:09 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:38:09 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:38:09 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:38:24 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:38:24 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:38:24 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:38:24 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:38:24 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:38:24 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:38:24 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:38:24 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:38:24 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:38:24 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:38:24 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:38:24 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:38:24 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:38:24 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:38:24 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:38:24 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:38:40 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:38:40 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:38:40 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:38:40 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:38:40 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:38:40 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:38:40 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:38:40 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:38:40 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:38:40 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:38:40 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:38:40 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:38:40 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:38:41 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:38:41 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:38:41 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:38:56 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:38:56 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:38:56 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:38:56 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:38:56 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:38:56 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:38:56 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:38:56 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:38:56 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:38:56 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:38:56 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:38:56 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:38:56 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:38:56 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:38:56 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:38:56 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:39:11 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:39:11 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:39:11 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:39:11 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:39:11 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:39:11 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:39:11 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:39:11 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:39:11 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:39:11 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:39:11 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:39:11 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:39:11 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:39:11 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:39:11 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:39:11 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:39:26 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:39:26 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:39:26 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:39:26 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:39:26 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:39:26 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:39:26 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:39:26 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:39:26 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:39:26 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:39:26 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:39:26 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:39:26 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:39:27 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:39:27 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:39:27 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:39:42 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:39:42 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:39:42 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:39:42 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:39:42 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:39:42 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:39:42 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:39:42 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:39:42 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:39:42 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:39:42 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:39:42 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:39:42 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:39:42 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:39:42 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:39:42 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:39:57 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:39:57 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:39:57 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:39:57 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:39:57 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:39:57 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:39:57 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:39:57 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:39:57 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:39:57 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:39:57 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:39:57 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:39:57 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:39:57 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:39:57 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:39:57 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:40:12 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:40:12 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:40:12 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:40:12 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:40:12 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:40:12 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:40:12 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:40:12 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:40:12 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:40:12 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:40:12 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:40:12 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:40:12 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:40:12 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:40:13 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:40:13 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:40:28 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:40:28 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:40:28 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:40:28 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:40:28 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:40:28 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:40:28 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:40:28 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:40:28 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:40:28 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:40:28 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:40:28 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:40:28 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:40:28 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:40:28 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:40:28 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:40:43 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:40:43 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:40:43 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:40:44 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:40:44 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:40:44 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:40:44 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:40:44 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:40:44 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:40:44 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:40:44 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:40:44 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:40:44 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:40:44 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:40:44 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:40:44 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:40:59 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:40:59 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:40:59 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:40:59 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:40:59 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:40:59 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:40:59 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:40:59 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:40:59 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:40:59 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:40:59 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:40:59 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:40:59 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:40:59 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:40:59 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:40:59 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:41:14 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:41:14 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:41:14 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:41:14 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:41:14 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:41:14 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:41:14 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:41:14 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:41:14 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:41:14 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:41:14 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:41:14 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:41:14 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:41:14 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:41:14 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:41:14 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:41:30 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:41:30 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:41:30 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:41:30 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:41:30 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:41:30 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:41:30 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:41:30 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:41:30 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:41:30 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:41:30 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:41:30 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:41:30 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:41:30 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:41:30 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:41:30 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:41:45 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:41:45 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:41:45 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:41:45 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:41:45 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:41:45 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:41:45 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:41:45 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:41:45 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:41:45 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:41:45 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:41:45 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:41:45 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:41:45 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:41:46 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:41:46 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:42:01 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:42:01 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:42:01 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:42:01 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:42:01 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:42:01 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:42:01 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:42:01 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:42:01 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:42:01 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:42:01 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:42:01 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:42:01 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:42:01 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:42:01 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:42:01 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:42:17 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:42:17 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:42:17 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:42:17 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:42:17 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:42:17 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:42:17 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:42:17 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:42:17 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:42:17 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:42:17 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:42:17 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:42:17 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:42:17 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:42:17 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:42:17 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:42:32 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:42:32 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:42:32 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:42:32 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:42:32 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:42:32 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:42:32 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:42:32 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:42:32 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:42:32 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:42:32 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:42:32 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:42:32 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:42:32 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:42:32 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:42:32 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:42:47 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:42:47 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:42:47 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:42:47 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:42:48 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:42:48 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:42:48 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:42:48 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:42:48 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:42:48 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:42:48 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:42:48 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:42:48 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:42:48 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:42:48 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:42:48 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:43:03 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:43:03 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:43:03 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:43:03 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:43:03 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:43:03 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:43:03 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:43:03 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:43:03 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:43:03 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:43:03 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:43:03 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:43:03 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:43:03 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:43:03 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:43:03 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:43:18 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:43:18 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:43:18 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:43:18 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:43:18 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:43:18 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:43:18 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:43:18 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:43:18 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:43:18 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:43:18 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:43:18 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:43:18 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:43:18 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:43:18 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:43:18 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:43:33 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:43:33 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:43:33 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:43:33 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:43:33 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:43:34 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:43:34 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:43:34 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:43:34 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:43:34 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:43:34 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:43:34 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:43:34 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:43:34 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:43:34 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:43:34 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:43:49 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:43:49 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:43:49 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:43:49 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:43:49 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:43:49 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:43:49 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:43:49 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:43:49 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:43:49 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:43:49 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:43:49 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:43:49 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:43:49 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:43:49 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:43:49 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:44:05 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:44:05 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:44:05 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:44:05 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:44:05 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:44:05 INFO None 5066066: status RUNNING/PENDING
2026-07-02 13:44:05 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:44:05 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:44:05 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:44:05 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:44:05 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:44:05 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:44:05 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:44:05 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:44:05 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:44:05 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:44:20 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:44:20 INFO None 5066060: status RUNNING/PENDING
2026-07-02 13:44:20 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:44:20 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:44:20 INFO None 5066064: status RUNNING/PENDING
2026-07-02 13:44:20 INFO None 5066066: status FINISHED
2026-07-02 13:44:20 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:44:20 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:44:20 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:44:20 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:44:20 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:44:20 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:44:20 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:44:20 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:44:20 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:44:20 INFO Jobs still running: ['5066059', '5066060', '5066061', '5066062', '5066064', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:44:35 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:44:35 INFO None 5066060: status FINISHED
2026-07-02 13:44:35 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:44:35 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:44:35 INFO None 5066064: status FINISHED
2026-07-02 13:44:35 INFO None 5066066: status FINISHED
2026-07-02 13:44:35 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:44:35 INFO None 5066070: status RUNNING/PENDING
2026-07-02 13:44:35 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:44:35 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:44:36 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:44:36 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:44:36 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:44:36 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:44:36 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:44:36 INFO Jobs still running: ['5066059', '5066061', '5066062', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:44:51 INFO None 5066059: status RUNNING/PENDING
2026-07-02 13:44:51 INFO None 5066060: status FINISHED
2026-07-02 13:44:51 INFO None 5066061: status RUNNING/PENDING
2026-07-02 13:44:51 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:44:51 INFO None 5066064: status FINISHED
2026-07-02 13:44:51 INFO None 5066066: status FINISHED
2026-07-02 13:44:51 INFO None 5066067: status RUNNING/PENDING
2026-07-02 13:44:51 INFO None 5066070: status FINISHED
2026-07-02 13:44:51 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:44:51 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:44:51 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:44:51 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:44:51 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:44:51 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:44:51 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:44:51 INFO Jobs still running: ['5066059', '5066061', '5066062', '5066067', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:45:06 INFO None 5066059: status FINISHED
2026-07-02 13:45:06 INFO None 5066060: status FINISHED
2026-07-02 13:45:06 INFO None 5066061: status FINISHED
2026-07-02 13:45:06 INFO None 5066062: status RUNNING/PENDING
2026-07-02 13:45:06 INFO None 5066064: status FINISHED
2026-07-02 13:45:06 INFO None 5066066: status FINISHED
2026-07-02 13:45:06 INFO None 5066067: status FINISHED
2026-07-02 13:45:06 INFO None 5066070: status FINISHED
2026-07-02 13:45:06 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:45:06 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:45:06 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:45:06 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:45:06 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:45:06 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:45:06 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:45:06 INFO Jobs still running: ['5066062', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:45:21 INFO None 5066059: status FINISHED
2026-07-02 13:45:21 INFO None 5066060: status FINISHED
2026-07-02 13:45:21 INFO None 5066061: status FINISHED
2026-07-02 13:45:21 INFO None 5066062: status FINISHED
2026-07-02 13:45:21 INFO None 5066064: status FINISHED
2026-07-02 13:45:21 INFO None 5066066: status FINISHED
2026-07-02 13:45:21 INFO None 5066067: status FINISHED
2026-07-02 13:45:21 INFO None 5066070: status FINISHED
2026-07-02 13:45:21 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:45:21 INFO None 5066073: status RUNNING/PENDING
2026-07-02 13:45:21 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:45:22 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:45:22 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:45:22 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:45:22 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:45:22 INFO Jobs still running: ['5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:45:37 INFO None 5066059: status FINISHED
2026-07-02 13:45:37 INFO None 5066060: status FINISHED
2026-07-02 13:45:37 INFO None 5066061: status FINISHED
2026-07-02 13:45:37 INFO None 5066062: status FINISHED
2026-07-02 13:45:37 INFO None 5066064: status FINISHED
2026-07-02 13:45:37 INFO None 5066066: status FINISHED
2026-07-02 13:45:37 INFO None 5066067: status FINISHED
2026-07-02 13:45:37 INFO None 5066070: status FINISHED
2026-07-02 13:45:37 INFO None 5066072: status RUNNING/PENDING
2026-07-02 13:45:37 INFO None 5066073: status FINISHED
2026-07-02 13:45:37 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:45:37 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:45:37 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:45:37 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:45:37 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:45:37 INFO Jobs still running: ['5066072', '5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:45:52 INFO None 5066059: status FINISHED
2026-07-02 13:45:52 INFO None 5066060: status FINISHED
2026-07-02 13:45:52 INFO None 5066061: status FINISHED
2026-07-02 13:45:52 INFO None 5066062: status FINISHED
2026-07-02 13:45:52 INFO None 5066064: status FINISHED
2026-07-02 13:45:52 INFO None 5066066: status FINISHED
2026-07-02 13:45:52 INFO None 5066067: status FINISHED
2026-07-02 13:45:52 INFO None 5066070: status FINISHED
2026-07-02 13:45:52 INFO None 5066072: status FINISHED
2026-07-02 13:45:52 INFO None 5066073: status FINISHED
2026-07-02 13:45:52 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:45:52 INFO None 5066075: status RUNNING/PENDING
2026-07-02 13:45:52 INFO None 5066076: status RUNNING/PENDING
2026-07-02 13:45:52 INFO None 5066078: status RUNNING/PENDING
2026-07-02 13:45:52 INFO None 5066079: status RUNNING/PENDING
2026-07-02 13:45:52 INFO Jobs still running: ['5066074', '5066075', '5066076', '5066078', '5066079']. Waiting...
2026-07-02 13:46:07 INFO None 5066059: status FINISHED
2026-07-02 13:46:07 INFO None 5066060: status FINISHED
2026-07-02 13:46:07 INFO None 5066061: status FINISHED
2026-07-02 13:46:07 INFO None 5066062: status FINISHED
2026-07-02 13:46:07 INFO None 5066064: status FINISHED
2026-07-02 13:46:07 INFO None 5066066: status FINISHED
2026-07-02 13:46:07 INFO None 5066067: status FINISHED
2026-07-02 13:46:07 INFO None 5066070: status FINISHED
2026-07-02 13:46:07 INFO None 5066072: status FINISHED
2026-07-02 13:46:07 INFO None 5066073: status FINISHED
2026-07-02 13:46:07 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:46:07 INFO None 5066075: status FINISHED
2026-07-02 13:46:07 INFO None 5066076: status FINISHED
2026-07-02 13:46:07 INFO None 5066078: status FINISHED
2026-07-02 13:46:08 INFO None 5066079: status FINISHED
2026-07-02 13:46:08 INFO Jobs still running: ['5066074']. Waiting...
2026-07-02 13:46:23 INFO None 5066059: status FINISHED
2026-07-02 13:46:23 INFO None 5066060: status FINISHED
2026-07-02 13:46:23 INFO None 5066061: status FINISHED
2026-07-02 13:46:23 INFO None 5066062: status FINISHED
2026-07-02 13:46:23 INFO None 5066064: status FINISHED
2026-07-02 13:46:23 INFO None 5066066: status FINISHED
2026-07-02 13:46:23 INFO None 5066067: status FINISHED
2026-07-02 13:46:23 INFO None 5066070: status FINISHED
2026-07-02 13:46:23 INFO None 5066072: status FINISHED
2026-07-02 13:46:23 INFO None 5066073: status FINISHED
2026-07-02 13:46:23 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:46:23 INFO None 5066075: status FINISHED
2026-07-02 13:46:23 INFO None 5066076: status FINISHED
2026-07-02 13:46:23 INFO None 5066078: status FINISHED
2026-07-02 13:46:23 INFO None 5066079: status FINISHED
2026-07-02 13:46:23 INFO Jobs still running: ['5066074']. Waiting...
2026-07-02 13:46:38 INFO None 5066059: status FINISHED
2026-07-02 13:46:38 INFO None 5066060: status FINISHED
2026-07-02 13:46:38 INFO None 5066061: status FINISHED
2026-07-02 13:46:38 INFO None 5066062: status FINISHED
2026-07-02 13:46:38 INFO None 5066064: status FINISHED
2026-07-02 13:46:38 INFO None 5066066: status FINISHED
2026-07-02 13:46:38 INFO None 5066067: status FINISHED
2026-07-02 13:46:38 INFO None 5066070: status FINISHED
2026-07-02 13:46:38 INFO None 5066072: status FINISHED
2026-07-02 13:46:38 INFO None 5066073: status FINISHED
2026-07-02 13:46:38 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:46:38 INFO None 5066075: status FINISHED
2026-07-02 13:46:38 INFO None 5066076: status FINISHED
2026-07-02 13:46:38 INFO None 5066078: status FINISHED
2026-07-02 13:46:38 INFO None 5066079: status FINISHED
2026-07-02 13:46:38 INFO Jobs still running: ['5066074']. Waiting...
2026-07-02 13:46:53 INFO None 5066059: status FINISHED
2026-07-02 13:46:53 INFO None 5066060: status FINISHED
2026-07-02 13:46:53 INFO None 5066061: status FINISHED
2026-07-02 13:46:53 INFO None 5066062: status FINISHED
2026-07-02 13:46:53 INFO None 5066064: status FINISHED
2026-07-02 13:46:53 INFO None 5066066: status FINISHED
2026-07-02 13:46:53 INFO None 5066067: status FINISHED
2026-07-02 13:46:54 INFO None 5066070: status FINISHED
2026-07-02 13:46:54 INFO None 5066072: status FINISHED
2026-07-02 13:46:54 INFO None 5066073: status FINISHED
2026-07-02 13:46:54 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:46:54 INFO None 5066075: status FINISHED
2026-07-02 13:46:54 INFO None 5066076: status FINISHED
2026-07-02 13:46:54 INFO None 5066078: status FINISHED
2026-07-02 13:46:54 INFO None 5066079: status FINISHED
2026-07-02 13:46:54 INFO Jobs still running: ['5066074']. Waiting...
2026-07-02 13:47:09 INFO None 5066059: status FINISHED
2026-07-02 13:47:09 INFO None 5066060: status FINISHED
2026-07-02 13:47:09 INFO None 5066061: status FINISHED
2026-07-02 13:47:09 INFO None 5066062: status FINISHED
2026-07-02 13:47:09 INFO None 5066064: status FINISHED
2026-07-02 13:47:09 INFO None 5066066: status FINISHED
2026-07-02 13:47:09 INFO None 5066067: status FINISHED
2026-07-02 13:47:09 INFO None 5066070: status FINISHED
2026-07-02 13:47:09 INFO None 5066072: status FINISHED
2026-07-02 13:47:09 INFO None 5066073: status FINISHED
2026-07-02 13:47:09 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:47:09 INFO None 5066075: status FINISHED
2026-07-02 13:47:09 INFO None 5066076: status FINISHED
2026-07-02 13:47:09 INFO None 5066078: status FINISHED
2026-07-02 13:47:09 INFO None 5066079: status FINISHED
2026-07-02 13:47:09 INFO Jobs still running: ['5066074']. Waiting...
2026-07-02 13:47:24 INFO None 5066059: status FINISHED
2026-07-02 13:47:24 INFO None 5066060: status FINISHED
2026-07-02 13:47:24 INFO None 5066061: status FINISHED
2026-07-02 13:47:24 INFO None 5066062: status FINISHED
2026-07-02 13:47:24 INFO None 5066064: status FINISHED
2026-07-02 13:47:24 INFO None 5066066: status FINISHED
2026-07-02 13:47:24 INFO None 5066067: status FINISHED
2026-07-02 13:47:24 INFO None 5066070: status FINISHED
2026-07-02 13:47:24 INFO None 5066072: status FINISHED
2026-07-02 13:47:24 INFO None 5066073: status FINISHED
2026-07-02 13:47:24 INFO None 5066074: status RUNNING/PENDING
2026-07-02 13:47:24 INFO None 5066075: status FINISHED
2026-07-02 13:47:24 INFO None 5066076: status FINISHED
2026-07-02 13:47:24 INFO None 5066078: status FINISHED
2026-07-02 13:47:24 INFO None 5066079: status FINISHED
2026-07-02 13:47:24 INFO Jobs still running: ['5066074']. Waiting...
2026-07-02 13:47:39 INFO None 5066059: status FINISHED
2026-07-02 13:47:39 INFO None 5066060: status FINISHED
2026-07-02 13:47:39 INFO None 5066061: status FINISHED
2026-07-02 13:47:39 INFO None 5066062: status FINISHED
2026-07-02 13:47:39 INFO None 5066064: status FINISHED
2026-07-02 13:47:39 INFO None 5066066: status FINISHED
2026-07-02 13:47:39 INFO None 5066067: status FINISHED
2026-07-02 13:47:39 INFO None 5066070: status FINISHED
2026-07-02 13:47:39 INFO None 5066072: status FINISHED
2026-07-02 13:47:39 INFO None 5066073: status FINISHED
2026-07-02 13:47:39 INFO None 5066074: status FINISHED
2026-07-02 13:47:39 INFO None 5066075: status FINISHED
2026-07-02 13:47:40 INFO None 5066076: status FINISHED
2026-07-02 13:47:40 INFO None 5066078: status FINISHED
2026-07-02 13:47:40 INFO None 5066079: status FINISHED
2026-07-02 13:47:40 INFO Jobs ['5066059', '5066060', '5066061', '5066062', '5066064', '5066066', '5066067', '5066070', '5066072', '5066073', '5066074', '5066075', '5066076', '5066078', '5066079'] have finished
2026-07-02 13:47:40 INFO Checking restart files were created ...
2026-07-02 13:47:40 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020614_10_ENS1.nc(3673513755 bytes)
2026-07-02 13:47:40 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020614_10_ENS2.nc(3673513755 bytes)
2026-07-02 13:47:40 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020614_10_ENS3.nc(3673513755 bytes)
2026-07-02 13:47:40 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020614_10_ENS4.nc(3673513755 bytes)
2026-07-02 13:47:40 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020614_10_ENS5.nc(3673513755 bytes)
2026-07-02 13:47:40 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020614_10_ENS6.nc(3673513755 bytes)
2026-07-02 13:47:40 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020614_10_ENS7.nc(3673513755 bytes)
2026-07-02 13:47:40 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020614_10_ENS8.nc(3673513755 bytes)
2026-07-02 13:47:40 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020614_10_ENS9.nc(3673513755 bytes)
2026-07-02 13:47:40 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020614_10_ENS10.nc(3673513755 bytes)
2026-07-02 13:47:40 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020614_10_ENS11.nc(3673513755 bytes)
2026-07-02 13:47:40 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020614_10_ENS12.nc(3673513755 bytes)
2026-07-02 13:47:40 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020614_10_ENS13.nc(3673513755 bytes)
2026-07-02 13:47:40 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020614_10_ENS14.nc(3673513755 bytes)
2026-07-02 13:47:40 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020614_10_ENS15.nc(3673513755 bytes)
2026-07-02 13:47:40 INFO  Run_model() completed successfully.
2026-07-02 13:47:40 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 14:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 13:47:40 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 00:00:00 days=153073 seconds=0
2026-07-02 13:47:40 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-02 13:47:40 INFO [TIME] increment current_time 2020-02-06 14:00:00 -> 2020-02-07 00:00:00
2026-07-02 13:47:40 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 00:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 13:47:40 INFO ---------->>> Running process_satellite_data()
2026-07-02 13:47:40 INFO [DART] No satellite data found, skipping assimilation
2026-07-02 13:47:40 INFO after_assimilation() skipped
2026-07-02 13:47:40 INFO Next run starts from 2020-02-07 00:00:00
2026-07-02 13:47:40 INFO Cycle is DONE; starting a new loop!
2026-07-02 13:47:40 INFO [TIME] step_end current_time=2020-02-07 00:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 13:47:40 INFO [TIME] step_start current_time=2020-02-07 00:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 13:47:40 INFO [TIME] window start=2020-02-07 00:00:00 end=2020-02-07 01:00:00 run_hours=1 has_assimilation=False
2026-07-02 13:47:40 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:47:41 INFO Hourly dataset computed and listing created
2026-07-02 13:47:44 INFO Hourly dataset computed
2026-07-02 13:47:44 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:47:45 INFO Hourly dataset computed and listing created
2026-07-02 13:47:46 INFO Hourly dataset computed
2026-07-02 13:47:46 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:47:47 INFO Hourly dataset computed and listing created
2026-07-02 13:47:47 INFO Hourly dataset computed
2026-07-02 13:47:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:47:48 INFO Hourly dataset computed and listing created
2026-07-02 13:47:48 INFO Hourly dataset computed
2026-07-02 13:47:49 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:47:49 INFO Hourly dataset computed and listing created
2026-07-02 13:47:50 INFO Hourly dataset computed
2026-07-02 13:47:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:47:51 INFO Hourly dataset computed and listing created
2026-07-02 13:47:51 INFO Hourly dataset computed
2026-07-02 13:47:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:47:52 INFO Hourly dataset computed and listing created
2026-07-02 13:47:53 INFO Hourly dataset computed
2026-07-02 13:47:53 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:47:54 INFO Hourly dataset computed and listing created
2026-07-02 13:47:54 INFO Hourly dataset computed
2026-07-02 13:47:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:47:55 INFO Hourly dataset computed and listing created
2026-07-02 13:47:56 INFO Hourly dataset computed
2026-07-02 13:47:56 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:47:57 INFO Hourly dataset computed and listing created
2026-07-02 13:47:57 INFO Hourly dataset computed
2026-07-02 13:47:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:47:58 INFO Hourly dataset computed and listing created
2026-07-02 13:47:59 INFO Hourly dataset computed
2026-07-02 13:47:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:48:00 INFO Hourly dataset computed and listing created
2026-07-02 13:48:00 INFO Hourly dataset computed
2026-07-02 13:48:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:48:01 INFO Hourly dataset computed and listing created
2026-07-02 13:48:02 INFO Hourly dataset computed
2026-07-02 13:48:02 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:48:03 INFO Hourly dataset computed and listing created
2026-07-02 13:48:04 INFO Hourly dataset computed
2026-07-02 13:48:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:48:05 INFO Hourly dataset computed and listing created
2026-07-02 13:48:05 INFO Hourly dataset computed
2026-07-02 13:48:05 INFO ---------->>> Running CHIMERE model from 2020-02-07 00:00:00 to 2020-02-07 01:00:00
2026-07-02 13:48:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:48:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1
2026-07-02 13:48:05 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020614_10_ENS1.nc
2026-07-02 13:48:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-02 13:48:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:48:05 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-02 13:48:05 INFO Queuing job for member 1...
2026-07-02 13:48:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:48:05 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-02 13:48:06 INFO Found: ['5066131']
2026-07-02 13:48:11 INFO [TGCC-IRENE] Submitted job with ID:['5066131']
2026-07-02 13:48:11 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:48:11 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2
2026-07-02 13:48:11 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020614_10_ENS2.nc
2026-07-02 13:48:11 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-02 13:48:11 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:48:11 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-02 13:48:11 INFO Queuing job for member 2...
2026-07-02 13:48:11 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:48:11 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-02 13:48:12 INFO Found: ['5066132']
2026-07-02 13:48:17 INFO [TGCC-IRENE] Submitted job with ID:['5066132']
2026-07-02 13:48:17 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:48:17 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3
2026-07-02 13:48:17 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020614_10_ENS3.nc
2026-07-02 13:48:17 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-02 13:48:17 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:48:17 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-02 13:48:17 INFO Queuing job for member 3...
2026-07-02 13:48:17 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:48:17 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-02 13:48:19 INFO Found: ['5066133']
2026-07-02 13:48:24 INFO [TGCC-IRENE] Submitted job with ID:['5066133']
2026-07-02 13:48:24 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:48:24 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4
2026-07-02 13:48:24 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020614_10_ENS4.nc
2026-07-02 13:48:24 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-02 13:48:24 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:48:24 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-02 13:48:24 INFO Queuing job for member 4...
2026-07-02 13:48:24 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:48:24 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-02 13:48:25 INFO Found: ['5066134']
2026-07-02 13:48:30 INFO [TGCC-IRENE] Submitted job with ID:['5066134']
2026-07-02 13:48:30 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:48:30 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5
2026-07-02 13:48:30 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020614_10_ENS5.nc
2026-07-02 13:48:30 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-02 13:48:30 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:48:30 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-02 13:48:30 INFO Queuing job for member 5...
2026-07-02 13:48:30 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:48:30 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-02 13:48:31 INFO Found: ['5066136']
2026-07-02 13:48:36 INFO [TGCC-IRENE] Submitted job with ID:['5066136']
2026-07-02 13:48:36 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:48:36 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6
2026-07-02 13:48:36 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020614_10_ENS6.nc
2026-07-02 13:48:36 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-02 13:48:36 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:48:36 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-02 13:48:36 INFO Queuing job for member 6...
2026-07-02 13:48:36 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:48:36 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-02 13:48:37 INFO Found: ['5066137']
2026-07-02 13:48:42 INFO [TGCC-IRENE] Submitted job with ID:['5066137']
2026-07-02 13:48:42 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:48:42 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7
2026-07-02 13:48:42 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020614_10_ENS7.nc
2026-07-02 13:48:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-02 13:48:42 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:48:42 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-02 13:48:42 INFO Queuing job for member 7...
2026-07-02 13:48:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:48:42 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-02 13:48:42 INFO Found: ['5066139']
2026-07-02 13:48:47 INFO [TGCC-IRENE] Submitted job with ID:['5066139']
2026-07-02 13:48:47 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:48:47 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8
2026-07-02 13:48:47 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020614_10_ENS8.nc
2026-07-02 13:48:47 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-02 13:48:47 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:48:47 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-02 13:48:47 INFO Queuing job for member 8...
2026-07-02 13:48:47 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:48:47 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-02 13:48:48 INFO Found: ['5066143']
2026-07-02 13:48:53 INFO [TGCC-IRENE] Submitted job with ID:['5066143']
2026-07-02 13:48:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:48:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9
2026-07-02 13:48:53 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020614_10_ENS9.nc
2026-07-02 13:48:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-02 13:48:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:48:53 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-02 13:48:53 INFO Queuing job for member 9...
2026-07-02 13:48:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:48:53 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-02 13:48:54 INFO Found: ['5066168']
2026-07-02 13:48:59 INFO [TGCC-IRENE] Submitted job with ID:['5066168']
2026-07-02 13:48:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:48:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10
2026-07-02 13:48:59 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020614_10_ENS10.nc
2026-07-02 13:48:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-02 13:48:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:48:59 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-02 13:48:59 INFO Queuing job for member 10...
2026-07-02 13:48:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:48:59 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-02 13:49:00 INFO Found: ['5066179']
2026-07-02 13:49:05 INFO [TGCC-IRENE] Submitted job with ID:['5066179']
2026-07-02 13:49:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:49:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11
2026-07-02 13:49:05 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020614_10_ENS11.nc
2026-07-02 13:49:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-02 13:49:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:49:05 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-02 13:49:05 INFO Queuing job for member 11...
2026-07-02 13:49:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:49:05 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-02 13:49:07 INFO Found: ['5066182']
2026-07-02 13:49:12 INFO [TGCC-IRENE] Submitted job with ID:['5066182']
2026-07-02 13:49:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:49:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12
2026-07-02 13:49:12 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020614_10_ENS12.nc
2026-07-02 13:49:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-02 13:49:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:49:12 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-02 13:49:13 INFO Queuing job for member 12...
2026-07-02 13:49:13 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:49:13 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-02 13:49:13 INFO Found: ['5066183']
2026-07-02 13:49:18 INFO [TGCC-IRENE] Submitted job with ID:['5066183']
2026-07-02 13:49:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:49:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13
2026-07-02 13:49:18 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020614_10_ENS13.nc
2026-07-02 13:49:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-02 13:49:18 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:49:18 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-02 13:49:18 INFO Queuing job for member 13...
2026-07-02 13:49:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:49:18 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-02 13:49:19 INFO Found: ['5066184']
2026-07-02 13:49:24 INFO [TGCC-IRENE] Submitted job with ID:['5066184']
2026-07-02 13:49:24 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:49:24 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14
2026-07-02 13:49:24 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020614_10_ENS14.nc
2026-07-02 13:49:24 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-02 13:49:24 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:49:24 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-02 13:49:24 INFO Queuing job for member 14...
2026-07-02 13:49:24 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:49:24 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-02 13:49:25 INFO Found: ['5066185']
2026-07-02 13:49:30 INFO [TGCC-IRENE] Submitted job with ID:['5066185']
2026-07-02 13:49:30 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:49:30 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15
2026-07-02 13:49:30 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020614_10_ENS15.nc
2026-07-02 13:49:30 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-02 13:49:30 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:49:30 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-02 13:49:30 INFO Queuing job for member 15...
2026-07-02 13:49:30 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:49:30 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-02 13:49:31 INFO Found: ['5066186']
2026-07-02 13:49:36 INFO [TGCC-IRENE] Submitted job with ID:['5066186']
2026-07-02 13:49:36 INFO Checking job status ...
2026-07-02 13:49:36 INFO None 5066131: status RUNNING/PENDING
2026-07-02 13:49:36 INFO None 5066132: status RUNNING/PENDING
2026-07-02 13:49:36 INFO None 5066133: status RUNNING/PENDING
2026-07-02 13:49:36 INFO None 5066134: status RUNNING/PENDING
2026-07-02 13:49:36 INFO None 5066136: status RUNNING/PENDING
2026-07-02 13:49:36 INFO None 5066137: status RUNNING/PENDING
2026-07-02 13:49:36 INFO None 5066139: status RUNNING/PENDING
2026-07-02 13:49:36 INFO None 5066143: status RUNNING/PENDING
2026-07-02 13:49:36 INFO None 5066168: status RUNNING/PENDING
2026-07-02 13:49:36 INFO None 5066179: status RUNNING/PENDING
2026-07-02 13:49:36 INFO None 5066182: status RUNNING/PENDING
2026-07-02 13:49:36 INFO None 5066183: status RUNNING/PENDING
2026-07-02 13:49:36 INFO None 5066184: status RUNNING/PENDING
2026-07-02 13:49:36 INFO None 5066185: status RUNNING/PENDING
2026-07-02 13:49:36 INFO None 5066186: status RUNNING/PENDING
2026-07-02 13:49:36 INFO Jobs still running: ['5066131', '5066132', '5066133', '5066134', '5066136', '5066137', '5066139', '5066143', '5066168', '5066179', '5066182', '5066183', '5066184', '5066185', '5066186']. Waiting...
2026-07-02 13:49:51 INFO None 5066131: status RUNNING/PENDING
2026-07-02 13:49:51 INFO None 5066132: status RUNNING/PENDING
2026-07-02 13:49:51 INFO None 5066133: status RUNNING/PENDING
2026-07-02 13:49:51 INFO None 5066134: status RUNNING/PENDING
2026-07-02 13:49:51 INFO None 5066136: status RUNNING/PENDING
2026-07-02 13:49:51 INFO None 5066137: status RUNNING/PENDING
2026-07-02 13:49:51 INFO None 5066139: status RUNNING/PENDING
2026-07-02 13:49:51 INFO None 5066143: status RUNNING/PENDING
2026-07-02 13:49:51 INFO None 5066168: status RUNNING/PENDING
2026-07-02 13:49:51 INFO None 5066179: status RUNNING/PENDING
2026-07-02 13:49:51 INFO None 5066182: status RUNNING/PENDING
2026-07-02 13:49:51 INFO None 5066183: status RUNNING/PENDING
2026-07-02 13:49:51 INFO None 5066184: status RUNNING/PENDING
2026-07-02 13:49:51 INFO None 5066185: status RUNNING/PENDING
2026-07-02 13:49:51 INFO None 5066186: status RUNNING/PENDING
2026-07-02 13:49:51 INFO Jobs still running: ['5066131', '5066132', '5066133', '5066134', '5066136', '5066137', '5066139', '5066143', '5066168', '5066179', '5066182', '5066183', '5066184', '5066185', '5066186']. Waiting...
2026-07-02 13:50:06 INFO None 5066131: status RUNNING/PENDING
2026-07-02 13:50:06 INFO None 5066132: status FINISHED
2026-07-02 13:50:06 INFO None 5066133: status RUNNING/PENDING
2026-07-02 13:50:06 INFO None 5066134: status RUNNING/PENDING
2026-07-02 13:50:06 INFO None 5066136: status RUNNING/PENDING
2026-07-02 13:50:06 INFO None 5066137: status RUNNING/PENDING
2026-07-02 13:50:07 INFO None 5066139: status RUNNING/PENDING
2026-07-02 13:50:07 INFO None 5066143: status RUNNING/PENDING
2026-07-02 13:50:07 INFO None 5066168: status RUNNING/PENDING
2026-07-02 13:50:07 INFO None 5066179: status RUNNING/PENDING
2026-07-02 13:50:07 INFO None 5066182: status RUNNING/PENDING
2026-07-02 13:50:07 INFO None 5066183: status RUNNING/PENDING
2026-07-02 13:50:07 INFO None 5066184: status RUNNING/PENDING
2026-07-02 13:50:07 INFO None 5066185: status RUNNING/PENDING
2026-07-02 13:50:07 INFO None 5066186: status RUNNING/PENDING
2026-07-02 13:50:07 INFO Jobs still running: ['5066131', '5066133', '5066134', '5066136', '5066137', '5066139', '5066143', '5066168', '5066179', '5066182', '5066183', '5066184', '5066185', '5066186']. Waiting...
2026-07-02 13:50:22 INFO None 5066131: status FINISHED
2026-07-02 13:50:22 INFO None 5066132: status FINISHED
2026-07-02 13:50:22 INFO None 5066133: status RUNNING/PENDING
2026-07-02 13:50:22 INFO None 5066134: status RUNNING/PENDING
2026-07-02 13:50:22 INFO None 5066136: status RUNNING/PENDING
2026-07-02 13:50:22 INFO None 5066137: status RUNNING/PENDING
2026-07-02 13:50:22 INFO None 5066139: status RUNNING/PENDING
2026-07-02 13:50:22 INFO None 5066143: status RUNNING/PENDING
2026-07-02 13:50:22 INFO None 5066168: status RUNNING/PENDING
2026-07-02 13:50:22 INFO None 5066179: status RUNNING/PENDING
2026-07-02 13:50:22 INFO None 5066182: status RUNNING/PENDING
2026-07-02 13:50:22 INFO None 5066183: status RUNNING/PENDING
2026-07-02 13:50:22 INFO None 5066184: status RUNNING/PENDING
2026-07-02 13:50:22 INFO None 5066185: status RUNNING/PENDING
2026-07-02 13:50:22 INFO None 5066186: status RUNNING/PENDING
2026-07-02 13:50:22 INFO Jobs still running: ['5066133', '5066134', '5066136', '5066137', '5066139', '5066143', '5066168', '5066179', '5066182', '5066183', '5066184', '5066185', '5066186']. Waiting...
2026-07-02 13:50:37 INFO None 5066131: status FINISHED
2026-07-02 13:50:37 INFO None 5066132: status FINISHED
2026-07-02 13:50:37 INFO None 5066133: status RUNNING/PENDING
2026-07-02 13:50:37 INFO None 5066134: status RUNNING/PENDING
2026-07-02 13:50:37 INFO None 5066136: status RUNNING/PENDING
2026-07-02 13:50:37 INFO None 5066137: status RUNNING/PENDING
2026-07-02 13:50:37 INFO None 5066139: status RUNNING/PENDING
2026-07-02 13:50:37 INFO None 5066143: status RUNNING/PENDING
2026-07-02 13:50:37 INFO None 5066168: status RUNNING/PENDING
2026-07-02 13:50:37 INFO None 5066179: status RUNNING/PENDING
2026-07-02 13:50:37 INFO None 5066182: status RUNNING/PENDING
2026-07-02 13:50:37 INFO None 5066183: status RUNNING/PENDING
2026-07-02 13:50:37 INFO None 5066184: status RUNNING/PENDING
2026-07-02 13:50:37 INFO None 5066185: status RUNNING/PENDING
2026-07-02 13:50:37 INFO None 5066186: status RUNNING/PENDING
2026-07-02 13:50:37 INFO Jobs still running: ['5066133', '5066134', '5066136', '5066137', '5066139', '5066143', '5066168', '5066179', '5066182', '5066183', '5066184', '5066185', '5066186']. Waiting...
2026-07-02 13:50:53 INFO None 5066131: status FINISHED
2026-07-02 13:50:53 INFO None 5066132: status FINISHED
2026-07-02 13:50:53 INFO None 5066133: status FINISHED
2026-07-02 13:50:53 INFO None 5066134: status RUNNING/PENDING
2026-07-02 13:50:53 INFO None 5066136: status RUNNING/PENDING
2026-07-02 13:50:53 INFO None 5066137: status RUNNING/PENDING
2026-07-02 13:50:53 INFO None 5066139: status FINISHED
2026-07-02 13:50:53 INFO None 5066143: status RUNNING/PENDING
2026-07-02 13:50:53 INFO None 5066168: status RUNNING/PENDING
2026-07-02 13:50:53 INFO None 5066179: status RUNNING/PENDING
2026-07-02 13:50:53 INFO None 5066182: status RUNNING/PENDING
2026-07-02 13:50:53 INFO None 5066183: status RUNNING/PENDING
2026-07-02 13:50:53 INFO None 5066184: status RUNNING/PENDING
2026-07-02 13:50:53 INFO None 5066185: status RUNNING/PENDING
2026-07-02 13:50:53 INFO None 5066186: status RUNNING/PENDING
2026-07-02 13:50:53 INFO Jobs still running: ['5066134', '5066136', '5066137', '5066143', '5066168', '5066179', '5066182', '5066183', '5066184', '5066185', '5066186']. Waiting...
2026-07-02 13:51:08 INFO None 5066131: status FINISHED
2026-07-02 13:51:08 INFO None 5066132: status FINISHED
2026-07-02 13:51:08 INFO None 5066133: status FINISHED
2026-07-02 13:51:08 INFO None 5066134: status RUNNING/PENDING
2026-07-02 13:51:08 INFO None 5066136: status RUNNING/PENDING
2026-07-02 13:51:08 INFO None 5066137: status RUNNING/PENDING
2026-07-02 13:51:08 INFO None 5066139: status FINISHED
2026-07-02 13:51:08 INFO None 5066143: status RUNNING/PENDING
2026-07-02 13:51:08 INFO None 5066168: status RUNNING/PENDING
2026-07-02 13:51:08 INFO None 5066179: status RUNNING/PENDING
2026-07-02 13:51:08 INFO None 5066182: status RUNNING/PENDING
2026-07-02 13:51:08 INFO None 5066183: status RUNNING/PENDING
2026-07-02 13:51:08 INFO None 5066184: status RUNNING/PENDING
2026-07-02 13:51:08 INFO None 5066185: status RUNNING/PENDING
2026-07-02 13:51:08 INFO None 5066186: status RUNNING/PENDING
2026-07-02 13:51:08 INFO Jobs still running: ['5066134', '5066136', '5066137', '5066143', '5066168', '5066179', '5066182', '5066183', '5066184', '5066185', '5066186']. Waiting...
2026-07-02 13:51:23 INFO None 5066131: status FINISHED
2026-07-02 13:51:23 INFO None 5066132: status FINISHED
2026-07-02 13:51:23 INFO None 5066133: status FINISHED
2026-07-02 13:51:23 INFO None 5066134: status FINISHED
2026-07-02 13:51:23 INFO None 5066136: status FINISHED
2026-07-02 13:51:23 INFO None 5066137: status RUNNING/PENDING
2026-07-02 13:51:23 INFO None 5066139: status FINISHED
2026-07-02 13:51:24 INFO None 5066143: status FINISHED
2026-07-02 13:51:24 INFO None 5066168: status FINISHED
2026-07-02 13:51:24 INFO None 5066179: status RUNNING/PENDING
2026-07-02 13:51:24 INFO None 5066182: status RUNNING/PENDING
2026-07-02 13:51:24 INFO None 5066183: status RUNNING/PENDING
2026-07-02 13:51:24 INFO None 5066184: status RUNNING/PENDING
2026-07-02 13:51:24 INFO None 5066185: status RUNNING/PENDING
2026-07-02 13:51:24 INFO None 5066186: status RUNNING/PENDING
2026-07-02 13:51:24 INFO Jobs still running: ['5066137', '5066179', '5066182', '5066183', '5066184', '5066185', '5066186']. Waiting...
2026-07-02 13:51:40 INFO None 5066131: status FINISHED
2026-07-02 13:51:40 INFO None 5066132: status FINISHED
2026-07-02 13:51:40 INFO None 5066133: status FINISHED
2026-07-02 13:51:40 INFO None 5066134: status FINISHED
2026-07-02 13:51:40 INFO None 5066136: status FINISHED
2026-07-02 13:51:40 INFO None 5066137: status RUNNING/PENDING
2026-07-02 13:51:40 INFO None 5066139: status FINISHED
2026-07-02 13:51:40 INFO None 5066143: status FINISHED
2026-07-02 13:51:40 INFO None 5066168: status FINISHED
2026-07-02 13:51:40 INFO None 5066179: status FINISHED
2026-07-02 13:51:40 INFO None 5066182: status RUNNING/PENDING
2026-07-02 13:51:40 INFO None 5066183: status RUNNING/PENDING
2026-07-02 13:51:40 INFO None 5066184: status RUNNING/PENDING
2026-07-02 13:51:40 INFO None 5066185: status RUNNING/PENDING
2026-07-02 13:51:40 INFO None 5066186: status RUNNING/PENDING
2026-07-02 13:51:40 INFO Jobs still running: ['5066137', '5066182', '5066183', '5066184', '5066185', '5066186']. Waiting...
2026-07-02 13:51:55 INFO None 5066131: status FINISHED
2026-07-02 13:51:55 INFO None 5066132: status FINISHED
2026-07-02 13:51:55 INFO None 5066133: status FINISHED
2026-07-02 13:51:55 INFO None 5066134: status FINISHED
2026-07-02 13:51:55 INFO None 5066136: status FINISHED
2026-07-02 13:51:55 INFO None 5066137: status RUNNING/PENDING
2026-07-02 13:51:55 INFO None 5066139: status FINISHED
2026-07-02 13:51:55 INFO None 5066143: status FINISHED
2026-07-02 13:51:55 INFO None 5066168: status FINISHED
2026-07-02 13:51:55 INFO None 5066179: status FINISHED
2026-07-02 13:51:55 INFO None 5066182: status RUNNING/PENDING
2026-07-02 13:51:55 INFO None 5066183: status FINISHED
2026-07-02 13:51:55 INFO None 5066184: status FINISHED
2026-07-02 13:51:55 INFO None 5066185: status FINISHED
2026-07-02 13:51:55 INFO None 5066186: status RUNNING/PENDING
2026-07-02 13:51:55 INFO Jobs still running: ['5066137', '5066182', '5066186']. Waiting...
2026-07-02 13:52:10 INFO None 5066131: status FINISHED
2026-07-02 13:52:10 INFO None 5066132: status FINISHED
2026-07-02 13:52:10 INFO None 5066133: status FINISHED
2026-07-02 13:52:10 INFO None 5066134: status FINISHED
2026-07-02 13:52:10 INFO None 5066136: status FINISHED
2026-07-02 13:52:10 INFO None 5066137: status FINISHED
2026-07-02 13:52:10 INFO None 5066139: status FINISHED
2026-07-02 13:52:11 INFO None 5066143: status FINISHED
2026-07-02 13:52:11 INFO None 5066168: status FINISHED
2026-07-02 13:52:11 INFO None 5066179: status FINISHED
2026-07-02 13:52:11 INFO None 5066182: status RUNNING/PENDING
2026-07-02 13:52:11 INFO None 5066183: status FINISHED
2026-07-02 13:52:11 INFO None 5066184: status FINISHED
2026-07-02 13:52:11 INFO None 5066185: status FINISHED
2026-07-02 13:52:11 INFO None 5066186: status RUNNING/PENDING
2026-07-02 13:52:11 INFO Jobs still running: ['5066182', '5066186']. Waiting...
2026-07-02 13:52:26 INFO None 5066131: status FINISHED
2026-07-02 13:52:26 INFO None 5066132: status FINISHED
2026-07-02 13:52:26 INFO None 5066133: status FINISHED
2026-07-02 13:52:26 INFO None 5066134: status FINISHED
2026-07-02 13:52:26 INFO None 5066136: status FINISHED
2026-07-02 13:52:26 INFO None 5066137: status FINISHED
2026-07-02 13:52:26 INFO None 5066139: status FINISHED
2026-07-02 13:52:26 INFO None 5066143: status FINISHED
2026-07-02 13:52:26 INFO None 5066168: status FINISHED
2026-07-02 13:52:26 INFO None 5066179: status FINISHED
2026-07-02 13:52:26 INFO None 5066182: status RUNNING/PENDING
2026-07-02 13:52:26 INFO None 5066183: status FINISHED
2026-07-02 13:52:26 INFO None 5066184: status FINISHED
2026-07-02 13:52:26 INFO None 5066185: status FINISHED
2026-07-02 13:52:26 INFO None 5066186: status RUNNING/PENDING
2026-07-02 13:52:26 INFO Jobs still running: ['5066182', '5066186']. Waiting...
2026-07-02 13:52:41 INFO None 5066131: status FINISHED
2026-07-02 13:52:41 INFO None 5066132: status FINISHED
2026-07-02 13:52:41 INFO None 5066133: status FINISHED
2026-07-02 13:52:41 INFO None 5066134: status FINISHED
2026-07-02 13:52:41 INFO None 5066136: status FINISHED
2026-07-02 13:52:41 INFO None 5066137: status FINISHED
2026-07-02 13:52:41 INFO None 5066139: status FINISHED
2026-07-02 13:52:41 INFO None 5066143: status FINISHED
2026-07-02 13:52:41 INFO None 5066168: status FINISHED
2026-07-02 13:52:41 INFO None 5066179: status FINISHED
2026-07-02 13:52:41 INFO None 5066182: status RUNNING/PENDING
2026-07-02 13:52:41 INFO None 5066183: status FINISHED
2026-07-02 13:52:41 INFO None 5066184: status FINISHED
2026-07-02 13:52:41 INFO None 5066185: status FINISHED
2026-07-02 13:52:41 INFO None 5066186: status FINISHED
2026-07-02 13:52:41 INFO Jobs still running: ['5066182']. Waiting...
2026-07-02 13:52:56 INFO None 5066131: status FINISHED
2026-07-02 13:52:56 INFO None 5066132: status FINISHED
2026-07-02 13:52:56 INFO None 5066133: status FINISHED
2026-07-02 13:52:56 INFO None 5066134: status FINISHED
2026-07-02 13:52:56 INFO None 5066136: status FINISHED
2026-07-02 13:52:56 INFO None 5066137: status FINISHED
2026-07-02 13:52:56 INFO None 5066139: status FINISHED
2026-07-02 13:52:56 INFO None 5066143: status FINISHED
2026-07-02 13:52:56 INFO None 5066168: status FINISHED
2026-07-02 13:52:56 INFO None 5066179: status FINISHED
2026-07-02 13:52:56 INFO None 5066182: status RUNNING/PENDING
2026-07-02 13:52:56 INFO None 5066183: status FINISHED
2026-07-02 13:52:56 INFO None 5066184: status FINISHED
2026-07-02 13:52:57 INFO None 5066185: status FINISHED
2026-07-02 13:52:57 INFO None 5066186: status FINISHED
2026-07-02 13:52:57 INFO Jobs still running: ['5066182']. Waiting...
2026-07-02 13:53:12 INFO None 5066131: status FINISHED
2026-07-02 13:53:12 INFO None 5066132: status FINISHED
2026-07-02 13:53:12 INFO None 5066133: status FINISHED
2026-07-02 13:53:12 INFO None 5066134: status FINISHED
2026-07-02 13:53:12 INFO None 5066136: status FINISHED
2026-07-02 13:53:12 INFO None 5066137: status FINISHED
2026-07-02 13:53:12 INFO None 5066139: status FINISHED
2026-07-02 13:53:12 INFO None 5066143: status FINISHED
2026-07-02 13:53:12 INFO None 5066168: status FINISHED
2026-07-02 13:53:12 INFO None 5066179: status FINISHED
2026-07-02 13:53:12 INFO None 5066182: status FINISHED
2026-07-02 13:53:12 INFO None 5066183: status FINISHED
2026-07-02 13:53:12 INFO None 5066184: status FINISHED
2026-07-02 13:53:12 INFO None 5066185: status FINISHED
2026-07-02 13:53:12 INFO None 5066186: status FINISHED
2026-07-02 13:53:12 INFO Jobs ['5066131', '5066132', '5066133', '5066134', '5066136', '5066137', '5066139', '5066143', '5066168', '5066179', '5066182', '5066183', '5066184', '5066185', '5066186'] have finished
2026-07-02 13:53:12 INFO Checking restart files were created ...
2026-07-02 13:53:12 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020700_1_ENS1.nc(668832435 bytes)
2026-07-02 13:53:12 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020700_1_ENS2.nc(668832435 bytes)
2026-07-02 13:53:12 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020700_1_ENS3.nc(668832435 bytes)
2026-07-02 13:53:12 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020700_1_ENS4.nc(668832435 bytes)
2026-07-02 13:53:12 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020700_1_ENS5.nc(668832435 bytes)
2026-07-02 13:53:12 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020700_1_ENS6.nc(668832435 bytes)
2026-07-02 13:53:12 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020700_1_ENS7.nc(668832435 bytes)
2026-07-02 13:53:12 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020700_1_ENS8.nc(668832435 bytes)
2026-07-02 13:53:12 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020700_1_ENS9.nc(668832435 bytes)
2026-07-02 13:53:12 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020700_1_ENS10.nc(668832435 bytes)
2026-07-02 13:53:12 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020700_1_ENS11.nc(668832435 bytes)
2026-07-02 13:53:12 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020700_1_ENS12.nc(668832435 bytes)
2026-07-02 13:53:12 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020700_1_ENS13.nc(668832435 bytes)
2026-07-02 13:53:12 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020700_1_ENS14.nc(668832435 bytes)
2026-07-02 13:53:12 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020700_1_ENS15.nc(668832435 bytes)
2026-07-02 13:53:12 INFO  Run_model() completed successfully.
2026-07-02 13:53:12 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 00:00:00 simulated_time=2020-02-07 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 13:53:12 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 01:00:00 days=153073 seconds=3600
2026-07-02 13:53:12 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-02 13:53:12 INFO [TIME] increment current_time 2020-02-07 00:00:00 -> 2020-02-07 01:00:00
2026-07-02 13:53:12 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 01:00:00 simulated_time=2020-02-07 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 13:53:12 INFO ---------->>> Running process_satellite_data()
2026-07-02 13:53:12 INFO [DART] No satellite data found, skipping assimilation
2026-07-02 13:53:12 INFO after_assimilation() skipped
2026-07-02 13:53:12 INFO Next run starts from 2020-02-07 01:00:00
2026-07-02 13:53:12 INFO Cycle is DONE; starting a new loop!
2026-07-02 13:53:12 INFO [TIME] step_end current_time=2020-02-07 01:00:00 simulated_time=2020-02-07 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 13:53:12 INFO [TIME] step_start current_time=2020-02-07 01:00:00 simulated_time=2020-02-07 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 13:53:12 INFO [TIME] window start=2020-02-07 01:00:00 end=2020-02-07 09:00:00 run_hours=8 has_assimilation=True
2026-07-02 13:53:12 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:53:13 INFO Hourly dataset computed and listing created
2026-07-02 13:53:26 INFO Hourly dataset computed
2026-07-02 13:53:26 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:53:27 INFO Hourly dataset computed and listing created
2026-07-02 13:53:29 INFO Hourly dataset computed
2026-07-02 13:53:29 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:53:30 INFO Hourly dataset computed and listing created
2026-07-02 13:53:32 INFO Hourly dataset computed
2026-07-02 13:53:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:53:33 INFO Hourly dataset computed and listing created
2026-07-02 13:53:35 INFO Hourly dataset computed
2026-07-02 13:53:35 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:53:36 INFO Hourly dataset computed and listing created
2026-07-02 13:53:38 INFO Hourly dataset computed
2026-07-02 13:53:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:53:40 INFO Hourly dataset computed and listing created
2026-07-02 13:53:42 INFO Hourly dataset computed
2026-07-02 13:53:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:53:43 INFO Hourly dataset computed and listing created
2026-07-02 13:53:45 INFO Hourly dataset computed
2026-07-02 13:53:45 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:53:46 INFO Hourly dataset computed and listing created
2026-07-02 13:53:48 INFO Hourly dataset computed
2026-07-02 13:53:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:53:49 INFO Hourly dataset computed and listing created
2026-07-02 13:53:51 INFO Hourly dataset computed
2026-07-02 13:53:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:53:52 INFO Hourly dataset computed and listing created
2026-07-02 13:53:54 INFO Hourly dataset computed
2026-07-02 13:53:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:53:56 INFO Hourly dataset computed and listing created
2026-07-02 13:53:58 INFO Hourly dataset computed
2026-07-02 13:53:58 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:53:59 INFO Hourly dataset computed and listing created
2026-07-02 13:54:01 INFO Hourly dataset computed
2026-07-02 13:54:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:54:02 INFO Hourly dataset computed and listing created
2026-07-02 13:54:04 INFO Hourly dataset computed
2026-07-02 13:54:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:54:05 INFO Hourly dataset computed and listing created
2026-07-02 13:54:08 INFO Hourly dataset computed
2026-07-02 13:54:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 13:54:09 INFO Hourly dataset computed and listing created
2026-07-02 13:54:11 INFO Hourly dataset computed
2026-07-02 13:54:11 INFO ---------->>> Running CHIMERE model from 2020-02-07 01:00:00 to 2020-02-07 09:00:00
2026-07-02 13:54:11 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:54:11 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1
2026-07-02 13:54:11 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020700_1_ENS1.nc
2026-07-02 13:54:11 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-02 13:54:11 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:54:11 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-02 13:54:11 INFO Queuing job for member 1...
2026-07-02 13:54:11 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:54:11 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-02 13:54:12 INFO Found: ['5066204']
2026-07-02 13:54:17 INFO [TGCC-IRENE] Submitted job with ID:['5066204']
2026-07-02 13:54:17 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:54:17 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2
2026-07-02 13:54:17 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020700_1_ENS2.nc
2026-07-02 13:54:17 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-02 13:54:17 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:54:17 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-02 13:54:17 INFO Queuing job for member 2...
2026-07-02 13:54:17 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:54:17 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-02 13:54:18 INFO Found: ['5066205']
2026-07-02 13:54:23 INFO [TGCC-IRENE] Submitted job with ID:['5066205']
2026-07-02 13:54:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:54:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3
2026-07-02 13:54:23 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020700_1_ENS3.nc
2026-07-02 13:54:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-02 13:54:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:54:23 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-02 13:54:23 INFO Queuing job for member 3...
2026-07-02 13:54:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:54:23 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-02 13:54:24 INFO Found: ['5066206']
2026-07-02 13:54:29 INFO [TGCC-IRENE] Submitted job with ID:['5066206']
2026-07-02 13:54:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:54:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4
2026-07-02 13:54:29 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020700_1_ENS4.nc
2026-07-02 13:54:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-02 13:54:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:54:29 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-02 13:54:29 INFO Queuing job for member 4...
2026-07-02 13:54:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:54:29 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-02 13:54:29 INFO Found: ['5066207']
2026-07-02 13:54:34 INFO [TGCC-IRENE] Submitted job with ID:['5066207']
2026-07-02 13:54:34 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:54:34 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5
2026-07-02 13:54:34 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020700_1_ENS5.nc
2026-07-02 13:54:34 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-02 13:54:34 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:54:34 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-02 13:54:34 INFO Queuing job for member 5...
2026-07-02 13:54:34 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:54:34 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-02 13:54:35 INFO Found: ['5066209']
2026-07-02 13:54:40 INFO [TGCC-IRENE] Submitted job with ID:['5066209']
2026-07-02 13:54:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:54:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6
2026-07-02 13:54:40 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020700_1_ENS6.nc
2026-07-02 13:54:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-02 13:54:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:54:40 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-02 13:54:40 INFO Queuing job for member 6...
2026-07-02 13:54:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:54:40 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-02 13:54:41 INFO Found: ['5066212']
2026-07-02 13:54:46 INFO [TGCC-IRENE] Submitted job with ID:['5066212']
2026-07-02 13:54:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:54:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7
2026-07-02 13:54:46 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020700_1_ENS7.nc
2026-07-02 13:54:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-02 13:54:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:54:46 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-02 13:54:46 INFO Queuing job for member 7...
2026-07-02 13:54:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:54:46 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-02 13:54:47 INFO Found: ['5066214']
2026-07-02 13:54:52 INFO [TGCC-IRENE] Submitted job with ID:['5066214']
2026-07-02 13:54:52 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:54:52 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8
2026-07-02 13:54:52 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020700_1_ENS8.nc
2026-07-02 13:54:52 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-02 13:54:52 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:54:52 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-02 13:54:52 INFO Queuing job for member 8...
2026-07-02 13:54:52 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:54:52 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-02 13:54:53 INFO Found: ['5066216']
2026-07-02 13:54:58 INFO [TGCC-IRENE] Submitted job with ID:['5066216']
2026-07-02 13:54:58 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:54:58 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9
2026-07-02 13:54:58 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020700_1_ENS9.nc
2026-07-02 13:54:58 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-02 13:54:58 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:54:58 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-02 13:54:58 INFO Queuing job for member 9...
2026-07-02 13:54:58 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:54:58 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-02 13:55:00 INFO Found: ['5066217']
2026-07-02 13:55:05 INFO [TGCC-IRENE] Submitted job with ID:['5066217']
2026-07-02 13:55:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:55:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10
2026-07-02 13:55:05 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020700_1_ENS10.nc
2026-07-02 13:55:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-02 13:55:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:55:05 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-02 13:55:05 INFO Queuing job for member 10...
2026-07-02 13:55:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:55:05 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-02 13:55:06 INFO Found: ['5066219']
2026-07-02 13:55:11 INFO [TGCC-IRENE] Submitted job with ID:['5066219']
2026-07-02 13:55:11 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:55:11 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11
2026-07-02 13:55:11 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020700_1_ENS11.nc
2026-07-02 13:55:11 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-02 13:55:11 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:55:11 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-02 13:55:11 INFO Queuing job for member 11...
2026-07-02 13:55:11 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:55:11 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-02 13:55:12 INFO Found: ['5066220']
2026-07-02 13:55:17 INFO [TGCC-IRENE] Submitted job with ID:['5066220']
2026-07-02 13:55:17 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:55:17 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12
2026-07-02 13:55:17 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020700_1_ENS12.nc
2026-07-02 13:55:17 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-02 13:55:17 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:55:17 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-02 13:55:17 INFO Queuing job for member 12...
2026-07-02 13:55:17 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:55:17 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-02 13:55:17 INFO Found: ['5066223']
2026-07-02 13:55:22 INFO [TGCC-IRENE] Submitted job with ID:['5066223']
2026-07-02 13:55:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:55:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13
2026-07-02 13:55:22 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020700_1_ENS13.nc
2026-07-02 13:55:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-02 13:55:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:55:23 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-02 13:55:23 INFO Queuing job for member 13...
2026-07-02 13:55:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:55:23 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-02 13:55:23 INFO Found: ['5066242']
2026-07-02 13:55:28 INFO [TGCC-IRENE] Submitted job with ID:['5066242']
2026-07-02 13:55:28 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:55:28 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14
2026-07-02 13:55:28 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020700_1_ENS14.nc
2026-07-02 13:55:28 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-02 13:55:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:55:28 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-02 13:55:28 INFO Queuing job for member 14...
2026-07-02 13:55:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:55:28 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-02 13:55:29 INFO Found: ['5066262']
2026-07-02 13:55:34 INFO [TGCC-IRENE] Submitted job with ID:['5066262']
2026-07-02 13:55:34 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 13:55:34 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15
2026-07-02 13:55:34 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020700_1_ENS15.nc
2026-07-02 13:55:34 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-02 13:55:34 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 13:55:34 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-02 13:55:34 INFO Queuing job for member 15...
2026-07-02 13:55:34 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 13:55:34 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-02 13:55:35 INFO Found: ['5066280']
2026-07-02 13:55:40 INFO [TGCC-IRENE] Submitted job with ID:['5066280']
2026-07-02 13:55:40 INFO Checking job status ...
2026-07-02 13:55:40 INFO None 5066204: status RUNNING/PENDING
2026-07-02 13:55:40 INFO None 5066205: status RUNNING/PENDING
2026-07-02 13:55:40 INFO None 5066206: status RUNNING/PENDING
2026-07-02 13:55:40 INFO None 5066207: status RUNNING/PENDING
2026-07-02 13:55:40 INFO None 5066209: status RUNNING/PENDING
2026-07-02 13:55:40 INFO None 5066212: status RUNNING/PENDING
2026-07-02 13:55:40 INFO None 5066214: status RUNNING/PENDING
2026-07-02 13:55:40 INFO None 5066216: status RUNNING/PENDING
2026-07-02 13:55:40 INFO None 5066217: status RUNNING/PENDING
2026-07-02 13:55:40 INFO None 5066219: status RUNNING/PENDING
2026-07-02 13:55:40 INFO None 5066220: status RUNNING/PENDING
2026-07-02 13:55:40 INFO None 5066223: status RUNNING/PENDING
2026-07-02 13:55:40 INFO None 5066242: status RUNNING/PENDING
2026-07-02 13:55:40 INFO None 5066262: status RUNNING/PENDING
2026-07-02 13:55:40 INFO None 5066280: status RUNNING/PENDING
2026-07-02 13:55:40 INFO Jobs still running: ['5066204', '5066205', '5066206', '5066207', '5066209', '5066212', '5066214', '5066216', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 13:55:55 INFO None 5066204: status RUNNING/PENDING
2026-07-02 13:55:55 INFO None 5066205: status RUNNING/PENDING
2026-07-02 13:55:55 INFO None 5066206: status RUNNING/PENDING
2026-07-02 13:55:55 INFO None 5066207: status RUNNING/PENDING
2026-07-02 13:55:55 INFO None 5066209: status RUNNING/PENDING
2026-07-02 13:55:55 INFO None 5066212: status RUNNING/PENDING
2026-07-02 13:55:55 INFO None 5066214: status RUNNING/PENDING
2026-07-02 13:55:55 INFO None 5066216: status RUNNING/PENDING
2026-07-02 13:55:55 INFO None 5066217: status RUNNING/PENDING
2026-07-02 13:55:55 INFO None 5066219: status RUNNING/PENDING
2026-07-02 13:55:55 INFO None 5066220: status RUNNING/PENDING
2026-07-02 13:55:55 INFO None 5066223: status RUNNING/PENDING
2026-07-02 13:55:55 INFO None 5066242: status RUNNING/PENDING
2026-07-02 13:55:55 INFO None 5066262: status RUNNING/PENDING
2026-07-02 13:55:55 INFO None 5066280: status RUNNING/PENDING
2026-07-02 13:55:55 INFO Jobs still running: ['5066204', '5066205', '5066206', '5066207', '5066209', '5066212', '5066214', '5066216', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 13:56:10 INFO None 5066204: status RUNNING/PENDING
2026-07-02 13:56:10 INFO None 5066205: status RUNNING/PENDING
2026-07-02 13:56:10 INFO None 5066206: status RUNNING/PENDING
2026-07-02 13:56:10 INFO None 5066207: status RUNNING/PENDING
2026-07-02 13:56:10 INFO None 5066209: status RUNNING/PENDING
2026-07-02 13:56:10 INFO None 5066212: status RUNNING/PENDING
2026-07-02 13:56:10 INFO None 5066214: status RUNNING/PENDING
2026-07-02 13:56:10 INFO None 5066216: status RUNNING/PENDING
2026-07-02 13:56:10 INFO None 5066217: status RUNNING/PENDING
2026-07-02 13:56:10 INFO None 5066219: status RUNNING/PENDING
2026-07-02 13:56:11 INFO None 5066220: status RUNNING/PENDING
2026-07-02 13:56:11 INFO None 5066223: status RUNNING/PENDING
2026-07-02 13:56:11 INFO None 5066242: status RUNNING/PENDING
2026-07-02 13:56:11 INFO None 5066262: status RUNNING/PENDING
2026-07-02 13:56:11 INFO None 5066280: status RUNNING/PENDING
2026-07-02 13:56:11 INFO Jobs still running: ['5066204', '5066205', '5066206', '5066207', '5066209', '5066212', '5066214', '5066216', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 13:56:26 INFO None 5066204: status RUNNING/PENDING
2026-07-02 13:56:26 INFO None 5066205: status RUNNING/PENDING
2026-07-02 13:56:26 INFO None 5066206: status RUNNING/PENDING
2026-07-02 13:56:26 INFO None 5066207: status RUNNING/PENDING
2026-07-02 13:56:26 INFO None 5066209: status RUNNING/PENDING
2026-07-02 13:56:26 INFO None 5066212: status RUNNING/PENDING
2026-07-02 13:56:26 INFO None 5066214: status RUNNING/PENDING
2026-07-02 13:56:26 INFO None 5066216: status RUNNING/PENDING
2026-07-02 13:56:26 INFO None 5066217: status RUNNING/PENDING
2026-07-02 13:56:26 INFO None 5066219: status RUNNING/PENDING
2026-07-02 13:56:26 INFO None 5066220: status RUNNING/PENDING
2026-07-02 13:56:26 INFO None 5066223: status RUNNING/PENDING
2026-07-02 13:56:26 INFO None 5066242: status RUNNING/PENDING
2026-07-02 13:56:26 INFO None 5066262: status RUNNING/PENDING
2026-07-02 13:56:26 INFO None 5066280: status RUNNING/PENDING
2026-07-02 13:56:26 INFO Jobs still running: ['5066204', '5066205', '5066206', '5066207', '5066209', '5066212', '5066214', '5066216', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 13:56:42 INFO None 5066204: status RUNNING/PENDING
2026-07-02 13:56:42 INFO None 5066205: status RUNNING/PENDING
2026-07-02 13:56:42 INFO None 5066206: status RUNNING/PENDING
2026-07-02 13:56:42 INFO None 5066207: status RUNNING/PENDING
2026-07-02 13:56:42 INFO None 5066209: status RUNNING/PENDING
2026-07-02 13:56:42 INFO None 5066212: status RUNNING/PENDING
2026-07-02 13:56:42 INFO None 5066214: status RUNNING/PENDING
2026-07-02 13:56:42 INFO None 5066216: status RUNNING/PENDING
2026-07-02 13:56:42 INFO None 5066217: status RUNNING/PENDING
2026-07-02 13:56:42 INFO None 5066219: status RUNNING/PENDING
2026-07-02 13:56:42 INFO None 5066220: status RUNNING/PENDING
2026-07-02 13:56:42 INFO None 5066223: status RUNNING/PENDING
2026-07-02 13:56:42 INFO None 5066242: status RUNNING/PENDING
2026-07-02 13:56:42 INFO None 5066262: status RUNNING/PENDING
2026-07-02 13:56:42 INFO None 5066280: status RUNNING/PENDING
2026-07-02 13:56:42 INFO Jobs still running: ['5066204', '5066205', '5066206', '5066207', '5066209', '5066212', '5066214', '5066216', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 13:56:57 INFO None 5066204: status RUNNING/PENDING
2026-07-02 13:56:57 INFO None 5066205: status RUNNING/PENDING
2026-07-02 13:56:57 INFO None 5066206: status RUNNING/PENDING
2026-07-02 13:56:57 INFO None 5066207: status RUNNING/PENDING
2026-07-02 13:56:57 INFO None 5066209: status RUNNING/PENDING
2026-07-02 13:56:57 INFO None 5066212: status RUNNING/PENDING
2026-07-02 13:56:57 INFO None 5066214: status RUNNING/PENDING
2026-07-02 13:56:57 INFO None 5066216: status RUNNING/PENDING
2026-07-02 13:56:57 INFO None 5066217: status RUNNING/PENDING
2026-07-02 13:56:57 INFO None 5066219: status RUNNING/PENDING
2026-07-02 13:56:57 INFO None 5066220: status RUNNING/PENDING
2026-07-02 13:56:57 INFO None 5066223: status RUNNING/PENDING
2026-07-02 13:56:57 INFO None 5066242: status RUNNING/PENDING
2026-07-02 13:56:57 INFO None 5066262: status RUNNING/PENDING
2026-07-02 13:56:57 INFO None 5066280: status RUNNING/PENDING
2026-07-02 13:56:57 INFO Jobs still running: ['5066204', '5066205', '5066206', '5066207', '5066209', '5066212', '5066214', '5066216', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 13:57:12 INFO None 5066204: status RUNNING/PENDING
2026-07-02 13:57:12 INFO None 5066205: status RUNNING/PENDING
2026-07-02 13:57:12 INFO None 5066206: status RUNNING/PENDING
2026-07-02 13:57:12 INFO None 5066207: status RUNNING/PENDING
2026-07-02 13:57:12 INFO None 5066209: status RUNNING/PENDING
2026-07-02 13:57:12 INFO None 5066212: status RUNNING/PENDING
2026-07-02 13:57:12 INFO None 5066214: status RUNNING/PENDING
2026-07-02 13:57:12 INFO None 5066216: status RUNNING/PENDING
2026-07-02 13:57:12 INFO None 5066217: status RUNNING/PENDING
2026-07-02 13:57:12 INFO None 5066219: status RUNNING/PENDING
2026-07-02 13:57:12 INFO None 5066220: status RUNNING/PENDING
2026-07-02 13:57:12 INFO None 5066223: status RUNNING/PENDING
2026-07-02 13:57:13 INFO None 5066242: status RUNNING/PENDING
2026-07-02 13:57:13 INFO None 5066262: status RUNNING/PENDING
2026-07-02 13:57:13 INFO None 5066280: status RUNNING/PENDING
2026-07-02 13:57:13 INFO Jobs still running: ['5066204', '5066205', '5066206', '5066207', '5066209', '5066212', '5066214', '5066216', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 13:57:28 INFO None 5066204: status RUNNING/PENDING
2026-07-02 13:57:28 INFO None 5066205: status RUNNING/PENDING
2026-07-02 13:57:28 INFO None 5066206: status RUNNING/PENDING
2026-07-02 13:57:28 INFO None 5066207: status RUNNING/PENDING
2026-07-02 13:57:28 INFO None 5066209: status RUNNING/PENDING
2026-07-02 13:57:28 INFO None 5066212: status RUNNING/PENDING
2026-07-02 13:57:28 INFO None 5066214: status RUNNING/PENDING
2026-07-02 13:57:28 INFO None 5066216: status RUNNING/PENDING
2026-07-02 13:57:28 INFO None 5066217: status RUNNING/PENDING
2026-07-02 13:57:28 INFO None 5066219: status RUNNING/PENDING
2026-07-02 13:57:28 INFO None 5066220: status RUNNING/PENDING
2026-07-02 13:57:28 INFO None 5066223: status RUNNING/PENDING
2026-07-02 13:57:28 INFO None 5066242: status RUNNING/PENDING
2026-07-02 13:57:28 INFO None 5066262: status RUNNING/PENDING
2026-07-02 13:57:28 INFO None 5066280: status RUNNING/PENDING
2026-07-02 13:57:28 INFO Jobs still running: ['5066204', '5066205', '5066206', '5066207', '5066209', '5066212', '5066214', '5066216', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 13:57:43 INFO None 5066204: status RUNNING/PENDING
2026-07-02 13:57:43 INFO None 5066205: status RUNNING/PENDING
2026-07-02 13:57:43 INFO None 5066206: status RUNNING/PENDING
2026-07-02 13:57:43 INFO None 5066207: status RUNNING/PENDING
2026-07-02 13:57:43 INFO None 5066209: status RUNNING/PENDING
2026-07-02 13:57:43 INFO None 5066212: status RUNNING/PENDING
2026-07-02 13:57:43 INFO None 5066214: status RUNNING/PENDING
2026-07-02 13:57:43 INFO None 5066216: status RUNNING/PENDING
2026-07-02 13:57:43 INFO None 5066217: status RUNNING/PENDING
2026-07-02 13:57:43 INFO None 5066219: status RUNNING/PENDING
2026-07-02 13:57:43 INFO None 5066220: status RUNNING/PENDING
2026-07-02 13:57:43 INFO None 5066223: status RUNNING/PENDING
2026-07-02 13:57:43 INFO None 5066242: status RUNNING/PENDING
2026-07-02 13:57:43 INFO None 5066262: status RUNNING/PENDING
2026-07-02 13:57:43 INFO None 5066280: status RUNNING/PENDING
2026-07-02 13:57:43 INFO Jobs still running: ['5066204', '5066205', '5066206', '5066207', '5066209', '5066212', '5066214', '5066216', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 13:57:58 INFO None 5066204: status RUNNING/PENDING
2026-07-02 13:57:58 INFO None 5066205: status RUNNING/PENDING
2026-07-02 13:57:58 INFO None 5066206: status RUNNING/PENDING
2026-07-02 13:57:58 INFO None 5066207: status RUNNING/PENDING
2026-07-02 13:57:58 INFO None 5066209: status RUNNING/PENDING
2026-07-02 13:57:58 INFO None 5066212: status RUNNING/PENDING
2026-07-02 13:57:58 INFO None 5066214: status RUNNING/PENDING
2026-07-02 13:57:58 INFO None 5066216: status RUNNING/PENDING
2026-07-02 13:57:58 INFO None 5066217: status RUNNING/PENDING
2026-07-02 13:57:58 INFO None 5066219: status RUNNING/PENDING
2026-07-02 13:57:58 INFO None 5066220: status RUNNING/PENDING
2026-07-02 13:57:58 INFO None 5066223: status RUNNING/PENDING
2026-07-02 13:57:58 INFO None 5066242: status RUNNING/PENDING
2026-07-02 13:57:58 INFO None 5066262: status RUNNING/PENDING
2026-07-02 13:57:58 INFO None 5066280: status RUNNING/PENDING
2026-07-02 13:57:58 INFO Jobs still running: ['5066204', '5066205', '5066206', '5066207', '5066209', '5066212', '5066214', '5066216', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 13:58:13 INFO None 5066204: status RUNNING/PENDING
2026-07-02 13:58:14 INFO None 5066205: status RUNNING/PENDING
2026-07-02 13:58:14 INFO None 5066206: status RUNNING/PENDING
2026-07-02 13:58:14 INFO None 5066207: status RUNNING/PENDING
2026-07-02 13:58:14 INFO None 5066209: status RUNNING/PENDING
2026-07-02 13:58:14 INFO None 5066212: status RUNNING/PENDING
2026-07-02 13:58:14 INFO None 5066214: status RUNNING/PENDING
2026-07-02 13:58:14 INFO None 5066216: status RUNNING/PENDING
2026-07-02 13:58:14 INFO None 5066217: status RUNNING/PENDING
2026-07-02 13:58:14 INFO None 5066219: status RUNNING/PENDING
2026-07-02 13:58:14 INFO None 5066220: status RUNNING/PENDING
2026-07-02 13:58:14 INFO None 5066223: status RUNNING/PENDING
2026-07-02 13:58:14 INFO None 5066242: status RUNNING/PENDING
2026-07-02 13:58:14 INFO None 5066262: status RUNNING/PENDING
2026-07-02 13:58:14 INFO None 5066280: status RUNNING/PENDING
2026-07-02 13:58:14 INFO Jobs still running: ['5066204', '5066205', '5066206', '5066207', '5066209', '5066212', '5066214', '5066216', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 13:58:29 INFO None 5066204: status RUNNING/PENDING
2026-07-02 13:58:29 INFO None 5066205: status RUNNING/PENDING
2026-07-02 13:58:29 INFO None 5066206: status RUNNING/PENDING
2026-07-02 13:58:29 INFO None 5066207: status RUNNING/PENDING
2026-07-02 13:58:29 INFO None 5066209: status RUNNING/PENDING
2026-07-02 13:58:29 INFO None 5066212: status RUNNING/PENDING
2026-07-02 13:58:29 INFO None 5066214: status RUNNING/PENDING
2026-07-02 13:58:29 INFO None 5066216: status RUNNING/PENDING
2026-07-02 13:58:29 INFO None 5066217: status RUNNING/PENDING
2026-07-02 13:58:29 INFO None 5066219: status RUNNING/PENDING
2026-07-02 13:58:29 INFO None 5066220: status RUNNING/PENDING
2026-07-02 13:58:29 INFO None 5066223: status RUNNING/PENDING
2026-07-02 13:58:29 INFO None 5066242: status RUNNING/PENDING
2026-07-02 13:58:29 INFO None 5066262: status RUNNING/PENDING
2026-07-02 13:58:29 INFO None 5066280: status RUNNING/PENDING
2026-07-02 13:58:29 INFO Jobs still running: ['5066204', '5066205', '5066206', '5066207', '5066209', '5066212', '5066214', '5066216', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 13:58:44 INFO None 5066204: status RUNNING/PENDING
2026-07-02 13:58:44 INFO None 5066205: status RUNNING/PENDING
2026-07-02 13:58:44 INFO None 5066206: status RUNNING/PENDING
2026-07-02 13:58:44 INFO None 5066207: status RUNNING/PENDING
2026-07-02 13:58:44 INFO None 5066209: status RUNNING/PENDING
2026-07-02 13:58:44 INFO None 5066212: status RUNNING/PENDING
2026-07-02 13:58:44 INFO None 5066214: status RUNNING/PENDING
2026-07-02 13:58:44 INFO None 5066216: status RUNNING/PENDING
2026-07-02 13:58:44 INFO None 5066217: status RUNNING/PENDING
2026-07-02 13:58:44 INFO None 5066219: status RUNNING/PENDING
2026-07-02 13:58:44 INFO None 5066220: status RUNNING/PENDING
2026-07-02 13:58:44 INFO None 5066223: status RUNNING/PENDING
2026-07-02 13:58:44 INFO None 5066242: status RUNNING/PENDING
2026-07-02 13:58:44 INFO None 5066262: status RUNNING/PENDING
2026-07-02 13:58:44 INFO None 5066280: status RUNNING/PENDING
2026-07-02 13:58:44 INFO Jobs still running: ['5066204', '5066205', '5066206', '5066207', '5066209', '5066212', '5066214', '5066216', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 13:58:59 INFO None 5066204: status RUNNING/PENDING
2026-07-02 13:58:59 INFO None 5066205: status RUNNING/PENDING
2026-07-02 13:58:59 INFO None 5066206: status RUNNING/PENDING
2026-07-02 13:59:00 INFO None 5066207: status RUNNING/PENDING
2026-07-02 13:59:00 INFO None 5066209: status RUNNING/PENDING
2026-07-02 13:59:00 INFO None 5066212: status RUNNING/PENDING
2026-07-02 13:59:00 INFO None 5066214: status RUNNING/PENDING
2026-07-02 13:59:00 INFO None 5066216: status RUNNING/PENDING
2026-07-02 13:59:00 INFO None 5066217: status RUNNING/PENDING
2026-07-02 13:59:00 INFO None 5066219: status RUNNING/PENDING
2026-07-02 13:59:00 INFO None 5066220: status RUNNING/PENDING
2026-07-02 13:59:00 INFO None 5066223: status RUNNING/PENDING
2026-07-02 13:59:00 INFO None 5066242: status RUNNING/PENDING
2026-07-02 13:59:00 INFO None 5066262: status RUNNING/PENDING
2026-07-02 13:59:00 INFO None 5066280: status RUNNING/PENDING
2026-07-02 13:59:00 INFO Jobs still running: ['5066204', '5066205', '5066206', '5066207', '5066209', '5066212', '5066214', '5066216', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 13:59:15 INFO None 5066204: status RUNNING/PENDING
2026-07-02 13:59:15 INFO None 5066205: status RUNNING/PENDING
2026-07-02 13:59:15 INFO None 5066206: status RUNNING/PENDING
2026-07-02 13:59:15 INFO None 5066207: status RUNNING/PENDING
2026-07-02 13:59:15 INFO None 5066209: status RUNNING/PENDING
2026-07-02 13:59:15 INFO None 5066212: status RUNNING/PENDING
2026-07-02 13:59:15 INFO None 5066214: status RUNNING/PENDING
2026-07-02 13:59:15 INFO None 5066216: status RUNNING/PENDING
2026-07-02 13:59:15 INFO None 5066217: status RUNNING/PENDING
2026-07-02 13:59:15 INFO None 5066219: status RUNNING/PENDING
2026-07-02 13:59:15 INFO None 5066220: status RUNNING/PENDING
2026-07-02 13:59:15 INFO None 5066223: status RUNNING/PENDING
2026-07-02 13:59:15 INFO None 5066242: status RUNNING/PENDING
2026-07-02 13:59:15 INFO None 5066262: status RUNNING/PENDING
2026-07-02 13:59:15 INFO None 5066280: status RUNNING/PENDING
2026-07-02 13:59:15 INFO Jobs still running: ['5066204', '5066205', '5066206', '5066207', '5066209', '5066212', '5066214', '5066216', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 13:59:30 INFO None 5066204: status RUNNING/PENDING
2026-07-02 13:59:30 INFO None 5066205: status RUNNING/PENDING
2026-07-02 13:59:30 INFO None 5066206: status RUNNING/PENDING
2026-07-02 13:59:30 INFO None 5066207: status RUNNING/PENDING
2026-07-02 13:59:30 INFO None 5066209: status RUNNING/PENDING
2026-07-02 13:59:30 INFO None 5066212: status RUNNING/PENDING
2026-07-02 13:59:30 INFO None 5066214: status RUNNING/PENDING
2026-07-02 13:59:30 INFO None 5066216: status RUNNING/PENDING
2026-07-02 13:59:30 INFO None 5066217: status RUNNING/PENDING
2026-07-02 13:59:30 INFO None 5066219: status RUNNING/PENDING
2026-07-02 13:59:30 INFO None 5066220: status RUNNING/PENDING
2026-07-02 13:59:30 INFO None 5066223: status RUNNING/PENDING
2026-07-02 13:59:30 INFO None 5066242: status RUNNING/PENDING
2026-07-02 13:59:30 INFO None 5066262: status RUNNING/PENDING
2026-07-02 13:59:30 INFO None 5066280: status RUNNING/PENDING
2026-07-02 13:59:30 INFO Jobs still running: ['5066204', '5066205', '5066206', '5066207', '5066209', '5066212', '5066214', '5066216', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 13:59:45 INFO None 5066204: status RUNNING/PENDING
2026-07-02 13:59:45 INFO None 5066205: status RUNNING/PENDING
2026-07-02 13:59:45 INFO None 5066206: status RUNNING/PENDING
2026-07-02 13:59:45 INFO None 5066207: status RUNNING/PENDING
2026-07-02 13:59:45 INFO None 5066209: status RUNNING/PENDING
2026-07-02 13:59:45 INFO None 5066212: status RUNNING/PENDING
2026-07-02 13:59:45 INFO None 5066214: status RUNNING/PENDING
2026-07-02 13:59:45 INFO None 5066216: status RUNNING/PENDING
2026-07-02 13:59:45 INFO None 5066217: status RUNNING/PENDING
2026-07-02 13:59:45 INFO None 5066219: status RUNNING/PENDING
2026-07-02 13:59:46 INFO None 5066220: status RUNNING/PENDING
2026-07-02 13:59:46 INFO None 5066223: status RUNNING/PENDING
2026-07-02 13:59:46 INFO None 5066242: status RUNNING/PENDING
2026-07-02 13:59:46 INFO None 5066262: status RUNNING/PENDING
2026-07-02 13:59:46 INFO None 5066280: status RUNNING/PENDING
2026-07-02 13:59:46 INFO Jobs still running: ['5066204', '5066205', '5066206', '5066207', '5066209', '5066212', '5066214', '5066216', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 14:00:01 INFO None 5066204: status RUNNING/PENDING
2026-07-02 14:00:01 INFO None 5066205: status RUNNING/PENDING
2026-07-02 14:00:01 INFO None 5066206: status RUNNING/PENDING
2026-07-02 14:00:01 INFO None 5066207: status RUNNING/PENDING
2026-07-02 14:00:01 INFO None 5066209: status RUNNING/PENDING
2026-07-02 14:00:01 INFO None 5066212: status RUNNING/PENDING
2026-07-02 14:00:01 INFO None 5066214: status RUNNING/PENDING
2026-07-02 14:00:01 INFO None 5066216: status RUNNING/PENDING
2026-07-02 14:00:01 INFO None 5066217: status RUNNING/PENDING
2026-07-02 14:00:01 INFO None 5066219: status RUNNING/PENDING
2026-07-02 14:00:01 INFO None 5066220: status RUNNING/PENDING
2026-07-02 14:00:01 INFO None 5066223: status RUNNING/PENDING
2026-07-02 14:00:01 INFO None 5066242: status RUNNING/PENDING
2026-07-02 14:00:01 INFO None 5066262: status RUNNING/PENDING
2026-07-02 14:00:01 INFO None 5066280: status RUNNING/PENDING
2026-07-02 14:00:01 INFO Jobs still running: ['5066204', '5066205', '5066206', '5066207', '5066209', '5066212', '5066214', '5066216', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 14:00:16 INFO None 5066204: status RUNNING/PENDING
2026-07-02 14:00:16 INFO None 5066205: status RUNNING/PENDING
2026-07-02 14:00:16 INFO None 5066206: status RUNNING/PENDING
2026-07-02 14:00:16 INFO None 5066207: status RUNNING/PENDING
2026-07-02 14:00:16 INFO None 5066209: status RUNNING/PENDING
2026-07-02 14:00:16 INFO None 5066212: status RUNNING/PENDING
2026-07-02 14:00:16 INFO None 5066214: status RUNNING/PENDING
2026-07-02 14:00:16 INFO None 5066216: status RUNNING/PENDING
2026-07-02 14:00:16 INFO None 5066217: status RUNNING/PENDING
2026-07-02 14:00:16 INFO None 5066219: status RUNNING/PENDING
2026-07-02 14:00:16 INFO None 5066220: status RUNNING/PENDING
2026-07-02 14:00:16 INFO None 5066223: status RUNNING/PENDING
2026-07-02 14:00:16 INFO None 5066242: status RUNNING/PENDING
2026-07-02 14:00:16 INFO None 5066262: status RUNNING/PENDING
2026-07-02 14:00:16 INFO None 5066280: status RUNNING/PENDING
2026-07-02 14:00:16 INFO Jobs still running: ['5066204', '5066205', '5066206', '5066207', '5066209', '5066212', '5066214', '5066216', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 14:00:31 INFO None 5066204: status RUNNING/PENDING
2026-07-02 14:00:31 INFO None 5066205: status RUNNING/PENDING
2026-07-02 14:00:31 INFO None 5066206: status RUNNING/PENDING
2026-07-02 14:00:31 INFO None 5066207: status RUNNING/PENDING
2026-07-02 14:00:32 INFO None 5066209: status RUNNING/PENDING
2026-07-02 14:00:32 INFO None 5066212: status RUNNING/PENDING
2026-07-02 14:00:32 INFO None 5066214: status RUNNING/PENDING
2026-07-02 14:00:32 INFO None 5066216: status RUNNING/PENDING
2026-07-02 14:00:32 INFO None 5066217: status RUNNING/PENDING
2026-07-02 14:00:32 INFO None 5066219: status RUNNING/PENDING
2026-07-02 14:00:32 INFO None 5066220: status RUNNING/PENDING
2026-07-02 14:00:32 INFO None 5066223: status RUNNING/PENDING
2026-07-02 14:00:32 INFO None 5066242: status RUNNING/PENDING
2026-07-02 14:00:32 INFO None 5066262: status RUNNING/PENDING
2026-07-02 14:00:32 INFO None 5066280: status RUNNING/PENDING
2026-07-02 14:00:32 INFO Jobs still running: ['5066204', '5066205', '5066206', '5066207', '5066209', '5066212', '5066214', '5066216', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 14:00:47 INFO None 5066204: status RUNNING/PENDING
2026-07-02 14:00:47 INFO None 5066205: status RUNNING/PENDING
2026-07-02 14:00:47 INFO None 5066206: status RUNNING/PENDING
2026-07-02 14:00:47 INFO None 5066207: status RUNNING/PENDING
2026-07-02 14:00:47 INFO None 5066209: status RUNNING/PENDING
2026-07-02 14:00:47 INFO None 5066212: status RUNNING/PENDING
2026-07-02 14:00:47 INFO None 5066214: status RUNNING/PENDING
2026-07-02 14:00:47 INFO None 5066216: status RUNNING/PENDING
2026-07-02 14:00:47 INFO None 5066217: status RUNNING/PENDING
2026-07-02 14:00:47 INFO None 5066219: status RUNNING/PENDING
2026-07-02 14:00:47 INFO None 5066220: status RUNNING/PENDING
2026-07-02 14:00:47 INFO None 5066223: status RUNNING/PENDING
2026-07-02 14:00:47 INFO None 5066242: status RUNNING/PENDING
2026-07-02 14:00:47 INFO None 5066262: status RUNNING/PENDING
2026-07-02 14:00:47 INFO None 5066280: status RUNNING/PENDING
2026-07-02 14:00:47 INFO Jobs still running: ['5066204', '5066205', '5066206', '5066207', '5066209', '5066212', '5066214', '5066216', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 14:01:03 INFO None 5066204: status RUNNING/PENDING
2026-07-02 14:01:03 INFO None 5066205: status RUNNING/PENDING
2026-07-02 14:01:03 INFO None 5066206: status RUNNING/PENDING
2026-07-02 14:01:03 INFO None 5066207: status RUNNING/PENDING
2026-07-02 14:01:03 INFO None 5066209: status RUNNING/PENDING
2026-07-02 14:01:03 INFO None 5066212: status RUNNING/PENDING
2026-07-02 14:01:03 INFO None 5066214: status RUNNING/PENDING
2026-07-02 14:01:03 INFO None 5066216: status RUNNING/PENDING
2026-07-02 14:01:03 INFO None 5066217: status RUNNING/PENDING
2026-07-02 14:01:03 INFO None 5066219: status RUNNING/PENDING
2026-07-02 14:01:03 INFO None 5066220: status RUNNING/PENDING
2026-07-02 14:01:03 INFO None 5066223: status RUNNING/PENDING
2026-07-02 14:01:03 INFO None 5066242: status RUNNING/PENDING
2026-07-02 14:01:03 INFO None 5066262: status RUNNING/PENDING
2026-07-02 14:01:03 INFO None 5066280: status RUNNING/PENDING
2026-07-02 14:01:03 INFO Jobs still running: ['5066204', '5066205', '5066206', '5066207', '5066209', '5066212', '5066214', '5066216', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 14:01:19 INFO None 5066204: status RUNNING/PENDING
2026-07-02 14:01:19 INFO None 5066205: status RUNNING/PENDING
2026-07-02 14:01:19 INFO None 5066206: status RUNNING/PENDING
2026-07-02 14:01:19 INFO None 5066207: status RUNNING/PENDING
2026-07-02 14:01:19 INFO None 5066209: status RUNNING/PENDING
2026-07-02 14:01:19 INFO None 5066212: status RUNNING/PENDING
2026-07-02 14:01:19 INFO None 5066214: status RUNNING/PENDING
2026-07-02 14:01:19 INFO None 5066216: status RUNNING/PENDING
2026-07-02 14:01:19 INFO None 5066217: status RUNNING/PENDING
2026-07-02 14:01:19 INFO None 5066219: status RUNNING/PENDING
2026-07-02 14:01:19 INFO None 5066220: status RUNNING/PENDING
2026-07-02 14:01:19 INFO None 5066223: status RUNNING/PENDING
2026-07-02 14:01:19 INFO None 5066242: status RUNNING/PENDING
2026-07-02 14:01:19 INFO None 5066262: status RUNNING/PENDING
2026-07-02 14:01:19 INFO None 5066280: status RUNNING/PENDING
2026-07-02 14:01:19 INFO Jobs still running: ['5066204', '5066205', '5066206', '5066207', '5066209', '5066212', '5066214', '5066216', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 14:01:34 INFO None 5066204: status RUNNING/PENDING
2026-07-02 14:01:34 INFO None 5066205: status RUNNING/PENDING
2026-07-02 14:01:34 INFO None 5066206: status RUNNING/PENDING
2026-07-02 14:01:34 INFO None 5066207: status RUNNING/PENDING
2026-07-02 14:01:34 INFO None 5066209: status RUNNING/PENDING
2026-07-02 14:01:34 INFO None 5066212: status RUNNING/PENDING
2026-07-02 14:01:34 INFO None 5066214: status RUNNING/PENDING
2026-07-02 14:01:34 INFO None 5066216: status RUNNING/PENDING
2026-07-02 14:01:34 INFO None 5066217: status RUNNING/PENDING
2026-07-02 14:01:34 INFO None 5066219: status RUNNING/PENDING
2026-07-02 14:01:34 INFO None 5066220: status RUNNING/PENDING
2026-07-02 14:01:34 INFO None 5066223: status RUNNING/PENDING
2026-07-02 14:01:34 INFO None 5066242: status RUNNING/PENDING
2026-07-02 14:01:34 INFO None 5066262: status RUNNING/PENDING
2026-07-02 14:01:34 INFO None 5066280: status RUNNING/PENDING
2026-07-02 14:01:34 INFO Jobs still running: ['5066204', '5066205', '5066206', '5066207', '5066209', '5066212', '5066214', '5066216', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 14:01:49 INFO None 5066204: status RUNNING/PENDING
2026-07-02 14:01:49 INFO None 5066205: status RUNNING/PENDING
2026-07-02 14:01:49 INFO None 5066206: status RUNNING/PENDING
2026-07-02 14:01:49 INFO None 5066207: status RUNNING/PENDING
2026-07-02 14:01:49 INFO None 5066209: status RUNNING/PENDING
2026-07-02 14:01:49 INFO None 5066212: status RUNNING/PENDING
2026-07-02 14:01:49 INFO None 5066214: status RUNNING/PENDING
2026-07-02 14:01:49 INFO None 5066216: status RUNNING/PENDING
2026-07-02 14:01:49 INFO None 5066217: status RUNNING/PENDING
2026-07-02 14:01:49 INFO None 5066219: status RUNNING/PENDING
2026-07-02 14:01:49 INFO None 5066220: status RUNNING/PENDING
2026-07-02 14:01:49 INFO None 5066223: status RUNNING/PENDING
2026-07-02 14:01:49 INFO None 5066242: status RUNNING/PENDING
2026-07-02 14:01:49 INFO None 5066262: status RUNNING/PENDING
2026-07-02 14:01:49 INFO None 5066280: status RUNNING/PENDING
2026-07-02 14:01:49 INFO Jobs still running: ['5066204', '5066205', '5066206', '5066207', '5066209', '5066212', '5066214', '5066216', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 14:02:04 INFO None 5066204: status RUNNING/PENDING
2026-07-02 14:02:04 INFO None 5066205: status RUNNING/PENDING
2026-07-02 14:02:04 INFO None 5066206: status RUNNING/PENDING
2026-07-02 14:02:04 INFO None 5066207: status RUNNING/PENDING
2026-07-02 14:02:04 INFO None 5066209: status RUNNING/PENDING
2026-07-02 14:02:04 INFO None 5066212: status RUNNING/PENDING
2026-07-02 14:02:04 INFO None 5066214: status RUNNING/PENDING
2026-07-02 14:02:05 INFO None 5066216: status RUNNING/PENDING
2026-07-02 14:02:05 INFO None 5066217: status RUNNING/PENDING
2026-07-02 14:02:05 INFO None 5066219: status RUNNING/PENDING
2026-07-02 14:02:05 INFO None 5066220: status RUNNING/PENDING
2026-07-02 14:02:05 INFO None 5066223: status RUNNING/PENDING
2026-07-02 14:02:05 INFO None 5066242: status RUNNING/PENDING
2026-07-02 14:02:05 INFO None 5066262: status RUNNING/PENDING
2026-07-02 14:02:05 INFO None 5066280: status RUNNING/PENDING
2026-07-02 14:02:05 INFO Jobs still running: ['5066204', '5066205', '5066206', '5066207', '5066209', '5066212', '5066214', '5066216', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 14:02:20 INFO None 5066204: status RUNNING/PENDING
2026-07-02 14:02:20 INFO None 5066205: status RUNNING/PENDING
2026-07-02 14:02:20 INFO None 5066206: status RUNNING/PENDING
2026-07-02 14:02:20 INFO None 5066207: status RUNNING/PENDING
2026-07-02 14:02:20 INFO None 5066209: status RUNNING/PENDING
2026-07-02 14:02:20 INFO None 5066212: status RUNNING/PENDING
2026-07-02 14:02:20 INFO None 5066214: status RUNNING/PENDING
2026-07-02 14:02:20 INFO None 5066216: status RUNNING/PENDING
2026-07-02 14:02:20 INFO None 5066217: status RUNNING/PENDING
2026-07-02 14:02:20 INFO None 5066219: status RUNNING/PENDING
2026-07-02 14:02:20 INFO None 5066220: status RUNNING/PENDING
2026-07-02 14:02:20 INFO None 5066223: status RUNNING/PENDING
2026-07-02 14:02:20 INFO None 5066242: status RUNNING/PENDING
2026-07-02 14:02:20 INFO None 5066262: status RUNNING/PENDING
2026-07-02 14:02:20 INFO None 5066280: status RUNNING/PENDING
2026-07-02 14:02:20 INFO Jobs still running: ['5066204', '5066205', '5066206', '5066207', '5066209', '5066212', '5066214', '5066216', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 14:02:35 INFO None 5066204: status RUNNING/PENDING
2026-07-02 14:02:35 INFO None 5066205: status RUNNING/PENDING
2026-07-02 14:02:35 INFO None 5066206: status RUNNING/PENDING
2026-07-02 14:02:35 INFO None 5066207: status RUNNING/PENDING
2026-07-02 14:02:35 INFO None 5066209: status RUNNING/PENDING
2026-07-02 14:02:35 INFO None 5066212: status RUNNING/PENDING
2026-07-02 14:02:35 INFO None 5066214: status FINISHED
2026-07-02 14:02:35 INFO None 5066216: status RUNNING/PENDING
2026-07-02 14:02:35 INFO None 5066217: status RUNNING/PENDING
2026-07-02 14:02:35 INFO None 5066219: status RUNNING/PENDING
2026-07-02 14:02:35 INFO None 5066220: status RUNNING/PENDING
2026-07-02 14:02:35 INFO None 5066223: status RUNNING/PENDING
2026-07-02 14:02:35 INFO None 5066242: status RUNNING/PENDING
2026-07-02 14:02:35 INFO None 5066262: status RUNNING/PENDING
2026-07-02 14:02:35 INFO None 5066280: status RUNNING/PENDING
2026-07-02 14:02:35 INFO Jobs still running: ['5066204', '5066205', '5066206', '5066207', '5066209', '5066212', '5066216', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 14:02:52 INFO None 5066204: status FINISHED
2026-07-02 14:02:52 INFO None 5066205: status RUNNING/PENDING
2026-07-02 14:02:52 INFO None 5066206: status RUNNING/PENDING
2026-07-02 14:02:52 INFO None 5066207: status FINISHED
2026-07-02 14:02:52 INFO None 5066209: status RUNNING/PENDING
2026-07-02 14:02:52 INFO None 5066212: status RUNNING/PENDING
2026-07-02 14:02:52 INFO None 5066214: status FINISHED
2026-07-02 14:02:52 INFO None 5066216: status FINISHED
2026-07-02 14:02:52 INFO None 5066217: status RUNNING/PENDING
2026-07-02 14:02:52 INFO None 5066219: status RUNNING/PENDING
2026-07-02 14:02:52 INFO None 5066220: status RUNNING/PENDING
2026-07-02 14:02:52 INFO None 5066223: status RUNNING/PENDING
2026-07-02 14:02:52 INFO None 5066242: status RUNNING/PENDING
2026-07-02 14:02:52 INFO None 5066262: status RUNNING/PENDING
2026-07-02 14:02:52 INFO None 5066280: status RUNNING/PENDING
2026-07-02 14:02:52 INFO Jobs still running: ['5066205', '5066206', '5066209', '5066212', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 14:03:07 INFO None 5066204: status FINISHED
2026-07-02 14:03:08 INFO None 5066205: status RUNNING/PENDING
2026-07-02 14:03:08 INFO None 5066206: status FINISHED
2026-07-02 14:03:08 INFO None 5066207: status FINISHED
2026-07-02 14:03:08 INFO None 5066209: status FINISHED
2026-07-02 14:03:08 INFO None 5066212: status RUNNING/PENDING
2026-07-02 14:03:08 INFO None 5066214: status FINISHED
2026-07-02 14:03:08 INFO None 5066216: status FINISHED
2026-07-02 14:03:08 INFO None 5066217: status RUNNING/PENDING
2026-07-02 14:03:08 INFO None 5066219: status RUNNING/PENDING
2026-07-02 14:03:08 INFO None 5066220: status RUNNING/PENDING
2026-07-02 14:03:08 INFO None 5066223: status RUNNING/PENDING
2026-07-02 14:03:08 INFO None 5066242: status RUNNING/PENDING
2026-07-02 14:03:08 INFO None 5066262: status RUNNING/PENDING
2026-07-02 14:03:08 INFO None 5066280: status RUNNING/PENDING
2026-07-02 14:03:08 INFO Jobs still running: ['5066205', '5066212', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 14:03:23 INFO None 5066204: status FINISHED
2026-07-02 14:03:23 INFO None 5066205: status RUNNING/PENDING
2026-07-02 14:03:23 INFO None 5066206: status FINISHED
2026-07-02 14:03:23 INFO None 5066207: status FINISHED
2026-07-02 14:03:23 INFO None 5066209: status FINISHED
2026-07-02 14:03:23 INFO None 5066212: status RUNNING/PENDING
2026-07-02 14:03:23 INFO None 5066214: status FINISHED
2026-07-02 14:03:23 INFO None 5066216: status FINISHED
2026-07-02 14:03:23 INFO None 5066217: status FINISHED
2026-07-02 14:03:23 INFO None 5066219: status RUNNING/PENDING
2026-07-02 14:03:23 INFO None 5066220: status RUNNING/PENDING
2026-07-02 14:03:23 INFO None 5066223: status RUNNING/PENDING
2026-07-02 14:03:23 INFO None 5066242: status RUNNING/PENDING
2026-07-02 14:03:23 INFO None 5066262: status RUNNING/PENDING
2026-07-02 14:03:23 INFO None 5066280: status RUNNING/PENDING
2026-07-02 14:03:23 INFO Jobs still running: ['5066205', '5066212', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 14:03:38 INFO None 5066204: status FINISHED
2026-07-02 14:03:38 INFO None 5066205: status FINISHED
2026-07-02 14:03:38 INFO None 5066206: status FINISHED
2026-07-02 14:03:38 INFO None 5066207: status FINISHED
2026-07-02 14:03:38 INFO None 5066209: status FINISHED
2026-07-02 14:03:38 INFO None 5066212: status RUNNING/PENDING
2026-07-02 14:03:38 INFO None 5066214: status FINISHED
2026-07-02 14:03:38 INFO None 5066216: status FINISHED
2026-07-02 14:03:38 INFO None 5066217: status FINISHED
2026-07-02 14:03:38 INFO None 5066219: status RUNNING/PENDING
2026-07-02 14:03:38 INFO None 5066220: status RUNNING/PENDING
2026-07-02 14:03:38 INFO None 5066223: status RUNNING/PENDING
2026-07-02 14:03:38 INFO None 5066242: status RUNNING/PENDING
2026-07-02 14:03:38 INFO None 5066262: status RUNNING/PENDING
2026-07-02 14:03:38 INFO None 5066280: status RUNNING/PENDING
2026-07-02 14:03:38 INFO Jobs still running: ['5066212', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 14:03:53 INFO None 5066204: status FINISHED
2026-07-02 14:03:53 INFO None 5066205: status FINISHED
2026-07-02 14:03:53 INFO None 5066206: status FINISHED
2026-07-02 14:03:53 INFO None 5066207: status FINISHED
2026-07-02 14:03:53 INFO None 5066209: status FINISHED
2026-07-02 14:03:53 INFO None 5066212: status FINISHED
2026-07-02 14:03:53 INFO None 5066214: status FINISHED
2026-07-02 14:03:53 INFO None 5066216: status FINISHED
2026-07-02 14:03:54 INFO None 5066217: status FINISHED
2026-07-02 14:03:54 INFO None 5066219: status FINISHED
2026-07-02 14:03:54 INFO None 5066220: status RUNNING/PENDING
2026-07-02 14:03:54 INFO None 5066223: status FINISHED
2026-07-02 14:03:54 INFO None 5066242: status RUNNING/PENDING
2026-07-02 14:03:54 INFO None 5066262: status RUNNING/PENDING
2026-07-02 14:03:54 INFO None 5066280: status RUNNING/PENDING
2026-07-02 14:03:54 INFO Jobs still running: ['5066220', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 14:04:09 INFO None 5066204: status FINISHED
2026-07-02 14:04:09 INFO None 5066205: status FINISHED
2026-07-02 14:04:09 INFO None 5066206: status FINISHED
2026-07-02 14:04:09 INFO None 5066207: status FINISHED
2026-07-02 14:04:09 INFO None 5066209: status FINISHED
2026-07-02 14:04:09 INFO None 5066212: status FINISHED
2026-07-02 14:04:09 INFO None 5066214: status FINISHED
2026-07-02 14:04:09 INFO None 5066216: status FINISHED
2026-07-02 14:04:09 INFO None 5066217: status FINISHED
2026-07-02 14:04:09 INFO None 5066219: status FINISHED
2026-07-02 14:04:09 INFO None 5066220: status RUNNING/PENDING
2026-07-02 14:04:09 INFO None 5066223: status FINISHED
2026-07-02 14:04:09 INFO None 5066242: status RUNNING/PENDING
2026-07-02 14:04:09 INFO None 5066262: status RUNNING/PENDING
2026-07-02 14:04:09 INFO None 5066280: status RUNNING/PENDING
2026-07-02 14:04:09 INFO Jobs still running: ['5066220', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 14:04:24 INFO None 5066204: status FINISHED
2026-07-02 14:04:24 INFO None 5066205: status FINISHED
2026-07-02 14:04:24 INFO None 5066206: status FINISHED
2026-07-02 14:04:24 INFO None 5066207: status FINISHED
2026-07-02 14:04:24 INFO None 5066209: status FINISHED
2026-07-02 14:04:24 INFO None 5066212: status FINISHED
2026-07-02 14:04:24 INFO None 5066214: status FINISHED
2026-07-02 14:04:24 INFO None 5066216: status FINISHED
2026-07-02 14:04:24 INFO None 5066217: status FINISHED
2026-07-02 14:04:24 INFO None 5066219: status FINISHED
2026-07-02 14:04:24 INFO None 5066220: status RUNNING/PENDING
2026-07-02 14:04:24 INFO None 5066223: status FINISHED
2026-07-02 14:04:24 INFO None 5066242: status RUNNING/PENDING
2026-07-02 14:04:24 INFO None 5066262: status RUNNING/PENDING
2026-07-02 14:04:24 INFO None 5066280: status RUNNING/PENDING
2026-07-02 14:04:24 INFO Jobs still running: ['5066220', '5066242', '5066262', '5066280']. Waiting...
2026-07-02 14:04:40 INFO None 5066204: status FINISHED
2026-07-02 14:04:41 INFO None 5066205: status FINISHED
2026-07-02 14:04:41 INFO None 5066206: status FINISHED
2026-07-02 14:04:41 INFO None 5066207: status FINISHED
2026-07-02 14:04:41 INFO None 5066209: status FINISHED
2026-07-02 14:04:41 INFO None 5066212: status FINISHED
2026-07-02 14:04:41 INFO None 5066214: status FINISHED
2026-07-02 14:04:41 INFO None 5066216: status FINISHED
2026-07-02 14:04:41 INFO None 5066217: status FINISHED
2026-07-02 14:04:41 INFO None 5066219: status FINISHED
2026-07-02 14:04:41 INFO None 5066220: status FINISHED
2026-07-02 14:04:41 INFO None 5066223: status FINISHED
2026-07-02 14:04:41 INFO None 5066242: status RUNNING/PENDING
2026-07-02 14:04:41 INFO None 5066262: status RUNNING/PENDING
2026-07-02 14:04:41 INFO None 5066280: status RUNNING/PENDING
2026-07-02 14:04:41 INFO Jobs still running: ['5066242', '5066262', '5066280']. Waiting...
2026-07-02 14:04:56 INFO None 5066204: status FINISHED
2026-07-02 14:04:56 INFO None 5066205: status FINISHED
2026-07-02 14:04:56 INFO None 5066206: status FINISHED
2026-07-02 14:04:56 INFO None 5066207: status FINISHED
2026-07-02 14:04:56 INFO None 5066209: status FINISHED
2026-07-02 14:04:56 INFO None 5066212: status FINISHED
2026-07-02 14:04:56 INFO None 5066214: status FINISHED
2026-07-02 14:04:56 INFO None 5066216: status FINISHED
2026-07-02 14:04:56 INFO None 5066217: status FINISHED
2026-07-02 14:04:56 INFO None 5066219: status FINISHED
2026-07-02 14:04:56 INFO None 5066220: status FINISHED
2026-07-02 14:04:56 INFO None 5066223: status FINISHED
2026-07-02 14:04:56 INFO None 5066242: status RUNNING/PENDING
2026-07-02 14:04:56 INFO None 5066262: status RUNNING/PENDING
2026-07-02 14:04:56 INFO None 5066280: status FINISHED
2026-07-02 14:04:56 INFO Jobs still running: ['5066242', '5066262']. Waiting...
2026-07-02 14:05:11 INFO None 5066204: status FINISHED
2026-07-02 14:05:11 INFO None 5066205: status FINISHED
2026-07-02 14:05:11 INFO None 5066206: status FINISHED
2026-07-02 14:05:11 INFO None 5066207: status FINISHED
2026-07-02 14:05:11 INFO None 5066209: status FINISHED
2026-07-02 14:05:11 INFO None 5066212: status FINISHED
2026-07-02 14:05:11 INFO None 5066214: status FINISHED
2026-07-02 14:05:11 INFO None 5066216: status FINISHED
2026-07-02 14:05:11 INFO None 5066217: status FINISHED
2026-07-02 14:05:11 INFO None 5066219: status FINISHED
2026-07-02 14:05:11 INFO None 5066220: status FINISHED
2026-07-02 14:05:11 INFO None 5066223: status FINISHED
2026-07-02 14:05:11 INFO None 5066242: status RUNNING/PENDING
2026-07-02 14:05:11 INFO None 5066262: status RUNNING/PENDING
2026-07-02 14:05:11 INFO None 5066280: status FINISHED
2026-07-02 14:05:11 INFO Jobs still running: ['5066242', '5066262']. Waiting...
2026-07-02 14:05:26 INFO None 5066204: status FINISHED
2026-07-02 14:05:26 INFO None 5066205: status FINISHED
2026-07-02 14:05:26 INFO None 5066206: status FINISHED
2026-07-02 14:05:26 INFO None 5066207: status FINISHED
2026-07-02 14:05:26 INFO None 5066209: status FINISHED
2026-07-02 14:05:26 INFO None 5066212: status FINISHED
2026-07-02 14:05:26 INFO None 5066214: status FINISHED
2026-07-02 14:05:26 INFO None 5066216: status FINISHED
2026-07-02 14:05:26 INFO None 5066217: status FINISHED
2026-07-02 14:05:26 INFO None 5066219: status FINISHED
2026-07-02 14:05:27 INFO None 5066220: status FINISHED
2026-07-02 14:05:27 INFO None 5066223: status FINISHED
2026-07-02 14:05:27 INFO None 5066242: status RUNNING/PENDING
2026-07-02 14:05:27 INFO None 5066262: status RUNNING/PENDING
2026-07-02 14:05:27 INFO None 5066280: status FINISHED
2026-07-02 14:05:27 INFO Jobs still running: ['5066242', '5066262']. Waiting...
2026-07-02 14:05:42 INFO None 5066204: status FINISHED
2026-07-02 14:05:42 INFO None 5066205: status FINISHED
2026-07-02 14:05:42 INFO None 5066206: status FINISHED
2026-07-02 14:05:42 INFO None 5066207: status FINISHED
2026-07-02 14:05:42 INFO None 5066209: status FINISHED
2026-07-02 14:05:42 INFO None 5066212: status FINISHED
2026-07-02 14:05:42 INFO None 5066214: status FINISHED
2026-07-02 14:05:42 INFO None 5066216: status FINISHED
2026-07-02 14:05:42 INFO None 5066217: status FINISHED
2026-07-02 14:05:42 INFO None 5066219: status FINISHED
2026-07-02 14:05:42 INFO None 5066220: status FINISHED
2026-07-02 14:05:42 INFO None 5066223: status FINISHED
2026-07-02 14:05:42 INFO None 5066242: status RUNNING/PENDING
2026-07-02 14:05:42 INFO None 5066262: status RUNNING/PENDING
2026-07-02 14:05:42 INFO None 5066280: status FINISHED
2026-07-02 14:05:42 INFO Jobs still running: ['5066242', '5066262']. Waiting...
2026-07-02 14:05:57 INFO None 5066204: status FINISHED
2026-07-02 14:05:57 INFO None 5066205: status FINISHED
2026-07-02 14:05:57 INFO None 5066206: status FINISHED
2026-07-02 14:05:57 INFO None 5066207: status FINISHED
2026-07-02 14:05:57 INFO None 5066209: status FINISHED
2026-07-02 14:05:57 INFO None 5066212: status FINISHED
2026-07-02 14:05:57 INFO None 5066214: status FINISHED
2026-07-02 14:05:57 INFO None 5066216: status FINISHED
2026-07-02 14:05:57 INFO None 5066217: status FINISHED
2026-07-02 14:05:57 INFO None 5066219: status FINISHED
2026-07-02 14:05:57 INFO None 5066220: status FINISHED
2026-07-02 14:05:57 INFO None 5066223: status FINISHED
2026-07-02 14:05:57 INFO None 5066242: status RUNNING/PENDING
2026-07-02 14:05:57 INFO None 5066262: status RUNNING/PENDING
2026-07-02 14:05:57 INFO None 5066280: status FINISHED
2026-07-02 14:05:57 INFO Jobs still running: ['5066242', '5066262']. Waiting...
2026-07-02 14:06:12 INFO None 5066204: status FINISHED
2026-07-02 14:06:12 INFO None 5066205: status FINISHED
2026-07-02 14:06:12 INFO None 5066206: status FINISHED
2026-07-02 14:06:12 INFO None 5066207: status FINISHED
2026-07-02 14:06:12 INFO None 5066209: status FINISHED
2026-07-02 14:06:12 INFO None 5066212: status FINISHED
2026-07-02 14:06:12 INFO None 5066214: status FINISHED
2026-07-02 14:06:12 INFO None 5066216: status FINISHED
2026-07-02 14:06:12 INFO None 5066217: status FINISHED
2026-07-02 14:06:12 INFO None 5066219: status FINISHED
2026-07-02 14:06:12 INFO None 5066220: status FINISHED
2026-07-02 14:06:12 INFO None 5066223: status FINISHED
2026-07-02 14:06:12 INFO None 5066242: status RUNNING/PENDING
2026-07-02 14:06:12 INFO None 5066262: status RUNNING/PENDING
2026-07-02 14:06:12 INFO None 5066280: status FINISHED
2026-07-02 14:06:12 INFO Jobs still running: ['5066242', '5066262']. Waiting...
2026-07-02 14:06:28 INFO None 5066204: status FINISHED
2026-07-02 14:06:28 INFO None 5066205: status FINISHED
2026-07-02 14:06:28 INFO None 5066206: status FINISHED
2026-07-02 14:06:28 INFO None 5066207: status FINISHED
2026-07-02 14:06:28 INFO None 5066209: status FINISHED
2026-07-02 14:06:28 INFO None 5066212: status FINISHED
2026-07-02 14:06:28 INFO None 5066214: status FINISHED
2026-07-02 14:06:28 INFO None 5066216: status FINISHED
2026-07-02 14:06:28 INFO None 5066217: status FINISHED
2026-07-02 14:06:28 INFO None 5066219: status FINISHED
2026-07-02 14:06:28 INFO None 5066220: status FINISHED
2026-07-02 14:06:28 INFO None 5066223: status FINISHED
2026-07-02 14:06:28 INFO None 5066242: status RUNNING/PENDING
2026-07-02 14:06:28 INFO None 5066262: status RUNNING/PENDING
2026-07-02 14:06:28 INFO None 5066280: status FINISHED
2026-07-02 14:06:28 INFO Jobs still running: ['5066242', '5066262']. Waiting...
2026-07-02 14:06:43 INFO None 5066204: status FINISHED
2026-07-02 14:06:43 INFO None 5066205: status FINISHED
2026-07-02 14:06:43 INFO None 5066206: status FINISHED
2026-07-02 14:06:43 INFO None 5066207: status FINISHED
2026-07-02 14:06:43 INFO None 5066209: status FINISHED
2026-07-02 14:06:43 INFO None 5066212: status FINISHED
2026-07-02 14:06:43 INFO None 5066214: status FINISHED
2026-07-02 14:06:43 INFO None 5066216: status FINISHED
2026-07-02 14:06:43 INFO None 5066217: status FINISHED
2026-07-02 14:06:43 INFO None 5066219: status FINISHED
2026-07-02 14:06:43 INFO None 5066220: status FINISHED
2026-07-02 14:06:43 INFO None 5066223: status FINISHED
2026-07-02 14:06:43 INFO None 5066242: status FINISHED
2026-07-02 14:06:43 INFO None 5066262: status FINISHED
2026-07-02 14:06:43 INFO None 5066280: status FINISHED
2026-07-02 14:06:43 INFO Jobs ['5066204', '5066205', '5066206', '5066207', '5066209', '5066212', '5066214', '5066216', '5066217', '5066219', '5066220', '5066223', '5066242', '5066262', '5066280'] have finished
2026-07-02 14:06:43 INFO Checking restart files were created ...
2026-07-02 14:06:43 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020701_8_ENS1.nc(3005806795 bytes)
2026-07-02 14:06:43 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020701_8_ENS2.nc(3005806795 bytes)
2026-07-02 14:06:43 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020701_8_ENS3.nc(3005806795 bytes)
2026-07-02 14:06:43 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020701_8_ENS4.nc(3005806795 bytes)
2026-07-02 14:06:43 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020701_8_ENS5.nc(3005806795 bytes)
2026-07-02 14:06:43 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020701_8_ENS6.nc(3005806795 bytes)
2026-07-02 14:06:43 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020701_8_ENS7.nc(3005806795 bytes)
2026-07-02 14:06:43 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020701_8_ENS8.nc(3005806795 bytes)
2026-07-02 14:06:43 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020701_8_ENS9.nc(3005806795 bytes)
2026-07-02 14:06:43 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020701_8_ENS10.nc(3005806795 bytes)
2026-07-02 14:06:43 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020701_8_ENS11.nc(3005806795 bytes)
2026-07-02 14:06:43 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020701_8_ENS12.nc(3005806795 bytes)
2026-07-02 14:06:43 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020701_8_ENS13.nc(3005806795 bytes)
2026-07-02 14:06:43 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020701_8_ENS14.nc(3005806795 bytes)
2026-07-02 14:06:43 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020701_8_ENS15.nc(3005806795 bytes)
2026-07-02 14:06:43 INFO  Run_model() completed successfully.
2026-07-02 14:06:43 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 01:00:00 simulated_time=2020-02-07 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 14:06:43 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 09:00:00 days=153073 seconds=32400
2026-07-02 14:06:43 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-02 14:06:43 INFO [TIME] increment current_time 2020-02-07 01:00:00 -> 2020-02-07 09:00:00
2026-07-02 14:06:43 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 09:00:00 simulated_time=2020-02-07 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 14:06:43 INFO ---------->>> Running process_satellite_data()
2026-07-02 14:06:43 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12016.nc
2026-07-02 14:06:43 INFO ---------->>> Running run_obs_converter()
2026-07-02 14:06:43 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_32090_153073.out
2026-07-02 14:06:43 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_32090_153073.out
2026-07-02 14:06:43 INFO ---------->>> Running DART
2026-07-02 14:06:43 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-02 14:06:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/chim_ENS1_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/chim_ENS1_2020020709_1_out_toDART.nc
2026-07-02 14:06:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/chim_ENS2_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/chim_ENS2_2020020709_1_out_toDART.nc
2026-07-02 14:06:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/chim_ENS3_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/chim_ENS3_2020020709_1_out_toDART.nc
2026-07-02 14:06:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/chim_ENS4_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/chim_ENS4_2020020709_1_out_toDART.nc
2026-07-02 14:06:45 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/chim_ENS5_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/chim_ENS5_2020020709_1_out_toDART.nc
2026-07-02 14:06:45 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/chim_ENS6_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/chim_ENS6_2020020709_1_out_toDART.nc
2026-07-02 14:06:45 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/chim_ENS7_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/chim_ENS7_2020020709_1_out_toDART.nc
2026-07-02 14:06:46 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/chim_ENS8_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/chim_ENS8_2020020709_1_out_toDART.nc
2026-07-02 14:06:46 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/chim_ENS9_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/chim_ENS9_2020020709_1_out_toDART.nc
2026-07-02 14:06:46 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/chim_ENS10_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/chim_ENS10_2020020709_1_out_toDART.nc
2026-07-02 14:06:47 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/chim_ENS11_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/chim_ENS11_2020020709_1_out_toDART.nc
2026-07-02 14:06:47 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/chim_ENS12_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/chim_ENS12_2020020709_1_out_toDART.nc
2026-07-02 14:06:47 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/chim_ENS13_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/chim_ENS13_2020020709_1_out_toDART.nc
2026-07-02 14:06:48 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/chim_ENS14_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/chim_ENS14_2020020709_1_out_toDART.nc
2026-07-02 14:06:48 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/chim_ENS15_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/chim_ENS15_2020020709_1_out_toDART.nc
2026-07-02 14:06:48 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-02 14:06:48 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-02 14:06:48 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-02 14:06:49 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-02 14:06:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-02 14:06:49 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-02 14:06:54 INFO Found: []
2026-07-02 14:06:54 INFO No job id returned by command ./run_filter.bsh
2026-07-02 14:06:54 INFO No monitoring will be performed
2026-07-02 14:06:54 INFO Moving DART output files to analysis and preassim directories for date 2020020709 if present ...
2026-07-02 14:06:54 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020709'
2026-07-02 14:06:54 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020709'
2026-07-02 14:06:54 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020709'
2026-07-02 14:06:54 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020709'
2026-07-02 14:06:54 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020709'
2026-07-02 14:06:54 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020709'
2026-07-02 14:06:54 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020709'
2026-07-02 14:06:54 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020709'
2026-07-02 14:06:54 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020709'
2026-07-02 14:06:54 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020709'
2026-07-02 14:06:54 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020709'
2026-07-02 14:06:54 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020709'
2026-07-02 14:06:54 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020709'
2026-07-02 14:06:54 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020709'
2026-07-02 14:06:54 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020709'
2026-07-02 14:06:54 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020709'
2026-07-02 14:06:54 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020709'
2026-07-02 14:06:54 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020709'
2026-07-02 14:06:54 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020709'
2026-07-02 14:06:54 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020709'
2026-07-02 14:06:55 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020709'
2026-07-02 14:06:55 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020709'
2026-07-02 14:06:55 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020709'
2026-07-02 14:06:55 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020709'
2026-07-02 14:06:55 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020709'
2026-07-02 14:06:55 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020709'
2026-07-02 14:06:55 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020709'
2026-07-02 14:06:55 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020709'
2026-07-02 14:06:55 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020709'
2026-07-02 14:06:55 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020709'
2026-07-02 14:06:55 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020709'
2026-07-02 14:06:55 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020709'
2026-07-02 14:06:55 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020709'
2026-07-02 14:06:55 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020709'
2026-07-02 14:06:55 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-02 14:06:55 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-02 14:06:55 INFO run_dart() is DONE.
2026-07-02 14:06:55 INFO ---------->>> Running update_pollutant_in_end()
2026-07-02 14:06:55 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:06:55 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:06:56 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:06:56 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:06:56 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:06:57 INFO [Posterior Diff ENS1] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020709/diff_posterior_ENS1_2020020709.nc
2026-07-02 14:06:57 INFO [Posterior Diff ENS1] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020709/diff_posterior_ENS1_2020020709.relative.nc
2026-07-02 14:06:57 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:06:57 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:06:58 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:06:58 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:06:58 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:06:59 INFO [Posterior Diff ENS2] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020709/diff_posterior_ENS2_2020020709.nc
2026-07-02 14:06:59 INFO [Posterior Diff ENS2] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020709/diff_posterior_ENS2_2020020709.relative.nc
2026-07-02 14:06:59 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:06:59 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:07:00 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:07:00 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:07:00 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:07:01 INFO [Posterior Diff ENS3] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020709/diff_posterior_ENS3_2020020709.nc
2026-07-02 14:07:01 INFO [Posterior Diff ENS3] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020709/diff_posterior_ENS3_2020020709.relative.nc
2026-07-02 14:07:01 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:07:01 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:07:02 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:07:02 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:07:02 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:07:03 INFO [Posterior Diff ENS4] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020709/diff_posterior_ENS4_2020020709.nc
2026-07-02 14:07:03 INFO [Posterior Diff ENS4] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020709/diff_posterior_ENS4_2020020709.relative.nc
2026-07-02 14:07:03 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:07:03 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:07:04 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:07:04 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:07:04 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:07:05 INFO [Posterior Diff ENS5] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020709/diff_posterior_ENS5_2020020709.nc
2026-07-02 14:07:05 INFO [Posterior Diff ENS5] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020709/diff_posterior_ENS5_2020020709.relative.nc
2026-07-02 14:07:05 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:07:05 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:07:06 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:07:06 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:07:06 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:07:07 INFO [Posterior Diff ENS6] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020709/diff_posterior_ENS6_2020020709.nc
2026-07-02 14:07:07 INFO [Posterior Diff ENS6] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020709/diff_posterior_ENS6_2020020709.relative.nc
2026-07-02 14:07:07 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:07:07 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:07:08 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:07:08 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:07:09 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:07:09 INFO [Posterior Diff ENS7] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020709/diff_posterior_ENS7_2020020709.nc
2026-07-02 14:07:09 INFO [Posterior Diff ENS7] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020709/diff_posterior_ENS7_2020020709.relative.nc
2026-07-02 14:07:09 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:07:09 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:07:10 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:07:10 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:07:11 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:07:11 INFO [Posterior Diff ENS8] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020709/diff_posterior_ENS8_2020020709.nc
2026-07-02 14:07:11 INFO [Posterior Diff ENS8] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020709/diff_posterior_ENS8_2020020709.relative.nc
2026-07-02 14:07:11 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:07:11 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:07:12 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:07:12 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:07:13 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:07:13 INFO [Posterior Diff ENS9] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020709/diff_posterior_ENS9_2020020709.nc
2026-07-02 14:07:13 INFO [Posterior Diff ENS9] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020709/diff_posterior_ENS9_2020020709.relative.nc
2026-07-02 14:07:13 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:07:13 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:07:14 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:07:14 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:07:14 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:07:15 INFO [Posterior Diff ENS10] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020709/diff_posterior_ENS10_2020020709.nc
2026-07-02 14:07:15 INFO [Posterior Diff ENS10] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020709/diff_posterior_ENS10_2020020709.relative.nc
2026-07-02 14:07:15 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:07:15 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:07:16 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:07:16 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:07:17 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:07:17 INFO [Posterior Diff ENS11] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020709/diff_posterior_ENS11_2020020709.nc
2026-07-02 14:07:17 INFO [Posterior Diff ENS11] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020709/diff_posterior_ENS11_2020020709.relative.nc
2026-07-02 14:07:17 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:07:17 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:07:18 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:07:18 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:07:19 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:07:19 INFO [Posterior Diff ENS12] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020709/diff_posterior_ENS12_2020020709.nc
2026-07-02 14:07:19 INFO [Posterior Diff ENS12] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020709/diff_posterior_ENS12_2020020709.relative.nc
2026-07-02 14:07:19 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:07:19 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:07:20 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:07:20 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:07:21 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:07:21 INFO [Posterior Diff ENS13] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020709/diff_posterior_ENS13_2020020709.nc
2026-07-02 14:07:21 INFO [Posterior Diff ENS13] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020709/diff_posterior_ENS13_2020020709.relative.nc
2026-07-02 14:07:21 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:07:21 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:07:22 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:07:22 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:07:23 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:07:23 INFO [Posterior Diff ENS14] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020709/diff_posterior_ENS14_2020020709.nc
2026-07-02 14:07:23 INFO [Posterior Diff ENS14] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020709/diff_posterior_ENS14_2020020709.relative.nc
2026-07-02 14:07:23 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:07:23 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:07:24 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:07:24 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:07:25 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:07:25 INFO [Posterior Diff ENS15] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020709/diff_posterior_ENS15_2020020709.nc
2026-07-02 14:07:25 INFO [Posterior Diff ENS15] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020709/diff_posterior_ENS15_2020020709.relative.nc
2026-07-02 14:07:25 INFO Next run starts from 2020-02-07 09:00:00
2026-07-02 14:07:25 INFO Cycle is DONE; starting a new loop!
2026-07-02 14:07:25 INFO [TIME] step_end current_time=2020-02-07 09:00:00 simulated_time=2020-02-07 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 14:07:25 INFO [TIME] step_start current_time=2020-02-07 09:00:00 simulated_time=2020-02-07 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 14:07:25 INFO [TIME] window start=2020-02-07 09:00:00 end=2020-02-07 11:00:00 run_hours=2 has_assimilation=True
2026-07-02 14:07:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:07:26 INFO Hourly dataset computed and listing created
2026-07-02 14:07:30 INFO Hourly dataset computed
2026-07-02 14:07:30 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:07:31 INFO Hourly dataset computed and listing created
2026-07-02 14:07:31 INFO Hourly dataset computed
2026-07-02 14:07:31 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:07:32 INFO Hourly dataset computed and listing created
2026-07-02 14:07:33 INFO Hourly dataset computed
2026-07-02 14:07:33 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:07:34 INFO Hourly dataset computed and listing created
2026-07-02 14:07:35 INFO Hourly dataset computed
2026-07-02 14:07:35 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:07:36 INFO Hourly dataset computed and listing created
2026-07-02 14:07:36 INFO Hourly dataset computed
2026-07-02 14:07:36 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:07:37 INFO Hourly dataset computed and listing created
2026-07-02 14:07:38 INFO Hourly dataset computed
2026-07-02 14:07:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:07:39 INFO Hourly dataset computed and listing created
2026-07-02 14:07:40 INFO Hourly dataset computed
2026-07-02 14:07:40 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:07:41 INFO Hourly dataset computed and listing created
2026-07-02 14:07:41 INFO Hourly dataset computed
2026-07-02 14:07:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:07:42 INFO Hourly dataset computed and listing created
2026-07-02 14:07:43 INFO Hourly dataset computed
2026-07-02 14:07:43 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:07:44 INFO Hourly dataset computed and listing created
2026-07-02 14:07:45 INFO Hourly dataset computed
2026-07-02 14:07:45 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:07:46 INFO Hourly dataset computed and listing created
2026-07-02 14:07:46 INFO Hourly dataset computed
2026-07-02 14:07:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:07:47 INFO Hourly dataset computed and listing created
2026-07-02 14:07:48 INFO Hourly dataset computed
2026-07-02 14:07:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:07:49 INFO Hourly dataset computed and listing created
2026-07-02 14:07:50 INFO Hourly dataset computed
2026-07-02 14:07:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:07:51 INFO Hourly dataset computed and listing created
2026-07-02 14:07:52 INFO Hourly dataset computed
2026-07-02 14:07:52 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:07:52 INFO Hourly dataset computed and listing created
2026-07-02 14:07:53 INFO Hourly dataset computed
2026-07-02 14:07:53 INFO ---------->>> Running CHIMERE model from 2020-02-07 09:00:00 to 2020-02-07 11:00:00
2026-07-02 14:07:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:07:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1
2026-07-02 14:07:53 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020701_8_ENS1.nc
2026-07-02 14:07:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-02 14:07:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:07:53 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-02 14:07:53 INFO Queuing job for member 1...
2026-07-02 14:07:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:07:53 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-02 14:07:54 INFO Found: ['5066357']
2026-07-02 14:07:59 INFO [TGCC-IRENE] Submitted job with ID:['5066357']
2026-07-02 14:07:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:07:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2
2026-07-02 14:07:59 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020701_8_ENS2.nc
2026-07-02 14:07:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-02 14:07:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:07:59 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-02 14:07:59 INFO Queuing job for member 2...
2026-07-02 14:07:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:07:59 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-02 14:08:00 INFO Found: ['5066358']
2026-07-02 14:08:05 INFO [TGCC-IRENE] Submitted job with ID:['5066358']
2026-07-02 14:08:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:08:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3
2026-07-02 14:08:05 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020701_8_ENS3.nc
2026-07-02 14:08:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-02 14:08:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:08:05 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-02 14:08:05 INFO Queuing job for member 3...
2026-07-02 14:08:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:08:05 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-02 14:08:07 INFO Found: ['5066360']
2026-07-02 14:08:12 INFO [TGCC-IRENE] Submitted job with ID:['5066360']
2026-07-02 14:08:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:08:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4
2026-07-02 14:08:12 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020701_8_ENS4.nc
2026-07-02 14:08:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-02 14:08:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:08:12 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-02 14:08:12 INFO Queuing job for member 4...
2026-07-02 14:08:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:08:12 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-02 14:08:13 INFO Found: ['5066361']
2026-07-02 14:08:18 INFO [TGCC-IRENE] Submitted job with ID:['5066361']
2026-07-02 14:08:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:08:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5
2026-07-02 14:08:18 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020701_8_ENS5.nc
2026-07-02 14:08:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-02 14:08:18 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:08:18 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-02 14:08:18 INFO Queuing job for member 5...
2026-07-02 14:08:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:08:18 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-02 14:08:19 INFO Found: ['5066362']
2026-07-02 14:08:24 INFO [TGCC-IRENE] Submitted job with ID:['5066362']
2026-07-02 14:08:24 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:08:24 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6
2026-07-02 14:08:24 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020701_8_ENS6.nc
2026-07-02 14:08:24 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-02 14:08:24 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:08:24 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-02 14:08:24 INFO Queuing job for member 6...
2026-07-02 14:08:24 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:08:24 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-02 14:08:25 INFO Found: ['5066363']
2026-07-02 14:08:30 INFO [TGCC-IRENE] Submitted job with ID:['5066363']
2026-07-02 14:08:30 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:08:30 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7
2026-07-02 14:08:30 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020701_8_ENS7.nc
2026-07-02 14:08:30 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-02 14:08:30 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:08:30 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-02 14:08:30 INFO Queuing job for member 7...
2026-07-02 14:08:30 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:08:30 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-02 14:08:30 INFO Found: ['5066365']
2026-07-02 14:08:35 INFO [TGCC-IRENE] Submitted job with ID:['5066365']
2026-07-02 14:08:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:08:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8
2026-07-02 14:08:35 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020701_8_ENS8.nc
2026-07-02 14:08:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-02 14:08:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:08:35 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-02 14:08:35 INFO Queuing job for member 8...
2026-07-02 14:08:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:08:35 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-02 14:08:36 INFO Found: ['5066366']
2026-07-02 14:08:41 INFO [TGCC-IRENE] Submitted job with ID:['5066366']
2026-07-02 14:08:41 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:08:41 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9
2026-07-02 14:08:41 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020701_8_ENS9.nc
2026-07-02 14:08:41 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-02 14:08:41 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:08:41 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-02 14:08:41 INFO Queuing job for member 9...
2026-07-02 14:08:41 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:08:41 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-02 14:08:42 INFO Found: ['5066368']
2026-07-02 14:08:47 INFO [TGCC-IRENE] Submitted job with ID:['5066368']
2026-07-02 14:08:47 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:08:47 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10
2026-07-02 14:08:47 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020701_8_ENS10.nc
2026-07-02 14:08:47 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-02 14:08:47 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:08:47 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-02 14:08:47 INFO Queuing job for member 10...
2026-07-02 14:08:47 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:08:47 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-02 14:08:48 INFO Found: ['5066369']
2026-07-02 14:08:53 INFO [TGCC-IRENE] Submitted job with ID:['5066369']
2026-07-02 14:08:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:08:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11
2026-07-02 14:08:53 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020701_8_ENS11.nc
2026-07-02 14:08:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-02 14:08:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:08:53 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-02 14:08:53 INFO Queuing job for member 11...
2026-07-02 14:08:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:08:53 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-02 14:08:53 INFO Found: ['5066371']
2026-07-02 14:08:58 INFO [TGCC-IRENE] Submitted job with ID:['5066371']
2026-07-02 14:08:58 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:08:58 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12
2026-07-02 14:08:58 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020701_8_ENS12.nc
2026-07-02 14:08:58 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-02 14:08:58 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:08:58 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-02 14:08:58 INFO Queuing job for member 12...
2026-07-02 14:08:58 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:08:58 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-02 14:09:00 INFO Found: ['5066372']
2026-07-02 14:09:05 INFO [TGCC-IRENE] Submitted job with ID:['5066372']
2026-07-02 14:09:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:09:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13
2026-07-02 14:09:05 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020701_8_ENS13.nc
2026-07-02 14:09:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-02 14:09:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:09:05 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-02 14:09:05 INFO Queuing job for member 13...
2026-07-02 14:09:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:09:05 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-02 14:09:06 INFO Found: ['5066375']
2026-07-02 14:09:11 INFO [TGCC-IRENE] Submitted job with ID:['5066375']
2026-07-02 14:09:11 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:09:11 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14
2026-07-02 14:09:11 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020701_8_ENS14.nc
2026-07-02 14:09:11 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-02 14:09:11 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:09:11 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-02 14:09:11 INFO Queuing job for member 14...
2026-07-02 14:09:11 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:09:11 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-02 14:09:12 INFO Found: ['5066376']
2026-07-02 14:09:17 INFO [TGCC-IRENE] Submitted job with ID:['5066376']
2026-07-02 14:09:17 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:09:17 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15
2026-07-02 14:09:17 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020701_8_ENS15.nc
2026-07-02 14:09:17 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-02 14:09:17 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:09:17 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-02 14:09:17 INFO Queuing job for member 15...
2026-07-02 14:09:17 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:09:17 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-02 14:09:18 INFO Found: ['5066377']
2026-07-02 14:09:23 INFO [TGCC-IRENE] Submitted job with ID:['5066377']
2026-07-02 14:09:23 INFO Checking job status ...
2026-07-02 14:09:23 INFO None 5066357: status RUNNING/PENDING
2026-07-02 14:09:23 INFO None 5066358: status RUNNING/PENDING
2026-07-02 14:09:23 INFO None 5066360: status RUNNING/PENDING
2026-07-02 14:09:23 INFO None 5066361: status RUNNING/PENDING
2026-07-02 14:09:23 INFO None 5066362: status RUNNING/PENDING
2026-07-02 14:09:23 INFO None 5066363: status RUNNING/PENDING
2026-07-02 14:09:23 INFO None 5066365: status RUNNING/PENDING
2026-07-02 14:09:23 INFO None 5066366: status RUNNING/PENDING
2026-07-02 14:09:23 INFO None 5066368: status RUNNING/PENDING
2026-07-02 14:09:23 INFO None 5066369: status RUNNING/PENDING
2026-07-02 14:09:23 INFO None 5066371: status RUNNING/PENDING
2026-07-02 14:09:23 INFO None 5066372: status RUNNING/PENDING
2026-07-02 14:09:23 INFO None 5066375: status RUNNING/PENDING
2026-07-02 14:09:23 INFO None 5066376: status RUNNING/PENDING
2026-07-02 14:09:23 INFO None 5066377: status RUNNING/PENDING
2026-07-02 14:09:23 INFO Jobs still running: ['5066357', '5066358', '5066360', '5066361', '5066362', '5066363', '5066365', '5066366', '5066368', '5066369', '5066371', '5066372', '5066375', '5066376', '5066377']. Waiting...
2026-07-02 14:09:38 INFO None 5066357: status RUNNING/PENDING
2026-07-02 14:09:38 INFO None 5066358: status RUNNING/PENDING
2026-07-02 14:09:38 INFO None 5066360: status RUNNING/PENDING
2026-07-02 14:09:38 INFO None 5066361: status RUNNING/PENDING
2026-07-02 14:09:38 INFO None 5066362: status RUNNING/PENDING
2026-07-02 14:09:38 INFO None 5066363: status RUNNING/PENDING
2026-07-02 14:09:38 INFO None 5066365: status RUNNING/PENDING
2026-07-02 14:09:38 INFO None 5066366: status RUNNING/PENDING
2026-07-02 14:09:38 INFO None 5066368: status RUNNING/PENDING
2026-07-02 14:09:38 INFO None 5066369: status RUNNING/PENDING
2026-07-02 14:09:38 INFO None 5066371: status RUNNING/PENDING
2026-07-02 14:09:38 INFO None 5066372: status RUNNING/PENDING
2026-07-02 14:09:38 INFO None 5066375: status RUNNING/PENDING
2026-07-02 14:09:38 INFO None 5066376: status RUNNING/PENDING
2026-07-02 14:09:38 INFO None 5066377: status RUNNING/PENDING
2026-07-02 14:09:38 INFO Jobs still running: ['5066357', '5066358', '5066360', '5066361', '5066362', '5066363', '5066365', '5066366', '5066368', '5066369', '5066371', '5066372', '5066375', '5066376', '5066377']. Waiting...
2026-07-02 14:09:54 INFO None 5066357: status RUNNING/PENDING
2026-07-02 14:09:54 INFO None 5066358: status RUNNING/PENDING
2026-07-02 14:09:55 INFO None 5066360: status RUNNING/PENDING
2026-07-02 14:09:55 INFO None 5066361: status RUNNING/PENDING
2026-07-02 14:09:55 INFO None 5066362: status RUNNING/PENDING
2026-07-02 14:09:55 INFO None 5066363: status RUNNING/PENDING
2026-07-02 14:09:55 INFO None 5066365: status RUNNING/PENDING
2026-07-02 14:09:55 INFO None 5066366: status RUNNING/PENDING
2026-07-02 14:09:55 INFO None 5066368: status RUNNING/PENDING
2026-07-02 14:09:55 INFO None 5066369: status RUNNING/PENDING
2026-07-02 14:09:55 INFO None 5066371: status RUNNING/PENDING
2026-07-02 14:09:55 INFO None 5066372: status RUNNING/PENDING
2026-07-02 14:09:55 INFO None 5066375: status RUNNING/PENDING
2026-07-02 14:09:55 INFO None 5066376: status RUNNING/PENDING
2026-07-02 14:09:55 INFO None 5066377: status RUNNING/PENDING
2026-07-02 14:09:55 INFO Jobs still running: ['5066357', '5066358', '5066360', '5066361', '5066362', '5066363', '5066365', '5066366', '5066368', '5066369', '5066371', '5066372', '5066375', '5066376', '5066377']. Waiting...
2026-07-02 14:10:10 INFO None 5066357: status RUNNING/PENDING
2026-07-02 14:10:10 INFO None 5066358: status RUNNING/PENDING
2026-07-02 14:10:10 INFO None 5066360: status RUNNING/PENDING
2026-07-02 14:10:10 INFO None 5066361: status RUNNING/PENDING
2026-07-02 14:10:10 INFO None 5066362: status RUNNING/PENDING
2026-07-02 14:10:10 INFO None 5066363: status RUNNING/PENDING
2026-07-02 14:10:10 INFO None 5066365: status RUNNING/PENDING
2026-07-02 14:10:10 INFO None 5066366: status RUNNING/PENDING
2026-07-02 14:10:10 INFO None 5066368: status RUNNING/PENDING
2026-07-02 14:10:10 INFO None 5066369: status RUNNING/PENDING
2026-07-02 14:10:10 INFO None 5066371: status RUNNING/PENDING
2026-07-02 14:10:10 INFO None 5066372: status RUNNING/PENDING
2026-07-02 14:10:10 INFO None 5066375: status RUNNING/PENDING
2026-07-02 14:10:10 INFO None 5066376: status RUNNING/PENDING
2026-07-02 14:10:10 INFO None 5066377: status RUNNING/PENDING
2026-07-02 14:10:10 INFO Jobs still running: ['5066357', '5066358', '5066360', '5066361', '5066362', '5066363', '5066365', '5066366', '5066368', '5066369', '5066371', '5066372', '5066375', '5066376', '5066377']. Waiting...
2026-07-02 14:10:25 INFO None 5066357: status RUNNING/PENDING
2026-07-02 14:10:25 INFO None 5066358: status RUNNING/PENDING
2026-07-02 14:10:25 INFO None 5066360: status RUNNING/PENDING
2026-07-02 14:10:25 INFO None 5066361: status RUNNING/PENDING
2026-07-02 14:10:25 INFO None 5066362: status RUNNING/PENDING
2026-07-02 14:10:25 INFO None 5066363: status RUNNING/PENDING
2026-07-02 14:10:25 INFO None 5066365: status RUNNING/PENDING
2026-07-02 14:10:25 INFO None 5066366: status RUNNING/PENDING
2026-07-02 14:10:25 INFO None 5066368: status RUNNING/PENDING
2026-07-02 14:10:25 INFO None 5066369: status RUNNING/PENDING
2026-07-02 14:10:25 INFO None 5066371: status RUNNING/PENDING
2026-07-02 14:10:25 INFO None 5066372: status RUNNING/PENDING
2026-07-02 14:10:25 INFO None 5066375: status RUNNING/PENDING
2026-07-02 14:10:25 INFO None 5066376: status RUNNING/PENDING
2026-07-02 14:10:25 INFO None 5066377: status RUNNING/PENDING
2026-07-02 14:10:25 INFO Jobs still running: ['5066357', '5066358', '5066360', '5066361', '5066362', '5066363', '5066365', '5066366', '5066368', '5066369', '5066371', '5066372', '5066375', '5066376', '5066377']. Waiting...
2026-07-02 14:10:40 INFO None 5066357: status RUNNING/PENDING
2026-07-02 14:10:40 INFO None 5066358: status RUNNING/PENDING
2026-07-02 14:10:40 INFO None 5066360: status RUNNING/PENDING
2026-07-02 14:10:40 INFO None 5066361: status RUNNING/PENDING
2026-07-02 14:10:40 INFO None 5066362: status RUNNING/PENDING
2026-07-02 14:10:40 INFO None 5066363: status RUNNING/PENDING
2026-07-02 14:10:40 INFO None 5066365: status RUNNING/PENDING
2026-07-02 14:10:40 INFO None 5066366: status RUNNING/PENDING
2026-07-02 14:10:40 INFO None 5066368: status RUNNING/PENDING
2026-07-02 14:10:41 INFO None 5066369: status RUNNING/PENDING
2026-07-02 14:10:41 INFO None 5066371: status RUNNING/PENDING
2026-07-02 14:10:41 INFO None 5066372: status RUNNING/PENDING
2026-07-02 14:10:41 INFO None 5066375: status RUNNING/PENDING
2026-07-02 14:10:41 INFO None 5066376: status RUNNING/PENDING
2026-07-02 14:10:41 INFO None 5066377: status RUNNING/PENDING
2026-07-02 14:10:41 INFO Jobs still running: ['5066357', '5066358', '5066360', '5066361', '5066362', '5066363', '5066365', '5066366', '5066368', '5066369', '5066371', '5066372', '5066375', '5066376', '5066377']. Waiting...
2026-07-02 14:10:56 INFO None 5066357: status RUNNING/PENDING
2026-07-02 14:10:56 INFO None 5066358: status FINISHED
2026-07-02 14:10:56 INFO None 5066360: status RUNNING/PENDING
2026-07-02 14:10:56 INFO None 5066361: status RUNNING/PENDING
2026-07-02 14:10:56 INFO None 5066362: status RUNNING/PENDING
2026-07-02 14:10:56 INFO None 5066363: status RUNNING/PENDING
2026-07-02 14:10:56 INFO None 5066365: status RUNNING/PENDING
2026-07-02 14:10:56 INFO None 5066366: status RUNNING/PENDING
2026-07-02 14:10:56 INFO None 5066368: status RUNNING/PENDING
2026-07-02 14:10:56 INFO None 5066369: status RUNNING/PENDING
2026-07-02 14:10:56 INFO None 5066371: status RUNNING/PENDING
2026-07-02 14:10:56 INFO None 5066372: status RUNNING/PENDING
2026-07-02 14:10:56 INFO None 5066375: status RUNNING/PENDING
2026-07-02 14:10:56 INFO None 5066376: status RUNNING/PENDING
2026-07-02 14:10:56 INFO None 5066377: status RUNNING/PENDING
2026-07-02 14:10:56 INFO Jobs still running: ['5066357', '5066360', '5066361', '5066362', '5066363', '5066365', '5066366', '5066368', '5066369', '5066371', '5066372', '5066375', '5066376', '5066377']. Waiting...
2026-07-02 14:11:11 INFO None 5066357: status FINISHED
2026-07-02 14:11:11 INFO None 5066358: status FINISHED
2026-07-02 14:11:11 INFO None 5066360: status RUNNING/PENDING
2026-07-02 14:11:11 INFO None 5066361: status RUNNING/PENDING
2026-07-02 14:11:11 INFO None 5066362: status RUNNING/PENDING
2026-07-02 14:11:11 INFO None 5066363: status RUNNING/PENDING
2026-07-02 14:11:11 INFO None 5066365: status RUNNING/PENDING
2026-07-02 14:11:11 INFO None 5066366: status RUNNING/PENDING
2026-07-02 14:11:11 INFO None 5066368: status RUNNING/PENDING
2026-07-02 14:11:11 INFO None 5066369: status RUNNING/PENDING
2026-07-02 14:11:11 INFO None 5066371: status RUNNING/PENDING
2026-07-02 14:11:11 INFO None 5066372: status RUNNING/PENDING
2026-07-02 14:11:11 INFO None 5066375: status RUNNING/PENDING
2026-07-02 14:11:11 INFO None 5066376: status RUNNING/PENDING
2026-07-02 14:11:11 INFO None 5066377: status RUNNING/PENDING
2026-07-02 14:11:11 INFO Jobs still running: ['5066360', '5066361', '5066362', '5066363', '5066365', '5066366', '5066368', '5066369', '5066371', '5066372', '5066375', '5066376', '5066377']. Waiting...
2026-07-02 14:11:26 INFO None 5066357: status FINISHED
2026-07-02 14:11:26 INFO None 5066358: status FINISHED
2026-07-02 14:11:26 INFO None 5066360: status RUNNING/PENDING
2026-07-02 14:11:26 INFO None 5066361: status RUNNING/PENDING
2026-07-02 14:11:26 INFO None 5066362: status RUNNING/PENDING
2026-07-02 14:11:26 INFO None 5066363: status RUNNING/PENDING
2026-07-02 14:11:26 INFO None 5066365: status RUNNING/PENDING
2026-07-02 14:11:26 INFO None 5066366: status RUNNING/PENDING
2026-07-02 14:11:26 INFO None 5066368: status RUNNING/PENDING
2026-07-02 14:11:26 INFO None 5066369: status RUNNING/PENDING
2026-07-02 14:11:26 INFO None 5066371: status RUNNING/PENDING
2026-07-02 14:11:26 INFO None 5066372: status RUNNING/PENDING
2026-07-02 14:11:26 INFO None 5066375: status RUNNING/PENDING
2026-07-02 14:11:26 INFO None 5066376: status RUNNING/PENDING
2026-07-02 14:11:26 INFO None 5066377: status RUNNING/PENDING
2026-07-02 14:11:26 INFO Jobs still running: ['5066360', '5066361', '5066362', '5066363', '5066365', '5066366', '5066368', '5066369', '5066371', '5066372', '5066375', '5066376', '5066377']. Waiting...
2026-07-02 14:11:43 INFO None 5066357: status FINISHED
2026-07-02 14:11:43 INFO None 5066358: status FINISHED
2026-07-02 14:11:43 INFO None 5066360: status RUNNING/PENDING
2026-07-02 14:11:43 INFO None 5066361: status RUNNING/PENDING
2026-07-02 14:11:43 INFO None 5066362: status RUNNING/PENDING
2026-07-02 14:11:43 INFO None 5066363: status FINISHED
2026-07-02 14:11:43 INFO None 5066365: status FINISHED
2026-07-02 14:11:44 INFO None 5066366: status RUNNING/PENDING
2026-07-02 14:11:44 INFO None 5066368: status RUNNING/PENDING
2026-07-02 14:11:44 INFO None 5066369: status RUNNING/PENDING
2026-07-02 14:11:44 INFO None 5066371: status RUNNING/PENDING
2026-07-02 14:11:44 INFO None 5066372: status RUNNING/PENDING
2026-07-02 14:11:44 INFO None 5066375: status RUNNING/PENDING
2026-07-02 14:11:44 INFO None 5066376: status RUNNING/PENDING
2026-07-02 14:11:44 INFO None 5066377: status RUNNING/PENDING
2026-07-02 14:11:44 INFO Jobs still running: ['5066360', '5066361', '5066362', '5066366', '5066368', '5066369', '5066371', '5066372', '5066375', '5066376', '5066377']. Waiting...
2026-07-02 14:11:59 INFO None 5066357: status FINISHED
2026-07-02 14:11:59 INFO None 5066358: status FINISHED
2026-07-02 14:11:59 INFO None 5066360: status RUNNING/PENDING
2026-07-02 14:11:59 INFO None 5066361: status FINISHED
2026-07-02 14:11:59 INFO None 5066362: status RUNNING/PENDING
2026-07-02 14:11:59 INFO None 5066363: status FINISHED
2026-07-02 14:11:59 INFO None 5066365: status FINISHED
2026-07-02 14:11:59 INFO None 5066366: status RUNNING/PENDING
2026-07-02 14:11:59 INFO None 5066368: status RUNNING/PENDING
2026-07-02 14:11:59 INFO None 5066369: status RUNNING/PENDING
2026-07-02 14:11:59 INFO None 5066371: status RUNNING/PENDING
2026-07-02 14:11:59 INFO None 5066372: status RUNNING/PENDING
2026-07-02 14:11:59 INFO None 5066375: status RUNNING/PENDING
2026-07-02 14:11:59 INFO None 5066376: status RUNNING/PENDING
2026-07-02 14:11:59 INFO None 5066377: status RUNNING/PENDING
2026-07-02 14:11:59 INFO Jobs still running: ['5066360', '5066362', '5066366', '5066368', '5066369', '5066371', '5066372', '5066375', '5066376', '5066377']. Waiting...
2026-07-02 14:12:14 INFO None 5066357: status FINISHED
2026-07-02 14:12:14 INFO None 5066358: status FINISHED
2026-07-02 14:12:14 INFO None 5066360: status RUNNING/PENDING
2026-07-02 14:12:14 INFO None 5066361: status FINISHED
2026-07-02 14:12:14 INFO None 5066362: status FINISHED
2026-07-02 14:12:14 INFO None 5066363: status FINISHED
2026-07-02 14:12:14 INFO None 5066365: status FINISHED
2026-07-02 14:12:14 INFO None 5066366: status RUNNING/PENDING
2026-07-02 14:12:14 INFO None 5066368: status RUNNING/PENDING
2026-07-02 14:12:14 INFO None 5066369: status FINISHED
2026-07-02 14:12:14 INFO None 5066371: status FINISHED
2026-07-02 14:12:14 INFO None 5066372: status RUNNING/PENDING
2026-07-02 14:12:14 INFO None 5066375: status RUNNING/PENDING
2026-07-02 14:12:14 INFO None 5066376: status RUNNING/PENDING
2026-07-02 14:12:14 INFO None 5066377: status RUNNING/PENDING
2026-07-02 14:12:14 INFO Jobs still running: ['5066360', '5066366', '5066368', '5066372', '5066375', '5066376', '5066377']. Waiting...
2026-07-02 14:12:31 INFO None 5066357: status FINISHED
2026-07-02 14:12:31 INFO None 5066358: status FINISHED
2026-07-02 14:12:31 INFO None 5066360: status RUNNING/PENDING
2026-07-02 14:12:31 INFO None 5066361: status FINISHED
2026-07-02 14:12:31 INFO None 5066362: status FINISHED
2026-07-02 14:12:31 INFO None 5066363: status FINISHED
2026-07-02 14:12:31 INFO None 5066365: status FINISHED
2026-07-02 14:12:31 INFO None 5066366: status RUNNING/PENDING
2026-07-02 14:12:31 INFO None 5066368: status RUNNING/PENDING
2026-07-02 14:12:31 INFO None 5066369: status FINISHED
2026-07-02 14:12:31 INFO None 5066371: status FINISHED
2026-07-02 14:12:31 INFO None 5066372: status RUNNING/PENDING
2026-07-02 14:12:31 INFO None 5066375: status RUNNING/PENDING
2026-07-02 14:12:31 INFO None 5066376: status RUNNING/PENDING
2026-07-02 14:12:31 INFO None 5066377: status FINISHED
2026-07-02 14:12:31 INFO Jobs still running: ['5066360', '5066366', '5066368', '5066372', '5066375', '5066376']. Waiting...
2026-07-02 14:12:46 INFO None 5066357: status FINISHED
2026-07-02 14:12:46 INFO None 5066358: status FINISHED
2026-07-02 14:12:46 INFO None 5066360: status RUNNING/PENDING
2026-07-02 14:12:46 INFO None 5066361: status FINISHED
2026-07-02 14:12:46 INFO None 5066362: status FINISHED
2026-07-02 14:12:46 INFO None 5066363: status FINISHED
2026-07-02 14:12:46 INFO None 5066365: status FINISHED
2026-07-02 14:12:46 INFO None 5066366: status FINISHED
2026-07-02 14:12:46 INFO None 5066368: status RUNNING/PENDING
2026-07-02 14:12:46 INFO None 5066369: status FINISHED
2026-07-02 14:12:46 INFO None 5066371: status FINISHED
2026-07-02 14:12:46 INFO None 5066372: status RUNNING/PENDING
2026-07-02 14:12:46 INFO None 5066375: status RUNNING/PENDING
2026-07-02 14:12:46 INFO None 5066376: status RUNNING/PENDING
2026-07-02 14:12:47 INFO None 5066377: status FINISHED
2026-07-02 14:12:47 INFO Jobs still running: ['5066360', '5066368', '5066372', '5066375', '5066376']. Waiting...
2026-07-02 14:13:02 INFO None 5066357: status FINISHED
2026-07-02 14:13:02 INFO None 5066358: status FINISHED
2026-07-02 14:13:02 INFO None 5066360: status FINISHED
2026-07-02 14:13:02 INFO None 5066361: status FINISHED
2026-07-02 14:13:02 INFO None 5066362: status FINISHED
2026-07-02 14:13:02 INFO None 5066363: status FINISHED
2026-07-02 14:13:02 INFO None 5066365: status FINISHED
2026-07-02 14:13:02 INFO None 5066366: status FINISHED
2026-07-02 14:13:02 INFO None 5066368: status FINISHED
2026-07-02 14:13:02 INFO None 5066369: status FINISHED
2026-07-02 14:13:02 INFO None 5066371: status FINISHED
2026-07-02 14:13:02 INFO None 5066372: status FINISHED
2026-07-02 14:13:02 INFO None 5066375: status FINISHED
2026-07-02 14:13:02 INFO None 5066376: status FINISHED
2026-07-02 14:13:02 INFO None 5066377: status FINISHED
2026-07-02 14:13:02 INFO Jobs ['5066357', '5066358', '5066360', '5066361', '5066362', '5066363', '5066365', '5066366', '5066368', '5066369', '5066371', '5066372', '5066375', '5066376', '5066377'] have finished
2026-07-02 14:13:02 INFO Checking restart files were created ...
2026-07-02 14:13:02 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020709_2_ENS1.nc(1002685915 bytes)
2026-07-02 14:13:02 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020709_2_ENS2.nc(1002685915 bytes)
2026-07-02 14:13:02 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020709_2_ENS3.nc(1002685915 bytes)
2026-07-02 14:13:02 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020709_2_ENS4.nc(1002685915 bytes)
2026-07-02 14:13:02 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020709_2_ENS5.nc(1002685915 bytes)
2026-07-02 14:13:02 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020709_2_ENS6.nc(1002685915 bytes)
2026-07-02 14:13:02 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020709_2_ENS7.nc(1002685915 bytes)
2026-07-02 14:13:02 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020709_2_ENS8.nc(1002685915 bytes)
2026-07-02 14:13:02 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020709_2_ENS9.nc(1002685915 bytes)
2026-07-02 14:13:02 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020709_2_ENS10.nc(1002685915 bytes)
2026-07-02 14:13:02 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020709_2_ENS11.nc(1002685915 bytes)
2026-07-02 14:13:02 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020709_2_ENS12.nc(1002685915 bytes)
2026-07-02 14:13:02 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020709_2_ENS13.nc(1002685915 bytes)
2026-07-02 14:13:02 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020709_2_ENS14.nc(1002685915 bytes)
2026-07-02 14:13:02 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020709_2_ENS15.nc(1002685915 bytes)
2026-07-02 14:13:02 INFO  Run_model() completed successfully.
2026-07-02 14:13:02 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 09:00:00 simulated_time=2020-02-07 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 14:13:02 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 11:00:00 days=153073 seconds=39600
2026-07-02 14:13:02 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-02 14:13:02 INFO [TIME] increment current_time 2020-02-07 09:00:00 -> 2020-02-07 11:00:00
2026-07-02 14:13:02 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 11:00:00 simulated_time=2020-02-07 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 14:13:02 INFO ---------->>> Running process_satellite_data()
2026-07-02 14:13:02 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12017.nc
2026-07-02 14:13:02 INFO ---------->>> Running run_obs_converter()
2026-07-02 14:13:02 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_37907_153073.out
2026-07-02 14:13:02 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_37907_153073.out
2026-07-02 14:13:02 INFO ---------->>> Running DART
2026-07-02 14:13:02 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-02 14:13:02 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/chim_ENS1_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/chim_ENS1_2020020711_1_out_toDART.nc
2026-07-02 14:13:02 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/chim_ENS2_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/chim_ENS2_2020020711_1_out_toDART.nc
2026-07-02 14:13:03 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/chim_ENS3_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/chim_ENS3_2020020711_1_out_toDART.nc
2026-07-02 14:13:03 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/chim_ENS4_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/chim_ENS4_2020020711_1_out_toDART.nc
2026-07-02 14:13:03 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/chim_ENS5_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/chim_ENS5_2020020711_1_out_toDART.nc
2026-07-02 14:13:04 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/chim_ENS6_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/chim_ENS6_2020020711_1_out_toDART.nc
2026-07-02 14:13:04 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/chim_ENS7_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/chim_ENS7_2020020711_1_out_toDART.nc
2026-07-02 14:13:04 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/chim_ENS8_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/chim_ENS8_2020020711_1_out_toDART.nc
2026-07-02 14:13:05 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/chim_ENS9_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/chim_ENS9_2020020711_1_out_toDART.nc
2026-07-02 14:13:05 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/chim_ENS10_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/chim_ENS10_2020020711_1_out_toDART.nc
2026-07-02 14:13:05 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/chim_ENS11_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/chim_ENS11_2020020711_1_out_toDART.nc
2026-07-02 14:13:05 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/chim_ENS12_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/chim_ENS12_2020020711_1_out_toDART.nc
2026-07-02 14:13:06 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/chim_ENS13_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/chim_ENS13_2020020711_1_out_toDART.nc
2026-07-02 14:13:06 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/chim_ENS14_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/chim_ENS14_2020020711_1_out_toDART.nc
2026-07-02 14:13:06 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/chim_ENS15_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/chim_ENS15_2020020711_1_out_toDART.nc
2026-07-02 14:13:07 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-02 14:13:07 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-02 14:13:07 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-02 14:13:07 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-02 14:13:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-02 14:13:07 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-02 14:13:17 INFO Found: []
2026-07-02 14:13:17 INFO No job id returned by command ./run_filter.bsh
2026-07-02 14:13:17 INFO No monitoring will be performed
2026-07-02 14:13:17 INFO Moving DART output files to analysis and preassim directories for date 2020020711 if present ...
2026-07-02 14:13:17 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020711'
2026-07-02 14:13:17 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020711'
2026-07-02 14:13:17 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020711'
2026-07-02 14:13:17 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020711'
2026-07-02 14:13:17 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020711'
2026-07-02 14:13:17 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020711'
2026-07-02 14:13:17 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020711'
2026-07-02 14:13:17 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020711'
2026-07-02 14:13:17 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020711'
2026-07-02 14:13:17 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020711'
2026-07-02 14:13:17 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020711'
2026-07-02 14:13:17 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020711'
2026-07-02 14:13:17 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020711'
2026-07-02 14:13:17 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020711'
2026-07-02 14:13:17 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020711'
2026-07-02 14:13:17 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020711'
2026-07-02 14:13:17 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020711'
2026-07-02 14:13:17 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020711'
2026-07-02 14:13:17 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020711'
2026-07-02 14:13:17 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020711'
2026-07-02 14:13:17 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020711'
2026-07-02 14:13:17 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020711'
2026-07-02 14:13:17 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020711'
2026-07-02 14:13:17 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020711'
2026-07-02 14:13:17 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020711'
2026-07-02 14:13:17 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020711'
2026-07-02 14:13:17 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020711'
2026-07-02 14:13:17 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020711'
2026-07-02 14:13:17 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020711'
2026-07-02 14:13:17 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020711'
2026-07-02 14:13:17 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020711'
2026-07-02 14:13:17 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020711'
2026-07-02 14:13:17 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020711'
2026-07-02 14:13:17 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020711'
2026-07-02 14:13:17 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-02 14:13:17 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-02 14:13:17 INFO run_dart() is DONE.
2026-07-02 14:13:17 INFO ---------->>> Running update_pollutant_in_end()
2026-07-02 14:13:18 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:13:18 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:13:18 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:13:18 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:13:19 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:13:19 INFO [Posterior Diff ENS1] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020711/diff_posterior_ENS1_2020020711.nc
2026-07-02 14:13:19 INFO [Posterior Diff ENS1] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020711/diff_posterior_ENS1_2020020711.relative.nc
2026-07-02 14:13:20 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:13:20 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:13:20 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:13:20 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:13:21 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:13:21 INFO [Posterior Diff ENS2] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020711/diff_posterior_ENS2_2020020711.nc
2026-07-02 14:13:21 INFO [Posterior Diff ENS2] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020711/diff_posterior_ENS2_2020020711.relative.nc
2026-07-02 14:13:22 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:13:22 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:13:22 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:13:22 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:13:23 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:13:23 INFO [Posterior Diff ENS3] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020711/diff_posterior_ENS3_2020020711.nc
2026-07-02 14:13:23 INFO [Posterior Diff ENS3] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020711/diff_posterior_ENS3_2020020711.relative.nc
2026-07-02 14:13:23 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:13:23 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:13:24 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:13:24 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:13:25 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:13:25 INFO [Posterior Diff ENS4] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020711/diff_posterior_ENS4_2020020711.nc
2026-07-02 14:13:25 INFO [Posterior Diff ENS4] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020711/diff_posterior_ENS4_2020020711.relative.nc
2026-07-02 14:13:26 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:13:26 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:13:26 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:13:26 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:13:27 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:13:27 INFO [Posterior Diff ENS5] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020711/diff_posterior_ENS5_2020020711.nc
2026-07-02 14:13:27 INFO [Posterior Diff ENS5] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020711/diff_posterior_ENS5_2020020711.relative.nc
2026-07-02 14:13:27 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:13:27 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:13:28 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:13:28 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:13:29 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:13:29 INFO [Posterior Diff ENS6] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020711/diff_posterior_ENS6_2020020711.nc
2026-07-02 14:13:29 INFO [Posterior Diff ENS6] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020711/diff_posterior_ENS6_2020020711.relative.nc
2026-07-02 14:13:30 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:13:30 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:13:30 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:13:30 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:13:31 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:13:31 INFO [Posterior Diff ENS7] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020711/diff_posterior_ENS7_2020020711.nc
2026-07-02 14:13:31 INFO [Posterior Diff ENS7] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020711/diff_posterior_ENS7_2020020711.relative.nc
2026-07-02 14:13:32 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:13:32 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:13:32 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:13:32 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:13:33 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:13:33 INFO [Posterior Diff ENS8] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020711/diff_posterior_ENS8_2020020711.nc
2026-07-02 14:13:33 INFO [Posterior Diff ENS8] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020711/diff_posterior_ENS8_2020020711.relative.nc
2026-07-02 14:13:34 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:13:34 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:13:34 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:13:34 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:13:35 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:13:35 INFO [Posterior Diff ENS9] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020711/diff_posterior_ENS9_2020020711.nc
2026-07-02 14:13:35 INFO [Posterior Diff ENS9] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020711/diff_posterior_ENS9_2020020711.relative.nc
2026-07-02 14:13:36 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:13:36 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:13:36 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:13:36 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:13:37 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:13:37 INFO [Posterior Diff ENS10] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020711/diff_posterior_ENS10_2020020711.nc
2026-07-02 14:13:37 INFO [Posterior Diff ENS10] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020711/diff_posterior_ENS10_2020020711.relative.nc
2026-07-02 14:13:38 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:13:38 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:13:38 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:13:38 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:13:39 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:13:39 INFO [Posterior Diff ENS11] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020711/diff_posterior_ENS11_2020020711.nc
2026-07-02 14:13:39 INFO [Posterior Diff ENS11] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020711/diff_posterior_ENS11_2020020711.relative.nc
2026-07-02 14:13:40 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:13:40 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:13:40 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:13:40 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:13:41 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:13:41 INFO [Posterior Diff ENS12] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020711/diff_posterior_ENS12_2020020711.nc
2026-07-02 14:13:41 INFO [Posterior Diff ENS12] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020711/diff_posterior_ENS12_2020020711.relative.nc
2026-07-02 14:13:42 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:13:42 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:13:42 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:13:42 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:13:43 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:13:43 INFO [Posterior Diff ENS13] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020711/diff_posterior_ENS13_2020020711.nc
2026-07-02 14:13:43 INFO [Posterior Diff ENS13] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020711/diff_posterior_ENS13_2020020711.relative.nc
2026-07-02 14:13:44 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:13:44 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:13:45 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:13:45 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:13:45 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:13:46 INFO [Posterior Diff ENS14] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020711/diff_posterior_ENS14_2020020711.nc
2026-07-02 14:13:46 INFO [Posterior Diff ENS14] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020711/diff_posterior_ENS14_2020020711.relative.nc
2026-07-02 14:13:46 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:13:46 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:13:47 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:13:47 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:13:47 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:13:48 INFO [Posterior Diff ENS15] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020711/diff_posterior_ENS15_2020020711.nc
2026-07-02 14:13:48 INFO [Posterior Diff ENS15] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020711/diff_posterior_ENS15_2020020711.relative.nc
2026-07-02 14:13:48 INFO Next run starts from 2020-02-07 11:00:00
2026-07-02 14:13:48 INFO Cycle is DONE; starting a new loop!
2026-07-02 14:13:48 INFO [TIME] step_end current_time=2020-02-07 11:00:00 simulated_time=2020-02-07 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 14:13:48 INFO [TIME] step_start current_time=2020-02-07 11:00:00 simulated_time=2020-02-07 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 14:13:48 INFO [TIME] window start=2020-02-07 11:00:00 end=2020-02-07 12:00:00 run_hours=1 has_assimilation=True
2026-07-02 14:13:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:13:49 INFO Hourly dataset computed and listing created
2026-07-02 14:13:50 INFO Hourly dataset computed
2026-07-02 14:13:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:13:51 INFO Hourly dataset computed and listing created
2026-07-02 14:13:52 INFO Hourly dataset computed
2026-07-02 14:13:52 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:13:53 INFO Hourly dataset computed and listing created
2026-07-02 14:13:53 INFO Hourly dataset computed
2026-07-02 14:13:53 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:13:54 INFO Hourly dataset computed and listing created
2026-07-02 14:13:55 INFO Hourly dataset computed
2026-07-02 14:13:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:13:56 INFO Hourly dataset computed and listing created
2026-07-02 14:13:56 INFO Hourly dataset computed
2026-07-02 14:13:56 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:13:57 INFO Hourly dataset computed and listing created
2026-07-02 14:13:58 INFO Hourly dataset computed
2026-07-02 14:13:58 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:13:59 INFO Hourly dataset computed and listing created
2026-07-02 14:13:59 INFO Hourly dataset computed
2026-07-02 14:13:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:14:00 INFO Hourly dataset computed and listing created
2026-07-02 14:14:01 INFO Hourly dataset computed
2026-07-02 14:14:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:14:01 INFO Hourly dataset computed and listing created
2026-07-02 14:14:02 INFO Hourly dataset computed
2026-07-02 14:14:02 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:14:03 INFO Hourly dataset computed and listing created
2026-07-02 14:14:04 INFO Hourly dataset computed
2026-07-02 14:14:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:14:04 INFO Hourly dataset computed and listing created
2026-07-02 14:14:05 INFO Hourly dataset computed
2026-07-02 14:14:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:14:06 INFO Hourly dataset computed and listing created
2026-07-02 14:14:06 INFO Hourly dataset computed
2026-07-02 14:14:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:14:07 INFO Hourly dataset computed and listing created
2026-07-02 14:14:08 INFO Hourly dataset computed
2026-07-02 14:14:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:14:09 INFO Hourly dataset computed and listing created
2026-07-02 14:14:09 INFO Hourly dataset computed
2026-07-02 14:14:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:14:10 INFO Hourly dataset computed and listing created
2026-07-02 14:14:11 INFO Hourly dataset computed
2026-07-02 14:14:11 INFO ---------->>> Running CHIMERE model from 2020-02-07 11:00:00 to 2020-02-07 12:00:00
2026-07-02 14:14:11 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:14:11 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1
2026-07-02 14:14:11 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020709_2_ENS1.nc
2026-07-02 14:14:11 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-02 14:14:11 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:14:11 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-02 14:14:11 INFO Queuing job for member 1...
2026-07-02 14:14:11 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:14:11 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-02 14:14:12 INFO Found: ['5066400']
2026-07-02 14:14:17 INFO [TGCC-IRENE] Submitted job with ID:['5066400']
2026-07-02 14:14:17 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:14:17 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2
2026-07-02 14:14:17 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020709_2_ENS2.nc
2026-07-02 14:14:17 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-02 14:14:17 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:14:17 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-02 14:14:17 INFO Queuing job for member 2...
2026-07-02 14:14:17 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:14:17 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-02 14:14:17 INFO Found: ['5066401']
2026-07-02 14:14:22 INFO [TGCC-IRENE] Submitted job with ID:['5066401']
2026-07-02 14:14:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:14:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3
2026-07-02 14:14:22 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020709_2_ENS3.nc
2026-07-02 14:14:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-02 14:14:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:14:22 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-02 14:14:22 INFO Queuing job for member 3...
2026-07-02 14:14:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:14:22 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-02 14:14:23 INFO Found: ['5066402']
2026-07-02 14:14:28 INFO [TGCC-IRENE] Submitted job with ID:['5066402']
2026-07-02 14:14:28 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:14:28 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4
2026-07-02 14:14:28 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020709_2_ENS4.nc
2026-07-02 14:14:28 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-02 14:14:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:14:28 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-02 14:14:28 INFO Queuing job for member 4...
2026-07-02 14:14:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:14:28 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-02 14:14:29 INFO Found: ['5066403']
2026-07-02 14:14:34 INFO [TGCC-IRENE] Submitted job with ID:['5066403']
2026-07-02 14:14:34 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:14:34 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5
2026-07-02 14:14:34 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020709_2_ENS5.nc
2026-07-02 14:14:34 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-02 14:14:34 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:14:34 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-02 14:14:34 INFO Queuing job for member 5...
2026-07-02 14:14:34 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:14:34 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-02 14:14:34 INFO Found: ['5066405']
2026-07-02 14:14:39 INFO [TGCC-IRENE] Submitted job with ID:['5066405']
2026-07-02 14:14:39 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:14:39 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6
2026-07-02 14:14:39 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020709_2_ENS6.nc
2026-07-02 14:14:39 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-02 14:14:39 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:14:39 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-02 14:14:39 INFO Queuing job for member 6...
2026-07-02 14:14:39 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:14:39 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-02 14:14:40 INFO Found: ['5066406']
2026-07-02 14:14:45 INFO [TGCC-IRENE] Submitted job with ID:['5066406']
2026-07-02 14:14:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:14:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7
2026-07-02 14:14:45 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020709_2_ENS7.nc
2026-07-02 14:14:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-02 14:14:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:14:45 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-02 14:14:45 INFO Queuing job for member 7...
2026-07-02 14:14:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:14:45 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-02 14:14:46 INFO Found: ['5066407']
2026-07-02 14:14:51 INFO [TGCC-IRENE] Submitted job with ID:['5066407']
2026-07-02 14:14:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:14:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8
2026-07-02 14:14:51 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020709_2_ENS8.nc
2026-07-02 14:14:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-02 14:14:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:14:51 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-02 14:14:51 INFO Queuing job for member 8...
2026-07-02 14:14:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:14:51 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-02 14:14:52 INFO Found: ['5066408']
2026-07-02 14:14:57 INFO [TGCC-IRENE] Submitted job with ID:['5066408']
2026-07-02 14:14:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:14:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9
2026-07-02 14:14:57 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020709_2_ENS9.nc
2026-07-02 14:14:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-02 14:14:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:14:57 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-02 14:14:57 INFO Queuing job for member 9...
2026-07-02 14:14:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:14:57 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-02 14:14:57 INFO Found: ['5066409']
2026-07-02 14:15:02 INFO [TGCC-IRENE] Submitted job with ID:['5066409']
2026-07-02 14:15:02 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:15:02 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10
2026-07-02 14:15:02 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020709_2_ENS10.nc
2026-07-02 14:15:02 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-02 14:15:02 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:15:02 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-02 14:15:02 INFO Queuing job for member 10...
2026-07-02 14:15:02 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:15:02 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-02 14:15:03 INFO Found: ['5066411']
2026-07-02 14:15:08 INFO [TGCC-IRENE] Submitted job with ID:['5066411']
2026-07-02 14:15:08 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:15:08 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11
2026-07-02 14:15:08 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020709_2_ENS11.nc
2026-07-02 14:15:08 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-02 14:15:08 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:15:08 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-02 14:15:08 INFO Queuing job for member 11...
2026-07-02 14:15:08 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:15:08 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-02 14:15:09 INFO Found: ['5066412']
2026-07-02 14:15:14 INFO [TGCC-IRENE] Submitted job with ID:['5066412']
2026-07-02 14:15:14 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:15:14 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12
2026-07-02 14:15:14 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020709_2_ENS12.nc
2026-07-02 14:15:14 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-02 14:15:14 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:15:14 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-02 14:15:14 INFO Queuing job for member 12...
2026-07-02 14:15:14 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:15:14 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-02 14:15:15 INFO Found: ['5066413']
2026-07-02 14:15:20 INFO [TGCC-IRENE] Submitted job with ID:['5066413']
2026-07-02 14:15:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:15:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13
2026-07-02 14:15:20 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020709_2_ENS13.nc
2026-07-02 14:15:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-02 14:15:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:15:20 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-02 14:15:20 INFO Queuing job for member 13...
2026-07-02 14:15:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:15:20 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-02 14:15:21 INFO Found: ['5066414']
2026-07-02 14:15:26 INFO [TGCC-IRENE] Submitted job with ID:['5066414']
2026-07-02 14:15:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:15:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14
2026-07-02 14:15:26 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020709_2_ENS14.nc
2026-07-02 14:15:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-02 14:15:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:15:26 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-02 14:15:26 INFO Queuing job for member 14...
2026-07-02 14:15:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:15:26 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-02 14:15:27 INFO Found: ['5066415']
2026-07-02 14:15:32 INFO [TGCC-IRENE] Submitted job with ID:['5066415']
2026-07-02 14:15:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:15:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15
2026-07-02 14:15:32 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020709_2_ENS15.nc
2026-07-02 14:15:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-02 14:15:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:15:32 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-02 14:15:32 INFO Queuing job for member 15...
2026-07-02 14:15:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:15:32 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-02 14:15:33 INFO Found: ['5066416']
2026-07-02 14:15:38 INFO [TGCC-IRENE] Submitted job with ID:['5066416']
2026-07-02 14:15:38 INFO Checking job status ...
2026-07-02 14:15:38 INFO None 5066400: status RUNNING/PENDING
2026-07-02 14:15:38 INFO None 5066401: status RUNNING/PENDING
2026-07-02 14:15:38 INFO None 5066402: status RUNNING/PENDING
2026-07-02 14:15:38 INFO None 5066403: status RUNNING/PENDING
2026-07-02 14:15:38 INFO None 5066405: status RUNNING/PENDING
2026-07-02 14:15:38 INFO None 5066406: status RUNNING/PENDING
2026-07-02 14:15:38 INFO None 5066407: status RUNNING/PENDING
2026-07-02 14:15:38 INFO None 5066408: status RUNNING/PENDING
2026-07-02 14:15:38 INFO None 5066409: status RUNNING/PENDING
2026-07-02 14:15:38 INFO None 5066411: status RUNNING/PENDING
2026-07-02 14:15:38 INFO None 5066412: status RUNNING/PENDING
2026-07-02 14:15:38 INFO None 5066413: status RUNNING/PENDING
2026-07-02 14:15:38 INFO None 5066414: status RUNNING/PENDING
2026-07-02 14:15:38 INFO None 5066415: status RUNNING/PENDING
2026-07-02 14:15:38 INFO None 5066416: status RUNNING/PENDING
2026-07-02 14:15:38 INFO Jobs still running: ['5066400', '5066401', '5066402', '5066403', '5066405', '5066406', '5066407', '5066408', '5066409', '5066411', '5066412', '5066413', '5066414', '5066415', '5066416']. Waiting...
2026-07-02 14:15:53 INFO None 5066400: status RUNNING/PENDING
2026-07-02 14:15:53 INFO None 5066401: status RUNNING/PENDING
2026-07-02 14:15:53 INFO None 5066402: status RUNNING/PENDING
2026-07-02 14:15:53 INFO None 5066403: status RUNNING/PENDING
2026-07-02 14:15:53 INFO None 5066405: status RUNNING/PENDING
2026-07-02 14:15:53 INFO None 5066406: status RUNNING/PENDING
2026-07-02 14:15:53 INFO None 5066407: status RUNNING/PENDING
2026-07-02 14:15:53 INFO None 5066408: status RUNNING/PENDING
2026-07-02 14:15:53 INFO None 5066409: status RUNNING/PENDING
2026-07-02 14:15:53 INFO None 5066411: status RUNNING/PENDING
2026-07-02 14:15:53 INFO None 5066412: status RUNNING/PENDING
2026-07-02 14:15:53 INFO None 5066413: status RUNNING/PENDING
2026-07-02 14:15:53 INFO None 5066414: status RUNNING/PENDING
2026-07-02 14:15:53 INFO None 5066415: status RUNNING/PENDING
2026-07-02 14:15:53 INFO None 5066416: status RUNNING/PENDING
2026-07-02 14:15:53 INFO Jobs still running: ['5066400', '5066401', '5066402', '5066403', '5066405', '5066406', '5066407', '5066408', '5066409', '5066411', '5066412', '5066413', '5066414', '5066415', '5066416']. Waiting...
2026-07-02 14:16:08 INFO None 5066400: status FINISHED
2026-07-02 14:16:08 INFO None 5066401: status RUNNING/PENDING
2026-07-02 14:16:08 INFO None 5066402: status RUNNING/PENDING
2026-07-02 14:16:08 INFO None 5066403: status RUNNING/PENDING
2026-07-02 14:16:08 INFO None 5066405: status RUNNING/PENDING
2026-07-02 14:16:08 INFO None 5066406: status RUNNING/PENDING
2026-07-02 14:16:08 INFO None 5066407: status RUNNING/PENDING
2026-07-02 14:16:08 INFO None 5066408: status RUNNING/PENDING
2026-07-02 14:16:08 INFO None 5066409: status RUNNING/PENDING
2026-07-02 14:16:08 INFO None 5066411: status RUNNING/PENDING
2026-07-02 14:16:08 INFO None 5066412: status RUNNING/PENDING
2026-07-02 14:16:08 INFO None 5066413: status RUNNING/PENDING
2026-07-02 14:16:08 INFO None 5066414: status RUNNING/PENDING
2026-07-02 14:16:08 INFO None 5066415: status RUNNING/PENDING
2026-07-02 14:16:08 INFO None 5066416: status RUNNING/PENDING
2026-07-02 14:16:08 INFO Jobs still running: ['5066401', '5066402', '5066403', '5066405', '5066406', '5066407', '5066408', '5066409', '5066411', '5066412', '5066413', '5066414', '5066415', '5066416']. Waiting...
2026-07-02 14:16:24 INFO None 5066400: status FINISHED
2026-07-02 14:16:24 INFO None 5066401: status RUNNING/PENDING
2026-07-02 14:16:24 INFO None 5066402: status RUNNING/PENDING
2026-07-02 14:16:24 INFO None 5066403: status RUNNING/PENDING
2026-07-02 14:16:24 INFO None 5066405: status RUNNING/PENDING
2026-07-02 14:16:24 INFO None 5066406: status RUNNING/PENDING
2026-07-02 14:16:24 INFO None 5066407: status RUNNING/PENDING
2026-07-02 14:16:24 INFO None 5066408: status RUNNING/PENDING
2026-07-02 14:16:24 INFO None 5066409: status RUNNING/PENDING
2026-07-02 14:16:24 INFO None 5066411: status RUNNING/PENDING
2026-07-02 14:16:24 INFO None 5066412: status RUNNING/PENDING
2026-07-02 14:16:24 INFO None 5066413: status RUNNING/PENDING
2026-07-02 14:16:24 INFO None 5066414: status RUNNING/PENDING
2026-07-02 14:16:24 INFO None 5066415: status RUNNING/PENDING
2026-07-02 14:16:24 INFO None 5066416: status RUNNING/PENDING
2026-07-02 14:16:24 INFO Jobs still running: ['5066401', '5066402', '5066403', '5066405', '5066406', '5066407', '5066408', '5066409', '5066411', '5066412', '5066413', '5066414', '5066415', '5066416']. Waiting...
2026-07-02 14:16:39 INFO None 5066400: status FINISHED
2026-07-02 14:16:39 INFO None 5066401: status RUNNING/PENDING
2026-07-02 14:16:39 INFO None 5066402: status RUNNING/PENDING
2026-07-02 14:16:39 INFO None 5066403: status RUNNING/PENDING
2026-07-02 14:16:39 INFO None 5066405: status RUNNING/PENDING
2026-07-02 14:16:39 INFO None 5066406: status RUNNING/PENDING
2026-07-02 14:16:39 INFO None 5066407: status RUNNING/PENDING
2026-07-02 14:16:39 INFO None 5066408: status RUNNING/PENDING
2026-07-02 14:16:39 INFO None 5066409: status RUNNING/PENDING
2026-07-02 14:16:39 INFO None 5066411: status RUNNING/PENDING
2026-07-02 14:16:39 INFO None 5066412: status RUNNING/PENDING
2026-07-02 14:16:39 INFO None 5066413: status RUNNING/PENDING
2026-07-02 14:16:39 INFO None 5066414: status RUNNING/PENDING
2026-07-02 14:16:39 INFO None 5066415: status RUNNING/PENDING
2026-07-02 14:16:39 INFO None 5066416: status RUNNING/PENDING
2026-07-02 14:16:39 INFO Jobs still running: ['5066401', '5066402', '5066403', '5066405', '5066406', '5066407', '5066408', '5066409', '5066411', '5066412', '5066413', '5066414', '5066415', '5066416']. Waiting...
2026-07-02 14:16:55 INFO None 5066400: status FINISHED
2026-07-02 14:16:55 INFO None 5066401: status RUNNING/PENDING
2026-07-02 14:16:55 INFO None 5066402: status FINISHED
2026-07-02 14:16:55 INFO None 5066403: status FINISHED
2026-07-02 14:16:55 INFO None 5066405: status RUNNING/PENDING
2026-07-02 14:16:55 INFO None 5066406: status FINISHED
2026-07-02 14:16:55 INFO None 5066407: status FINISHED
2026-07-02 14:16:55 INFO None 5066408: status RUNNING/PENDING
2026-07-02 14:16:55 INFO None 5066409: status RUNNING/PENDING
2026-07-02 14:16:55 INFO None 5066411: status RUNNING/PENDING
2026-07-02 14:16:55 INFO None 5066412: status RUNNING/PENDING
2026-07-02 14:16:55 INFO None 5066413: status RUNNING/PENDING
2026-07-02 14:16:55 INFO None 5066414: status RUNNING/PENDING
2026-07-02 14:16:55 INFO None 5066415: status RUNNING/PENDING
2026-07-02 14:16:55 INFO None 5066416: status RUNNING/PENDING
2026-07-02 14:16:55 INFO Jobs still running: ['5066401', '5066405', '5066408', '5066409', '5066411', '5066412', '5066413', '5066414', '5066415', '5066416']. Waiting...
2026-07-02 14:17:10 INFO None 5066400: status FINISHED
2026-07-02 14:17:10 INFO None 5066401: status FINISHED
2026-07-02 14:17:10 INFO None 5066402: status FINISHED
2026-07-02 14:17:10 INFO None 5066403: status FINISHED
2026-07-02 14:17:10 INFO None 5066405: status RUNNING/PENDING
2026-07-02 14:17:10 INFO None 5066406: status FINISHED
2026-07-02 14:17:10 INFO None 5066407: status FINISHED
2026-07-02 14:17:10 INFO None 5066408: status RUNNING/PENDING
2026-07-02 14:17:10 INFO None 5066409: status RUNNING/PENDING
2026-07-02 14:17:10 INFO None 5066411: status RUNNING/PENDING
2026-07-02 14:17:10 INFO None 5066412: status RUNNING/PENDING
2026-07-02 14:17:10 INFO None 5066413: status RUNNING/PENDING
2026-07-02 14:17:10 INFO None 5066414: status RUNNING/PENDING
2026-07-02 14:17:10 INFO None 5066415: status RUNNING/PENDING
2026-07-02 14:17:10 INFO None 5066416: status RUNNING/PENDING
2026-07-02 14:17:10 INFO Jobs still running: ['5066405', '5066408', '5066409', '5066411', '5066412', '5066413', '5066414', '5066415', '5066416']. Waiting...
2026-07-02 14:17:25 INFO None 5066400: status FINISHED
2026-07-02 14:17:26 INFO None 5066401: status FINISHED
2026-07-02 14:17:26 INFO None 5066402: status FINISHED
2026-07-02 14:17:26 INFO None 5066403: status FINISHED
2026-07-02 14:17:26 INFO None 5066405: status RUNNING/PENDING
2026-07-02 14:17:26 INFO None 5066406: status FINISHED
2026-07-02 14:17:26 INFO None 5066407: status FINISHED
2026-07-02 14:17:26 INFO None 5066408: status FINISHED
2026-07-02 14:17:26 INFO None 5066409: status FINISHED
2026-07-02 14:17:26 INFO None 5066411: status FINISHED
2026-07-02 14:17:26 INFO None 5066412: status RUNNING/PENDING
2026-07-02 14:17:26 INFO None 5066413: status RUNNING/PENDING
2026-07-02 14:17:26 INFO None 5066414: status RUNNING/PENDING
2026-07-02 14:17:26 INFO None 5066415: status RUNNING/PENDING
2026-07-02 14:17:26 INFO None 5066416: status RUNNING/PENDING
2026-07-02 14:17:26 INFO Jobs still running: ['5066405', '5066412', '5066413', '5066414', '5066415', '5066416']. Waiting...
2026-07-02 14:17:41 INFO None 5066400: status FINISHED
2026-07-02 14:17:41 INFO None 5066401: status FINISHED
2026-07-02 14:17:41 INFO None 5066402: status FINISHED
2026-07-02 14:17:41 INFO None 5066403: status FINISHED
2026-07-02 14:17:41 INFO None 5066405: status RUNNING/PENDING
2026-07-02 14:17:41 INFO None 5066406: status FINISHED
2026-07-02 14:17:41 INFO None 5066407: status FINISHED
2026-07-02 14:17:41 INFO None 5066408: status FINISHED
2026-07-02 14:17:41 INFO None 5066409: status FINISHED
2026-07-02 14:17:41 INFO None 5066411: status FINISHED
2026-07-02 14:17:41 INFO None 5066412: status RUNNING/PENDING
2026-07-02 14:17:41 INFO None 5066413: status RUNNING/PENDING
2026-07-02 14:17:41 INFO None 5066414: status RUNNING/PENDING
2026-07-02 14:17:41 INFO None 5066415: status RUNNING/PENDING
2026-07-02 14:17:41 INFO None 5066416: status RUNNING/PENDING
2026-07-02 14:17:41 INFO Jobs still running: ['5066405', '5066412', '5066413', '5066414', '5066415', '5066416']. Waiting...
2026-07-02 14:17:56 INFO None 5066400: status FINISHED
2026-07-02 14:17:56 INFO None 5066401: status FINISHED
2026-07-02 14:17:56 INFO None 5066402: status FINISHED
2026-07-02 14:17:56 INFO None 5066403: status FINISHED
2026-07-02 14:17:56 INFO None 5066405: status RUNNING/PENDING
2026-07-02 14:17:56 INFO None 5066406: status FINISHED
2026-07-02 14:17:56 INFO None 5066407: status FINISHED
2026-07-02 14:17:56 INFO None 5066408: status FINISHED
2026-07-02 14:17:56 INFO None 5066409: status FINISHED
2026-07-02 14:17:56 INFO None 5066411: status FINISHED
2026-07-02 14:17:56 INFO None 5066412: status RUNNING/PENDING
2026-07-02 14:17:56 INFO None 5066413: status RUNNING/PENDING
2026-07-02 14:17:56 INFO None 5066414: status RUNNING/PENDING
2026-07-02 14:17:56 INFO None 5066415: status RUNNING/PENDING
2026-07-02 14:17:56 INFO None 5066416: status RUNNING/PENDING
2026-07-02 14:17:56 INFO Jobs still running: ['5066405', '5066412', '5066413', '5066414', '5066415', '5066416']. Waiting...
2026-07-02 14:18:11 INFO None 5066400: status FINISHED
2026-07-02 14:18:11 INFO None 5066401: status FINISHED
2026-07-02 14:18:11 INFO None 5066402: status FINISHED
2026-07-02 14:18:11 INFO None 5066403: status FINISHED
2026-07-02 14:18:11 INFO None 5066405: status FINISHED
2026-07-02 14:18:11 INFO None 5066406: status FINISHED
2026-07-02 14:18:11 INFO None 5066407: status FINISHED
2026-07-02 14:18:11 INFO None 5066408: status FINISHED
2026-07-02 14:18:11 INFO None 5066409: status FINISHED
2026-07-02 14:18:12 INFO None 5066411: status FINISHED
2026-07-02 14:18:12 INFO None 5066412: status RUNNING/PENDING
2026-07-02 14:18:12 INFO None 5066413: status RUNNING/PENDING
2026-07-02 14:18:12 INFO None 5066414: status RUNNING/PENDING
2026-07-02 14:18:12 INFO None 5066415: status FINISHED
2026-07-02 14:18:12 INFO None 5066416: status FINISHED
2026-07-02 14:18:12 INFO Jobs still running: ['5066412', '5066413', '5066414']. Waiting...
2026-07-02 14:18:27 INFO None 5066400: status FINISHED
2026-07-02 14:18:27 INFO None 5066401: status FINISHED
2026-07-02 14:18:27 INFO None 5066402: status FINISHED
2026-07-02 14:18:27 INFO None 5066403: status FINISHED
2026-07-02 14:18:27 INFO None 5066405: status FINISHED
2026-07-02 14:18:27 INFO None 5066406: status FINISHED
2026-07-02 14:18:27 INFO None 5066407: status FINISHED
2026-07-02 14:18:27 INFO None 5066408: status FINISHED
2026-07-02 14:18:27 INFO None 5066409: status FINISHED
2026-07-02 14:18:27 INFO None 5066411: status FINISHED
2026-07-02 14:18:27 INFO None 5066412: status FINISHED
2026-07-02 14:18:27 INFO None 5066413: status RUNNING/PENDING
2026-07-02 14:18:27 INFO None 5066414: status RUNNING/PENDING
2026-07-02 14:18:27 INFO None 5066415: status FINISHED
2026-07-02 14:18:27 INFO None 5066416: status FINISHED
2026-07-02 14:18:27 INFO Jobs still running: ['5066413', '5066414']. Waiting...
2026-07-02 14:18:42 INFO None 5066400: status FINISHED
2026-07-02 14:18:42 INFO None 5066401: status FINISHED
2026-07-02 14:18:42 INFO None 5066402: status FINISHED
2026-07-02 14:18:42 INFO None 5066403: status FINISHED
2026-07-02 14:18:42 INFO None 5066405: status FINISHED
2026-07-02 14:18:42 INFO None 5066406: status FINISHED
2026-07-02 14:18:42 INFO None 5066407: status FINISHED
2026-07-02 14:18:42 INFO None 5066408: status FINISHED
2026-07-02 14:18:43 INFO None 5066409: status FINISHED
2026-07-02 14:18:43 INFO None 5066411: status FINISHED
2026-07-02 14:18:43 INFO None 5066412: status FINISHED
2026-07-02 14:18:43 INFO None 5066413: status RUNNING/PENDING
2026-07-02 14:18:43 INFO None 5066414: status RUNNING/PENDING
2026-07-02 14:18:43 INFO None 5066415: status FINISHED
2026-07-02 14:18:43 INFO None 5066416: status FINISHED
2026-07-02 14:18:43 INFO Jobs still running: ['5066413', '5066414']. Waiting...
2026-07-02 14:18:58 INFO None 5066400: status FINISHED
2026-07-02 14:18:58 INFO None 5066401: status FINISHED
2026-07-02 14:18:58 INFO None 5066402: status FINISHED
2026-07-02 14:18:58 INFO None 5066403: status FINISHED
2026-07-02 14:18:58 INFO None 5066405: status FINISHED
2026-07-02 14:18:58 INFO None 5066406: status FINISHED
2026-07-02 14:18:58 INFO None 5066407: status FINISHED
2026-07-02 14:18:58 INFO None 5066408: status FINISHED
2026-07-02 14:18:58 INFO None 5066409: status FINISHED
2026-07-02 14:18:58 INFO None 5066411: status FINISHED
2026-07-02 14:18:58 INFO None 5066412: status FINISHED
2026-07-02 14:18:58 INFO None 5066413: status FINISHED
2026-07-02 14:18:58 INFO None 5066414: status RUNNING/PENDING
2026-07-02 14:18:58 INFO None 5066415: status FINISHED
2026-07-02 14:18:58 INFO None 5066416: status FINISHED
2026-07-02 14:18:58 INFO Jobs still running: ['5066414']. Waiting...
2026-07-02 14:19:13 INFO None 5066400: status FINISHED
2026-07-02 14:19:13 INFO None 5066401: status FINISHED
2026-07-02 14:19:13 INFO None 5066402: status FINISHED
2026-07-02 14:19:13 INFO None 5066403: status FINISHED
2026-07-02 14:19:13 INFO None 5066405: status FINISHED
2026-07-02 14:19:13 INFO None 5066406: status FINISHED
2026-07-02 14:19:13 INFO None 5066407: status FINISHED
2026-07-02 14:19:13 INFO None 5066408: status FINISHED
2026-07-02 14:19:13 INFO None 5066409: status FINISHED
2026-07-02 14:19:13 INFO None 5066411: status FINISHED
2026-07-02 14:19:13 INFO None 5066412: status FINISHED
2026-07-02 14:19:13 INFO None 5066413: status FINISHED
2026-07-02 14:19:13 INFO None 5066414: status RUNNING/PENDING
2026-07-02 14:19:13 INFO None 5066415: status FINISHED
2026-07-02 14:19:13 INFO None 5066416: status FINISHED
2026-07-02 14:19:13 INFO Jobs still running: ['5066414']. Waiting...
2026-07-02 14:19:28 INFO None 5066400: status FINISHED
2026-07-02 14:19:28 INFO None 5066401: status FINISHED
2026-07-02 14:19:28 INFO None 5066402: status FINISHED
2026-07-02 14:19:28 INFO None 5066403: status FINISHED
2026-07-02 14:19:28 INFO None 5066405: status FINISHED
2026-07-02 14:19:28 INFO None 5066406: status FINISHED
2026-07-02 14:19:28 INFO None 5066407: status FINISHED
2026-07-02 14:19:28 INFO None 5066408: status FINISHED
2026-07-02 14:19:28 INFO None 5066409: status FINISHED
2026-07-02 14:19:28 INFO None 5066411: status FINISHED
2026-07-02 14:19:28 INFO None 5066412: status FINISHED
2026-07-02 14:19:28 INFO None 5066413: status FINISHED
2026-07-02 14:19:28 INFO None 5066414: status FINISHED
2026-07-02 14:19:28 INFO None 5066415: status FINISHED
2026-07-02 14:19:28 INFO None 5066416: status FINISHED
2026-07-02 14:19:28 INFO Jobs ['5066400', '5066401', '5066402', '5066403', '5066405', '5066406', '5066407', '5066408', '5066409', '5066411', '5066412', '5066413', '5066414', '5066415', '5066416'] have finished
2026-07-02 14:19:28 INFO Checking restart files were created ...
2026-07-02 14:19:28 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020711_1_ENS1.nc(668832435 bytes)
2026-07-02 14:19:28 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020711_1_ENS2.nc(668832435 bytes)
2026-07-02 14:19:28 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020711_1_ENS3.nc(668832435 bytes)
2026-07-02 14:19:28 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020711_1_ENS4.nc(668832435 bytes)
2026-07-02 14:19:28 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020711_1_ENS5.nc(668832435 bytes)
2026-07-02 14:19:28 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020711_1_ENS6.nc(668832435 bytes)
2026-07-02 14:19:28 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020711_1_ENS7.nc(668832435 bytes)
2026-07-02 14:19:28 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020711_1_ENS8.nc(668832435 bytes)
2026-07-02 14:19:28 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020711_1_ENS9.nc(668832435 bytes)
2026-07-02 14:19:28 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020711_1_ENS10.nc(668832435 bytes)
2026-07-02 14:19:28 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020711_1_ENS11.nc(668832435 bytes)
2026-07-02 14:19:28 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020711_1_ENS12.nc(668832435 bytes)
2026-07-02 14:19:28 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020711_1_ENS13.nc(668832435 bytes)
2026-07-02 14:19:28 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020711_1_ENS14.nc(668832435 bytes)
2026-07-02 14:19:28 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020711_1_ENS15.nc(668832435 bytes)
2026-07-02 14:19:28 INFO  Run_model() completed successfully.
2026-07-02 14:19:28 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 11:00:00 simulated_time=2020-02-07 12:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 14:19:28 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 12:00:00 days=153073 seconds=43200
2026-07-02 14:19:28 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-02 14:19:28 INFO [TIME] increment current_time 2020-02-07 11:00:00 -> 2020-02-07 12:00:00
2026-07-02 14:19:28 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 12:00:00 simulated_time=2020-02-07 12:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 14:19:28 INFO ---------->>> Running process_satellite_data()
2026-07-02 14:19:29 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12018.nc
2026-07-02 14:19:29 INFO ---------->>> Running run_obs_converter()
2026-07-02 14:19:29 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_43993_153073.out
2026-07-02 14:19:29 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_43993_153073.out
2026-07-02 14:19:29 INFO ---------->>> Running DART
2026-07-02 14:19:29 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-02 14:19:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/chim_ENS1_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/chim_ENS1_2020020712_1_out_toDART.nc
2026-07-02 14:19:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/chim_ENS2_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/chim_ENS2_2020020712_1_out_toDART.nc
2026-07-02 14:19:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/chim_ENS3_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/chim_ENS3_2020020712_1_out_toDART.nc
2026-07-02 14:19:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/chim_ENS4_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/chim_ENS4_2020020712_1_out_toDART.nc
2026-07-02 14:19:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/chim_ENS5_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/chim_ENS5_2020020712_1_out_toDART.nc
2026-07-02 14:19:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/chim_ENS6_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/chim_ENS6_2020020712_1_out_toDART.nc
2026-07-02 14:19:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/chim_ENS7_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/chim_ENS7_2020020712_1_out_toDART.nc
2026-07-02 14:19:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/chim_ENS8_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/chim_ENS8_2020020712_1_out_toDART.nc
2026-07-02 14:19:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/chim_ENS9_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/chim_ENS9_2020020712_1_out_toDART.nc
2026-07-02 14:19:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/chim_ENS10_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/chim_ENS10_2020020712_1_out_toDART.nc
2026-07-02 14:19:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/chim_ENS11_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/chim_ENS11_2020020712_1_out_toDART.nc
2026-07-02 14:19:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/chim_ENS12_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/chim_ENS12_2020020712_1_out_toDART.nc
2026-07-02 14:19:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/chim_ENS13_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/chim_ENS13_2020020712_1_out_toDART.nc
2026-07-02 14:19:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/chim_ENS14_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/chim_ENS14_2020020712_1_out_toDART.nc
2026-07-02 14:19:33 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/chim_ENS15_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/chim_ENS15_2020020712_1_out_toDART.nc
2026-07-02 14:19:33 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-02 14:19:33 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-02 14:19:33 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-02 14:19:33 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-02 14:19:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-02 14:19:33 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-02 14:19:47 INFO Found: []
2026-07-02 14:19:47 INFO No job id returned by command ./run_filter.bsh
2026-07-02 14:19:47 INFO No monitoring will be performed
2026-07-02 14:19:47 INFO Moving DART output files to analysis and preassim directories for date 2020020712 if present ...
2026-07-02 14:19:47 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020712'
2026-07-02 14:19:47 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020712'
2026-07-02 14:19:47 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020712'
2026-07-02 14:19:47 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020712'
2026-07-02 14:19:47 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020712'
2026-07-02 14:19:47 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020712'
2026-07-02 14:19:47 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020712'
2026-07-02 14:19:47 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020712'
2026-07-02 14:19:47 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020712'
2026-07-02 14:19:47 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020712'
2026-07-02 14:19:47 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020712'
2026-07-02 14:19:47 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020712'
2026-07-02 14:19:47 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020712'
2026-07-02 14:19:47 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020712'
2026-07-02 14:19:47 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020712'
2026-07-02 14:19:47 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020712'
2026-07-02 14:19:47 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020712'
2026-07-02 14:19:47 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020712'
2026-07-02 14:19:47 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020712'
2026-07-02 14:19:47 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020712'
2026-07-02 14:19:47 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020712'
2026-07-02 14:19:47 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020712'
2026-07-02 14:19:47 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020712'
2026-07-02 14:19:47 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020712'
2026-07-02 14:19:47 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020712'
2026-07-02 14:19:47 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020712'
2026-07-02 14:19:47 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020712'
2026-07-02 14:19:47 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020712'
2026-07-02 14:19:47 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020712'
2026-07-02 14:19:47 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020712'
2026-07-02 14:19:47 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020712'
2026-07-02 14:19:47 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020712'
2026-07-02 14:19:47 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020712'
2026-07-02 14:19:48 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020712'
2026-07-02 14:19:48 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-02 14:19:48 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-02 14:19:48 INFO run_dart() is DONE.
2026-07-02 14:19:48 INFO ---------->>> Running update_pollutant_in_end()
2026-07-02 14:19:48 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:19:48 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:19:49 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:19:49 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:19:49 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:19:50 INFO [Posterior Diff ENS1] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020712/diff_posterior_ENS1_2020020712.nc
2026-07-02 14:19:50 INFO [Posterior Diff ENS1] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020712/diff_posterior_ENS1_2020020712.relative.nc
2026-07-02 14:19:50 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:19:50 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:19:51 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:19:51 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:19:51 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:19:51 INFO [Posterior Diff ENS2] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020712/diff_posterior_ENS2_2020020712.nc
2026-07-02 14:19:51 INFO [Posterior Diff ENS2] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020712/diff_posterior_ENS2_2020020712.relative.nc
2026-07-02 14:19:52 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:19:52 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:19:53 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:19:53 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:19:53 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:19:53 INFO [Posterior Diff ENS3] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020712/diff_posterior_ENS3_2020020712.nc
2026-07-02 14:19:53 INFO [Posterior Diff ENS3] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020712/diff_posterior_ENS3_2020020712.relative.nc
2026-07-02 14:19:54 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:19:54 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:19:54 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:19:54 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:19:55 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:19:55 INFO [Posterior Diff ENS4] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020712/diff_posterior_ENS4_2020020712.nc
2026-07-02 14:19:55 INFO [Posterior Diff ENS4] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020712/diff_posterior_ENS4_2020020712.relative.nc
2026-07-02 14:19:56 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:19:56 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:19:57 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:19:57 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:19:57 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:19:57 INFO [Posterior Diff ENS5] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020712/diff_posterior_ENS5_2020020712.nc
2026-07-02 14:19:57 INFO [Posterior Diff ENS5] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020712/diff_posterior_ENS5_2020020712.relative.nc
2026-07-02 14:19:58 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:19:58 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:19:58 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:19:58 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:19:59 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:19:59 INFO [Posterior Diff ENS6] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020712/diff_posterior_ENS6_2020020712.nc
2026-07-02 14:19:59 INFO [Posterior Diff ENS6] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020712/diff_posterior_ENS6_2020020712.relative.nc
2026-07-02 14:20:00 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:20:00 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:20:01 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:20:01 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:20:01 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:20:01 INFO [Posterior Diff ENS7] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020712/diff_posterior_ENS7_2020020712.nc
2026-07-02 14:20:01 INFO [Posterior Diff ENS7] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020712/diff_posterior_ENS7_2020020712.relative.nc
2026-07-02 14:20:02 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:20:02 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:20:03 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:20:03 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:20:03 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:20:03 INFO [Posterior Diff ENS8] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020712/diff_posterior_ENS8_2020020712.nc
2026-07-02 14:20:03 INFO [Posterior Diff ENS8] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020712/diff_posterior_ENS8_2020020712.relative.nc
2026-07-02 14:20:04 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:20:04 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:20:05 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:20:05 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:20:05 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:20:05 INFO [Posterior Diff ENS9] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020712/diff_posterior_ENS9_2020020712.nc
2026-07-02 14:20:05 INFO [Posterior Diff ENS9] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020712/diff_posterior_ENS9_2020020712.relative.nc
2026-07-02 14:20:06 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:20:06 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:20:07 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:20:07 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:20:07 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:20:07 INFO [Posterior Diff ENS10] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020712/diff_posterior_ENS10_2020020712.nc
2026-07-02 14:20:07 INFO [Posterior Diff ENS10] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020712/diff_posterior_ENS10_2020020712.relative.nc
2026-07-02 14:20:08 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:20:08 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:20:09 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:20:09 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:20:09 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:20:09 INFO [Posterior Diff ENS11] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020712/diff_posterior_ENS11_2020020712.nc
2026-07-02 14:20:09 INFO [Posterior Diff ENS11] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020712/diff_posterior_ENS11_2020020712.relative.nc
2026-07-02 14:20:10 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:20:10 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:20:11 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:20:11 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:20:11 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:20:11 INFO [Posterior Diff ENS12] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020712/diff_posterior_ENS12_2020020712.nc
2026-07-02 14:20:11 INFO [Posterior Diff ENS12] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020712/diff_posterior_ENS12_2020020712.relative.nc
2026-07-02 14:20:12 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:20:12 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:20:13 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:20:13 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:20:13 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:20:13 INFO [Posterior Diff ENS13] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020712/diff_posterior_ENS13_2020020712.nc
2026-07-02 14:20:13 INFO [Posterior Diff ENS13] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020712/diff_posterior_ENS13_2020020712.relative.nc
2026-07-02 14:20:14 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:20:14 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:20:14 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:20:14 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:20:15 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:20:15 INFO [Posterior Diff ENS14] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020712/diff_posterior_ENS14_2020020712.nc
2026-07-02 14:20:15 INFO [Posterior Diff ENS14] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020712/diff_posterior_ENS14_2020020712.relative.nc
2026-07-02 14:20:16 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:20:16 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:20:17 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:20:17 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:20:17 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:20:17 INFO [Posterior Diff ENS15] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020712/diff_posterior_ENS15_2020020712.nc
2026-07-02 14:20:17 INFO [Posterior Diff ENS15] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020712/diff_posterior_ENS15_2020020712.relative.nc
2026-07-02 14:20:17 INFO Next run starts from 2020-02-07 12:00:00
2026-07-02 14:20:17 INFO Cycle is DONE; starting a new loop!
2026-07-02 14:20:17 INFO [TIME] step_end current_time=2020-02-07 12:00:00 simulated_time=2020-02-07 12:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 14:20:17 INFO [TIME] step_start current_time=2020-02-07 12:00:00 simulated_time=2020-02-07 12:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 14:20:17 INFO [TIME] window start=2020-02-07 12:00:00 end=2020-02-07 14:00:00 run_hours=2 has_assimilation=True
2026-07-02 14:20:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:20:19 INFO Hourly dataset computed and listing created
2026-07-02 14:20:25 INFO Hourly dataset computed
2026-07-02 14:20:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:20:26 INFO Hourly dataset computed and listing created
2026-07-02 14:20:27 INFO Hourly dataset computed
2026-07-02 14:20:27 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:20:28 INFO Hourly dataset computed and listing created
2026-07-02 14:20:28 INFO Hourly dataset computed
2026-07-02 14:20:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:20:29 INFO Hourly dataset computed and listing created
2026-07-02 14:20:30 INFO Hourly dataset computed
2026-07-02 14:20:30 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:20:31 INFO Hourly dataset computed and listing created
2026-07-02 14:20:32 INFO Hourly dataset computed
2026-07-02 14:20:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:20:33 INFO Hourly dataset computed and listing created
2026-07-02 14:20:33 INFO Hourly dataset computed
2026-07-02 14:20:33 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:20:34 INFO Hourly dataset computed and listing created
2026-07-02 14:20:35 INFO Hourly dataset computed
2026-07-02 14:20:35 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:20:36 INFO Hourly dataset computed and listing created
2026-07-02 14:20:37 INFO Hourly dataset computed
2026-07-02 14:20:37 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:20:38 INFO Hourly dataset computed and listing created
2026-07-02 14:20:38 INFO Hourly dataset computed
2026-07-02 14:20:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:20:39 INFO Hourly dataset computed and listing created
2026-07-02 14:20:40 INFO Hourly dataset computed
2026-07-02 14:20:40 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:20:41 INFO Hourly dataset computed and listing created
2026-07-02 14:20:42 INFO Hourly dataset computed
2026-07-02 14:20:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:20:43 INFO Hourly dataset computed and listing created
2026-07-02 14:20:44 INFO Hourly dataset computed
2026-07-02 14:20:44 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:20:44 INFO Hourly dataset computed and listing created
2026-07-02 14:20:45 INFO Hourly dataset computed
2026-07-02 14:20:45 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:20:46 INFO Hourly dataset computed and listing created
2026-07-02 14:20:47 INFO Hourly dataset computed
2026-07-02 14:20:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:20:48 INFO Hourly dataset computed and listing created
2026-07-02 14:20:49 INFO Hourly dataset computed
2026-07-02 14:20:49 INFO ---------->>> Running CHIMERE model from 2020-02-07 12:00:00 to 2020-02-07 14:00:00
2026-07-02 14:20:49 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:20:49 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1
2026-07-02 14:20:49 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020711_1_ENS1.nc
2026-07-02 14:20:49 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-02 14:20:49 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:20:49 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-02 14:20:49 INFO Queuing job for member 1...
2026-07-02 14:20:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:20:49 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-02 14:20:49 INFO Found: ['5066455']
2026-07-02 14:20:54 INFO [TGCC-IRENE] Submitted job with ID:['5066455']
2026-07-02 14:20:54 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:20:54 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2
2026-07-02 14:20:54 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020711_1_ENS2.nc
2026-07-02 14:20:54 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-02 14:20:54 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:20:54 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-02 14:20:54 INFO Queuing job for member 2...
2026-07-02 14:20:54 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:20:54 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-02 14:20:55 INFO Found: ['5066457']
2026-07-02 14:21:00 INFO [TGCC-IRENE] Submitted job with ID:['5066457']
2026-07-02 14:21:00 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:21:00 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3
2026-07-02 14:21:00 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020711_1_ENS3.nc
2026-07-02 14:21:00 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-02 14:21:00 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:21:00 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-02 14:21:00 INFO Queuing job for member 3...
2026-07-02 14:21:00 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:21:00 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-02 14:21:01 INFO Found: ['5066459']
2026-07-02 14:21:06 INFO [TGCC-IRENE] Submitted job with ID:['5066459']
2026-07-02 14:21:06 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:21:06 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4
2026-07-02 14:21:06 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020711_1_ENS4.nc
2026-07-02 14:21:06 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-02 14:21:06 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:21:06 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-02 14:21:06 INFO Queuing job for member 4...
2026-07-02 14:21:06 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:21:06 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-02 14:21:07 INFO Found: ['5066460']
2026-07-02 14:21:12 INFO [TGCC-IRENE] Submitted job with ID:['5066460']
2026-07-02 14:21:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:21:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5
2026-07-02 14:21:12 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020711_1_ENS5.nc
2026-07-02 14:21:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-02 14:21:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:21:12 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-02 14:21:12 INFO Queuing job for member 5...
2026-07-02 14:21:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:21:12 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-02 14:21:12 INFO Found: ['5066461']
2026-07-02 14:21:17 INFO [TGCC-IRENE] Submitted job with ID:['5066461']
2026-07-02 14:21:17 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:21:17 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6
2026-07-02 14:21:17 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020711_1_ENS6.nc
2026-07-02 14:21:17 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-02 14:21:17 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:21:17 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-02 14:21:17 INFO Queuing job for member 6...
2026-07-02 14:21:17 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:21:17 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-02 14:21:18 INFO Found: ['5066462']
2026-07-02 14:21:23 INFO [TGCC-IRENE] Submitted job with ID:['5066462']
2026-07-02 14:21:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:21:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7
2026-07-02 14:21:23 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020711_1_ENS7.nc
2026-07-02 14:21:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-02 14:21:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:21:23 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-02 14:21:23 INFO Queuing job for member 7...
2026-07-02 14:21:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:21:23 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-02 14:21:24 INFO Found: ['5066463']
2026-07-02 14:21:29 INFO [TGCC-IRENE] Submitted job with ID:['5066463']
2026-07-02 14:21:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:21:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8
2026-07-02 14:21:29 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020711_1_ENS8.nc
2026-07-02 14:21:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-02 14:21:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:21:29 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-02 14:21:29 INFO Queuing job for member 8...
2026-07-02 14:21:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:21:29 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-02 14:21:29 INFO Found: ['5066464']
2026-07-02 14:21:34 INFO [TGCC-IRENE] Submitted job with ID:['5066464']
2026-07-02 14:21:34 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:21:34 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9
2026-07-02 14:21:34 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020711_1_ENS9.nc
2026-07-02 14:21:34 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-02 14:21:34 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:21:34 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-02 14:21:34 INFO Queuing job for member 9...
2026-07-02 14:21:34 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:21:34 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-02 14:21:35 INFO Found: ['5066465']
2026-07-02 14:21:40 INFO [TGCC-IRENE] Submitted job with ID:['5066465']
2026-07-02 14:21:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:21:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10
2026-07-02 14:21:40 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020711_1_ENS10.nc
2026-07-02 14:21:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-02 14:21:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:21:40 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-02 14:21:40 INFO Queuing job for member 10...
2026-07-02 14:21:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:21:40 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-02 14:21:43 INFO Found: ['5066469']
2026-07-02 14:21:48 INFO [TGCC-IRENE] Submitted job with ID:['5066469']
2026-07-02 14:21:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:21:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11
2026-07-02 14:21:48 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020711_1_ENS11.nc
2026-07-02 14:21:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-02 14:21:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:21:48 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-02 14:21:48 INFO Queuing job for member 11...
2026-07-02 14:21:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:21:48 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-02 14:21:49 INFO Found: ['5066481']
2026-07-02 14:21:54 INFO [TGCC-IRENE] Submitted job with ID:['5066481']
2026-07-02 14:21:54 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:21:54 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12
2026-07-02 14:21:54 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020711_1_ENS12.nc
2026-07-02 14:21:54 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-02 14:21:54 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:21:54 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-02 14:21:54 INFO Queuing job for member 12...
2026-07-02 14:21:54 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:21:54 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-02 14:21:55 INFO Found: ['5066490']
2026-07-02 14:22:00 INFO [TGCC-IRENE] Submitted job with ID:['5066490']
2026-07-02 14:22:00 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:22:00 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13
2026-07-02 14:22:00 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020711_1_ENS13.nc
2026-07-02 14:22:00 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-02 14:22:00 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:22:00 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-02 14:22:00 INFO Queuing job for member 13...
2026-07-02 14:22:00 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:22:00 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-02 14:22:00 INFO Found: ['5066499']
2026-07-02 14:22:05 INFO [TGCC-IRENE] Submitted job with ID:['5066499']
2026-07-02 14:22:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:22:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14
2026-07-02 14:22:05 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020711_1_ENS14.nc
2026-07-02 14:22:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-02 14:22:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:22:05 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-02 14:22:05 INFO Queuing job for member 14...
2026-07-02 14:22:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:22:05 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-02 14:22:06 INFO Found: ['5066508']
2026-07-02 14:22:11 INFO [TGCC-IRENE] Submitted job with ID:['5066508']
2026-07-02 14:22:11 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:22:11 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15
2026-07-02 14:22:11 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020711_1_ENS15.nc
2026-07-02 14:22:11 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-02 14:22:11 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:22:11 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-02 14:22:11 INFO Queuing job for member 15...
2026-07-02 14:22:11 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:22:11 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-02 14:22:12 INFO Found: ['5066510']
2026-07-02 14:22:17 INFO [TGCC-IRENE] Submitted job with ID:['5066510']
2026-07-02 14:22:17 INFO Checking job status ...
2026-07-02 14:22:17 INFO None 5066455: status RUNNING/PENDING
2026-07-02 14:22:17 INFO None 5066457: status RUNNING/PENDING
2026-07-02 14:22:17 INFO None 5066459: status RUNNING/PENDING
2026-07-02 14:22:17 INFO None 5066460: status RUNNING/PENDING
2026-07-02 14:22:17 INFO None 5066461: status RUNNING/PENDING
2026-07-02 14:22:17 INFO None 5066462: status RUNNING/PENDING
2026-07-02 14:22:17 INFO None 5066463: status RUNNING/PENDING
2026-07-02 14:22:17 INFO None 5066464: status RUNNING/PENDING
2026-07-02 14:22:17 INFO None 5066465: status RUNNING/PENDING
2026-07-02 14:22:17 INFO None 5066469: status RUNNING/PENDING
2026-07-02 14:22:17 INFO None 5066481: status RUNNING/PENDING
2026-07-02 14:22:17 INFO None 5066490: status RUNNING/PENDING
2026-07-02 14:22:17 INFO None 5066499: status RUNNING/PENDING
2026-07-02 14:22:17 INFO None 5066508: status RUNNING/PENDING
2026-07-02 14:22:17 INFO None 5066510: status RUNNING/PENDING
2026-07-02 14:22:17 INFO Jobs still running: ['5066455', '5066457', '5066459', '5066460', '5066461', '5066462', '5066463', '5066464', '5066465', '5066469', '5066481', '5066490', '5066499', '5066508', '5066510']. Waiting...
2026-07-02 14:22:33 INFO None 5066455: status RUNNING/PENDING
2026-07-02 14:22:33 INFO None 5066457: status RUNNING/PENDING
2026-07-02 14:22:33 INFO None 5066459: status RUNNING/PENDING
2026-07-02 14:22:33 INFO None 5066460: status RUNNING/PENDING
2026-07-02 14:22:33 INFO None 5066461: status RUNNING/PENDING
2026-07-02 14:22:33 INFO None 5066462: status RUNNING/PENDING
2026-07-02 14:22:33 INFO None 5066463: status RUNNING/PENDING
2026-07-02 14:22:33 INFO None 5066464: status RUNNING/PENDING
2026-07-02 14:22:33 INFO None 5066465: status RUNNING/PENDING
2026-07-02 14:22:33 INFO None 5066469: status RUNNING/PENDING
2026-07-02 14:22:33 INFO None 5066481: status RUNNING/PENDING
2026-07-02 14:22:33 INFO None 5066490: status RUNNING/PENDING
2026-07-02 14:22:33 INFO None 5066499: status RUNNING/PENDING
2026-07-02 14:22:33 INFO None 5066508: status RUNNING/PENDING
2026-07-02 14:22:33 INFO None 5066510: status RUNNING/PENDING
2026-07-02 14:22:33 INFO Jobs still running: ['5066455', '5066457', '5066459', '5066460', '5066461', '5066462', '5066463', '5066464', '5066465', '5066469', '5066481', '5066490', '5066499', '5066508', '5066510']. Waiting...
2026-07-02 14:22:48 INFO None 5066455: status RUNNING/PENDING
2026-07-02 14:22:48 INFO None 5066457: status RUNNING/PENDING
2026-07-02 14:22:48 INFO None 5066459: status RUNNING/PENDING
2026-07-02 14:22:48 INFO None 5066460: status RUNNING/PENDING
2026-07-02 14:22:48 INFO None 5066461: status RUNNING/PENDING
2026-07-02 14:22:48 INFO None 5066462: status RUNNING/PENDING
2026-07-02 14:22:48 INFO None 5066463: status RUNNING/PENDING
2026-07-02 14:22:48 INFO None 5066464: status RUNNING/PENDING
2026-07-02 14:22:48 INFO None 5066465: status RUNNING/PENDING
2026-07-02 14:22:49 INFO None 5066469: status RUNNING/PENDING
2026-07-02 14:22:49 INFO None 5066481: status RUNNING/PENDING
2026-07-02 14:22:49 INFO None 5066490: status RUNNING/PENDING
2026-07-02 14:22:49 INFO None 5066499: status RUNNING/PENDING
2026-07-02 14:22:49 INFO None 5066508: status RUNNING/PENDING
2026-07-02 14:22:49 INFO None 5066510: status RUNNING/PENDING
2026-07-02 14:22:49 INFO Jobs still running: ['5066455', '5066457', '5066459', '5066460', '5066461', '5066462', '5066463', '5066464', '5066465', '5066469', '5066481', '5066490', '5066499', '5066508', '5066510']. Waiting...
2026-07-02 14:23:04 INFO None 5066455: status RUNNING/PENDING
2026-07-02 14:23:04 INFO None 5066457: status RUNNING/PENDING
2026-07-02 14:23:04 INFO None 5066459: status RUNNING/PENDING
2026-07-02 14:23:04 INFO None 5066460: status RUNNING/PENDING
2026-07-02 14:23:04 INFO None 5066461: status RUNNING/PENDING
2026-07-02 14:23:04 INFO None 5066462: status RUNNING/PENDING
2026-07-02 14:23:04 INFO None 5066463: status RUNNING/PENDING
2026-07-02 14:23:04 INFO None 5066464: status RUNNING/PENDING
2026-07-02 14:23:04 INFO None 5066465: status RUNNING/PENDING
2026-07-02 14:23:04 INFO None 5066469: status RUNNING/PENDING
2026-07-02 14:23:04 INFO None 5066481: status RUNNING/PENDING
2026-07-02 14:23:04 INFO None 5066490: status RUNNING/PENDING
2026-07-02 14:23:04 INFO None 5066499: status RUNNING/PENDING
2026-07-02 14:23:04 INFO None 5066508: status RUNNING/PENDING
2026-07-02 14:23:04 INFO None 5066510: status RUNNING/PENDING
2026-07-02 14:23:04 INFO Jobs still running: ['5066455', '5066457', '5066459', '5066460', '5066461', '5066462', '5066463', '5066464', '5066465', '5066469', '5066481', '5066490', '5066499', '5066508', '5066510']. Waiting...
2026-07-02 14:23:19 INFO None 5066455: status RUNNING/PENDING
2026-07-02 14:23:19 INFO None 5066457: status RUNNING/PENDING
2026-07-02 14:23:19 INFO None 5066459: status RUNNING/PENDING
2026-07-02 14:23:19 INFO None 5066460: status RUNNING/PENDING
2026-07-02 14:23:19 INFO None 5066461: status RUNNING/PENDING
2026-07-02 14:23:19 INFO None 5066462: status RUNNING/PENDING
2026-07-02 14:23:19 INFO None 5066463: status RUNNING/PENDING
2026-07-02 14:23:19 INFO None 5066464: status RUNNING/PENDING
2026-07-02 14:23:19 INFO None 5066465: status RUNNING/PENDING
2026-07-02 14:23:19 INFO None 5066469: status RUNNING/PENDING
2026-07-02 14:23:19 INFO None 5066481: status RUNNING/PENDING
2026-07-02 14:23:19 INFO None 5066490: status RUNNING/PENDING
2026-07-02 14:23:19 INFO None 5066499: status RUNNING/PENDING
2026-07-02 14:23:19 INFO None 5066508: status RUNNING/PENDING
2026-07-02 14:23:19 INFO None 5066510: status RUNNING/PENDING
2026-07-02 14:23:19 INFO Jobs still running: ['5066455', '5066457', '5066459', '5066460', '5066461', '5066462', '5066463', '5066464', '5066465', '5066469', '5066481', '5066490', '5066499', '5066508', '5066510']. Waiting...
2026-07-02 14:23:36 INFO None 5066455: status RUNNING/PENDING
2026-07-02 14:23:36 INFO None 5066457: status RUNNING/PENDING
2026-07-02 14:23:36 INFO None 5066459: status RUNNING/PENDING
2026-07-02 14:23:36 INFO None 5066460: status RUNNING/PENDING
2026-07-02 14:23:36 INFO None 5066461: status RUNNING/PENDING
2026-07-02 14:23:36 INFO None 5066462: status RUNNING/PENDING
2026-07-02 14:23:36 INFO None 5066463: status RUNNING/PENDING
2026-07-02 14:23:36 INFO None 5066464: status RUNNING/PENDING
2026-07-02 14:23:36 INFO None 5066465: status RUNNING/PENDING
2026-07-02 14:23:36 INFO None 5066469: status RUNNING/PENDING
2026-07-02 14:23:36 INFO None 5066481: status RUNNING/PENDING
2026-07-02 14:23:36 INFO None 5066490: status RUNNING/PENDING
2026-07-02 14:23:36 INFO None 5066499: status RUNNING/PENDING
2026-07-02 14:23:36 INFO None 5066508: status RUNNING/PENDING
2026-07-02 14:23:36 INFO None 5066510: status RUNNING/PENDING
2026-07-02 14:23:36 INFO Jobs still running: ['5066455', '5066457', '5066459', '5066460', '5066461', '5066462', '5066463', '5066464', '5066465', '5066469', '5066481', '5066490', '5066499', '5066508', '5066510']. Waiting...
2026-07-02 14:23:51 INFO None 5066455: status RUNNING/PENDING
2026-07-02 14:23:51 INFO None 5066457: status RUNNING/PENDING
2026-07-02 14:23:51 INFO None 5066459: status RUNNING/PENDING
2026-07-02 14:23:51 INFO None 5066460: status RUNNING/PENDING
2026-07-02 14:23:51 INFO None 5066461: status RUNNING/PENDING
2026-07-02 14:23:51 INFO None 5066462: status RUNNING/PENDING
2026-07-02 14:23:51 INFO None 5066463: status RUNNING/PENDING
2026-07-02 14:23:51 INFO None 5066464: status RUNNING/PENDING
2026-07-02 14:23:51 INFO None 5066465: status RUNNING/PENDING
2026-07-02 14:23:51 INFO None 5066469: status RUNNING/PENDING
2026-07-02 14:23:51 INFO None 5066481: status RUNNING/PENDING
2026-07-02 14:23:52 INFO None 5066490: status RUNNING/PENDING
2026-07-02 14:23:52 INFO None 5066499: status RUNNING/PENDING
2026-07-02 14:23:52 INFO None 5066508: status RUNNING/PENDING
2026-07-02 14:23:52 INFO None 5066510: status RUNNING/PENDING
2026-07-02 14:23:52 INFO Jobs still running: ['5066455', '5066457', '5066459', '5066460', '5066461', '5066462', '5066463', '5066464', '5066465', '5066469', '5066481', '5066490', '5066499', '5066508', '5066510']. Waiting...
2026-07-02 14:24:07 INFO None 5066455: status RUNNING/PENDING
2026-07-02 14:24:07 INFO None 5066457: status RUNNING/PENDING
2026-07-02 14:24:07 INFO None 5066459: status FINISHED
2026-07-02 14:24:07 INFO None 5066460: status RUNNING/PENDING
2026-07-02 14:24:07 INFO None 5066461: status RUNNING/PENDING
2026-07-02 14:24:07 INFO None 5066462: status RUNNING/PENDING
2026-07-02 14:24:07 INFO None 5066463: status RUNNING/PENDING
2026-07-02 14:24:07 INFO None 5066464: status RUNNING/PENDING
2026-07-02 14:24:07 INFO None 5066465: status RUNNING/PENDING
2026-07-02 14:24:07 INFO None 5066469: status RUNNING/PENDING
2026-07-02 14:24:07 INFO None 5066481: status RUNNING/PENDING
2026-07-02 14:24:07 INFO None 5066490: status RUNNING/PENDING
2026-07-02 14:24:07 INFO None 5066499: status RUNNING/PENDING
2026-07-02 14:24:07 INFO None 5066508: status RUNNING/PENDING
2026-07-02 14:24:07 INFO None 5066510: status RUNNING/PENDING
2026-07-02 14:24:07 INFO Jobs still running: ['5066455', '5066457', '5066460', '5066461', '5066462', '5066463', '5066464', '5066465', '5066469', '5066481', '5066490', '5066499', '5066508', '5066510']. Waiting...
2026-07-02 14:24:22 INFO None 5066455: status FINISHED
2026-07-02 14:24:22 INFO None 5066457: status RUNNING/PENDING
2026-07-02 14:24:22 INFO None 5066459: status FINISHED
2026-07-02 14:24:22 INFO None 5066460: status RUNNING/PENDING
2026-07-02 14:24:22 INFO None 5066461: status RUNNING/PENDING
2026-07-02 14:24:22 INFO None 5066462: status RUNNING/PENDING
2026-07-02 14:24:22 INFO None 5066463: status FINISHED
2026-07-02 14:24:22 INFO None 5066464: status RUNNING/PENDING
2026-07-02 14:24:22 INFO None 5066465: status RUNNING/PENDING
2026-07-02 14:24:22 INFO None 5066469: status RUNNING/PENDING
2026-07-02 14:24:22 INFO None 5066481: status RUNNING/PENDING
2026-07-02 14:24:22 INFO None 5066490: status RUNNING/PENDING
2026-07-02 14:24:22 INFO None 5066499: status RUNNING/PENDING
2026-07-02 14:24:22 INFO None 5066508: status RUNNING/PENDING
2026-07-02 14:24:22 INFO None 5066510: status RUNNING/PENDING
2026-07-02 14:24:22 INFO Jobs still running: ['5066457', '5066460', '5066461', '5066462', '5066464', '5066465', '5066469', '5066481', '5066490', '5066499', '5066508', '5066510']. Waiting...
2026-07-02 14:24:37 INFO None 5066455: status FINISHED
2026-07-02 14:24:37 INFO None 5066457: status RUNNING/PENDING
2026-07-02 14:24:37 INFO None 5066459: status FINISHED
2026-07-02 14:24:37 INFO None 5066460: status FINISHED
2026-07-02 14:24:37 INFO None 5066461: status RUNNING/PENDING
2026-07-02 14:24:37 INFO None 5066462: status RUNNING/PENDING
2026-07-02 14:24:37 INFO None 5066463: status FINISHED
2026-07-02 14:24:37 INFO None 5066464: status FINISHED
2026-07-02 14:24:37 INFO None 5066465: status RUNNING/PENDING
2026-07-02 14:24:37 INFO None 5066469: status RUNNING/PENDING
2026-07-02 14:24:37 INFO None 5066481: status RUNNING/PENDING
2026-07-02 14:24:37 INFO None 5066490: status RUNNING/PENDING
2026-07-02 14:24:37 INFO None 5066499: status RUNNING/PENDING
2026-07-02 14:24:37 INFO None 5066508: status RUNNING/PENDING
2026-07-02 14:24:37 INFO None 5066510: status RUNNING/PENDING
2026-07-02 14:24:37 INFO Jobs still running: ['5066457', '5066461', '5066462', '5066465', '5066469', '5066481', '5066490', '5066499', '5066508', '5066510']. Waiting...
2026-07-02 14:24:52 INFO None 5066455: status FINISHED
2026-07-02 14:24:52 INFO None 5066457: status RUNNING/PENDING
2026-07-02 14:24:52 INFO None 5066459: status FINISHED
2026-07-02 14:24:52 INFO None 5066460: status FINISHED
2026-07-02 14:24:53 INFO None 5066461: status RUNNING/PENDING
2026-07-02 14:24:53 INFO None 5066462: status RUNNING/PENDING
2026-07-02 14:24:53 INFO None 5066463: status FINISHED
2026-07-02 14:24:53 INFO None 5066464: status FINISHED
2026-07-02 14:24:53 INFO None 5066465: status RUNNING/PENDING
2026-07-02 14:24:53 INFO None 5066469: status RUNNING/PENDING
2026-07-02 14:24:53 INFO None 5066481: status RUNNING/PENDING
2026-07-02 14:24:53 INFO None 5066490: status RUNNING/PENDING
2026-07-02 14:24:53 INFO None 5066499: status RUNNING/PENDING
2026-07-02 14:24:53 INFO None 5066508: status RUNNING/PENDING
2026-07-02 14:24:53 INFO None 5066510: status RUNNING/PENDING
2026-07-02 14:24:53 INFO Jobs still running: ['5066457', '5066461', '5066462', '5066465', '5066469', '5066481', '5066490', '5066499', '5066508', '5066510']. Waiting...
2026-07-02 14:25:08 INFO None 5066455: status FINISHED
2026-07-02 14:25:08 INFO None 5066457: status FINISHED
2026-07-02 14:25:08 INFO None 5066459: status FINISHED
2026-07-02 14:25:08 INFO None 5066460: status FINISHED
2026-07-02 14:25:08 INFO None 5066461: status RUNNING/PENDING
2026-07-02 14:25:08 INFO None 5066462: status FINISHED
2026-07-02 14:25:08 INFO None 5066463: status FINISHED
2026-07-02 14:25:08 INFO None 5066464: status FINISHED
2026-07-02 14:25:08 INFO None 5066465: status RUNNING/PENDING
2026-07-02 14:25:08 INFO None 5066469: status RUNNING/PENDING
2026-07-02 14:25:08 INFO None 5066481: status RUNNING/PENDING
2026-07-02 14:25:08 INFO None 5066490: status RUNNING/PENDING
2026-07-02 14:25:08 INFO None 5066499: status RUNNING/PENDING
2026-07-02 14:25:08 INFO None 5066508: status RUNNING/PENDING
2026-07-02 14:25:08 INFO None 5066510: status FINISHED
2026-07-02 14:25:08 INFO Jobs still running: ['5066461', '5066465', '5066469', '5066481', '5066490', '5066499', '5066508']. Waiting...
2026-07-02 14:25:25 INFO None 5066455: status FINISHED
2026-07-02 14:25:25 INFO None 5066457: status FINISHED
2026-07-02 14:25:25 INFO None 5066459: status FINISHED
2026-07-02 14:25:25 INFO None 5066460: status FINISHED
2026-07-02 14:25:25 INFO None 5066461: status FINISHED
2026-07-02 14:25:25 INFO None 5066462: status FINISHED
2026-07-02 14:25:25 INFO None 5066463: status FINISHED
2026-07-02 14:25:25 INFO None 5066464: status FINISHED
2026-07-02 14:25:25 INFO None 5066465: status RUNNING/PENDING
2026-07-02 14:25:25 INFO None 5066469: status RUNNING/PENDING
2026-07-02 14:25:25 INFO None 5066481: status RUNNING/PENDING
2026-07-02 14:25:25 INFO None 5066490: status RUNNING/PENDING
2026-07-02 14:25:25 INFO None 5066499: status RUNNING/PENDING
2026-07-02 14:25:25 INFO None 5066508: status RUNNING/PENDING
2026-07-02 14:25:25 INFO None 5066510: status FINISHED
2026-07-02 14:25:25 INFO Jobs still running: ['5066465', '5066469', '5066481', '5066490', '5066499', '5066508']. Waiting...
2026-07-02 14:25:40 INFO None 5066455: status FINISHED
2026-07-02 14:25:40 INFO None 5066457: status FINISHED
2026-07-02 14:25:40 INFO None 5066459: status FINISHED
2026-07-02 14:25:40 INFO None 5066460: status FINISHED
2026-07-02 14:25:40 INFO None 5066461: status FINISHED
2026-07-02 14:25:40 INFO None 5066462: status FINISHED
2026-07-02 14:25:40 INFO None 5066463: status FINISHED
2026-07-02 14:25:40 INFO None 5066464: status FINISHED
2026-07-02 14:25:40 INFO None 5066465: status RUNNING/PENDING
2026-07-02 14:25:40 INFO None 5066469: status RUNNING/PENDING
2026-07-02 14:25:40 INFO None 5066481: status RUNNING/PENDING
2026-07-02 14:25:40 INFO None 5066490: status FINISHED
2026-07-02 14:25:40 INFO None 5066499: status RUNNING/PENDING
2026-07-02 14:25:40 INFO None 5066508: status RUNNING/PENDING
2026-07-02 14:25:40 INFO None 5066510: status FINISHED
2026-07-02 14:25:40 INFO Jobs still running: ['5066465', '5066469', '5066481', '5066499', '5066508']. Waiting...
2026-07-02 14:25:55 INFO None 5066455: status FINISHED
2026-07-02 14:25:55 INFO None 5066457: status FINISHED
2026-07-02 14:25:55 INFO None 5066459: status FINISHED
2026-07-02 14:25:55 INFO None 5066460: status FINISHED
2026-07-02 14:25:55 INFO None 5066461: status FINISHED
2026-07-02 14:25:55 INFO None 5066462: status FINISHED
2026-07-02 14:25:55 INFO None 5066463: status FINISHED
2026-07-02 14:25:55 INFO None 5066464: status FINISHED
2026-07-02 14:25:55 INFO None 5066465: status RUNNING/PENDING
2026-07-02 14:25:55 INFO None 5066469: status FINISHED
2026-07-02 14:25:55 INFO None 5066481: status FINISHED
2026-07-02 14:25:55 INFO None 5066490: status FINISHED
2026-07-02 14:25:55 INFO None 5066499: status FINISHED
2026-07-02 14:25:55 INFO None 5066508: status FINISHED
2026-07-02 14:25:55 INFO None 5066510: status FINISHED
2026-07-02 14:25:55 INFO Jobs still running: ['5066465']. Waiting...
2026-07-02 14:26:10 INFO None 5066455: status FINISHED
2026-07-02 14:26:10 INFO None 5066457: status FINISHED
2026-07-02 14:26:11 INFO None 5066459: status FINISHED
2026-07-02 14:26:11 INFO None 5066460: status FINISHED
2026-07-02 14:26:11 INFO None 5066461: status FINISHED
2026-07-02 14:26:11 INFO None 5066462: status FINISHED
2026-07-02 14:26:11 INFO None 5066463: status FINISHED
2026-07-02 14:26:11 INFO None 5066464: status FINISHED
2026-07-02 14:26:11 INFO None 5066465: status RUNNING/PENDING
2026-07-02 14:26:11 INFO None 5066469: status FINISHED
2026-07-02 14:26:11 INFO None 5066481: status FINISHED
2026-07-02 14:26:11 INFO None 5066490: status FINISHED
2026-07-02 14:26:11 INFO None 5066499: status FINISHED
2026-07-02 14:26:11 INFO None 5066508: status FINISHED
2026-07-02 14:26:11 INFO None 5066510: status FINISHED
2026-07-02 14:26:11 INFO Jobs still running: ['5066465']. Waiting...
2026-07-02 14:26:26 INFO None 5066455: status FINISHED
2026-07-02 14:26:26 INFO None 5066457: status FINISHED
2026-07-02 14:26:26 INFO None 5066459: status FINISHED
2026-07-02 14:26:26 INFO None 5066460: status FINISHED
2026-07-02 14:26:26 INFO None 5066461: status FINISHED
2026-07-02 14:26:26 INFO None 5066462: status FINISHED
2026-07-02 14:26:26 INFO None 5066463: status FINISHED
2026-07-02 14:26:26 INFO None 5066464: status FINISHED
2026-07-02 14:26:26 INFO None 5066465: status FINISHED
2026-07-02 14:26:26 INFO None 5066469: status FINISHED
2026-07-02 14:26:26 INFO None 5066481: status FINISHED
2026-07-02 14:26:26 INFO None 5066490: status FINISHED
2026-07-02 14:26:26 INFO None 5066499: status FINISHED
2026-07-02 14:26:26 INFO None 5066508: status FINISHED
2026-07-02 14:26:26 INFO None 5066510: status FINISHED
2026-07-02 14:26:26 INFO Jobs ['5066455', '5066457', '5066459', '5066460', '5066461', '5066462', '5066463', '5066464', '5066465', '5066469', '5066481', '5066490', '5066499', '5066508', '5066510'] have finished
2026-07-02 14:26:26 INFO Checking restart files were created ...
2026-07-02 14:26:26 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020712_2_ENS1.nc(1002685915 bytes)
2026-07-02 14:26:26 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020712_2_ENS2.nc(1002685915 bytes)
2026-07-02 14:26:26 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020712_2_ENS3.nc(1002685915 bytes)
2026-07-02 14:26:26 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020712_2_ENS4.nc(1002685915 bytes)
2026-07-02 14:26:26 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020712_2_ENS5.nc(1002685915 bytes)
2026-07-02 14:26:26 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020712_2_ENS6.nc(1002685915 bytes)
2026-07-02 14:26:26 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020712_2_ENS7.nc(1002685915 bytes)
2026-07-02 14:26:26 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020712_2_ENS8.nc(1002685915 bytes)
2026-07-02 14:26:26 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020712_2_ENS9.nc(1002685915 bytes)
2026-07-02 14:26:26 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020712_2_ENS10.nc(1002685915 bytes)
2026-07-02 14:26:26 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020712_2_ENS11.nc(1002685915 bytes)
2026-07-02 14:26:26 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020712_2_ENS12.nc(1002685915 bytes)
2026-07-02 14:26:26 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020712_2_ENS13.nc(1002685915 bytes)
2026-07-02 14:26:26 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020712_2_ENS14.nc(1002685915 bytes)
2026-07-02 14:26:26 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020712_2_ENS15.nc(1002685915 bytes)
2026-07-02 14:26:26 INFO  Run_model() completed successfully.
2026-07-02 14:26:26 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 12:00:00 simulated_time=2020-02-07 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 14:26:26 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 14:00:00 days=153073 seconds=50400
2026-07-02 14:26:26 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-02 14:26:26 INFO [TIME] increment current_time 2020-02-07 12:00:00 -> 2020-02-07 14:00:00
2026-07-02 14:26:26 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 14:00:00 simulated_time=2020-02-07 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 14:26:26 INFO ---------->>> Running process_satellite_data()
2026-07-02 14:26:26 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12019.nc
2026-07-02 14:26:26 INFO ---------->>> Running run_obs_converter()
2026-07-02 14:26:27 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_50083_153073.out
2026-07-02 14:26:27 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_50083_153073.out
2026-07-02 14:26:27 INFO ---------->>> Running DART
2026-07-02 14:26:27 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-02 14:26:27 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/chim_ENS1_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/chim_ENS1_2020020714_1_out_toDART.nc
2026-07-02 14:26:27 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/chim_ENS2_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/chim_ENS2_2020020714_1_out_toDART.nc
2026-07-02 14:26:27 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/chim_ENS3_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/chim_ENS3_2020020714_1_out_toDART.nc
2026-07-02 14:26:27 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/chim_ENS4_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/chim_ENS4_2020020714_1_out_toDART.nc
2026-07-02 14:26:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/chim_ENS5_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/chim_ENS5_2020020714_1_out_toDART.nc
2026-07-02 14:26:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/chim_ENS6_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/chim_ENS6_2020020714_1_out_toDART.nc
2026-07-02 14:26:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/chim_ENS7_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/chim_ENS7_2020020714_1_out_toDART.nc
2026-07-02 14:26:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/chim_ENS8_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/chim_ENS8_2020020714_1_out_toDART.nc
2026-07-02 14:26:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/chim_ENS9_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/chim_ENS9_2020020714_1_out_toDART.nc
2026-07-02 14:26:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/chim_ENS10_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/chim_ENS10_2020020714_1_out_toDART.nc
2026-07-02 14:26:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/chim_ENS11_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/chim_ENS11_2020020714_1_out_toDART.nc
2026-07-02 14:26:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/chim_ENS12_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/chim_ENS12_2020020714_1_out_toDART.nc
2026-07-02 14:26:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/chim_ENS13_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/chim_ENS13_2020020714_1_out_toDART.nc
2026-07-02 14:26:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/chim_ENS14_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/chim_ENS14_2020020714_1_out_toDART.nc
2026-07-02 14:26:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/chim_ENS15_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/chim_ENS15_2020020714_1_out_toDART.nc
2026-07-02 14:26:31 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-02 14:26:31 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-02 14:26:31 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-02 14:26:31 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-02 14:26:31 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-02 14:26:31 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-02 14:26:40 INFO Found: []
2026-07-02 14:26:40 INFO No job id returned by command ./run_filter.bsh
2026-07-02 14:26:40 INFO No monitoring will be performed
2026-07-02 14:26:40 INFO Moving DART output files to analysis and preassim directories for date 2020020714 if present ...
2026-07-02 14:26:40 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020714'
2026-07-02 14:26:40 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020714'
2026-07-02 14:26:40 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020714'
2026-07-02 14:26:40 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020714'
2026-07-02 14:26:40 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020714'
2026-07-02 14:26:40 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020714'
2026-07-02 14:26:40 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020714'
2026-07-02 14:26:40 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020714'
2026-07-02 14:26:40 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020714'
2026-07-02 14:26:40 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020714'
2026-07-02 14:26:40 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020714'
2026-07-02 14:26:40 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020714'
2026-07-02 14:26:40 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020714'
2026-07-02 14:26:40 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020714'
2026-07-02 14:26:40 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020714'
2026-07-02 14:26:40 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020714'
2026-07-02 14:26:40 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020714'
2026-07-02 14:26:40 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020714'
2026-07-02 14:26:40 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020714'
2026-07-02 14:26:40 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020714'
2026-07-02 14:26:40 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020714'
2026-07-02 14:26:40 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020714'
2026-07-02 14:26:40 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020714'
2026-07-02 14:26:40 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020714'
2026-07-02 14:26:40 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020714'
2026-07-02 14:26:40 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020714'
2026-07-02 14:26:40 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020714'
2026-07-02 14:26:40 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020714'
2026-07-02 14:26:40 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020714'
2026-07-02 14:26:40 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020714'
2026-07-02 14:26:40 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020714'
2026-07-02 14:26:40 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/analysis/2020020714'
2026-07-02 14:26:40 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020714'
2026-07-02 14:26:40 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/preassim/2020020714'
2026-07-02 14:26:40 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-02 14:26:40 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-02 14:26:40 INFO run_dart() is DONE.
2026-07-02 14:26:40 INFO ---------->>> Running update_pollutant_in_end()
2026-07-02 14:26:40 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:26:40 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:26:41 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:26:41 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:26:42 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:26:42 INFO [Posterior Diff ENS1] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020714/diff_posterior_ENS1_2020020714.nc
2026-07-02 14:26:42 INFO [Posterior Diff ENS1] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020714/diff_posterior_ENS1_2020020714.relative.nc
2026-07-02 14:26:42 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:26:42 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:26:43 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:26:43 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:26:44 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:26:44 INFO [Posterior Diff ENS2] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020714/diff_posterior_ENS2_2020020714.nc
2026-07-02 14:26:44 INFO [Posterior Diff ENS2] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020714/diff_posterior_ENS2_2020020714.relative.nc
2026-07-02 14:26:44 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:26:44 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:26:45 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:26:45 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:26:46 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:26:46 INFO [Posterior Diff ENS3] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020714/diff_posterior_ENS3_2020020714.nc
2026-07-02 14:26:46 INFO [Posterior Diff ENS3] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020714/diff_posterior_ENS3_2020020714.relative.nc
2026-07-02 14:26:46 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:26:46 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:26:47 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:26:47 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:26:48 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:26:48 INFO [Posterior Diff ENS4] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020714/diff_posterior_ENS4_2020020714.nc
2026-07-02 14:26:48 INFO [Posterior Diff ENS4] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020714/diff_posterior_ENS4_2020020714.relative.nc
2026-07-02 14:26:48 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:26:48 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:26:49 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:26:49 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:26:50 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:26:50 INFO [Posterior Diff ENS5] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020714/diff_posterior_ENS5_2020020714.nc
2026-07-02 14:26:50 INFO [Posterior Diff ENS5] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020714/diff_posterior_ENS5_2020020714.relative.nc
2026-07-02 14:26:50 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:26:50 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:26:51 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:26:51 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:26:52 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:26:52 INFO [Posterior Diff ENS6] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020714/diff_posterior_ENS6_2020020714.nc
2026-07-02 14:26:52 INFO [Posterior Diff ENS6] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020714/diff_posterior_ENS6_2020020714.relative.nc
2026-07-02 14:26:52 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:26:52 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:26:53 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:26:53 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:26:54 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:26:54 INFO [Posterior Diff ENS7] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020714/diff_posterior_ENS7_2020020714.nc
2026-07-02 14:26:54 INFO [Posterior Diff ENS7] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020714/diff_posterior_ENS7_2020020714.relative.nc
2026-07-02 14:26:54 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:26:54 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:26:55 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:26:55 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:26:56 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:26:56 INFO [Posterior Diff ENS8] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020714/diff_posterior_ENS8_2020020714.nc
2026-07-02 14:26:56 INFO [Posterior Diff ENS8] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020714/diff_posterior_ENS8_2020020714.relative.nc
2026-07-02 14:26:56 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:26:56 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:26:57 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:26:57 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:26:58 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:26:58 INFO [Posterior Diff ENS9] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020714/diff_posterior_ENS9_2020020714.nc
2026-07-02 14:26:58 INFO [Posterior Diff ENS9] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020714/diff_posterior_ENS9_2020020714.relative.nc
2026-07-02 14:26:58 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:26:58 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:26:59 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:26:59 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:27:00 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:27:00 INFO [Posterior Diff ENS10] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020714/diff_posterior_ENS10_2020020714.nc
2026-07-02 14:27:00 INFO [Posterior Diff ENS10] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020714/diff_posterior_ENS10_2020020714.relative.nc
2026-07-02 14:27:00 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:27:00 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:27:01 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:27:01 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:27:02 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:27:02 INFO [Posterior Diff ENS11] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020714/diff_posterior_ENS11_2020020714.nc
2026-07-02 14:27:02 INFO [Posterior Diff ENS11] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020714/diff_posterior_ENS11_2020020714.relative.nc
2026-07-02 14:27:02 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:27:02 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:27:03 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:27:03 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:27:04 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:27:04 INFO [Posterior Diff ENS12] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020714/diff_posterior_ENS12_2020020714.nc
2026-07-02 14:27:04 INFO [Posterior Diff ENS12] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020714/diff_posterior_ENS12_2020020714.relative.nc
2026-07-02 14:27:05 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:27:05 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:27:05 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:27:05 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:27:06 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:27:06 INFO [Posterior Diff ENS13] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020714/diff_posterior_ENS13_2020020714.nc
2026-07-02 14:27:06 INFO [Posterior Diff ENS13] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020714/diff_posterior_ENS13_2020020714.relative.nc
2026-07-02 14:27:06 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:27:06 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:27:07 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:27:07 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:27:08 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:27:08 INFO [Posterior Diff ENS14] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020714/diff_posterior_ENS14_2020020714.nc
2026-07-02 14:27:08 INFO [Posterior Diff ENS14] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020714/diff_posterior_ENS14_2020020714.relative.nc
2026-07-02 14:27:09 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:27:09 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:27:10 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-02 14:27:10 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc using posterior/prior ratio.
2026-07-02 14:27:10 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-02 14:27:10 INFO [Posterior Diff ENS15] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020714/diff_posterior_ENS15_2020020714.nc
2026-07-02 14:27:10 INFO [Posterior Diff ENS15] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmpemis_0615_15m_low_v2/posteriors/2020020714/diff_posterior_ENS15_2020020714.relative.nc
2026-07-02 14:27:10 INFO Next run starts from 2020-02-07 14:00:00
2026-07-02 14:27:10 INFO Cycle is DONE; starting a new loop!
2026-07-02 14:27:10 INFO [TIME] step_end current_time=2020-02-07 14:00:00 simulated_time=2020-02-07 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 14:27:10 INFO [TIME] step_start current_time=2020-02-07 14:00:00 simulated_time=2020-02-07 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 14:27:10 INFO [TIME] window start=2020-02-07 14:00:00 end=2020-02-08 00:00:00 run_hours=10 has_assimilation=False
2026-07-02 14:27:10 INFO Copying EMIS of next day ...
2026-07-02 14:27:11 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens01.nc
2026-07-02 14:27:11 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:27:12 INFO Hourly dataset computed and listing created
2026-07-02 14:27:27 INFO Hourly dataset computed
2026-07-02 14:27:27 INFO Copying EMIS of next day ...
2026-07-02 14:27:27 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens02.nc
2026-07-02 14:27:27 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:27:28 INFO Hourly dataset computed and listing created
2026-07-02 14:27:31 INFO Hourly dataset computed
2026-07-02 14:27:31 INFO Copying EMIS of next day ...
2026-07-02 14:27:32 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens03.nc
2026-07-02 14:27:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:27:33 INFO Hourly dataset computed and listing created
2026-07-02 14:27:35 INFO Hourly dataset computed
2026-07-02 14:27:35 INFO Copying EMIS of next day ...
2026-07-02 14:27:36 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens04.nc
2026-07-02 14:27:36 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:27:37 INFO Hourly dataset computed and listing created
2026-07-02 14:27:40 INFO Hourly dataset computed
2026-07-02 14:27:40 INFO Copying EMIS of next day ...
2026-07-02 14:27:40 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens05.nc
2026-07-02 14:27:40 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:27:41 INFO Hourly dataset computed and listing created
2026-07-02 14:27:43 INFO Hourly dataset computed
2026-07-02 14:27:43 INFO Copying EMIS of next day ...
2026-07-02 14:27:44 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens06.nc
2026-07-02 14:27:44 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:27:45 INFO Hourly dataset computed and listing created
2026-07-02 14:27:48 INFO Hourly dataset computed
2026-07-02 14:27:48 INFO Copying EMIS of next day ...
2026-07-02 14:27:49 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens07.nc
2026-07-02 14:27:49 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:27:50 INFO Hourly dataset computed and listing created
2026-07-02 14:27:53 INFO Hourly dataset computed
2026-07-02 14:27:53 INFO Copying EMIS of next day ...
2026-07-02 14:27:53 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens08.nc
2026-07-02 14:27:53 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:27:54 INFO Hourly dataset computed and listing created
2026-07-02 14:27:57 INFO Hourly dataset computed
2026-07-02 14:27:57 INFO Copying EMIS of next day ...
2026-07-02 14:27:57 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens09.nc
2026-07-02 14:27:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:27:58 INFO Hourly dataset computed and listing created
2026-07-02 14:28:01 INFO Hourly dataset computed
2026-07-02 14:28:01 INFO Copying EMIS of next day ...
2026-07-02 14:28:01 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens10.nc
2026-07-02 14:28:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:28:02 INFO Hourly dataset computed and listing created
2026-07-02 14:28:05 INFO Hourly dataset computed
2026-07-02 14:28:05 INFO Copying EMIS of next day ...
2026-07-02 14:28:05 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens11.nc
2026-07-02 14:28:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:28:07 INFO Hourly dataset computed and listing created
2026-07-02 14:28:09 INFO Hourly dataset computed
2026-07-02 14:28:09 INFO Copying EMIS of next day ...
2026-07-02 14:28:10 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens12.nc
2026-07-02 14:28:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:28:11 INFO Hourly dataset computed and listing created
2026-07-02 14:28:14 INFO Hourly dataset computed
2026-07-02 14:28:14 INFO Copying EMIS of next day ...
2026-07-02 14:28:14 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens13.nc
2026-07-02 14:28:14 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:28:15 INFO Hourly dataset computed and listing created
2026-07-02 14:28:17 INFO Hourly dataset computed
2026-07-02 14:28:18 INFO Copying EMIS of next day ...
2026-07-02 14:28:18 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens14.nc
2026-07-02 14:28:18 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:28:19 INFO Hourly dataset computed and listing created
2026-07-02 14:28:23 INFO Hourly dataset computed
2026-07-02 14:28:23 INFO Copying EMIS of next day ...
2026-07-02 14:28:23 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens15.nc
2026-07-02 14:28:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-02 14:28:25 INFO Hourly dataset computed and listing created
2026-07-02 14:28:33 INFO Hourly dataset computed
2026-07-02 14:28:33 INFO ---------->>> Running CHIMERE model from 2020-02-07 14:00:00 to 2020-02-08 00:00:00
2026-07-02 14:28:33 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:28:33 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1
2026-07-02 14:28:33 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020712_2_ENS1.nc
2026-07-02 14:28:33 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-02 14:28:33 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:28:33 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-02 14:28:33 INFO Queuing job for member 1...
2026-07-02 14:28:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:28:33 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-02 14:28:34 INFO Found: ['5066551']
2026-07-02 14:28:39 INFO [TGCC-IRENE] Submitted job with ID:['5066551']
2026-07-02 14:28:39 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:28:39 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2
2026-07-02 14:28:39 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020712_2_ENS2.nc
2026-07-02 14:28:39 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-02 14:28:39 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:28:39 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-02 14:28:39 INFO Queuing job for member 2...
2026-07-02 14:28:39 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:28:39 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-02 14:28:40 INFO Found: ['5066554']
2026-07-02 14:28:45 INFO [TGCC-IRENE] Submitted job with ID:['5066554']
2026-07-02 14:28:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:28:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3
2026-07-02 14:28:45 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020712_2_ENS3.nc
2026-07-02 14:28:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-02 14:28:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:28:45 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-02 14:28:45 INFO Queuing job for member 3...
2026-07-02 14:28:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:28:45 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-02 14:28:46 INFO Found: ['5066555']
2026-07-02 14:28:51 INFO [TGCC-IRENE] Submitted job with ID:['5066555']
2026-07-02 14:28:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:28:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4
2026-07-02 14:28:51 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020712_2_ENS4.nc
2026-07-02 14:28:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-02 14:28:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:28:51 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-02 14:28:51 INFO Queuing job for member 4...
2026-07-02 14:28:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:28:51 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-02 14:28:51 INFO Found: ['5066556']
2026-07-02 14:28:56 INFO [TGCC-IRENE] Submitted job with ID:['5066556']
2026-07-02 14:28:56 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:28:56 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5
2026-07-02 14:28:56 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020712_2_ENS5.nc
2026-07-02 14:28:56 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-02 14:28:56 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:28:56 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-02 14:28:56 INFO Queuing job for member 5...
2026-07-02 14:28:56 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:28:56 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-02 14:28:57 INFO Found: ['5066557']
2026-07-02 14:29:02 INFO [TGCC-IRENE] Submitted job with ID:['5066557']
2026-07-02 14:29:02 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:29:02 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6
2026-07-02 14:29:02 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020712_2_ENS6.nc
2026-07-02 14:29:02 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-02 14:29:02 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:29:02 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-02 14:29:02 INFO Queuing job for member 6...
2026-07-02 14:29:02 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:29:02 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-02 14:29:03 INFO Found: ['5066559']
2026-07-02 14:29:08 INFO [TGCC-IRENE] Submitted job with ID:['5066559']
2026-07-02 14:29:08 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:29:08 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7
2026-07-02 14:29:08 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020712_2_ENS7.nc
2026-07-02 14:29:08 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-02 14:29:08 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:29:08 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-02 14:29:08 INFO Queuing job for member 7...
2026-07-02 14:29:08 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:29:08 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-02 14:29:10 INFO Found: ['5066560']
2026-07-02 14:29:15 INFO [TGCC-IRENE] Submitted job with ID:['5066560']
2026-07-02 14:29:15 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:29:15 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8
2026-07-02 14:29:15 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020712_2_ENS8.nc
2026-07-02 14:29:15 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-02 14:29:15 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:29:15 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-02 14:29:15 INFO Queuing job for member 8...
2026-07-02 14:29:15 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:29:15 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-02 14:29:16 INFO Found: ['5066561']
2026-07-02 14:29:21 INFO [TGCC-IRENE] Submitted job with ID:['5066561']
2026-07-02 14:29:21 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:29:21 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9
2026-07-02 14:29:21 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020712_2_ENS9.nc
2026-07-02 14:29:21 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-02 14:29:21 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:29:21 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-02 14:29:21 INFO Queuing job for member 9...
2026-07-02 14:29:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:29:21 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-02 14:29:22 INFO Found: ['5066562']
2026-07-02 14:29:27 INFO [TGCC-IRENE] Submitted job with ID:['5066562']
2026-07-02 14:29:27 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:29:27 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10
2026-07-02 14:29:27 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020712_2_ENS10.nc
2026-07-02 14:29:27 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-02 14:29:27 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:29:27 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-02 14:29:27 INFO Queuing job for member 10...
2026-07-02 14:29:27 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:29:27 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-02 14:29:28 INFO Found: ['5066563']
2026-07-02 14:29:33 INFO [TGCC-IRENE] Submitted job with ID:['5066563']
2026-07-02 14:29:33 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:29:33 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11
2026-07-02 14:29:33 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020712_2_ENS11.nc
2026-07-02 14:29:33 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-02 14:29:33 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:29:33 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-02 14:29:33 INFO Queuing job for member 11...
2026-07-02 14:29:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:29:33 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-02 14:29:33 INFO Found: ['5066564']
2026-07-02 14:29:38 INFO [TGCC-IRENE] Submitted job with ID:['5066564']
2026-07-02 14:29:38 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:29:38 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12
2026-07-02 14:29:38 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020712_2_ENS12.nc
2026-07-02 14:29:38 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-02 14:29:38 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:29:38 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-02 14:29:38 INFO Queuing job for member 12...
2026-07-02 14:29:38 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:29:38 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-02 14:29:39 INFO Found: ['5066565']
2026-07-02 14:29:44 INFO [TGCC-IRENE] Submitted job with ID:['5066565']
2026-07-02 14:29:44 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:29:44 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13
2026-07-02 14:29:44 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020712_2_ENS13.nc
2026-07-02 14:29:44 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-02 14:29:44 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:29:44 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-02 14:29:44 INFO Queuing job for member 13...
2026-07-02 14:29:44 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:29:44 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-02 14:29:45 INFO Found: ['5066566']
2026-07-02 14:29:50 INFO [TGCC-IRENE] Submitted job with ID:['5066566']
2026-07-02 14:29:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:29:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14
2026-07-02 14:29:50 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020712_2_ENS14.nc
2026-07-02 14:29:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-02 14:29:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:29:50 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-02 14:29:50 INFO Queuing job for member 14...
2026-07-02 14:29:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:29:50 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-02 14:29:51 INFO Found: ['5066567']
2026-07-02 14:29:56 INFO [TGCC-IRENE] Submitted job with ID:['5066567']
2026-07-02 14:29:56 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-02 14:29:56 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15
2026-07-02 14:29:56 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020712_2_ENS15.nc
2026-07-02 14:29:56 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-02 14:29:56 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-02 14:29:56 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-02 14:29:56 INFO Queuing job for member 15...
2026-07-02 14:29:56 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-02 14:29:56 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-02 14:29:56 INFO Found: ['5066568']
2026-07-02 14:30:01 INFO [TGCC-IRENE] Submitted job with ID:['5066568']
2026-07-02 14:30:01 INFO Checking job status ...
2026-07-02 14:30:03 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:30:03 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:30:03 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:30:03 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:30:03 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:30:03 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:30:03 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:30:03 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:30:03 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:30:03 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:30:03 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:30:03 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:30:03 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:30:03 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:30:03 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:30:03 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:30:18 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:30:18 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:30:18 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:30:18 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:30:18 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:30:18 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:30:18 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:30:18 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:30:18 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:30:18 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:30:18 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:30:18 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:30:18 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:30:18 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:30:18 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:30:18 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:30:34 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:30:34 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:30:34 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:30:34 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:30:34 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:30:34 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:30:34 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:30:34 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:30:34 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:30:34 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:30:34 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:30:34 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:30:34 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:30:34 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:30:34 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:30:34 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:30:49 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:30:49 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:30:49 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:30:49 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:30:49 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:30:49 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:30:49 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:30:49 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:30:49 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:30:49 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:30:49 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:30:49 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:30:49 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:30:49 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:30:49 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:30:49 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:31:04 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:31:04 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:31:04 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:31:04 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:31:04 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:31:04 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:31:04 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:31:04 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:31:04 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:31:04 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:31:04 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:31:04 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:31:04 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:31:04 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:31:05 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:31:05 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:31:20 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:31:20 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:31:20 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:31:20 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:31:20 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:31:20 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:31:20 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:31:20 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:31:20 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:31:20 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:31:20 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:31:20 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:31:20 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:31:20 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:31:20 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:31:20 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:31:35 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:31:35 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:31:35 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:31:35 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:31:35 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:31:35 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:31:35 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:31:35 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:31:35 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:31:35 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:31:35 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:31:35 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:31:35 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:31:35 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:31:35 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:31:35 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:31:50 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:31:50 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:31:50 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:31:50 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:31:50 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:31:50 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:31:50 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:31:50 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:31:50 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:31:50 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:31:50 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:31:50 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:31:50 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:31:50 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:31:50 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:31:50 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:32:06 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:32:06 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:32:06 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:32:06 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:32:06 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:32:06 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:32:06 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:32:06 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:32:06 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:32:06 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:32:06 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:32:06 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:32:06 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:32:06 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:32:06 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:32:06 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:32:21 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:32:21 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:32:21 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:32:21 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:32:21 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:32:21 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:32:21 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:32:21 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:32:21 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:32:21 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:32:21 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:32:21 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:32:21 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:32:21 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:32:21 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:32:21 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:32:36 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:32:36 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:32:36 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:32:36 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:32:36 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:32:36 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:32:36 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:32:36 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:32:36 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:32:36 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:32:36 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:32:36 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:32:36 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:32:36 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:32:36 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:32:36 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:32:53 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:32:53 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:32:53 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:32:53 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:32:53 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:32:53 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:32:53 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:32:53 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:32:53 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:32:53 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:32:53 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:32:53 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:32:53 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:32:53 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:32:53 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:32:53 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:33:08 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:33:08 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:33:08 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:33:08 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:33:08 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:33:08 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:33:08 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:33:08 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:33:08 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:33:08 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:33:08 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:33:08 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:33:08 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:33:08 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:33:08 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:33:08 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:33:23 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:33:23 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:33:23 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:33:23 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:33:23 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:33:23 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:33:23 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:33:23 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:33:23 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:33:23 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:33:23 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:33:23 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:33:23 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:33:23 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:33:23 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:33:23 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:33:38 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:33:39 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:33:39 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:33:39 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:33:39 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:33:39 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:33:39 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:33:39 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:33:39 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:33:39 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:33:39 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:33:39 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:33:39 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:33:39 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:33:39 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:33:39 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:33:55 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:33:55 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:33:55 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:33:55 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:33:55 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:33:55 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:33:55 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:33:55 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:33:55 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:33:55 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:33:55 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:33:55 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:33:55 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:33:55 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:33:55 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:33:55 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:34:10 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:34:10 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:34:10 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:34:10 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:34:10 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:34:10 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:34:10 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:34:10 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:34:10 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:34:10 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:34:10 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:34:10 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:34:10 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:34:10 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:34:10 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:34:10 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:34:25 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:34:25 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:34:25 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:34:25 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:34:25 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:34:25 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:34:25 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:34:25 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:34:25 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:34:25 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:34:25 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:34:25 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:34:25 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:34:25 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:34:25 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:34:25 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:34:40 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:34:40 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:34:41 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:34:41 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:34:41 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:34:41 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:34:41 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:34:41 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:34:41 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:34:41 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:34:41 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:34:41 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:34:41 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:34:41 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:34:41 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:34:41 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:34:56 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:34:56 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:34:56 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:34:56 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:34:56 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:34:56 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:34:56 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:34:56 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:34:56 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:34:56 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:34:56 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:34:56 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:34:56 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:34:56 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:34:56 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:34:56 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:35:11 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:35:11 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:35:11 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:35:11 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:35:11 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:35:11 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:35:11 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:35:11 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:35:11 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:35:11 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:35:11 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:35:11 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:35:11 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:35:11 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:35:11 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:35:11 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:35:26 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:35:26 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:35:26 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:35:27 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:35:27 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:35:27 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:35:27 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:35:27 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:35:27 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:35:27 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:35:27 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:35:27 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:35:27 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:35:27 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:35:27 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:35:27 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:35:42 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:35:42 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:35:42 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:35:42 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:35:42 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:35:42 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:35:42 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:35:42 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:35:42 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:35:42 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:35:42 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:35:42 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:35:42 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:35:42 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:35:42 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:35:42 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:35:57 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:35:57 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:35:57 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:35:57 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:35:57 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:35:57 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:35:57 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:35:57 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:35:57 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:35:57 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:35:57 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:35:57 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:35:57 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:35:57 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:35:57 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:35:57 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:36:12 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:36:12 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:36:12 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:36:12 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:36:12 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:36:12 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:36:12 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:36:12 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:36:12 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:36:12 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:36:12 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:36:12 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:36:13 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:36:13 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:36:13 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:36:13 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:36:28 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:36:28 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:36:28 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:36:28 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:36:28 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:36:28 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:36:28 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:36:28 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:36:28 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:36:28 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:36:28 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:36:28 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:36:28 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:36:28 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:36:28 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:36:28 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:36:43 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:36:43 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:36:43 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:36:43 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:36:43 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:36:43 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:36:43 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:36:44 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:36:44 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:36:44 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:36:44 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:36:44 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:36:44 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:36:44 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:36:44 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:36:44 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:36:59 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:36:59 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:36:59 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:36:59 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:36:59 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:36:59 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:36:59 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:36:59 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:36:59 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:36:59 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:36:59 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:36:59 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:36:59 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:36:59 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:36:59 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:36:59 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:37:14 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:37:14 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:37:14 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:37:14 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:37:14 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:37:14 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:37:14 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:37:14 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:37:14 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:37:14 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:37:14 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:37:14 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:37:14 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:37:14 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:37:14 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:37:14 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:37:29 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:37:29 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:37:29 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:37:29 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:37:29 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:37:29 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:37:29 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:37:29 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:37:29 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:37:29 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:37:29 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:37:29 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:37:29 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:37:29 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:37:29 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:37:29 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:37:45 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:37:45 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:37:45 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:37:45 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:37:45 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:37:45 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:37:45 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:37:45 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:37:45 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:37:45 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:37:45 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:37:45 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:37:45 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:37:45 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:37:45 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:37:45 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:38:00 INFO None 5066551: status RUNNING/PENDING
2026-07-02 14:38:00 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:38:00 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:38:00 INFO None 5066556: status RUNNING/PENDING
2026-07-02 14:38:00 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:38:00 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:38:00 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:38:00 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:38:00 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:38:00 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:38:00 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:38:00 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:38:01 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:38:01 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:38:01 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:38:01 INFO Jobs still running: ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:38:16 INFO None 5066551: status FINISHED
2026-07-02 14:38:16 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:38:16 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:38:16 INFO None 5066556: status FINISHED
2026-07-02 14:38:16 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:38:16 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:38:16 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:38:16 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:38:16 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:38:16 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:38:16 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:38:16 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:38:16 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:38:16 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:38:16 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:38:16 INFO Jobs still running: ['5066554', '5066555', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:38:31 INFO None 5066551: status FINISHED
2026-07-02 14:38:31 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:38:31 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:38:31 INFO None 5066556: status FINISHED
2026-07-02 14:38:31 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:38:31 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:38:31 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:38:31 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:38:31 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:38:31 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:38:31 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:38:31 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:38:31 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:38:31 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:38:31 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:38:31 INFO Jobs still running: ['5066554', '5066555', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:38:48 INFO None 5066551: status FINISHED
2026-07-02 14:38:48 INFO None 5066554: status RUNNING/PENDING
2026-07-02 14:38:48 INFO None 5066555: status RUNNING/PENDING
2026-07-02 14:38:48 INFO None 5066556: status FINISHED
2026-07-02 14:38:48 INFO None 5066557: status RUNNING/PENDING
2026-07-02 14:38:48 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:38:48 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:38:48 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:38:48 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:38:48 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:38:48 INFO None 5066564: status RUNNING/PENDING
2026-07-02 14:38:48 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:38:48 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:38:48 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:38:48 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:38:48 INFO Jobs still running: ['5066554', '5066555', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:39:03 INFO None 5066551: status FINISHED
2026-07-02 14:39:03 INFO None 5066554: status FINISHED
2026-07-02 14:39:03 INFO None 5066555: status FINISHED
2026-07-02 14:39:03 INFO None 5066556: status FINISHED
2026-07-02 14:39:03 INFO None 5066557: status FINISHED
2026-07-02 14:39:03 INFO None 5066559: status RUNNING/PENDING
2026-07-02 14:39:03 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:39:03 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:39:03 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:39:03 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:39:03 INFO None 5066564: status FINISHED
2026-07-02 14:39:04 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:39:04 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:39:04 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:39:04 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:39:04 INFO Jobs still running: ['5066559', '5066560', '5066561', '5066562', '5066563', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:39:19 INFO None 5066551: status FINISHED
2026-07-02 14:39:19 INFO None 5066554: status FINISHED
2026-07-02 14:39:19 INFO None 5066555: status FINISHED
2026-07-02 14:39:19 INFO None 5066556: status FINISHED
2026-07-02 14:39:19 INFO None 5066557: status FINISHED
2026-07-02 14:39:19 INFO None 5066559: status FINISHED
2026-07-02 14:39:19 INFO None 5066560: status RUNNING/PENDING
2026-07-02 14:39:19 INFO None 5066561: status RUNNING/PENDING
2026-07-02 14:39:19 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:39:19 INFO None 5066563: status RUNNING/PENDING
2026-07-02 14:39:19 INFO None 5066564: status FINISHED
2026-07-02 14:39:19 INFO None 5066565: status RUNNING/PENDING
2026-07-02 14:39:19 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:39:19 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:39:19 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:39:19 INFO Jobs still running: ['5066560', '5066561', '5066562', '5066563', '5066565', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:39:34 INFO None 5066551: status FINISHED
2026-07-02 14:39:34 INFO None 5066554: status FINISHED
2026-07-02 14:39:34 INFO None 5066555: status FINISHED
2026-07-02 14:39:34 INFO None 5066556: status FINISHED
2026-07-02 14:39:34 INFO None 5066557: status FINISHED
2026-07-02 14:39:34 INFO None 5066559: status FINISHED
2026-07-02 14:39:34 INFO None 5066560: status FINISHED
2026-07-02 14:39:34 INFO None 5066561: status FINISHED
2026-07-02 14:39:34 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:39:34 INFO None 5066563: status FINISHED
2026-07-02 14:39:34 INFO None 5066564: status FINISHED
2026-07-02 14:39:34 INFO None 5066565: status FINISHED
2026-07-02 14:39:34 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:39:34 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:39:34 INFO None 5066568: status RUNNING/PENDING
2026-07-02 14:39:34 INFO Jobs still running: ['5066562', '5066566', '5066567', '5066568']. Waiting...
2026-07-02 14:39:50 INFO None 5066551: status FINISHED
2026-07-02 14:39:50 INFO None 5066554: status FINISHED
2026-07-02 14:39:50 INFO None 5066555: status FINISHED
2026-07-02 14:39:50 INFO None 5066556: status FINISHED
2026-07-02 14:39:50 INFO None 5066557: status FINISHED
2026-07-02 14:39:50 INFO None 5066559: status FINISHED
2026-07-02 14:39:50 INFO None 5066560: status FINISHED
2026-07-02 14:39:50 INFO None 5066561: status FINISHED
2026-07-02 14:39:50 INFO None 5066562: status RUNNING/PENDING
2026-07-02 14:39:50 INFO None 5066563: status FINISHED
2026-07-02 14:39:50 INFO None 5066564: status FINISHED
2026-07-02 14:39:50 INFO None 5066565: status FINISHED
2026-07-02 14:39:50 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:39:50 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:39:50 INFO None 5066568: status FINISHED
2026-07-02 14:39:50 INFO Jobs still running: ['5066562', '5066566', '5066567']. Waiting...
2026-07-02 14:40:05 INFO None 5066551: status FINISHED
2026-07-02 14:40:05 INFO None 5066554: status FINISHED
2026-07-02 14:40:05 INFO None 5066555: status FINISHED
2026-07-02 14:40:05 INFO None 5066556: status FINISHED
2026-07-02 14:40:05 INFO None 5066557: status FINISHED
2026-07-02 14:40:05 INFO None 5066559: status FINISHED
2026-07-02 14:40:05 INFO None 5066560: status FINISHED
2026-07-02 14:40:05 INFO None 5066561: status FINISHED
2026-07-02 14:40:05 INFO None 5066562: status FINISHED
2026-07-02 14:40:05 INFO None 5066563: status FINISHED
2026-07-02 14:40:05 INFO None 5066564: status FINISHED
2026-07-02 14:40:05 INFO None 5066565: status FINISHED
2026-07-02 14:40:06 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:40:06 INFO None 5066567: status RUNNING/PENDING
2026-07-02 14:40:06 INFO None 5066568: status FINISHED
2026-07-02 14:40:06 INFO Jobs still running: ['5066566', '5066567']. Waiting...
2026-07-02 14:40:21 INFO None 5066551: status FINISHED
2026-07-02 14:40:21 INFO None 5066554: status FINISHED
2026-07-02 14:40:21 INFO None 5066555: status FINISHED
2026-07-02 14:40:21 INFO None 5066556: status FINISHED
2026-07-02 14:40:21 INFO None 5066557: status FINISHED
2026-07-02 14:40:21 INFO None 5066559: status FINISHED
2026-07-02 14:40:21 INFO None 5066560: status FINISHED
2026-07-02 14:40:21 INFO None 5066561: status FINISHED
2026-07-02 14:40:21 INFO None 5066562: status FINISHED
2026-07-02 14:40:21 INFO None 5066563: status FINISHED
2026-07-02 14:40:21 INFO None 5066564: status FINISHED
2026-07-02 14:40:21 INFO None 5066565: status FINISHED
2026-07-02 14:40:21 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:40:21 INFO None 5066567: status FINISHED
2026-07-02 14:40:21 INFO None 5066568: status FINISHED
2026-07-02 14:40:21 INFO Jobs still running: ['5066566']. Waiting...
2026-07-02 14:40:36 INFO None 5066551: status FINISHED
2026-07-02 14:40:36 INFO None 5066554: status FINISHED
2026-07-02 14:40:36 INFO None 5066555: status FINISHED
2026-07-02 14:40:36 INFO None 5066556: status FINISHED
2026-07-02 14:40:36 INFO None 5066557: status FINISHED
2026-07-02 14:40:36 INFO None 5066559: status FINISHED
2026-07-02 14:40:36 INFO None 5066560: status FINISHED
2026-07-02 14:40:36 INFO None 5066561: status FINISHED
2026-07-02 14:40:36 INFO None 5066562: status FINISHED
2026-07-02 14:40:36 INFO None 5066563: status FINISHED
2026-07-02 14:40:36 INFO None 5066564: status FINISHED
2026-07-02 14:40:36 INFO None 5066565: status FINISHED
2026-07-02 14:40:36 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:40:36 INFO None 5066567: status FINISHED
2026-07-02 14:40:36 INFO None 5066568: status FINISHED
2026-07-02 14:40:36 INFO Jobs still running: ['5066566']. Waiting...
2026-07-02 14:40:51 INFO None 5066551: status FINISHED
2026-07-02 14:40:51 INFO None 5066554: status FINISHED
2026-07-02 14:40:51 INFO None 5066555: status FINISHED
2026-07-02 14:40:51 INFO None 5066556: status FINISHED
2026-07-02 14:40:51 INFO None 5066557: status FINISHED
2026-07-02 14:40:51 INFO None 5066559: status FINISHED
2026-07-02 14:40:52 INFO None 5066560: status FINISHED
2026-07-02 14:40:52 INFO None 5066561: status FINISHED
2026-07-02 14:40:52 INFO None 5066562: status FINISHED
2026-07-02 14:40:52 INFO None 5066563: status FINISHED
2026-07-02 14:40:52 INFO None 5066564: status FINISHED
2026-07-02 14:40:52 INFO None 5066565: status FINISHED
2026-07-02 14:40:52 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:40:52 INFO None 5066567: status FINISHED
2026-07-02 14:40:52 INFO None 5066568: status FINISHED
2026-07-02 14:40:52 INFO Jobs still running: ['5066566']. Waiting...
2026-07-02 14:41:07 INFO None 5066551: status FINISHED
2026-07-02 14:41:07 INFO None 5066554: status FINISHED
2026-07-02 14:41:07 INFO None 5066555: status FINISHED
2026-07-02 14:41:07 INFO None 5066556: status FINISHED
2026-07-02 14:41:07 INFO None 5066557: status FINISHED
2026-07-02 14:41:07 INFO None 5066559: status FINISHED
2026-07-02 14:41:07 INFO None 5066560: status FINISHED
2026-07-02 14:41:07 INFO None 5066561: status FINISHED
2026-07-02 14:41:07 INFO None 5066562: status FINISHED
2026-07-02 14:41:07 INFO None 5066563: status FINISHED
2026-07-02 14:41:07 INFO None 5066564: status FINISHED
2026-07-02 14:41:07 INFO None 5066565: status FINISHED
2026-07-02 14:41:07 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:41:07 INFO None 5066567: status FINISHED
2026-07-02 14:41:07 INFO None 5066568: status FINISHED
2026-07-02 14:41:07 INFO Jobs still running: ['5066566']. Waiting...
2026-07-02 14:41:22 INFO None 5066551: status FINISHED
2026-07-02 14:41:22 INFO None 5066554: status FINISHED
2026-07-02 14:41:22 INFO None 5066555: status FINISHED
2026-07-02 14:41:22 INFO None 5066556: status FINISHED
2026-07-02 14:41:22 INFO None 5066557: status FINISHED
2026-07-02 14:41:22 INFO None 5066559: status FINISHED
2026-07-02 14:41:22 INFO None 5066560: status FINISHED
2026-07-02 14:41:22 INFO None 5066561: status FINISHED
2026-07-02 14:41:22 INFO None 5066562: status FINISHED
2026-07-02 14:41:22 INFO None 5066563: status FINISHED
2026-07-02 14:41:22 INFO None 5066564: status FINISHED
2026-07-02 14:41:22 INFO None 5066565: status FINISHED
2026-07-02 14:41:22 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:41:22 INFO None 5066567: status FINISHED
2026-07-02 14:41:22 INFO None 5066568: status FINISHED
2026-07-02 14:41:22 INFO Jobs still running: ['5066566']. Waiting...
2026-07-02 14:41:37 INFO None 5066551: status FINISHED
2026-07-02 14:41:37 INFO None 5066554: status FINISHED
2026-07-02 14:41:37 INFO None 5066555: status FINISHED
2026-07-02 14:41:37 INFO None 5066556: status FINISHED
2026-07-02 14:41:37 INFO None 5066557: status FINISHED
2026-07-02 14:41:37 INFO None 5066559: status FINISHED
2026-07-02 14:41:37 INFO None 5066560: status FINISHED
2026-07-02 14:41:37 INFO None 5066561: status FINISHED
2026-07-02 14:41:37 INFO None 5066562: status FINISHED
2026-07-02 14:41:37 INFO None 5066563: status FINISHED
2026-07-02 14:41:37 INFO None 5066564: status FINISHED
2026-07-02 14:41:37 INFO None 5066565: status FINISHED
2026-07-02 14:41:37 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:41:37 INFO None 5066567: status FINISHED
2026-07-02 14:41:37 INFO None 5066568: status FINISHED
2026-07-02 14:41:37 INFO Jobs still running: ['5066566']. Waiting...
2026-07-02 14:41:54 INFO None 5066551: status FINISHED
2026-07-02 14:41:54 INFO None 5066554: status FINISHED
2026-07-02 14:41:54 INFO None 5066555: status FINISHED
2026-07-02 14:41:54 INFO None 5066556: status FINISHED
2026-07-02 14:41:54 INFO None 5066557: status FINISHED
2026-07-02 14:41:54 INFO None 5066559: status FINISHED
2026-07-02 14:41:54 INFO None 5066560: status FINISHED
2026-07-02 14:41:54 INFO None 5066561: status FINISHED
2026-07-02 14:41:54 INFO None 5066562: status FINISHED
2026-07-02 14:41:54 INFO None 5066563: status FINISHED
2026-07-02 14:41:54 INFO None 5066564: status FINISHED
2026-07-02 14:41:54 INFO None 5066565: status FINISHED
2026-07-02 14:41:54 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:41:54 INFO None 5066567: status FINISHED
2026-07-02 14:41:54 INFO None 5066568: status FINISHED
2026-07-02 14:41:54 INFO Jobs still running: ['5066566']. Waiting...
2026-07-02 14:42:09 INFO None 5066551: status FINISHED
2026-07-02 14:42:09 INFO None 5066554: status FINISHED
2026-07-02 14:42:09 INFO None 5066555: status FINISHED
2026-07-02 14:42:09 INFO None 5066556: status FINISHED
2026-07-02 14:42:09 INFO None 5066557: status FINISHED
2026-07-02 14:42:09 INFO None 5066559: status FINISHED
2026-07-02 14:42:09 INFO None 5066560: status FINISHED
2026-07-02 14:42:09 INFO None 5066561: status FINISHED
2026-07-02 14:42:09 INFO None 5066562: status FINISHED
2026-07-02 14:42:09 INFO None 5066563: status FINISHED
2026-07-02 14:42:09 INFO None 5066564: status FINISHED
2026-07-02 14:42:09 INFO None 5066565: status FINISHED
2026-07-02 14:42:10 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:42:10 INFO None 5066567: status FINISHED
2026-07-02 14:42:10 INFO None 5066568: status FINISHED
2026-07-02 14:42:10 INFO Jobs still running: ['5066566']. Waiting...
2026-07-02 14:42:25 INFO None 5066551: status FINISHED
2026-07-02 14:42:25 INFO None 5066554: status FINISHED
2026-07-02 14:42:25 INFO None 5066555: status FINISHED
2026-07-02 14:42:25 INFO None 5066556: status FINISHED
2026-07-02 14:42:25 INFO None 5066557: status FINISHED
2026-07-02 14:42:25 INFO None 5066559: status FINISHED
2026-07-02 14:42:25 INFO None 5066560: status FINISHED
2026-07-02 14:42:25 INFO None 5066561: status FINISHED
2026-07-02 14:42:25 INFO None 5066562: status FINISHED
2026-07-02 14:42:25 INFO None 5066563: status FINISHED
2026-07-02 14:42:25 INFO None 5066564: status FINISHED
2026-07-02 14:42:25 INFO None 5066565: status FINISHED
2026-07-02 14:42:25 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:42:25 INFO None 5066567: status FINISHED
2026-07-02 14:42:25 INFO None 5066568: status FINISHED
2026-07-02 14:42:25 INFO Jobs still running: ['5066566']. Waiting...
2026-07-02 14:42:40 INFO None 5066551: status FINISHED
2026-07-02 14:42:42 INFO None 5066554: status FINISHED
2026-07-02 14:42:42 INFO None 5066555: status FINISHED
2026-07-02 14:42:42 INFO None 5066556: status FINISHED
2026-07-02 14:42:42 INFO None 5066557: status FINISHED
2026-07-02 14:42:42 INFO None 5066559: status FINISHED
2026-07-02 14:42:42 INFO None 5066560: status FINISHED
2026-07-02 14:42:42 INFO None 5066561: status FINISHED
2026-07-02 14:42:42 INFO None 5066562: status FINISHED
2026-07-02 14:42:42 INFO None 5066563: status FINISHED
2026-07-02 14:42:42 INFO None 5066564: status FINISHED
2026-07-02 14:42:42 INFO None 5066565: status FINISHED
2026-07-02 14:42:42 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:42:42 INFO None 5066567: status FINISHED
2026-07-02 14:42:42 INFO None 5066568: status FINISHED
2026-07-02 14:42:42 INFO Jobs still running: ['5066566']. Waiting...
2026-07-02 14:42:57 INFO None 5066551: status FINISHED
2026-07-02 14:42:57 INFO None 5066554: status FINISHED
2026-07-02 14:42:57 INFO None 5066555: status FINISHED
2026-07-02 14:42:57 INFO None 5066556: status FINISHED
2026-07-02 14:42:57 INFO None 5066557: status FINISHED
2026-07-02 14:42:57 INFO None 5066559: status FINISHED
2026-07-02 14:42:57 INFO None 5066560: status FINISHED
2026-07-02 14:42:57 INFO None 5066561: status FINISHED
2026-07-02 14:42:57 INFO None 5066562: status FINISHED
2026-07-02 14:42:57 INFO None 5066563: status FINISHED
2026-07-02 14:42:57 INFO None 5066564: status FINISHED
2026-07-02 14:42:57 INFO None 5066565: status FINISHED
2026-07-02 14:42:57 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:42:57 INFO None 5066567: status FINISHED
2026-07-02 14:42:57 INFO None 5066568: status FINISHED
2026-07-02 14:42:57 INFO Jobs still running: ['5066566']. Waiting...
2026-07-02 14:43:12 INFO None 5066551: status FINISHED
2026-07-02 14:43:12 INFO None 5066554: status FINISHED
2026-07-02 14:43:12 INFO None 5066555: status FINISHED
2026-07-02 14:43:12 INFO None 5066556: status FINISHED
2026-07-02 14:43:12 INFO None 5066557: status FINISHED
2026-07-02 14:43:13 INFO None 5066559: status FINISHED
2026-07-02 14:43:13 INFO None 5066560: status FINISHED
2026-07-02 14:43:13 INFO None 5066561: status FINISHED
2026-07-02 14:43:13 INFO None 5066562: status FINISHED
2026-07-02 14:43:13 INFO None 5066563: status FINISHED
2026-07-02 14:43:13 INFO None 5066564: status FINISHED
2026-07-02 14:43:13 INFO None 5066565: status FINISHED
2026-07-02 14:43:13 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:43:13 INFO None 5066567: status FINISHED
2026-07-02 14:43:13 INFO None 5066568: status FINISHED
2026-07-02 14:43:13 INFO Jobs still running: ['5066566']. Waiting...
2026-07-02 14:43:28 INFO None 5066551: status FINISHED
2026-07-02 14:43:28 INFO None 5066554: status FINISHED
2026-07-02 14:43:28 INFO None 5066555: status FINISHED
2026-07-02 14:43:28 INFO None 5066556: status FINISHED
2026-07-02 14:43:28 INFO None 5066557: status FINISHED
2026-07-02 14:43:28 INFO None 5066559: status FINISHED
2026-07-02 14:43:28 INFO None 5066560: status FINISHED
2026-07-02 14:43:28 INFO None 5066561: status FINISHED
2026-07-02 14:43:28 INFO None 5066562: status FINISHED
2026-07-02 14:43:28 INFO None 5066563: status FINISHED
2026-07-02 14:43:28 INFO None 5066564: status FINISHED
2026-07-02 14:43:28 INFO None 5066565: status FINISHED
2026-07-02 14:43:28 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:43:28 INFO None 5066567: status FINISHED
2026-07-02 14:43:28 INFO None 5066568: status FINISHED
2026-07-02 14:43:28 INFO Jobs still running: ['5066566']. Waiting...
2026-07-02 14:43:44 INFO None 5066551: status FINISHED
2026-07-02 14:43:44 INFO None 5066554: status FINISHED
2026-07-02 14:43:44 INFO None 5066555: status FINISHED
2026-07-02 14:43:44 INFO None 5066556: status FINISHED
2026-07-02 14:43:44 INFO None 5066557: status FINISHED
2026-07-02 14:43:44 INFO None 5066559: status FINISHED
2026-07-02 14:43:44 INFO None 5066560: status FINISHED
2026-07-02 14:43:44 INFO None 5066561: status FINISHED
2026-07-02 14:43:44 INFO None 5066562: status FINISHED
2026-07-02 14:43:44 INFO None 5066563: status FINISHED
2026-07-02 14:43:44 INFO None 5066564: status FINISHED
2026-07-02 14:43:44 INFO None 5066565: status FINISHED
2026-07-02 14:43:44 INFO None 5066566: status RUNNING/PENDING
2026-07-02 14:43:44 INFO None 5066567: status FINISHED
2026-07-02 14:43:44 INFO None 5066568: status FINISHED
2026-07-02 14:43:44 INFO Jobs still running: ['5066566']. Waiting...
2026-07-02 14:43:59 INFO None 5066551: status FINISHED
2026-07-02 14:43:59 INFO None 5066554: status FINISHED
2026-07-02 14:43:59 INFO None 5066555: status FINISHED
2026-07-02 14:43:59 INFO None 5066556: status FINISHED
2026-07-02 14:43:59 INFO None 5066557: status FINISHED
2026-07-02 14:43:59 INFO None 5066559: status FINISHED
2026-07-02 14:43:59 INFO None 5066560: status FINISHED
2026-07-02 14:43:59 INFO None 5066561: status FINISHED
2026-07-02 14:43:59 INFO None 5066562: status FINISHED
2026-07-02 14:43:59 INFO None 5066563: status FINISHED
2026-07-02 14:43:59 INFO None 5066564: status FINISHED
2026-07-02 14:43:59 INFO None 5066565: status FINISHED
2026-07-02 14:43:59 INFO None 5066566: status FINISHED
2026-07-02 14:43:59 INFO None 5066567: status FINISHED
2026-07-02 14:43:59 INFO None 5066568: status FINISHED
2026-07-02 14:43:59 INFO Jobs ['5066551', '5066554', '5066555', '5066556', '5066557', '5066559', '5066560', '5066561', '5066562', '5066563', '5066564', '5066565', '5066566', '5066567', '5066568'] have finished
2026-07-02 14:43:59 INFO Checking restart files were created ...
2026-07-02 14:43:59 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020714_10_ENS1.nc(3673513755 bytes)
2026-07-02 14:43:59 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020714_10_ENS2.nc(3673513755 bytes)
2026-07-02 14:43:59 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020714_10_ENS3.nc(3673513755 bytes)
2026-07-02 14:43:59 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020714_10_ENS4.nc(3673513755 bytes)
2026-07-02 14:43:59 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020714_10_ENS5.nc(3673513755 bytes)
2026-07-02 14:43:59 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020714_10_ENS6.nc(3673513755 bytes)
2026-07-02 14:43:59 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020714_10_ENS7.nc(3673513755 bytes)
2026-07-02 14:43:59 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020714_10_ENS8.nc(3673513755 bytes)
2026-07-02 14:43:59 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020714_10_ENS9.nc(3673513755 bytes)
2026-07-02 14:43:59 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020714_10_ENS10.nc(3673513755 bytes)
2026-07-02 14:43:59 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020714_10_ENS11.nc(3673513755 bytes)
2026-07-02 14:43:59 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020714_10_ENS12.nc(3673513755 bytes)
2026-07-02 14:43:59 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020714_10_ENS13.nc(3673513755 bytes)
2026-07-02 14:43:59 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020714_10_ENS14.nc(3673513755 bytes)
2026-07-02 14:43:59 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020714_10_ENS15.nc(3673513755 bytes)
2026-07-02 14:43:59 INFO  Run_model() completed successfully.
2026-07-02 14:43:59 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 14:00:00 simulated_time=2020-02-08 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 14:43:59 INFO [TIME] gregorian_conversion simulated_time=2020-02-08 00:00:00 days=153074 seconds=0
2026-07-02 14:43:59 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-02 14:43:59 INFO [TIME] increment current_time 2020-02-07 14:00:00 -> 2020-02-08 00:00:00
2026-07-02 14:43:59 INFO [TIME] after_increment_before_assimilation current_time=2020-02-08 00:00:00 simulated_time=2020-02-08 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 14:43:59 INFO ---------->>> Running process_satellite_data()
2026-07-02 14:43:59 INFO [DART] No satellite data found, skipping assimilation
2026-07-02 14:43:59 INFO after_assimilation() skipped
2026-07-02 14:43:59 INFO Next run starts from 2020-02-08 00:00:00
2026-07-02 14:43:59 INFO Cycle is DONE; starting a new loop!
2026-07-02 14:43:59 INFO [TIME] step_end current_time=2020-02-08 00:00:00 simulated_time=2020-02-08 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-02 14:43:59 INFO [PIPELINE] ---- TIME LOOP END ----
+ exit 0
