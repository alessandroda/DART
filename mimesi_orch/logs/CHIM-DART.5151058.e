+ SCRIPT_PID=129283
+ /bin/bash -x /tmp/tmp.joJhZkw1PD
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
2026-07-13 18:28:28 INFO 
███    ███ ██ ███    ███ ███████ ███████ ██
████  ████ ██ ████  ████ ██      ██      ██
██ ████ ██ ██ ██ ████ ██ █████   ███████ ██
██  ██  ██ ██ ██  ██  ██ ██           ██ ██
██      ██ ██ ██      ██ ███████ ███████ ██



2026-07-13 18:28:28 INFO [PIPELINE] =======================================
2026-07-13 18:28:28 INFO [PIPELINE] Starting chimere–DART orchestrator
2026-07-13 18:28:28 INFO [PIPELINE] Config file: config/config_irene_IM.yaml
2026-07-13 18:28:28 INFO [PIPELINE] Run dir: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart
2026-07-13 18:28:28 INFO [PIPELINE] Log file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/mimesi_orchestrator_logs/chimere_DART_20260713_182828.log
2026-07-13 18:28:28 INFO [PIPELINE] =======================================
2026-07-13 18:28:28 INFO Running assimilation with model_type=ModelType.CHIMERE
2026-07-13 18:28:28 INFO Using scheduler=Scheduler.SLURM, queue=rome
2026-07-13 18:28:28 INFO [STEP] ---- TIME LOOP START ----
2026-07-13 18:28:28 INFO [TIME] step_start current_time=2020-02-06 00:00:00 simulated_time=None dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-13 18:28:28 INFO [TIME] window start=2020-02-06 00:00:00 end=2020-02-06 01:00:00 run_hours=1 has_assimilation=False
2026-07-13 18:28:28 INFO Creating directories and links for ENS1 to run chimere's parallel part
2026-07-13 18:28:28 INFO Copying EMIS ...
2026-07-13 18:28:29 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens01.nc
2026-07-13 18:28:29 INFO Linking END ...
2026-07-13 18:28:29 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS1/end.2020020500_24_ENS1.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-13 18:28:29 INFO >> Checking links...
2026-07-13 18:28:29 INFO >> All links are good for ENS1  ...
2026-07-13 18:28:29 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 18:28:35 INFO Hourly dataset computed and listing created
2026-07-13 18:28:39 INFO Hourly dataset computed
2026-07-13 18:28:39 INFO Creating directories and links for ENS2 to run chimere's parallel part
2026-07-13 18:28:39 INFO Copying EMIS ...
2026-07-13 18:28:39 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens02.nc
2026-07-13 18:28:39 INFO Linking END ...
2026-07-13 18:28:39 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS2/end.2020020500_24_ENS2.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-13 18:28:39 INFO >> Checking links...
2026-07-13 18:28:39 INFO >> All links are good for ENS2  ...
2026-07-13 18:28:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 18:28:40 INFO Hourly dataset computed and listing created
2026-07-13 18:28:41 INFO Hourly dataset computed
2026-07-13 18:28:41 INFO Creating directories and links for ENS3 to run chimere's parallel part
2026-07-13 18:28:41 INFO Copying EMIS ...
2026-07-13 18:28:41 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens03.nc
2026-07-13 18:28:41 INFO Linking END ...
2026-07-13 18:28:41 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS3/end.2020020500_24_ENS3.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-13 18:28:41 INFO >> Checking links...
2026-07-13 18:28:41 INFO >> All links are good for ENS3  ...
2026-07-13 18:28:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 18:28:42 INFO Hourly dataset computed and listing created
2026-07-13 18:28:42 INFO Hourly dataset computed
2026-07-13 18:28:42 INFO Creating directories and links for ENS4 to run chimere's parallel part
2026-07-13 18:28:42 INFO Copying EMIS ...
2026-07-13 18:28:43 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens04.nc
2026-07-13 18:28:43 INFO Linking END ...
2026-07-13 18:28:43 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS4/end.2020020500_24_ENS4.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-13 18:28:43 INFO >> Checking links...
2026-07-13 18:28:43 INFO >> All links are good for ENS4  ...
2026-07-13 18:28:43 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 18:28:44 INFO Hourly dataset computed and listing created
2026-07-13 18:28:44 INFO Hourly dataset computed
2026-07-13 18:28:44 INFO Creating directories and links for ENS5 to run chimere's parallel part
2026-07-13 18:28:44 INFO Copying EMIS ...
2026-07-13 18:28:45 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens05.nc
2026-07-13 18:28:45 INFO Linking END ...
2026-07-13 18:28:45 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS5/end.2020020500_24_ENS5.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-13 18:28:45 INFO >> Checking links...
2026-07-13 18:28:45 INFO >> All links are good for ENS5  ...
2026-07-13 18:28:45 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 18:28:45 INFO Hourly dataset computed and listing created
2026-07-13 18:28:46 INFO Hourly dataset computed
2026-07-13 18:28:46 INFO Creating directories and links for ENS6 to run chimere's parallel part
2026-07-13 18:28:46 INFO Copying EMIS ...
2026-07-13 18:28:46 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens06.nc
2026-07-13 18:28:46 INFO Linking END ...
2026-07-13 18:28:46 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS6/end.2020020500_24_ENS6.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-13 18:28:46 INFO >> Checking links...
2026-07-13 18:28:46 INFO >> All links are good for ENS6  ...
2026-07-13 18:28:46 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 18:28:47 INFO Hourly dataset computed and listing created
2026-07-13 18:28:48 INFO Hourly dataset computed
2026-07-13 18:28:48 INFO Creating directories and links for ENS7 to run chimere's parallel part
2026-07-13 18:28:48 INFO Copying EMIS ...
2026-07-13 18:28:48 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens07.nc
2026-07-13 18:28:48 INFO Linking END ...
2026-07-13 18:28:48 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS7/end.2020020500_24_ENS7.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-13 18:28:48 INFO >> Checking links...
2026-07-13 18:28:48 INFO >> All links are good for ENS7  ...
2026-07-13 18:28:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 18:28:49 INFO Hourly dataset computed and listing created
2026-07-13 18:28:50 INFO Hourly dataset computed
2026-07-13 18:28:50 INFO Creating directories and links for ENS8 to run chimere's parallel part
2026-07-13 18:28:50 INFO Copying EMIS ...
2026-07-13 18:28:50 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens08.nc
2026-07-13 18:28:50 INFO Linking END ...
2026-07-13 18:28:50 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS8/end.2020020500_24_ENS8.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-13 18:28:50 INFO >> Checking links...
2026-07-13 18:28:50 INFO >> All links are good for ENS8  ...
2026-07-13 18:28:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 18:28:51 INFO Hourly dataset computed and listing created
2026-07-13 18:28:55 INFO Hourly dataset computed
2026-07-13 18:28:55 INFO Creating directories and links for ENS9 to run chimere's parallel part
2026-07-13 18:28:55 INFO Copying EMIS ...
2026-07-13 18:28:55 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens09.nc
2026-07-13 18:28:55 INFO Linking END ...
2026-07-13 18:28:55 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS9/end.2020020500_24_ENS9.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-13 18:28:55 INFO >> Checking links...
2026-07-13 18:28:55 INFO >> All links are good for ENS9  ...
2026-07-13 18:28:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 18:28:56 INFO Hourly dataset computed and listing created
2026-07-13 18:28:56 INFO Hourly dataset computed
2026-07-13 18:28:56 INFO Creating directories and links for ENS10 to run chimere's parallel part
2026-07-13 18:28:56 INFO Copying EMIS ...
2026-07-13 18:28:57 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens10.nc
2026-07-13 18:28:57 INFO Linking END ...
2026-07-13 18:28:57 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS10/end.2020020500_24_ENS10.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-13 18:28:57 INFO >> Checking links...
2026-07-13 18:28:57 INFO >> All links are good for ENS10  ...
2026-07-13 18:28:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 18:28:58 INFO Hourly dataset computed and listing created
2026-07-13 18:28:58 INFO Hourly dataset computed
2026-07-13 18:28:58 INFO Creating directories and links for ENS11 to run chimere's parallel part
2026-07-13 18:28:58 INFO Copying EMIS ...
2026-07-13 18:28:59 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens11.nc
2026-07-13 18:28:59 INFO Linking END ...
2026-07-13 18:28:59 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS11/end.2020020500_24_ENS11.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-13 18:28:59 INFO >> Checking links...
2026-07-13 18:28:59 INFO >> All links are good for ENS11  ...
2026-07-13 18:28:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 18:29:00 INFO Hourly dataset computed and listing created
2026-07-13 18:29:00 INFO Hourly dataset computed
2026-07-13 18:29:00 INFO Creating directories and links for ENS12 to run chimere's parallel part
2026-07-13 18:29:00 INFO Copying EMIS ...
2026-07-13 18:29:01 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens12.nc
2026-07-13 18:29:01 INFO Linking END ...
2026-07-13 18:29:01 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS12/end.2020020500_24_ENS12.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-13 18:29:01 INFO >> Checking links...
2026-07-13 18:29:01 INFO >> All links are good for ENS12  ...
2026-07-13 18:29:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 18:29:01 INFO Hourly dataset computed and listing created
2026-07-13 18:29:02 INFO Hourly dataset computed
2026-07-13 18:29:02 INFO Creating directories and links for ENS13 to run chimere's parallel part
2026-07-13 18:29:02 INFO Copying EMIS ...
2026-07-13 18:29:02 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens13.nc
2026-07-13 18:29:02 INFO Linking END ...
2026-07-13 18:29:03 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS13/end.2020020500_24_ENS13.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-13 18:29:03 INFO >> Checking links...
2026-07-13 18:29:03 INFO >> All links are good for ENS13  ...
2026-07-13 18:29:03 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 18:29:03 INFO Hourly dataset computed and listing created
2026-07-13 18:29:04 INFO Hourly dataset computed
2026-07-13 18:29:04 INFO Creating directories and links for ENS14 to run chimere's parallel part
2026-07-13 18:29:04 INFO Copying EMIS ...
2026-07-13 18:29:04 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens14.nc
2026-07-13 18:29:04 INFO Linking END ...
2026-07-13 18:29:04 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS14/end.2020020500_24_ENS14.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-13 18:29:04 INFO >> Checking links...
2026-07-13 18:29:04 INFO >> All links are good for ENS14  ...
2026-07-13 18:29:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 18:29:05 INFO Hourly dataset computed and listing created
2026-07-13 18:29:06 INFO Hourly dataset computed
2026-07-13 18:29:06 INFO Creating directories and links for ENS15 to run chimere's parallel part
2026-07-13 18:29:06 INFO Copying EMIS ...
2026-07-13 18:29:06 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens15.nc
2026-07-13 18:29:06 INFO Linking END ...
2026-07-13 18:29:06 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS15/end.2020020500_24_ENS15.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-13 18:29:06 INFO >> Checking links...
2026-07-13 18:29:06 INFO >> All links are good for ENS15  ...
2026-07-13 18:29:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 18:29:07 INFO Hourly dataset computed and listing created
2026-07-13 18:29:08 INFO Hourly dataset computed
2026-07-13 18:29:08 INFO ---------->>> Running CHIMERE model from 2020-02-06 00:00:00 to 2020-02-06 01:00:00
2026-07-13 18:29:08 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 18:29:08 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS1
2026-07-13 18:29:08 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS1/end.2020020500_24_ENS1.nc
2026-07-13 18:29:08 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-13 18:29:08 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 18:29:08 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-13 18:29:08 INFO Queuing job for member 1...
2026-07-13 18:29:08 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 18:29:08 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-13 18:29:09 INFO Found: ['5151145']
2026-07-13 18:29:14 INFO [TGCC-IRENE] Submitted job with ID:['5151145']
2026-07-13 18:29:14 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 18:29:14 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS2
2026-07-13 18:29:14 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS2/end.2020020500_24_ENS2.nc
2026-07-13 18:29:14 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-13 18:29:14 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 18:29:14 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-13 18:29:14 INFO Queuing job for member 2...
2026-07-13 18:29:14 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 18:29:14 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-13 18:29:15 INFO Found: ['5151151']
2026-07-13 18:29:20 INFO [TGCC-IRENE] Submitted job with ID:['5151151']
2026-07-13 18:29:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 18:29:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS3
2026-07-13 18:29:20 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS3/end.2020020500_24_ENS3.nc
2026-07-13 18:29:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-13 18:29:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 18:29:20 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-13 18:29:20 INFO Queuing job for member 3...
2026-07-13 18:29:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 18:29:20 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-13 18:29:21 INFO Found: ['5151165']
2026-07-13 18:29:26 INFO [TGCC-IRENE] Submitted job with ID:['5151165']
2026-07-13 18:29:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 18:29:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS4
2026-07-13 18:29:26 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS4/end.2020020500_24_ENS4.nc
2026-07-13 18:29:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-13 18:29:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 18:29:26 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-13 18:29:26 INFO Queuing job for member 4...
2026-07-13 18:29:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 18:29:26 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-13 18:29:27 INFO Found: ['5151167']
2026-07-13 18:29:32 INFO [TGCC-IRENE] Submitted job with ID:['5151167']
2026-07-13 18:29:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 18:29:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS5
2026-07-13 18:29:32 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS5/end.2020020500_24_ENS5.nc
2026-07-13 18:29:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-13 18:29:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 18:29:32 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-13 18:29:32 INFO Queuing job for member 5...
2026-07-13 18:29:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 18:29:32 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-13 18:29:32 INFO Found: ['5151170']
2026-07-13 18:29:37 INFO [TGCC-IRENE] Submitted job with ID:['5151170']
2026-07-13 18:29:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 18:29:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS6
2026-07-13 18:29:37 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS6/end.2020020500_24_ENS6.nc
2026-07-13 18:29:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-13 18:29:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 18:29:37 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-13 18:29:38 INFO Queuing job for member 6...
2026-07-13 18:29:38 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 18:29:38 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-13 18:29:38 INFO Found: ['5151171']
2026-07-13 18:29:43 INFO [TGCC-IRENE] Submitted job with ID:['5151171']
2026-07-13 18:29:43 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 18:29:43 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS7
2026-07-13 18:29:43 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS7/end.2020020500_24_ENS7.nc
2026-07-13 18:29:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-13 18:29:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 18:29:43 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-13 18:29:43 INFO Queuing job for member 7...
2026-07-13 18:29:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 18:29:43 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-13 18:29:44 INFO Found: ['5151174']
2026-07-13 18:29:49 INFO [TGCC-IRENE] Submitted job with ID:['5151174']
2026-07-13 18:29:49 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 18:29:49 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS8
2026-07-13 18:29:49 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS8/end.2020020500_24_ENS8.nc
2026-07-13 18:29:49 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-13 18:29:49 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 18:29:49 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-13 18:29:49 INFO Queuing job for member 8...
2026-07-13 18:29:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 18:29:49 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-13 18:29:50 INFO Found: ['5151176']
2026-07-13 18:29:55 INFO [TGCC-IRENE] Submitted job with ID:['5151176']
2026-07-13 18:29:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 18:29:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS9
2026-07-13 18:29:55 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS9/end.2020020500_24_ENS9.nc
2026-07-13 18:29:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-13 18:29:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 18:29:55 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-13 18:29:55 INFO Queuing job for member 9...
2026-07-13 18:29:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 18:29:55 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-13 18:29:56 INFO Found: ['5151189']
2026-07-13 18:30:01 INFO [TGCC-IRENE] Submitted job with ID:['5151189']
2026-07-13 18:30:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 18:30:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS10
2026-07-13 18:30:01 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS10/end.2020020500_24_ENS10.nc
2026-07-13 18:30:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-13 18:30:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 18:30:01 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-13 18:30:01 INFO Queuing job for member 10...
2026-07-13 18:30:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 18:30:01 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-13 18:30:01 INFO Found: ['5151196']
2026-07-13 18:30:06 INFO [TGCC-IRENE] Submitted job with ID:['5151196']
2026-07-13 18:30:06 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 18:30:06 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS11
2026-07-13 18:30:06 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS11/end.2020020500_24_ENS11.nc
2026-07-13 18:30:06 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-13 18:30:06 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 18:30:07 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-13 18:30:07 INFO Queuing job for member 11...
2026-07-13 18:30:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 18:30:07 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-13 18:30:07 INFO Found: ['5151199']
2026-07-13 18:30:12 INFO [TGCC-IRENE] Submitted job with ID:['5151199']
2026-07-13 18:30:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 18:30:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS12
2026-07-13 18:30:12 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS12/end.2020020500_24_ENS12.nc
2026-07-13 18:30:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-13 18:30:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 18:30:12 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-13 18:30:12 INFO Queuing job for member 12...
2026-07-13 18:30:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 18:30:12 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-13 18:30:13 INFO Found: ['5151204']
2026-07-13 18:30:18 INFO [TGCC-IRENE] Submitted job with ID:['5151204']
2026-07-13 18:30:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 18:30:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS13
2026-07-13 18:30:18 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS13/end.2020020500_24_ENS13.nc
2026-07-13 18:30:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-13 18:30:18 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 18:30:18 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-13 18:30:18 INFO Queuing job for member 13...
2026-07-13 18:30:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 18:30:18 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-13 18:30:19 INFO Found: ['5151211']
2026-07-13 18:30:24 INFO [TGCC-IRENE] Submitted job with ID:['5151211']
2026-07-13 18:30:24 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 18:30:24 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS14
2026-07-13 18:30:24 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS14/end.2020020500_24_ENS14.nc
2026-07-13 18:30:24 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-13 18:30:24 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 18:30:24 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-13 18:30:24 INFO Queuing job for member 14...
2026-07-13 18:30:24 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 18:30:24 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-13 18:30:25 INFO Found: ['5151216']
2026-07-13 18:30:30 INFO [TGCC-IRENE] Submitted job with ID:['5151216']
2026-07-13 18:30:30 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 18:30:30 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS15
2026-07-13 18:30:30 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS15/end.2020020500_24_ENS15.nc
2026-07-13 18:30:30 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-13 18:30:30 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 18:30:30 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-13 18:30:30 INFO Queuing job for member 15...
2026-07-13 18:30:30 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 18:30:30 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-13 18:30:30 INFO Found: ['5151217']
2026-07-13 18:30:35 INFO [TGCC-IRENE] Submitted job with ID:['5151217']
2026-07-13 18:30:35 INFO Checking job status ...
2026-07-13 18:30:35 INFO None 5151145: status RUNNING/PENDING
2026-07-13 18:30:35 INFO None 5151151: status RUNNING/PENDING
2026-07-13 18:30:35 INFO None 5151165: status RUNNING/PENDING
2026-07-13 18:30:35 INFO None 5151167: status RUNNING/PENDING
2026-07-13 18:30:35 INFO None 5151170: status RUNNING/PENDING
2026-07-13 18:30:35 INFO None 5151171: status RUNNING/PENDING
2026-07-13 18:30:36 INFO None 5151174: status RUNNING/PENDING
2026-07-13 18:30:36 INFO None 5151176: status RUNNING/PENDING
2026-07-13 18:30:36 INFO None 5151189: status RUNNING/PENDING
2026-07-13 18:30:36 INFO None 5151196: status RUNNING/PENDING
2026-07-13 18:30:36 INFO None 5151199: status RUNNING/PENDING
2026-07-13 18:30:36 INFO None 5151204: status RUNNING/PENDING
2026-07-13 18:30:36 INFO None 5151211: status RUNNING/PENDING
2026-07-13 18:30:36 INFO None 5151216: status RUNNING/PENDING
2026-07-13 18:30:36 INFO None 5151217: status RUNNING/PENDING
2026-07-13 18:30:36 INFO Jobs still running: ['5151145', '5151151', '5151165', '5151167', '5151170', '5151171', '5151174', '5151176', '5151189', '5151196', '5151199', '5151204', '5151211', '5151216', '5151217']. Waiting...
2026-07-13 18:30:51 INFO None 5151145: status RUNNING/PENDING
2026-07-13 18:30:51 INFO None 5151151: status RUNNING/PENDING
2026-07-13 18:30:51 INFO None 5151165: status RUNNING/PENDING
2026-07-13 18:30:51 INFO None 5151167: status RUNNING/PENDING
2026-07-13 18:30:51 INFO None 5151170: status RUNNING/PENDING
2026-07-13 18:30:51 INFO None 5151171: status RUNNING/PENDING
2026-07-13 18:30:51 INFO None 5151174: status RUNNING/PENDING
2026-07-13 18:30:51 INFO None 5151176: status RUNNING/PENDING
2026-07-13 18:30:51 INFO None 5151189: status RUNNING/PENDING
2026-07-13 18:30:51 INFO None 5151196: status RUNNING/PENDING
2026-07-13 18:30:51 INFO None 5151199: status RUNNING/PENDING
2026-07-13 18:30:51 INFO None 5151204: status RUNNING/PENDING
2026-07-13 18:30:51 INFO None 5151211: status RUNNING/PENDING
2026-07-13 18:30:51 INFO None 5151216: status RUNNING/PENDING
2026-07-13 18:30:51 INFO None 5151217: status RUNNING/PENDING
2026-07-13 18:30:51 INFO Jobs still running: ['5151145', '5151151', '5151165', '5151167', '5151170', '5151171', '5151174', '5151176', '5151189', '5151196', '5151199', '5151204', '5151211', '5151216', '5151217']. Waiting...
2026-07-13 18:31:06 INFO None 5151145: status RUNNING/PENDING
2026-07-13 18:31:06 INFO None 5151151: status RUNNING/PENDING
2026-07-13 18:31:06 INFO None 5151165: status RUNNING/PENDING
2026-07-13 18:31:06 INFO None 5151167: status RUNNING/PENDING
2026-07-13 18:31:06 INFO None 5151170: status RUNNING/PENDING
2026-07-13 18:31:06 INFO None 5151171: status RUNNING/PENDING
2026-07-13 18:31:06 INFO None 5151174: status RUNNING/PENDING
2026-07-13 18:31:06 INFO None 5151176: status RUNNING/PENDING
2026-07-13 18:31:06 INFO None 5151189: status RUNNING/PENDING
2026-07-13 18:31:06 INFO None 5151196: status RUNNING/PENDING
2026-07-13 18:31:06 INFO None 5151199: status RUNNING/PENDING
2026-07-13 18:31:06 INFO None 5151204: status RUNNING/PENDING
2026-07-13 18:31:06 INFO None 5151211: status RUNNING/PENDING
2026-07-13 18:31:06 INFO None 5151216: status RUNNING/PENDING
2026-07-13 18:31:06 INFO None 5151217: status RUNNING/PENDING
2026-07-13 18:31:06 INFO Jobs still running: ['5151145', '5151151', '5151165', '5151167', '5151170', '5151171', '5151174', '5151176', '5151189', '5151196', '5151199', '5151204', '5151211', '5151216', '5151217']. Waiting...
2026-07-13 18:31:21 INFO None 5151145: status RUNNING/PENDING
2026-07-13 18:33:19 INFO None 5151151: status FINISHED
2026-07-13 18:33:19 INFO None 5151165: status FINISHED
2026-07-13 18:33:19 INFO None 5151167: status FINISHED
2026-07-13 18:33:19 INFO None 5151170: status FINISHED
2026-07-13 18:33:19 INFO None 5151171: status FINISHED
2026-07-13 18:33:19 INFO None 5151174: status FINISHED
2026-07-13 18:33:19 INFO None 5151176: status RUNNING/PENDING
2026-07-13 18:33:19 INFO None 5151189: status RUNNING/PENDING
2026-07-13 18:33:19 INFO None 5151196: status FINISHED
2026-07-13 18:33:19 INFO None 5151199: status FINISHED
2026-07-13 18:33:19 INFO None 5151204: status FINISHED
2026-07-13 18:33:19 INFO None 5151211: status FINISHED
2026-07-13 18:33:19 INFO None 5151216: status RUNNING/PENDING
2026-07-13 18:33:19 INFO None 5151217: status RUNNING/PENDING
2026-07-13 18:33:19 INFO Jobs still running: ['5151145', '5151176', '5151189', '5151216', '5151217']. Waiting...
2026-07-13 18:33:34 INFO None 5151145: status RUNNING/PENDING
2026-07-13 18:33:35 INFO None 5151151: status FINISHED
2026-07-13 18:33:35 INFO None 5151165: status FINISHED
2026-07-13 18:33:35 INFO None 5151167: status FINISHED
2026-07-13 18:33:35 INFO None 5151170: status FINISHED
2026-07-13 18:33:35 INFO None 5151171: status FINISHED
2026-07-13 18:33:35 INFO None 5151174: status FINISHED
2026-07-13 18:33:35 INFO None 5151176: status RUNNING/PENDING
2026-07-13 18:33:35 INFO None 5151189: status RUNNING/PENDING
2026-07-13 18:33:35 INFO None 5151196: status FINISHED
2026-07-13 18:33:35 INFO None 5151199: status FINISHED
2026-07-13 18:33:35 INFO None 5151204: status FINISHED
2026-07-13 18:33:35 INFO None 5151211: status FINISHED
2026-07-13 18:33:35 INFO None 5151216: status RUNNING/PENDING
2026-07-13 18:33:35 INFO None 5151217: status RUNNING/PENDING
2026-07-13 18:33:35 INFO Jobs still running: ['5151145', '5151176', '5151189', '5151216', '5151217']. Waiting...
2026-07-13 18:33:50 INFO None 5151145: status RUNNING/PENDING
2026-07-13 18:33:50 INFO None 5151151: status FINISHED
2026-07-13 18:33:50 INFO None 5151165: status FINISHED
2026-07-13 18:33:50 INFO None 5151167: status FINISHED
2026-07-13 18:33:50 INFO None 5151170: status FINISHED
2026-07-13 18:33:50 INFO None 5151171: status FINISHED
2026-07-13 18:33:50 INFO None 5151174: status FINISHED
2026-07-13 18:33:50 INFO None 5151176: status RUNNING/PENDING
2026-07-13 18:33:50 INFO None 5151189: status RUNNING/PENDING
2026-07-13 18:33:50 INFO None 5151196: status FINISHED
2026-07-13 18:33:50 INFO None 5151199: status FINISHED
2026-07-13 18:33:50 INFO None 5151204: status FINISHED
2026-07-13 18:33:50 INFO None 5151211: status FINISHED
2026-07-13 18:33:50 INFO None 5151216: status RUNNING/PENDING
2026-07-13 18:33:50 INFO None 5151217: status RUNNING/PENDING
2026-07-13 18:33:50 INFO Jobs still running: ['5151145', '5151176', '5151189', '5151216', '5151217']. Waiting...
2026-07-13 18:34:05 INFO None 5151145: status RUNNING/PENDING
2026-07-13 18:34:05 INFO None 5151151: status FINISHED
2026-07-13 18:34:05 INFO None 5151165: status FINISHED
2026-07-13 18:34:05 INFO None 5151167: status FINISHED
2026-07-13 18:34:05 INFO None 5151170: status FINISHED
2026-07-13 18:34:05 INFO None 5151171: status FINISHED
2026-07-13 18:34:05 INFO None 5151174: status FINISHED
2026-07-13 18:34:05 INFO None 5151176: status FINISHED
2026-07-13 18:34:05 INFO None 5151189: status RUNNING/PENDING
2026-07-13 18:34:05 INFO None 5151196: status FINISHED
2026-07-13 18:34:05 INFO None 5151199: status FINISHED
2026-07-13 18:34:05 INFO None 5151204: status FINISHED
2026-07-13 18:34:05 INFO None 5151211: status FINISHED
2026-07-13 18:34:05 INFO None 5151216: status RUNNING/PENDING
2026-07-13 18:34:05 INFO None 5151217: status RUNNING/PENDING
2026-07-13 18:34:05 INFO Jobs still running: ['5151145', '5151189', '5151216', '5151217']. Waiting...
2026-07-13 18:34:20 INFO None 5151145: status RUNNING/PENDING
2026-07-13 18:34:20 INFO None 5151151: status FINISHED
2026-07-13 18:34:20 INFO None 5151165: status FINISHED
2026-07-13 18:34:20 INFO None 5151167: status FINISHED
2026-07-13 18:34:20 INFO None 5151170: status FINISHED
2026-07-13 18:34:20 INFO None 5151171: status FINISHED
2026-07-13 18:34:20 INFO None 5151174: status FINISHED
2026-07-13 18:34:21 INFO None 5151176: status FINISHED
2026-07-13 18:34:21 INFO None 5151189: status RUNNING/PENDING
2026-07-13 18:34:21 INFO None 5151196: status FINISHED
2026-07-13 18:34:21 INFO None 5151199: status FINISHED
2026-07-13 18:34:21 INFO None 5151204: status FINISHED
2026-07-13 18:34:21 INFO None 5151211: status FINISHED
2026-07-13 18:34:21 INFO None 5151216: status RUNNING/PENDING
2026-07-13 18:34:21 INFO None 5151217: status RUNNING/PENDING
2026-07-13 18:34:21 INFO Jobs still running: ['5151145', '5151189', '5151216', '5151217']. Waiting...
2026-07-13 18:34:36 INFO None 5151145: status RUNNING/PENDING
2026-07-13 18:34:36 INFO None 5151151: status FINISHED
2026-07-13 18:34:36 INFO None 5151165: status FINISHED
2026-07-13 18:34:36 INFO None 5151167: status FINISHED
2026-07-13 18:34:36 INFO None 5151170: status FINISHED
2026-07-13 18:34:36 INFO None 5151171: status FINISHED
2026-07-13 18:34:36 INFO None 5151174: status FINISHED
2026-07-13 18:34:36 INFO None 5151176: status FINISHED
2026-07-13 18:34:36 INFO None 5151189: status RUNNING/PENDING
2026-07-13 18:34:36 INFO None 5151196: status FINISHED
2026-07-13 18:34:36 INFO None 5151199: status FINISHED
2026-07-13 18:34:36 INFO None 5151204: status FINISHED
2026-07-13 18:34:36 INFO None 5151211: status FINISHED
2026-07-13 18:34:36 INFO None 5151216: status RUNNING/PENDING
2026-07-13 18:34:36 INFO None 5151217: status RUNNING/PENDING
2026-07-13 18:34:36 INFO Jobs still running: ['5151145', '5151189', '5151216', '5151217']. Waiting...
2026-07-13 18:34:51 INFO None 5151145: status FINISHED
2026-07-13 18:34:51 INFO None 5151151: status FINISHED
2026-07-13 18:34:51 INFO None 5151165: status FINISHED
2026-07-13 18:34:51 INFO None 5151167: status FINISHED
2026-07-13 18:34:51 INFO None 5151170: status FINISHED
2026-07-13 18:34:51 INFO None 5151171: status FINISHED
2026-07-13 18:34:51 INFO None 5151174: status FINISHED
2026-07-13 18:34:51 INFO None 5151176: status FINISHED
2026-07-13 18:34:51 INFO None 5151189: status RUNNING/PENDING
2026-07-13 18:34:51 INFO None 5151196: status FINISHED
2026-07-13 18:34:51 INFO None 5151199: status FINISHED
2026-07-13 18:34:51 INFO None 5151204: status FINISHED
2026-07-13 18:34:51 INFO None 5151211: status FINISHED
2026-07-13 18:34:51 INFO None 5151216: status FINISHED
2026-07-13 18:34:51 INFO None 5151217: status RUNNING/PENDING
2026-07-13 18:34:51 INFO Jobs still running: ['5151189', '5151217']. Waiting...
2026-07-13 18:35:06 INFO None 5151145: status FINISHED
2026-07-13 18:35:06 INFO None 5151151: status FINISHED
2026-07-13 18:35:06 INFO None 5151165: status FINISHED
2026-07-13 18:35:06 INFO None 5151167: status FINISHED
2026-07-13 18:35:06 INFO None 5151170: status FINISHED
2026-07-13 18:35:06 INFO None 5151171: status FINISHED
2026-07-13 18:35:06 INFO None 5151174: status FINISHED
2026-07-13 18:35:06 INFO None 5151176: status FINISHED
2026-07-13 18:35:06 INFO None 5151189: status RUNNING/PENDING
2026-07-13 18:35:06 INFO None 5151196: status FINISHED
2026-07-13 18:35:06 INFO None 5151199: status FINISHED
2026-07-13 18:35:06 INFO None 5151204: status FINISHED
2026-07-13 18:35:06 INFO None 5151211: status FINISHED
2026-07-13 18:35:06 INFO None 5151216: status FINISHED
2026-07-13 18:35:06 INFO None 5151217: status RUNNING/PENDING
2026-07-13 18:35:06 INFO Jobs still running: ['5151189', '5151217']. Waiting...
2026-07-13 18:35:21 INFO None 5151145: status FINISHED
2026-07-13 18:35:21 INFO None 5151151: status FINISHED
2026-07-13 18:35:22 INFO None 5151165: status FINISHED
2026-07-13 18:35:22 INFO None 5151167: status FINISHED
2026-07-13 18:35:22 INFO None 5151170: status FINISHED
2026-07-13 18:35:22 INFO None 5151171: status FINISHED
2026-07-13 18:35:22 INFO None 5151174: status FINISHED
2026-07-13 18:35:22 INFO None 5151176: status FINISHED
2026-07-13 18:35:22 INFO None 5151189: status FINISHED
2026-07-13 18:35:22 INFO None 5151196: status FINISHED
2026-07-13 18:35:22 INFO None 5151199: status FINISHED
2026-07-13 18:35:22 INFO None 5151204: status FINISHED
2026-07-13 18:35:22 INFO None 5151211: status FINISHED
2026-07-13 18:35:22 INFO None 5151216: status FINISHED
2026-07-13 18:35:22 INFO None 5151217: status FINISHED
2026-07-13 18:35:22 INFO Jobs ['5151145', '5151151', '5151165', '5151167', '5151170', '5151171', '5151174', '5151176', '5151189', '5151196', '5151199', '5151204', '5151211', '5151216', '5151217'] have finished
2026-07-13 18:35:22 INFO Checking restart files were created ...
2026-07-13 18:35:22 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS1/end.2020020600_1_ENS1.nc(668832435 bytes)
2026-07-13 18:35:22 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS2/end.2020020600_1_ENS2.nc(668832435 bytes)
2026-07-13 18:35:22 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS3/end.2020020600_1_ENS3.nc(668832435 bytes)
2026-07-13 18:35:22 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS4/end.2020020600_1_ENS4.nc(668832435 bytes)
2026-07-13 18:35:22 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS5/end.2020020600_1_ENS5.nc(668832435 bytes)
2026-07-13 18:35:22 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS6/end.2020020600_1_ENS6.nc(668832435 bytes)
2026-07-13 18:35:22 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS7/end.2020020600_1_ENS7.nc(668832435 bytes)
2026-07-13 18:35:22 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS8/end.2020020600_1_ENS8.nc(668832435 bytes)
2026-07-13 18:35:22 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS9/end.2020020600_1_ENS9.nc(668832435 bytes)
2026-07-13 18:35:22 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS10/end.2020020600_1_ENS10.nc(668832435 bytes)
2026-07-13 18:35:22 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS11/end.2020020600_1_ENS11.nc(668832435 bytes)
2026-07-13 18:35:22 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS12/end.2020020600_1_ENS12.nc(668832435 bytes)
2026-07-13 18:35:22 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS13/end.2020020600_1_ENS13.nc(668832435 bytes)
2026-07-13 18:35:22 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS14/end.2020020600_1_ENS14.nc(668832435 bytes)
2026-07-13 18:35:22 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS15/end.2020020600_1_ENS15.nc(668832435 bytes)
2026-07-13 18:35:22 INFO  Run_model() completed successfully.
2026-07-13 18:35:22 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 00:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-13 18:35:22 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 01:00:00 days=153072 seconds=3600
2026-07-13 18:35:22 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-13 18:35:22 INFO [TIME] increment current_time 2020-02-06 00:00:00 -> 2020-02-06 01:00:00
2026-07-13 18:35:22 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-13 18:35:22 INFO ---------->>> Running process_satellite_data()
2026-07-13 18:35:22 INFO [DART] No satellite data found, skipping assimilation
2026-07-13 18:35:22 INFO after_assimilation() skipped
2026-07-13 18:35:22 INFO Next run starts from 2020-02-06 01:00:00
2026-07-13 18:35:22 INFO Cycle is DONE; starting a new loop!
2026-07-13 18:35:22 INFO [TIME] step_end current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-13 18:35:22 INFO [TIME] step_start current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-13 18:35:22 INFO [TIME] window start=2020-02-06 01:00:00 end=2020-02-06 09:00:00 run_hours=8 has_assimilation=True
2026-07-13 18:35:22 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 18:35:23 INFO Hourly dataset computed and listing created
2026-07-13 18:35:38 INFO Hourly dataset computed
2026-07-13 18:35:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 18:35:39 INFO Hourly dataset computed and listing created
2026-07-13 18:35:41 INFO Hourly dataset computed
2026-07-13 18:35:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 18:35:42 INFO Hourly dataset computed and listing created
2026-07-13 18:35:45 INFO Hourly dataset computed
2026-07-13 18:35:45 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 18:35:46 INFO Hourly dataset computed and listing created
2026-07-13 18:35:48 INFO Hourly dataset computed
2026-07-13 18:35:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 18:35:49 INFO Hourly dataset computed and listing created
2026-07-13 18:35:51 INFO Hourly dataset computed
2026-07-13 18:35:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 18:35:52 INFO Hourly dataset computed and listing created
2026-07-13 18:35:54 INFO Hourly dataset computed
2026-07-13 18:35:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 18:35:55 INFO Hourly dataset computed and listing created
2026-07-13 18:35:57 INFO Hourly dataset computed
2026-07-13 18:35:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 18:35:58 INFO Hourly dataset computed and listing created
2026-07-13 18:36:00 INFO Hourly dataset computed
2026-07-13 18:36:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 18:38:42 INFO Hourly dataset computed and listing created
2026-07-13 18:38:44 INFO Hourly dataset computed
2026-07-13 18:38:44 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 18:38:45 INFO Hourly dataset computed and listing created
2026-07-13 18:38:47 INFO Hourly dataset computed
2026-07-13 18:38:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 18:38:48 INFO Hourly dataset computed and listing created
2026-07-13 18:38:50 INFO Hourly dataset computed
2026-07-13 18:38:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 18:38:52 INFO Hourly dataset computed and listing created
2026-07-13 18:38:54 INFO Hourly dataset computed
2026-07-13 18:38:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 18:38:55 INFO Hourly dataset computed and listing created
2026-07-13 18:38:56 INFO Hourly dataset computed
2026-07-13 18:38:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 18:38:58 INFO Hourly dataset computed and listing created
2026-07-13 18:39:00 INFO Hourly dataset computed
2026-07-13 18:39:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 18:39:01 INFO Hourly dataset computed and listing created
2026-07-13 18:39:03 INFO Hourly dataset computed
2026-07-13 18:39:03 INFO ---------->>> Running CHIMERE model from 2020-02-06 01:00:00 to 2020-02-06 09:00:00
2026-07-13 18:39:03 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 18:39:03 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS1
2026-07-13 18:39:03 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS1/end.2020020600_1_ENS1.nc
2026-07-13 18:39:03 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-13 18:39:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 18:39:03 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-13 18:39:03 INFO Queuing job for member 1...
2026-07-13 18:39:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 18:39:03 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-13 18:39:04 INFO Found: ['5151443']
2026-07-13 18:39:09 INFO [TGCC-IRENE] Submitted job with ID:['5151443']
2026-07-13 18:39:09 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 18:39:09 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS2
2026-07-13 18:39:09 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS2/end.2020020600_1_ENS2.nc
2026-07-13 18:39:09 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-13 18:39:09 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 18:39:09 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-13 18:39:09 INFO Queuing job for member 2...
2026-07-13 18:39:09 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 18:39:09 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-13 18:39:10 INFO Found: ['5151452']
2026-07-13 18:39:15 INFO [TGCC-IRENE] Submitted job with ID:['5151452']
2026-07-13 18:39:15 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 18:39:15 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS3
2026-07-13 18:39:15 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS3/end.2020020600_1_ENS3.nc
2026-07-13 18:39:15 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-13 18:39:15 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 18:39:15 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-13 18:39:15 INFO Queuing job for member 3...
2026-07-13 18:39:15 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 18:39:15 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-13 18:39:15 INFO Found: ['5151453']
2026-07-13 18:39:20 INFO [TGCC-IRENE] Submitted job with ID:['5151453']
2026-07-13 18:39:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 18:39:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS4
2026-07-13 18:39:20 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS4/end.2020020600_1_ENS4.nc
2026-07-13 18:39:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-13 18:39:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 18:39:20 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-13 18:39:20 INFO Queuing job for member 4...
2026-07-13 18:39:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 18:39:20 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-13 18:39:21 INFO Found: ['5151458']
2026-07-13 18:39:26 INFO [TGCC-IRENE] Submitted job with ID:['5151458']
2026-07-13 18:39:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 18:39:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS5
2026-07-13 18:39:26 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS5/end.2020020600_1_ENS5.nc
2026-07-13 18:39:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-13 18:39:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 18:39:26 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-13 18:39:26 INFO Queuing job for member 5...
2026-07-13 18:39:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 18:39:26 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-13 18:39:27 INFO Found: ['5151472']
2026-07-13 18:39:32 INFO [TGCC-IRENE] Submitted job with ID:['5151472']
2026-07-13 18:39:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 18:39:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS6
2026-07-13 18:39:32 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS6/end.2020020600_1_ENS6.nc
2026-07-13 18:39:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-13 18:39:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 18:39:32 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-13 18:39:32 INFO Queuing job for member 6...
2026-07-13 18:39:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 18:39:32 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-13 18:39:33 INFO Found: ['5151483']
2026-07-13 18:39:38 INFO [TGCC-IRENE] Submitted job with ID:['5151483']
2026-07-13 18:39:38 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 18:39:38 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS7
2026-07-13 18:39:38 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS7/end.2020020600_1_ENS7.nc
2026-07-13 18:39:38 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-13 18:39:38 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 18:39:38 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-13 18:39:38 INFO Queuing job for member 7...
2026-07-13 18:39:38 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 18:39:38 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-13 18:39:38 INFO Found: ['5151485']
2026-07-13 18:39:43 INFO [TGCC-IRENE] Submitted job with ID:['5151485']
2026-07-13 18:39:43 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 18:39:43 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS8
2026-07-13 18:39:43 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS8/end.2020020600_1_ENS8.nc
2026-07-13 18:39:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-13 18:39:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 18:39:43 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-13 18:39:43 INFO Queuing job for member 8...
2026-07-13 18:39:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 18:39:43 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-13 18:39:44 INFO Found: ['5151487']
2026-07-13 18:39:49 INFO [TGCC-IRENE] Submitted job with ID:['5151487']
2026-07-13 18:39:49 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 18:39:49 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS9
2026-07-13 18:39:49 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS9/end.2020020600_1_ENS9.nc
2026-07-13 18:39:49 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-13 18:39:49 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 18:39:49 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-13 18:39:49 INFO Queuing job for member 9...
2026-07-13 18:39:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 18:39:49 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-13 18:39:50 INFO Found: ['5151488']
2026-07-13 18:39:55 INFO [TGCC-IRENE] Submitted job with ID:['5151488']
2026-07-13 18:39:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 18:39:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS10
2026-07-13 18:39:55 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS10/end.2020020600_1_ENS10.nc
2026-07-13 18:39:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-13 18:39:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 18:39:55 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-13 18:39:55 INFO Queuing job for member 10...
2026-07-13 18:39:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 18:39:55 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-13 18:39:56 INFO Found: ['5151490']
2026-07-13 18:40:01 INFO [TGCC-IRENE] Submitted job with ID:['5151490']
2026-07-13 18:40:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 18:40:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS11
2026-07-13 18:40:01 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS11/end.2020020600_1_ENS11.nc
2026-07-13 18:40:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-13 18:40:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 18:40:01 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-13 18:40:01 INFO Queuing job for member 11...
2026-07-13 18:40:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 18:40:01 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-13 18:40:01 INFO Found: ['5151491']
2026-07-13 18:40:06 INFO [TGCC-IRENE] Submitted job with ID:['5151491']
2026-07-13 18:40:06 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 18:40:06 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS12
2026-07-13 18:40:06 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS12/end.2020020600_1_ENS12.nc
2026-07-13 18:40:06 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-13 18:40:06 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 18:40:06 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-13 18:40:06 INFO Queuing job for member 12...
2026-07-13 18:40:06 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 18:40:06 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-13 18:40:07 INFO Found: ['5151493']
2026-07-13 18:40:12 INFO [TGCC-IRENE] Submitted job with ID:['5151493']
2026-07-13 18:40:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 18:40:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS13
2026-07-13 18:40:12 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS13/end.2020020600_1_ENS13.nc
2026-07-13 18:40:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-13 18:40:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 18:40:12 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-13 18:40:12 INFO Queuing job for member 13...
2026-07-13 18:40:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 18:40:12 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-13 18:40:13 INFO Found: ['5151494']
2026-07-13 18:40:18 INFO [TGCC-IRENE] Submitted job with ID:['5151494']
2026-07-13 18:40:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 18:40:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS14
2026-07-13 18:40:18 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS14/end.2020020600_1_ENS14.nc
2026-07-13 18:40:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-13 18:40:18 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 18:40:18 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-13 18:40:18 INFO Queuing job for member 14...
2026-07-13 18:40:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 18:40:18 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-13 18:40:18 INFO Found: ['5151495']
2026-07-13 18:40:23 INFO [TGCC-IRENE] Submitted job with ID:['5151495']
2026-07-13 18:40:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 18:40:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS15
2026-07-13 18:40:23 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmpemis_0607_15m_low_v2/ENS15/end.2020020600_1_ENS15.nc
2026-07-13 18:40:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-13 18:40:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 18:40:23 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-13 18:40:23 INFO Queuing job for member 15...
2026-07-13 18:40:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 18:40:23 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-13 18:40:24 INFO Found: ['5151496']
2026-07-13 18:40:29 INFO [TGCC-IRENE] Submitted job with ID:['5151496']
2026-07-13 18:40:29 INFO Checking job status ...
2026-07-13 18:40:29 INFO None 5151443: status RUNNING/PENDING
2026-07-13 18:40:29 INFO None 5151452: status RUNNING/PENDING
2026-07-13 18:40:29 INFO None 5151453: status RUNNING/PENDING
2026-07-13 18:40:29 INFO None 5151458: status RUNNING/PENDING
2026-07-13 18:40:29 INFO None 5151472: status RUNNING/PENDING
2026-07-13 18:40:29 INFO None 5151483: status RUNNING/PENDING
2026-07-13 18:40:29 INFO None 5151485: status RUNNING/PENDING
2026-07-13 18:40:29 INFO None 5151487: status RUNNING/PENDING
2026-07-13 18:40:29 INFO None 5151488: status RUNNING/PENDING
2026-07-13 18:40:29 INFO None 5151490: status RUNNING/PENDING
2026-07-13 18:40:29 INFO None 5151491: status RUNNING/PENDING
2026-07-13 18:40:29 INFO None 5151493: status RUNNING/PENDING
2026-07-13 18:40:29 INFO None 5151494: status RUNNING/PENDING
2026-07-13 18:40:29 INFO None 5151495: status RUNNING/PENDING
2026-07-13 18:40:29 INFO None 5151496: status RUNNING/PENDING
2026-07-13 18:40:29 INFO Jobs still running: ['5151443', '5151452', '5151453', '5151458', '5151472', '5151483', '5151485', '5151487', '5151488', '5151490', '5151491', '5151493', '5151494', '5151495', '5151496']. Waiting...
2026-07-13 18:40:44 INFO None 5151443: status RUNNING/PENDING
2026-07-13 18:40:44 INFO None 5151452: status RUNNING/PENDING
2026-07-13 18:40:44 INFO None 5151453: status RUNNING/PENDING
2026-07-13 18:40:44 INFO None 5151458: status RUNNING/PENDING
2026-07-13 18:40:44 INFO None 5151472: status RUNNING/PENDING
2026-07-13 18:40:45 INFO None 5151483: status RUNNING/PENDING
2026-07-13 18:40:45 INFO None 5151485: status RUNNING/PENDING
2026-07-13 18:40:45 INFO None 5151487: status RUNNING/PENDING
2026-07-13 18:40:45 INFO None 5151488: status RUNNING/PENDING
2026-07-13 18:40:45 INFO None 5151490: status RUNNING/PENDING
2026-07-13 18:40:45 INFO None 5151491: status RUNNING/PENDING
2026-07-13 18:40:45 INFO None 5151493: status RUNNING/PENDING
2026-07-13 18:40:45 INFO None 5151494: status RUNNING/PENDING
2026-07-13 18:40:45 INFO None 5151495: status RUNNING/PENDING
2026-07-13 18:40:45 INFO None 5151496: status RUNNING/PENDING
2026-07-13 18:40:45 INFO Jobs still running: ['5151443', '5151452', '5151453', '5151458', '5151472', '5151483', '5151485', '5151487', '5151488', '5151490', '5151491', '5151493', '5151494', '5151495', '5151496']. Waiting...
2026-07-13 18:41:00 INFO None 5151443: status RUNNING/PENDING
2026-07-13 18:41:00 INFO None 5151452: status RUNNING/PENDING
2026-07-13 18:41:00 INFO None 5151453: status RUNNING/PENDING
2026-07-13 18:41:00 INFO None 5151458: status RUNNING/PENDING
2026-07-13 18:41:00 INFO None 5151472: status RUNNING/PENDING
2026-07-13 18:41:00 INFO None 5151483: status RUNNING/PENDING
2026-07-13 18:41:00 INFO None 5151485: status RUNNING/PENDING
2026-07-13 18:41:00 INFO None 5151487: status RUNNING/PENDING
2026-07-13 18:41:00 INFO None 5151488: status RUNNING/PENDING
2026-07-13 18:41:00 INFO None 5151490: status RUNNING/PENDING
2026-07-13 18:41:00 INFO None 5151491: status RUNNING/PENDING
2026-07-13 18:41:00 INFO None 5151493: status RUNNING/PENDING
2026-07-13 18:41:00 INFO None 5151494: status RUNNING/PENDING
2026-07-13 18:41:00 INFO None 5151495: status RUNNING/PENDING
2026-07-13 18:41:00 INFO None 5151496: status RUNNING/PENDING
2026-07-13 18:41:00 INFO Jobs still running: ['5151443', '5151452', '5151453', '5151458', '5151472', '5151483', '5151485', '5151487', '5151488', '5151490', '5151491', '5151493', '5151494', '5151495', '5151496']. Waiting...
2026-07-13 18:41:15 INFO None 5151443: status RUNNING/PENDING
2026-07-13 18:41:15 INFO None 5151452: status RUNNING/PENDING
2026-07-13 18:41:15 INFO None 5151453: status RUNNING/PENDING
2026-07-13 18:41:15 INFO None 5151458: status RUNNING/PENDING
2026-07-13 18:41:15 INFO None 5151472: status RUNNING/PENDING
2026-07-13 18:41:15 INFO None 5151483: status RUNNING/PENDING
2026-07-13 18:41:15 INFO None 5151485: status RUNNING/PENDING
2026-07-13 18:41:15 INFO None 5151487: status RUNNING/PENDING
2026-07-13 18:41:15 INFO None 5151488: status RUNNING/PENDING
2026-07-13 18:41:15 INFO None 5151490: status RUNNING/PENDING
2026-07-13 18:41:15 INFO None 5151491: status RUNNING/PENDING
2026-07-13 18:41:15 INFO None 5151493: status RUNNING/PENDING
2026-07-13 18:41:15 INFO None 5151494: status RUNNING/PENDING
2026-07-13 18:41:15 INFO None 5151495: status RUNNING/PENDING
2026-07-13 18:41:15 INFO None 5151496: status RUNNING/PENDING
2026-07-13 18:41:15 INFO Jobs still running: ['5151443', '5151452', '5151453', '5151458', '5151472', '5151483', '5151485', '5151487', '5151488', '5151490', '5151491', '5151493', '5151494', '5151495', '5151496']. Waiting...
2026-07-13 18:41:30 INFO None 5151443: status RUNNING/PENDING
2026-07-13 18:41:30 INFO None 5151452: status RUNNING/PENDING
2026-07-13 18:41:30 INFO None 5151453: status RUNNING/PENDING
2026-07-13 18:41:30 INFO None 5151458: status RUNNING/PENDING
2026-07-13 18:41:30 INFO None 5151472: status RUNNING/PENDING
2026-07-13 18:41:30 INFO None 5151483: status RUNNING/PENDING
2026-07-13 18:41:30 INFO None 5151485: status RUNNING/PENDING
2026-07-13 18:41:30 INFO None 5151487: status RUNNING/PENDING
2026-07-13 18:41:30 INFO None 5151488: status RUNNING/PENDING
2026-07-13 18:41:30 INFO None 5151490: status RUNNING/PENDING
2026-07-13 18:41:30 INFO None 5151491: status RUNNING/PENDING
2026-07-13 18:41:30 INFO None 5151493: status RUNNING/PENDING
2026-07-13 18:41:30 INFO None 5151494: status RUNNING/PENDING
2026-07-13 18:41:31 INFO None 5151495: status RUNNING/PENDING
2026-07-13 18:41:31 INFO None 5151496: status RUNNING/PENDING
2026-07-13 18:41:31 INFO Jobs still running: ['5151443', '5151452', '5151453', '5151458', '5151472', '5151483', '5151485', '5151487', '5151488', '5151490', '5151491', '5151493', '5151494', '5151495', '5151496']. Waiting...
[2026-07-13T18:42:28.765] error: *** JOB 5151058 ON irene4337 CANCELLED AT 2026-07-13T18:42:28 DUE to SIGNAL Terminated ***
[2026-07-13T18:42:35.262] error: _handle_signal_container: failed signal 15 pid 129187 StepId=5151058.batch No such process
[2026-07-13T18:43:36.007] error: _handle_signal_container: failed signal 15 pid 129187 StepId=5151058.batch No such process
[2026-07-13T18:43:36.007] error: _handle_signal_container: failed signal 15 pid 129187 StepId=5151058.batch No such process
[2026-07-13T18:43:36.008] error: _handle_signal_container: failed signal 15 pid 129187 StepId=5151058.batch No such process
[2026-07-13T18:43:36.008] error: _handle_signal_container: failed signal 15 pid 129187 StepId=5151058.batch No such process
[2026-07-13T18:43:36.008] error: _handle_signal_container: failed signal 15 pid 129187 StepId=5151058.batch No such process
2026-07-13 18:41:46 INFO None 5151443: status RUNNING/PENDING
