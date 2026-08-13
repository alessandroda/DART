+ SCRIPT_PID=1127796
+ /bin/bash -x /tmp/tmp.hVJbum2wdS
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
2026-06-13 00:41:27 INFO 
███    ███ ██ ███    ███ ███████ ███████ ██
████  ████ ██ ████  ████ ██      ██      ██
██ ████ ██ ██ ██ ████ ██ █████   ███████ ██
██  ██  ██ ██ ██  ██  ██ ██           ██ ██
██      ██ ██ ██      ██ ███████ ███████ ██



2026-06-13 00:41:27 INFO [PIPELINE] =======================================
2026-06-13 00:41:27 INFO [PIPELINE] Starting chimere–DART orchestrator
2026-06-13 00:41:27 INFO [PIPELINE] Config file: config/config_irene_IM.yaml
2026-06-13 00:41:27 INFO [PIPELINE] Run dir: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart
2026-06-13 00:41:27 INFO [PIPELINE] Log file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/mimesi_orchestrator_logs/chimere_DART_20260613_004127.log
2026-06-13 00:41:27 INFO [PIPELINE] =======================================
2026-06-13 00:41:27 INFO Running assimilation with model_type=ModelType.CHIMERE
2026-06-13 00:41:27 INFO Using scheduler=Scheduler.SLURM, queue=rome
2026-06-13 00:41:27 INFO [STEP] ---- TIME LOOP START ----
2026-06-13 00:41:27 INFO [TIME] step_start current_time=2020-02-06 00:00:00 simulated_time=None dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-13 00:41:27 INFO [TIME] window start=2020-02-06 00:00:00 end=2020-02-06 01:00:00 run_hours=1 has_assimilation=False
2026-06-13 00:41:27 INFO Creating directories and links for ENS1 to run chimere's parallel part
2026-06-13 00:41:27 INFO Linking EMIS ...
2026-06-13 00:41:27 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc -> /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens01.nc
2026-06-13 00:41:27 INFO Linking END ...
2026-06-13 00:41:27 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS1/end.2020020500_24_ENS1.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-13 00:41:27 INFO >> Checking links...
2026-06-13 00:41:27 INFO >> All links are good for ENS1  ...
2026-06-13 00:41:27 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-13 00:41:34 INFO Hourly dataset computed and listing created
2026-06-13 00:41:39 INFO Hourly dataset computed
2026-06-13 00:41:39 INFO Creating directories and links for ENS2 to run chimere's parallel part
2026-06-13 00:41:39 INFO Linking EMIS ...
2026-06-13 00:41:39 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc -> /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens02.nc
2026-06-13 00:41:39 INFO Linking END ...
2026-06-13 00:41:39 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS2/end.2020020500_24_ENS2.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-13 00:41:39 INFO >> Checking links...
2026-06-13 00:41:39 INFO >> All links are good for ENS2  ...
2026-06-13 00:41:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-13 00:41:40 INFO Hourly dataset computed and listing created
2026-06-13 00:41:41 INFO Hourly dataset computed
2026-06-13 00:41:41 INFO Creating directories and links for ENS3 to run chimere's parallel part
2026-06-13 00:41:41 INFO Linking EMIS ...
2026-06-13 00:41:41 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc -> /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens03.nc
2026-06-13 00:41:41 INFO Linking END ...
2026-06-13 00:41:41 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS3/end.2020020500_24_ENS3.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-13 00:41:41 INFO >> Checking links...
2026-06-13 00:41:41 INFO >> All links are good for ENS3  ...
2026-06-13 00:41:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-13 00:41:42 INFO Hourly dataset computed and listing created
2026-06-13 00:41:42 INFO Hourly dataset computed
2026-06-13 00:41:42 INFO Creating directories and links for ENS4 to run chimere's parallel part
2026-06-13 00:41:42 INFO Linking EMIS ...
2026-06-13 00:41:42 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc -> /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens04.nc
2026-06-13 00:41:42 INFO Linking END ...
2026-06-13 00:41:42 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS4/end.2020020500_24_ENS4.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-13 00:41:42 INFO >> Checking links...
2026-06-13 00:41:42 INFO >> All links are good for ENS4  ...
2026-06-13 00:41:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-13 00:41:43 INFO Hourly dataset computed and listing created
2026-06-13 00:41:44 INFO Hourly dataset computed
2026-06-13 00:41:44 INFO Creating directories and links for ENS5 to run chimere's parallel part
2026-06-13 00:41:44 INFO Linking EMIS ...
2026-06-13 00:41:44 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc -> /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens05.nc
2026-06-13 00:41:44 INFO Linking END ...
2026-06-13 00:41:44 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS5/end.2020020500_24_ENS5.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-13 00:41:44 INFO >> Checking links...
2026-06-13 00:41:44 INFO >> All links are good for ENS5  ...
2026-06-13 00:41:44 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-13 00:41:44 INFO Hourly dataset computed and listing created
2026-06-13 00:41:45 INFO Hourly dataset computed
2026-06-13 00:41:45 INFO Creating directories and links for ENS6 to run chimere's parallel part
2026-06-13 00:41:45 INFO Linking EMIS ...
2026-06-13 00:41:45 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc -> /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens06.nc
2026-06-13 00:41:45 INFO Linking END ...
2026-06-13 00:41:45 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS6/end.2020020500_24_ENS6.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-13 00:41:45 INFO >> Checking links...
2026-06-13 00:41:45 INFO >> All links are good for ENS6  ...
2026-06-13 00:41:45 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-13 00:41:46 INFO Hourly dataset computed and listing created
2026-06-13 00:41:47 INFO Hourly dataset computed
2026-06-13 00:41:47 INFO Creating directories and links for ENS7 to run chimere's parallel part
2026-06-13 00:41:47 INFO Linking EMIS ...
2026-06-13 00:41:47 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc -> /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens07.nc
2026-06-13 00:41:47 INFO Linking END ...
2026-06-13 00:41:47 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS7/end.2020020500_24_ENS7.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-13 00:41:47 INFO >> Checking links...
2026-06-13 00:41:47 INFO >> All links are good for ENS7  ...
2026-06-13 00:41:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-13 00:41:47 INFO Hourly dataset computed and listing created
2026-06-13 00:41:48 INFO Hourly dataset computed
2026-06-13 00:41:48 INFO Creating directories and links for ENS8 to run chimere's parallel part
2026-06-13 00:41:48 INFO Linking EMIS ...
2026-06-13 00:41:48 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc -> /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens08.nc
2026-06-13 00:41:48 INFO Linking END ...
2026-06-13 00:41:48 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS8/end.2020020500_24_ENS8.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-13 00:41:48 INFO >> Checking links...
2026-06-13 00:41:48 INFO >> All links are good for ENS8  ...
2026-06-13 00:41:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-13 00:41:49 INFO Hourly dataset computed and listing created
2026-06-13 00:41:49 INFO Hourly dataset computed
2026-06-13 00:41:49 INFO Creating directories and links for ENS9 to run chimere's parallel part
2026-06-13 00:41:49 INFO Linking EMIS ...
2026-06-13 00:41:49 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc -> /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens09.nc
2026-06-13 00:41:49 INFO Linking END ...
2026-06-13 00:41:49 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS9/end.2020020500_24_ENS9.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-13 00:41:49 INFO >> Checking links...
2026-06-13 00:41:49 INFO >> All links are good for ENS9  ...
2026-06-13 00:41:49 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-13 00:41:50 INFO Hourly dataset computed and listing created
2026-06-13 00:41:51 INFO Hourly dataset computed
2026-06-13 00:41:51 INFO Creating directories and links for ENS10 to run chimere's parallel part
2026-06-13 00:41:51 INFO Linking EMIS ...
2026-06-13 00:41:51 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc -> /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens10.nc
2026-06-13 00:41:51 INFO Linking END ...
2026-06-13 00:41:51 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS10/end.2020020500_24_ENS10.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-13 00:41:51 INFO >> Checking links...
2026-06-13 00:41:51 INFO >> All links are good for ENS10  ...
2026-06-13 00:41:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-13 00:41:52 INFO Hourly dataset computed and listing created
2026-06-13 00:41:52 INFO Hourly dataset computed
2026-06-13 00:41:52 INFO Creating directories and links for ENS11 to run chimere's parallel part
2026-06-13 00:41:52 INFO Linking EMIS ...
2026-06-13 00:41:52 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc -> /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens11.nc
2026-06-13 00:41:52 INFO Linking END ...
2026-06-13 00:41:52 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS11/end.2020020500_24_ENS11.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-13 00:41:52 INFO >> Checking links...
2026-06-13 00:41:52 INFO >> All links are good for ENS11  ...
2026-06-13 00:41:52 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-13 00:41:53 INFO Hourly dataset computed and listing created
2026-06-13 00:41:54 INFO Hourly dataset computed
2026-06-13 00:41:54 INFO Creating directories and links for ENS12 to run chimere's parallel part
2026-06-13 00:41:54 INFO Linking EMIS ...
2026-06-13 00:41:54 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc -> /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens12.nc
2026-06-13 00:41:54 INFO Linking END ...
2026-06-13 00:41:54 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS12/end.2020020500_24_ENS12.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-13 00:41:54 INFO >> Checking links...
2026-06-13 00:41:54 INFO >> All links are good for ENS12  ...
2026-06-13 00:41:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-13 00:41:54 INFO Hourly dataset computed and listing created
2026-06-13 00:41:55 INFO Hourly dataset computed
2026-06-13 00:41:55 INFO Creating directories and links for ENS13 to run chimere's parallel part
2026-06-13 00:41:55 INFO Linking EMIS ...
2026-06-13 00:41:55 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc -> /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens13.nc
2026-06-13 00:41:55 INFO Linking END ...
2026-06-13 00:41:55 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS13/end.2020020500_24_ENS13.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-13 00:41:55 INFO >> Checking links...
2026-06-13 00:41:55 INFO >> All links are good for ENS13  ...
2026-06-13 00:41:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-13 00:41:56 INFO Hourly dataset computed and listing created
2026-06-13 00:41:56 INFO Hourly dataset computed
2026-06-13 00:41:56 INFO Creating directories and links for ENS14 to run chimere's parallel part
2026-06-13 00:41:56 INFO Linking EMIS ...
2026-06-13 00:41:57 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc -> /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens14.nc
2026-06-13 00:41:57 INFO Linking END ...
2026-06-13 00:41:57 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS14/end.2020020500_24_ENS14.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-13 00:41:57 INFO >> Checking links...
2026-06-13 00:41:57 INFO >> All links are good for ENS14  ...
2026-06-13 00:41:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-13 00:41:57 INFO Hourly dataset computed and listing created
2026-06-13 00:41:58 INFO Hourly dataset computed
2026-06-13 00:41:58 INFO Creating directories and links for ENS15 to run chimere's parallel part
2026-06-13 00:41:58 INFO Linking EMIS ...
2026-06-13 00:41:58 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc -> /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens15.nc
2026-06-13 00:41:58 INFO Linking END ...
2026-06-13 00:41:58 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS15/end.2020020500_24_ENS15.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-13 00:41:58 INFO >> Checking links...
2026-06-13 00:41:58 INFO >> All links are good for ENS15  ...
2026-06-13 00:41:58 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-13 00:41:59 INFO Hourly dataset computed and listing created
2026-06-13 00:41:59 INFO Hourly dataset computed
2026-06-13 00:41:59 INFO ---------->>> Running CHIMERE model from 2020-02-06 00:00:00 to 2020-02-06 01:00:00
2026-06-13 00:41:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-13 00:41:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS1
2026-06-13 00:41:59 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS1/end.2020020500_24_ENS1.nc
2026-06-13 00:41:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-06-13 00:41:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-13 00:41:59 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-06-13 00:41:59 INFO Queuing job for member 1...
2026-06-13 00:41:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-13 00:41:59 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-06-13 00:42:02 INFO Found: ['4828449']
2026-06-13 00:42:07 INFO [TGCC-IRENE] Submitted job with ID:['4828449']
2026-06-13 00:42:07 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-13 00:42:07 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS2
2026-06-13 00:42:07 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS2/end.2020020500_24_ENS2.nc
2026-06-13 00:42:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-06-13 00:42:07 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-13 00:42:07 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-06-13 00:42:07 INFO Queuing job for member 2...
2026-06-13 00:42:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-13 00:42:07 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-06-13 00:42:09 INFO Found: ['4828451']
2026-06-13 00:42:14 INFO [TGCC-IRENE] Submitted job with ID:['4828451']
2026-06-13 00:42:14 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-13 00:42:14 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS3
2026-06-13 00:42:14 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS3/end.2020020500_24_ENS3.nc
2026-06-13 00:42:14 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-06-13 00:42:14 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-13 00:42:14 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-06-13 00:42:14 INFO Queuing job for member 3...
2026-06-13 00:42:14 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-13 00:42:14 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-06-13 00:42:17 INFO Found: ['4828452']
2026-06-13 00:42:22 INFO [TGCC-IRENE] Submitted job with ID:['4828452']
2026-06-13 00:42:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-13 00:42:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS4
2026-06-13 00:42:22 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS4/end.2020020500_24_ENS4.nc
2026-06-13 00:42:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-06-13 00:42:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-13 00:42:22 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-06-13 00:42:22 INFO Queuing job for member 4...
2026-06-13 00:42:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-13 00:42:22 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-06-13 00:42:24 INFO Found: ['4828453']
2026-06-13 00:42:29 INFO [TGCC-IRENE] Submitted job with ID:['4828453']
2026-06-13 00:42:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-13 00:42:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS5
2026-06-13 00:42:29 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS5/end.2020020500_24_ENS5.nc
2026-06-13 00:42:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-06-13 00:42:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-13 00:42:29 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-06-13 00:42:29 INFO Queuing job for member 5...
2026-06-13 00:42:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-13 00:42:29 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-06-13 00:42:32 INFO Found: ['4828456']
2026-06-13 00:42:37 INFO [TGCC-IRENE] Submitted job with ID:['4828456']
2026-06-13 00:42:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-13 00:42:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS6
2026-06-13 00:42:37 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS6/end.2020020500_24_ENS6.nc
2026-06-13 00:42:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-06-13 00:42:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-13 00:42:37 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-06-13 00:42:37 INFO Queuing job for member 6...
2026-06-13 00:42:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-13 00:42:37 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-06-13 00:42:38 INFO Found: ['4828457']
2026-06-13 00:42:43 INFO [TGCC-IRENE] Submitted job with ID:['4828457']
2026-06-13 00:42:43 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-13 00:42:43 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS7
2026-06-13 00:42:43 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS7/end.2020020500_24_ENS7.nc
2026-06-13 00:42:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-06-13 00:42:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-13 00:42:43 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-06-13 00:42:43 INFO Queuing job for member 7...
2026-06-13 00:42:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-13 00:42:43 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-06-13 00:42:44 INFO Found: ['4828459']
2026-06-13 00:42:49 INFO [TGCC-IRENE] Submitted job with ID:['4828459']
2026-06-13 00:42:49 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-13 00:42:49 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS8
2026-06-13 00:42:49 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS8/end.2020020500_24_ENS8.nc
2026-06-13 00:42:49 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-06-13 00:42:49 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-13 00:42:49 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-06-13 00:42:49 INFO Queuing job for member 8...
2026-06-13 00:42:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-13 00:42:49 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-06-13 00:42:49 INFO Found: ['4828461']
2026-06-13 00:42:54 INFO [TGCC-IRENE] Submitted job with ID:['4828461']
2026-06-13 00:42:54 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-13 00:42:54 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS9
2026-06-13 00:42:54 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS9/end.2020020500_24_ENS9.nc
2026-06-13 00:42:54 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-06-13 00:42:54 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-13 00:42:54 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-06-13 00:42:54 INFO Queuing job for member 9...
2026-06-13 00:42:54 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-13 00:42:54 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-06-13 00:42:55 INFO Found: ['4828463']
2026-06-13 00:43:00 INFO [TGCC-IRENE] Submitted job with ID:['4828463']
2026-06-13 00:43:00 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-13 00:43:00 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS10
2026-06-13 00:43:00 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS10/end.2020020500_24_ENS10.nc
2026-06-13 00:43:00 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-06-13 00:43:00 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-13 00:43:00 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-06-13 00:43:00 INFO Queuing job for member 10...
2026-06-13 00:43:00 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-13 00:43:00 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-06-13 00:43:01 INFO Found: ['4828464']
2026-06-13 00:43:06 INFO [TGCC-IRENE] Submitted job with ID:['4828464']
2026-06-13 00:43:06 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-13 00:43:06 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS11
2026-06-13 00:43:06 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS11/end.2020020500_24_ENS11.nc
2026-06-13 00:43:06 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-06-13 00:43:06 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-13 00:43:06 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-06-13 00:43:06 INFO Queuing job for member 11...
2026-06-13 00:43:06 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-13 00:43:06 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-06-13 00:43:06 INFO Found: ['4828466']
2026-06-13 00:43:11 INFO [TGCC-IRENE] Submitted job with ID:['4828466']
2026-06-13 00:43:11 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-13 00:43:11 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS12
2026-06-13 00:43:11 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS12/end.2020020500_24_ENS12.nc
2026-06-13 00:43:11 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-06-13 00:43:11 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-13 00:43:11 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-06-13 00:43:11 INFO Queuing job for member 12...
2026-06-13 00:43:11 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-13 00:43:11 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-06-13 00:43:12 INFO Found: ['4828467']
2026-06-13 00:43:17 INFO [TGCC-IRENE] Submitted job with ID:['4828467']
2026-06-13 00:43:17 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-13 00:43:17 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS13
2026-06-13 00:43:17 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS13/end.2020020500_24_ENS13.nc
2026-06-13 00:43:17 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-06-13 00:43:17 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-13 00:43:17 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-06-13 00:43:17 INFO Queuing job for member 13...
2026-06-13 00:43:17 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-13 00:43:17 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-06-13 00:43:18 INFO Found: ['4828468']
2026-06-13 00:43:23 INFO [TGCC-IRENE] Submitted job with ID:['4828468']
2026-06-13 00:43:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-13 00:43:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS14
2026-06-13 00:43:23 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS14/end.2020020500_24_ENS14.nc
2026-06-13 00:43:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-06-13 00:43:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-13 00:43:23 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-06-13 00:43:23 INFO Queuing job for member 14...
2026-06-13 00:43:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-13 00:43:23 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-06-13 00:43:25 INFO Found: ['4828469']
2026-06-13 00:43:30 INFO [TGCC-IRENE] Submitted job with ID:['4828469']
2026-06-13 00:43:30 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-13 00:43:30 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS15
2026-06-13 00:43:30 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS15/end.2020020500_24_ENS15.nc
2026-06-13 00:43:30 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-06-13 00:43:30 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-13 00:43:30 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-06-13 00:43:30 INFO Queuing job for member 15...
2026-06-13 00:43:30 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-13 00:43:30 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-06-13 00:43:32 INFO Found: ['4828470']
2026-06-13 00:43:37 INFO [TGCC-IRENE] Submitted job with ID:['4828470']
2026-06-13 00:43:37 INFO Checking job status ...
2026-06-13 00:43:37 INFO None 4828449: status RUNNING/PENDING
2026-06-13 00:43:37 INFO None 4828451: status RUNNING/PENDING
2026-06-13 00:43:37 INFO None 4828452: status RUNNING/PENDING
2026-06-13 00:43:37 INFO None 4828453: status RUNNING/PENDING
2026-06-13 00:43:37 INFO None 4828456: status RUNNING/PENDING
2026-06-13 00:43:37 INFO None 4828457: status RUNNING/PENDING
2026-06-13 00:43:37 INFO None 4828459: status RUNNING/PENDING
2026-06-13 00:43:37 INFO None 4828461: status RUNNING/PENDING
2026-06-13 00:43:37 INFO None 4828463: status RUNNING/PENDING
2026-06-13 00:43:37 INFO None 4828464: status RUNNING/PENDING
2026-06-13 00:43:37 INFO None 4828466: status RUNNING/PENDING
2026-06-13 00:43:37 INFO None 4828467: status RUNNING/PENDING
2026-06-13 00:43:37 INFO None 4828468: status RUNNING/PENDING
2026-06-13 00:43:37 INFO None 4828469: status RUNNING/PENDING
2026-06-13 00:43:37 INFO None 4828470: status RUNNING/PENDING
2026-06-13 00:43:37 INFO Jobs still running: ['4828449', '4828451', '4828452', '4828453', '4828456', '4828457', '4828459', '4828461', '4828463', '4828464', '4828466', '4828467', '4828468', '4828469', '4828470']. Waiting...
2026-06-13 00:43:52 INFO None 4828449: status RUNNING/PENDING
2026-06-13 00:43:52 INFO None 4828451: status RUNNING/PENDING
2026-06-13 00:43:52 INFO None 4828452: status RUNNING/PENDING
2026-06-13 00:43:52 INFO None 4828453: status RUNNING/PENDING
2026-06-13 00:43:52 INFO None 4828456: status RUNNING/PENDING
2026-06-13 00:43:52 INFO None 4828457: status RUNNING/PENDING
2026-06-13 00:43:52 INFO None 4828459: status RUNNING/PENDING
2026-06-13 00:43:52 INFO None 4828461: status RUNNING/PENDING
2026-06-13 00:43:52 INFO None 4828463: status RUNNING/PENDING
2026-06-13 00:43:53 INFO None 4828464: status RUNNING/PENDING
2026-06-13 00:43:53 INFO None 4828466: status RUNNING/PENDING
2026-06-13 00:43:55 INFO None 4828467: status RUNNING/PENDING
2026-06-13 00:43:55 INFO None 4828468: status RUNNING/PENDING
2026-06-13 00:43:55 INFO None 4828469: status RUNNING/PENDING
2026-06-13 00:43:55 INFO None 4828470: status RUNNING/PENDING
2026-06-13 00:43:55 INFO Jobs still running: ['4828449', '4828451', '4828452', '4828453', '4828456', '4828457', '4828459', '4828461', '4828463', '4828464', '4828466', '4828467', '4828468', '4828469', '4828470']. Waiting...
2026-06-13 00:44:10 INFO None 4828449: status RUNNING/PENDING
2026-06-13 00:44:10 INFO None 4828451: status RUNNING/PENDING
2026-06-13 00:44:10 INFO None 4828452: status RUNNING/PENDING
2026-06-13 00:44:10 INFO None 4828453: status RUNNING/PENDING
2026-06-13 00:44:10 INFO None 4828456: status RUNNING/PENDING
2026-06-13 00:44:10 INFO None 4828457: status RUNNING/PENDING
2026-06-13 00:44:10 INFO None 4828459: status RUNNING/PENDING
2026-06-13 00:44:10 INFO None 4828461: status RUNNING/PENDING
2026-06-13 00:44:10 INFO None 4828463: status RUNNING/PENDING
2026-06-13 00:44:10 INFO None 4828464: status RUNNING/PENDING
2026-06-13 00:44:10 INFO None 4828466: status RUNNING/PENDING
2026-06-13 00:44:10 INFO None 4828467: status RUNNING/PENDING
2026-06-13 00:44:10 INFO None 4828468: status RUNNING/PENDING
2026-06-13 00:44:10 INFO None 4828469: status RUNNING/PENDING
2026-06-13 00:44:10 INFO None 4828470: status RUNNING/PENDING
2026-06-13 00:44:10 INFO Jobs still running: ['4828449', '4828451', '4828452', '4828453', '4828456', '4828457', '4828459', '4828461', '4828463', '4828464', '4828466', '4828467', '4828468', '4828469', '4828470']. Waiting...
2026-06-13 00:44:25 INFO None 4828449: status FINISHED
2026-06-13 00:44:25 INFO None 4828451: status FINISHED
2026-06-13 00:44:25 INFO None 4828452: status FINISHED
2026-06-13 00:44:25 INFO None 4828453: status RUNNING/PENDING
2026-06-13 00:44:25 INFO None 4828456: status RUNNING/PENDING
2026-06-13 00:44:25 INFO None 4828457: status RUNNING/PENDING
2026-06-13 00:44:25 INFO None 4828459: status RUNNING/PENDING
2026-06-13 00:44:25 INFO None 4828461: status RUNNING/PENDING
2026-06-13 00:44:25 INFO None 4828463: status RUNNING/PENDING
2026-06-13 00:44:25 INFO None 4828464: status RUNNING/PENDING
2026-06-13 00:44:25 INFO None 4828466: status RUNNING/PENDING
2026-06-13 00:44:25 INFO None 4828467: status RUNNING/PENDING
2026-06-13 00:44:25 INFO None 4828468: status RUNNING/PENDING
2026-06-13 00:44:25 INFO None 4828469: status RUNNING/PENDING
2026-06-13 00:44:27 INFO None 4828470: status RUNNING/PENDING
2026-06-13 00:44:27 INFO Jobs still running: ['4828453', '4828456', '4828457', '4828459', '4828461', '4828463', '4828464', '4828466', '4828467', '4828468', '4828469', '4828470']. Waiting...
2026-06-13 00:44:42 INFO None 4828449: status FINISHED
2026-06-13 00:44:42 INFO None 4828451: status FINISHED
2026-06-13 00:44:42 INFO None 4828452: status FINISHED
2026-06-13 00:44:42 INFO None 4828453: status RUNNING/PENDING
2026-06-13 00:44:42 INFO None 4828456: status RUNNING/PENDING
2026-06-13 00:44:42 INFO None 4828457: status RUNNING/PENDING
2026-06-13 00:44:42 INFO None 4828459: status RUNNING/PENDING
2026-06-13 00:44:42 INFO None 4828461: status RUNNING/PENDING
2026-06-13 00:44:42 INFO None 4828463: status RUNNING/PENDING
2026-06-13 00:44:43 INFO None 4828464: status RUNNING/PENDING
2026-06-13 00:44:43 INFO None 4828466: status RUNNING/PENDING
2026-06-13 00:44:43 INFO None 4828467: status RUNNING/PENDING
2026-06-13 00:44:43 INFO None 4828468: status RUNNING/PENDING
2026-06-13 00:44:43 INFO None 4828469: status RUNNING/PENDING
2026-06-13 00:44:43 INFO None 4828470: status RUNNING/PENDING
2026-06-13 00:44:43 INFO Jobs still running: ['4828453', '4828456', '4828457', '4828459', '4828461', '4828463', '4828464', '4828466', '4828467', '4828468', '4828469', '4828470']. Waiting...
2026-06-13 00:44:58 INFO None 4828449: status FINISHED
2026-06-13 00:44:58 INFO None 4828451: status FINISHED
2026-06-13 00:44:58 INFO None 4828452: status FINISHED
2026-06-13 00:44:58 INFO None 4828453: status FINISHED
2026-06-13 00:44:58 INFO None 4828456: status FINISHED
2026-06-13 00:44:58 INFO None 4828457: status FINISHED
2026-06-13 00:44:58 INFO None 4828459: status FINISHED
2026-06-13 00:44:58 INFO None 4828461: status FINISHED
2026-06-13 00:44:58 INFO None 4828463: status FINISHED
2026-06-13 00:44:58 INFO None 4828464: status RUNNING/PENDING
2026-06-13 00:44:58 INFO None 4828466: status RUNNING/PENDING
2026-06-13 00:44:58 INFO None 4828467: status RUNNING/PENDING
2026-06-13 00:44:58 INFO None 4828468: status RUNNING/PENDING
2026-06-13 00:44:58 INFO None 4828469: status RUNNING/PENDING
2026-06-13 00:44:58 INFO None 4828470: status RUNNING/PENDING
2026-06-13 00:44:58 INFO Jobs still running: ['4828464', '4828466', '4828467', '4828468', '4828469', '4828470']. Waiting...
2026-06-13 00:45:13 INFO None 4828449: status FINISHED
2026-06-13 00:45:13 INFO None 4828451: status FINISHED
2026-06-13 00:45:13 INFO None 4828452: status FINISHED
2026-06-13 00:45:13 INFO None 4828453: status FINISHED
2026-06-13 00:45:13 INFO None 4828456: status FINISHED
2026-06-13 00:45:13 INFO None 4828457: status FINISHED
2026-06-13 00:45:13 INFO None 4828459: status FINISHED
2026-06-13 00:45:13 INFO None 4828461: status FINISHED
2026-06-13 00:45:13 INFO None 4828463: status FINISHED
2026-06-13 00:45:13 INFO None 4828464: status FINISHED
2026-06-13 00:45:13 INFO None 4828466: status FINISHED
2026-06-13 00:45:13 INFO None 4828467: status FINISHED
2026-06-13 00:45:13 INFO None 4828468: status RUNNING/PENDING
2026-06-13 00:45:13 INFO None 4828469: status RUNNING/PENDING
2026-06-13 00:45:13 INFO None 4828470: status RUNNING/PENDING
2026-06-13 00:45:13 INFO Jobs still running: ['4828468', '4828469', '4828470']. Waiting...
2026-06-13 00:45:28 INFO None 4828449: status FINISHED
2026-06-13 00:45:28 INFO None 4828451: status FINISHED
2026-06-13 00:45:28 INFO None 4828452: status FINISHED
2026-06-13 00:45:28 INFO None 4828453: status FINISHED
2026-06-13 00:45:28 INFO None 4828456: status FINISHED
2026-06-13 00:45:28 INFO None 4828457: status FINISHED
2026-06-13 00:45:28 INFO None 4828459: status FINISHED
2026-06-13 00:45:28 INFO None 4828461: status FINISHED
2026-06-13 00:45:28 INFO None 4828463: status FINISHED
2026-06-13 00:45:28 INFO None 4828464: status FINISHED
2026-06-13 00:45:28 INFO None 4828466: status FINISHED
2026-06-13 00:45:28 INFO None 4828467: status FINISHED
2026-06-13 00:45:28 INFO None 4828468: status RUNNING/PENDING
2026-06-13 00:45:28 INFO None 4828469: status RUNNING/PENDING
2026-06-13 00:45:28 INFO None 4828470: status RUNNING/PENDING
2026-06-13 00:45:28 INFO Jobs still running: ['4828468', '4828469', '4828470']. Waiting...
2026-06-13 00:45:44 INFO None 4828449: status FINISHED
2026-06-13 00:45:44 INFO None 4828451: status FINISHED
2026-06-13 00:45:44 INFO None 4828452: status FINISHED
2026-06-13 00:45:44 INFO None 4828453: status FINISHED
2026-06-13 00:45:44 INFO None 4828456: status FINISHED
2026-06-13 00:45:44 INFO None 4828457: status FINISHED
2026-06-13 00:45:44 INFO None 4828459: status FINISHED
2026-06-13 00:45:44 INFO None 4828461: status FINISHED
2026-06-13 00:45:44 INFO None 4828463: status FINISHED
2026-06-13 00:45:44 INFO None 4828464: status FINISHED
2026-06-13 00:45:44 INFO None 4828466: status FINISHED
2026-06-13 00:45:44 INFO None 4828467: status FINISHED
2026-06-13 00:45:44 INFO None 4828468: status RUNNING/PENDING
2026-06-13 00:45:44 INFO None 4828469: status RUNNING/PENDING
2026-06-13 00:45:44 INFO None 4828470: status RUNNING/PENDING
2026-06-13 00:45:44 INFO Jobs still running: ['4828468', '4828469', '4828470']. Waiting...
2026-06-13 00:45:59 INFO None 4828449: status FINISHED
2026-06-13 00:46:00 INFO None 4828451: status FINISHED
2026-06-13 00:46:00 INFO None 4828452: status FINISHED
2026-06-13 00:46:00 INFO None 4828453: status FINISHED
2026-06-13 00:46:00 INFO None 4828456: status FINISHED
2026-06-13 00:46:00 INFO None 4828457: status FINISHED
2026-06-13 00:46:00 INFO None 4828459: status FINISHED
2026-06-13 00:46:00 INFO None 4828461: status FINISHED
2026-06-13 00:46:00 INFO None 4828463: status FINISHED
2026-06-13 00:46:00 INFO None 4828464: status FINISHED
2026-06-13 00:46:02 INFO None 4828466: status FINISHED
2026-06-13 00:46:02 INFO None 4828467: status FINISHED
2026-06-13 00:46:02 INFO None 4828468: status RUNNING/PENDING
2026-06-13 00:46:02 INFO None 4828469: status FINISHED
2026-06-13 00:46:02 INFO None 4828470: status RUNNING/PENDING
2026-06-13 00:46:02 INFO Jobs still running: ['4828468', '4828470']. Waiting...
2026-06-13 00:46:17 INFO None 4828449: status FINISHED
2026-06-13 00:46:17 INFO None 4828451: status FINISHED
2026-06-13 00:46:17 INFO None 4828452: status FINISHED
2026-06-13 00:46:17 INFO None 4828453: status FINISHED
2026-06-13 00:46:17 INFO None 4828456: status FINISHED
2026-06-13 00:46:17 INFO None 4828457: status FINISHED
2026-06-13 00:46:17 INFO None 4828459: status FINISHED
2026-06-13 00:46:17 INFO None 4828461: status FINISHED
2026-06-13 00:46:17 INFO None 4828463: status FINISHED
2026-06-13 00:46:17 INFO None 4828464: status FINISHED
2026-06-13 00:46:17 INFO None 4828466: status FINISHED
2026-06-13 00:46:17 INFO None 4828467: status FINISHED
2026-06-13 00:46:17 INFO None 4828468: status RUNNING/PENDING
2026-06-13 00:46:17 INFO None 4828469: status FINISHED
2026-06-13 00:46:17 INFO None 4828470: status RUNNING/PENDING
2026-06-13 00:46:17 INFO Jobs still running: ['4828468', '4828470']. Waiting...
2026-06-13 00:46:32 INFO None 4828449: status FINISHED
2026-06-13 00:46:32 INFO None 4828451: status FINISHED
2026-06-13 00:46:32 INFO None 4828452: status FINISHED
2026-06-13 00:46:32 INFO None 4828453: status FINISHED
2026-06-13 00:46:32 INFO None 4828456: status FINISHED
2026-06-13 00:46:32 INFO None 4828457: status FINISHED
2026-06-13 00:46:32 INFO None 4828459: status FINISHED
2026-06-13 00:46:32 INFO None 4828461: status FINISHED
2026-06-13 00:46:32 INFO None 4828463: status FINISHED
2026-06-13 00:46:32 INFO None 4828464: status FINISHED
2026-06-13 00:46:32 INFO None 4828466: status FINISHED
2026-06-13 00:46:32 INFO None 4828467: status FINISHED
2026-06-13 00:46:32 INFO None 4828468: status FINISHED
2026-06-13 00:46:32 INFO None 4828469: status FINISHED
2026-06-13 00:46:32 INFO None 4828470: status FINISHED
2026-06-13 00:46:32 INFO Jobs ['4828449', '4828451', '4828452', '4828453', '4828456', '4828457', '4828459', '4828461', '4828463', '4828464', '4828466', '4828467', '4828468', '4828469', '4828470'] have finished
2026-06-13 00:46:32 INFO Checking restart files were created ...
2026-06-13 00:46:32 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS1/end.2020020600_1_ENS1.nc(668832435 bytes)
2026-06-13 00:46:32 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS2/end.2020020600_1_ENS2.nc(668832435 bytes)
2026-06-13 00:46:32 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS3/end.2020020600_1_ENS3.nc(668832435 bytes)
2026-06-13 00:46:32 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS4/end.2020020600_1_ENS4.nc(668832435 bytes)
2026-06-13 00:46:32 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS5/end.2020020600_1_ENS5.nc(668832435 bytes)
2026-06-13 00:46:32 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS6/end.2020020600_1_ENS6.nc(668832435 bytes)
2026-06-13 00:46:32 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS7/end.2020020600_1_ENS7.nc(668832435 bytes)
2026-06-13 00:46:32 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS8/end.2020020600_1_ENS8.nc(668832435 bytes)
2026-06-13 00:46:32 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS9/end.2020020600_1_ENS9.nc(668832435 bytes)
2026-06-13 00:46:32 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS10/end.2020020600_1_ENS10.nc(668832435 bytes)
2026-06-13 00:46:32 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS11/end.2020020600_1_ENS11.nc(668832435 bytes)
2026-06-13 00:46:32 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS12/end.2020020600_1_ENS12.nc(668832435 bytes)
2026-06-13 00:46:32 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS13/end.2020020600_1_ENS13.nc(668832435 bytes)
2026-06-13 00:46:32 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS14/end.2020020600_1_ENS14.nc(668832435 bytes)
2026-06-13 00:46:32 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS15/end.2020020600_1_ENS15.nc(668832435 bytes)
2026-06-13 00:46:32 INFO  Run_model() completed successfully.
2026-06-13 00:46:32 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 00:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-13 00:46:32 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 01:00:00 days=153072 seconds=3600
2026-06-13 00:46:32 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-06-13 00:46:32 INFO [TIME] increment current_time 2020-02-06 00:00:00 -> 2020-02-06 01:00:00
2026-06-13 00:46:32 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-13 00:46:32 INFO ---------->>> Running process_satellite_data()
2026-06-13 00:46:32 INFO [DART] No satellite data found, skipping assimilation
2026-06-13 00:46:32 INFO after_assimilation() skipped
2026-06-13 00:46:32 INFO Next run starts from 2020-02-06 01:00:00
2026-06-13 00:46:32 INFO Cycle is DONE; starting a new loop!
2026-06-13 00:46:32 INFO [TIME] step_end current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-13 00:46:32 INFO [TIME] step_start current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-13 00:46:32 INFO [TIME] window start=2020-02-06 01:00:00 end=2020-02-06 09:00:00 run_hours=8 has_assimilation=True
2026-06-13 00:46:32 INFO Creating directories and links for ENS1 to run chimere's parallel part
2026-06-13 00:46:32 INFO Linking EMIS ...
2026-06-13 00:46:32 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc
2026-06-13 00:46:32 INFO >> Checking links...
2026-06-13 00:46:33 INFO >> All links are good for ENS1  ...
2026-06-13 00:46:33 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-13 00:46:34 INFO Hourly dataset computed and listing created
2026-06-13 00:46:45 INFO Hourly dataset computed
2026-06-13 00:46:45 INFO Creating directories and links for ENS2 to run chimere's parallel part
2026-06-13 00:46:45 INFO Linking EMIS ...
2026-06-13 00:46:45 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc
2026-06-13 00:46:45 INFO >> Checking links...
2026-06-13 00:46:46 INFO >> All links are good for ENS2  ...
2026-06-13 00:46:46 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-13 00:46:47 INFO Hourly dataset computed and listing created
2026-06-13 00:46:49 INFO Hourly dataset computed
2026-06-13 00:46:49 INFO Creating directories and links for ENS3 to run chimere's parallel part
2026-06-13 00:46:49 INFO Linking EMIS ...
2026-06-13 00:46:49 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc
2026-06-13 00:46:49 INFO >> Checking links...
2026-06-13 00:46:49 INFO >> All links are good for ENS3  ...
2026-06-13 00:46:49 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-13 00:46:50 INFO Hourly dataset computed and listing created
2026-06-13 00:46:52 INFO Hourly dataset computed
2026-06-13 00:46:52 INFO Creating directories and links for ENS4 to run chimere's parallel part
2026-06-13 00:46:52 INFO Linking EMIS ...
2026-06-13 00:46:52 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc
2026-06-13 00:46:52 INFO >> Checking links...
2026-06-13 00:46:52 INFO >> All links are good for ENS4  ...
2026-06-13 00:46:52 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-13 00:46:53 INFO Hourly dataset computed and listing created
2026-06-13 00:46:56 INFO Hourly dataset computed
2026-06-13 00:46:56 INFO Creating directories and links for ENS5 to run chimere's parallel part
2026-06-13 00:46:56 INFO Linking EMIS ...
2026-06-13 00:46:56 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc
2026-06-13 00:46:56 INFO >> Checking links...
2026-06-13 00:46:56 INFO >> All links are good for ENS5  ...
2026-06-13 00:46:56 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-13 00:46:57 INFO Hourly dataset computed and listing created
2026-06-13 00:46:59 INFO Hourly dataset computed
2026-06-13 00:46:59 INFO Creating directories and links for ENS6 to run chimere's parallel part
2026-06-13 00:46:59 INFO Linking EMIS ...
2026-06-13 00:46:59 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc
2026-06-13 00:46:59 INFO >> Checking links...
2026-06-13 00:46:59 INFO >> All links are good for ENS6  ...
2026-06-13 00:46:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-13 00:47:00 INFO Hourly dataset computed and listing created
2026-06-13 00:47:02 INFO Hourly dataset computed
2026-06-13 00:47:02 INFO Creating directories and links for ENS7 to run chimere's parallel part
2026-06-13 00:47:02 INFO Linking EMIS ...
2026-06-13 00:47:02 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc
2026-06-13 00:47:02 INFO >> Checking links...
2026-06-13 00:47:03 INFO >> All links are good for ENS7  ...
2026-06-13 00:47:03 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-13 00:47:04 INFO Hourly dataset computed and listing created
2026-06-13 00:47:06 INFO Hourly dataset computed
2026-06-13 00:47:06 INFO Creating directories and links for ENS8 to run chimere's parallel part
2026-06-13 00:47:06 INFO Linking EMIS ...
2026-06-13 00:47:06 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc
2026-06-13 00:47:06 INFO >> Checking links...
2026-06-13 00:47:06 INFO >> All links are good for ENS8  ...
2026-06-13 00:47:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-13 00:47:07 INFO Hourly dataset computed and listing created
2026-06-13 00:47:09 INFO Hourly dataset computed
2026-06-13 00:47:09 INFO Creating directories and links for ENS9 to run chimere's parallel part
2026-06-13 00:47:09 INFO Linking EMIS ...
2026-06-13 00:47:09 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc
2026-06-13 00:47:09 INFO >> Checking links...
2026-06-13 00:47:09 INFO >> All links are good for ENS9  ...
2026-06-13 00:47:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-13 00:47:10 INFO Hourly dataset computed and listing created
2026-06-13 00:47:12 INFO Hourly dataset computed
2026-06-13 00:47:12 INFO Creating directories and links for ENS10 to run chimere's parallel part
2026-06-13 00:47:12 INFO Linking EMIS ...
2026-06-13 00:47:12 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc
2026-06-13 00:47:12 INFO >> Checking links...
2026-06-13 00:47:12 INFO >> All links are good for ENS10  ...
2026-06-13 00:47:12 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-13 00:47:13 INFO Hourly dataset computed and listing created
2026-06-13 00:47:15 INFO Hourly dataset computed
2026-06-13 00:47:15 INFO Creating directories and links for ENS11 to run chimere's parallel part
2026-06-13 00:47:15 INFO Linking EMIS ...
2026-06-13 00:47:15 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc
2026-06-13 00:47:15 INFO >> Checking links...
2026-06-13 00:47:16 INFO >> All links are good for ENS11  ...
2026-06-13 00:47:16 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-13 00:47:17 INFO Hourly dataset computed and listing created
2026-06-13 00:47:19 INFO Hourly dataset computed
2026-06-13 00:47:19 INFO Creating directories and links for ENS12 to run chimere's parallel part
2026-06-13 00:47:19 INFO Linking EMIS ...
2026-06-13 00:47:19 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc
2026-06-13 00:47:19 INFO >> Checking links...
2026-06-13 00:47:19 INFO >> All links are good for ENS12  ...
2026-06-13 00:47:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-13 00:47:20 INFO Hourly dataset computed and listing created
2026-06-13 00:47:22 INFO Hourly dataset computed
2026-06-13 00:47:22 INFO Creating directories and links for ENS13 to run chimere's parallel part
2026-06-13 00:47:22 INFO Linking EMIS ...
2026-06-13 00:47:22 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc
2026-06-13 00:47:22 INFO >> Checking links...
2026-06-13 00:47:22 INFO >> All links are good for ENS13  ...
2026-06-13 00:47:22 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-13 00:47:23 INFO Hourly dataset computed and listing created
2026-06-13 00:47:26 INFO Hourly dataset computed
2026-06-13 00:47:26 INFO Creating directories and links for ENS14 to run chimere's parallel part
2026-06-13 00:47:26 INFO Linking EMIS ...
2026-06-13 00:47:26 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc
2026-06-13 00:47:26 INFO >> Checking links...
2026-06-13 00:47:26 INFO >> All links are good for ENS14  ...
2026-06-13 00:47:26 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-13 00:47:27 INFO Hourly dataset computed and listing created
2026-06-13 00:47:29 INFO Hourly dataset computed
2026-06-13 00:47:29 INFO Creating directories and links for ENS15 to run chimere's parallel part
2026-06-13 00:47:29 INFO Linking EMIS ...
2026-06-13 00:47:29 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc
2026-06-13 00:47:29 INFO >> Checking links...
2026-06-13 00:47:30 INFO >> All links are good for ENS15  ...
2026-06-13 00:47:30 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-13 00:47:31 INFO Hourly dataset computed and listing created
2026-06-13 00:47:33 INFO Hourly dataset computed
2026-06-13 00:47:33 INFO ---------->>> Running CHIMERE model from 2020-02-06 01:00:00 to 2020-02-06 09:00:00
2026-06-13 00:47:33 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-13 00:47:33 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS1
2026-06-13 00:47:33 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS1/end.2020020600_1_ENS1.nc
2026-06-13 00:47:33 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-06-13 00:47:33 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-13 00:47:33 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-06-13 00:47:33 INFO Queuing job for member 1...
2026-06-13 00:47:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-13 00:47:33 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-06-13 00:47:33 INFO Found: ['4828484']
2026-06-13 00:47:38 INFO [TGCC-IRENE] Submitted job with ID:['4828484']
2026-06-13 00:47:38 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-13 00:47:38 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS2
2026-06-13 00:47:38 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS2/end.2020020600_1_ENS2.nc
2026-06-13 00:47:38 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-06-13 00:47:38 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-13 00:47:38 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-06-13 00:47:38 INFO Queuing job for member 2...
2026-06-13 00:47:38 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-13 00:47:38 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-06-13 00:47:39 INFO Found: ['4828486']
2026-06-13 00:47:44 INFO [TGCC-IRENE] Submitted job with ID:['4828486']
2026-06-13 00:47:44 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-13 00:47:44 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS3
2026-06-13 00:47:44 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS3/end.2020020600_1_ENS3.nc
2026-06-13 00:47:44 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-06-13 00:47:44 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-13 00:47:44 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-06-13 00:47:44 INFO Queuing job for member 3...
2026-06-13 00:47:44 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-13 00:47:44 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-06-13 00:47:45 INFO Found: ['4828487']
2026-06-13 00:47:50 INFO [TGCC-IRENE] Submitted job with ID:['4828487']
2026-06-13 00:47:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-13 00:47:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS4
2026-06-13 00:47:50 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS4/end.2020020600_1_ENS4.nc
2026-06-13 00:47:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-06-13 00:47:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-13 00:47:50 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-06-13 00:47:50 INFO Queuing job for member 4...
2026-06-13 00:47:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-13 00:47:50 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-06-13 00:47:52 INFO Found: ['4828488']
2026-06-13 00:47:57 INFO [TGCC-IRENE] Submitted job with ID:['4828488']
2026-06-13 00:47:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-13 00:47:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS5
2026-06-13 00:47:57 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS5/end.2020020600_1_ENS5.nc
2026-06-13 00:47:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-06-13 00:47:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-13 00:47:57 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-06-13 00:47:57 INFO Queuing job for member 5...
2026-06-13 00:47:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-13 00:47:57 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-06-13 00:48:00 INFO Found: ['4828489']
2026-06-13 00:48:05 INFO [TGCC-IRENE] Submitted job with ID:['4828489']
2026-06-13 00:48:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-13 00:48:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS6
2026-06-13 00:48:05 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS6/end.2020020600_1_ENS6.nc
2026-06-13 00:48:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-06-13 00:48:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-13 00:48:05 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-06-13 00:48:05 INFO Queuing job for member 6...
2026-06-13 00:48:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-13 00:48:05 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-06-13 00:48:07 INFO Found: ['4828491']
2026-06-13 00:48:12 INFO [TGCC-IRENE] Submitted job with ID:['4828491']
2026-06-13 00:48:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-13 00:48:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS7
2026-06-13 00:48:12 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS7/end.2020020600_1_ENS7.nc
2026-06-13 00:48:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-06-13 00:48:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-13 00:48:12 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-06-13 00:48:12 INFO Queuing job for member 7...
2026-06-13 00:48:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-13 00:48:12 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-06-13 00:48:15 INFO Found: ['4828492']
2026-06-13 00:48:20 INFO [TGCC-IRENE] Submitted job with ID:['4828492']
2026-06-13 00:48:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-13 00:48:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS8
2026-06-13 00:48:20 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS8/end.2020020600_1_ENS8.nc
2026-06-13 00:48:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-06-13 00:48:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-13 00:48:20 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-06-13 00:48:20 INFO Queuing job for member 8...
2026-06-13 00:48:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-13 00:48:20 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-06-13 00:48:22 INFO Found: ['4828493']
2026-06-13 00:48:27 INFO [TGCC-IRENE] Submitted job with ID:['4828493']
2026-06-13 00:48:27 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-13 00:48:27 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS9
2026-06-13 00:48:27 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS9/end.2020020600_1_ENS9.nc
2026-06-13 00:48:27 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-06-13 00:48:27 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-13 00:48:27 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-06-13 00:48:27 INFO Queuing job for member 9...
2026-06-13 00:48:27 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-13 00:48:27 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-06-13 00:48:30 INFO Found: ['4828494']
2026-06-13 00:48:35 INFO [TGCC-IRENE] Submitted job with ID:['4828494']
2026-06-13 00:48:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-13 00:48:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS10
2026-06-13 00:48:35 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS10/end.2020020600_1_ENS10.nc
2026-06-13 00:48:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-06-13 00:48:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-13 00:48:35 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-06-13 00:48:35 INFO Queuing job for member 10...
2026-06-13 00:48:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-13 00:48:35 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-06-13 00:48:37 INFO Found: ['4828495']
2026-06-13 00:48:42 INFO [TGCC-IRENE] Submitted job with ID:['4828495']
2026-06-13 00:48:42 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-13 00:48:42 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS11
2026-06-13 00:48:42 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS11/end.2020020600_1_ENS11.nc
2026-06-13 00:48:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-06-13 00:48:42 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-13 00:48:42 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-06-13 00:48:42 INFO Queuing job for member 11...
2026-06-13 00:48:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-13 00:48:42 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-06-13 00:48:45 INFO Found: ['4828498']
2026-06-13 00:48:50 INFO [TGCC-IRENE] Submitted job with ID:['4828498']
2026-06-13 00:48:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-13 00:48:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS12
2026-06-13 00:48:50 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS12/end.2020020600_1_ENS12.nc
2026-06-13 00:48:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-06-13 00:48:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-13 00:48:50 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-06-13 00:48:50 INFO Queuing job for member 12...
2026-06-13 00:48:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-13 00:48:50 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-06-13 00:48:52 INFO Found: ['4828500']
2026-06-13 00:48:57 INFO [TGCC-IRENE] Submitted job with ID:['4828500']
2026-06-13 00:48:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-13 00:48:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS13
2026-06-13 00:48:57 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS13/end.2020020600_1_ENS13.nc
2026-06-13 00:48:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-06-13 00:48:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-13 00:48:57 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-06-13 00:48:57 INFO Queuing job for member 13...
2026-06-13 00:48:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-13 00:48:57 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-06-13 00:49:00 INFO Found: ['4828503']
2026-06-13 00:49:05 INFO [TGCC-IRENE] Submitted job with ID:['4828503']
2026-06-13 00:49:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-13 00:49:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS14
2026-06-13 00:49:05 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS14/end.2020020600_1_ENS14.nc
2026-06-13 00:49:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-06-13 00:49:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-13 00:49:05 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-06-13 00:49:05 INFO Queuing job for member 14...
2026-06-13 00:49:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-13 00:49:05 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-06-13 00:49:07 INFO Found: ['4828508']
2026-06-13 00:49:12 INFO [TGCC-IRENE] Submitted job with ID:['4828508']
2026-06-13 00:49:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-13 00:49:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS15
2026-06-13 00:49:12 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS15/end.2020020600_1_ENS15.nc
2026-06-13 00:49:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-06-13 00:49:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-13 00:49:12 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-06-13 00:49:12 INFO Queuing job for member 15...
2026-06-13 00:49:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-13 00:49:12 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-06-13 00:49:15 INFO Found: ['4828513']
2026-06-13 00:49:20 INFO [TGCC-IRENE] Submitted job with ID:['4828513']
2026-06-13 00:49:20 INFO Checking job status ...
2026-06-13 00:49:20 INFO None 4828484: status RUNNING/PENDING
2026-06-13 00:49:20 INFO None 4828486: status RUNNING/PENDING
2026-06-13 00:49:20 INFO None 4828487: status RUNNING/PENDING
2026-06-13 00:49:20 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:49:20 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:49:20 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:49:20 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:49:20 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:49:20 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:49:20 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:49:20 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:49:20 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:49:20 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:49:20 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:49:20 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:49:20 INFO Jobs still running: ['4828484', '4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:49:35 INFO None 4828484: status RUNNING/PENDING
2026-06-13 00:49:35 INFO None 4828486: status RUNNING/PENDING
2026-06-13 00:49:36 INFO None 4828487: status RUNNING/PENDING
2026-06-13 00:49:36 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:49:36 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:49:36 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:49:36 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:49:36 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:49:36 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:49:36 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:49:36 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:49:36 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:49:36 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:49:36 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:49:36 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:49:36 INFO Jobs still running: ['4828484', '4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:49:51 INFO None 4828484: status RUNNING/PENDING
2026-06-13 00:49:51 INFO None 4828486: status RUNNING/PENDING
2026-06-13 00:49:51 INFO None 4828487: status RUNNING/PENDING
2026-06-13 00:49:51 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:49:51 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:49:51 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:49:52 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:49:52 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:49:52 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:49:52 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:49:52 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:49:52 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:49:52 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:49:52 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:49:52 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:49:52 INFO Jobs still running: ['4828484', '4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:50:09 INFO None 4828484: status RUNNING/PENDING
2026-06-13 00:50:09 INFO None 4828486: status RUNNING/PENDING
2026-06-13 00:50:09 INFO None 4828487: status RUNNING/PENDING
2026-06-13 00:50:09 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:50:09 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:50:09 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:50:09 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:50:09 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:50:09 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:50:09 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:50:09 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:50:09 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:50:09 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:50:09 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:50:09 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:50:09 INFO Jobs still running: ['4828484', '4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:50:24 INFO None 4828484: status RUNNING/PENDING
2026-06-13 00:50:24 INFO None 4828486: status RUNNING/PENDING
2026-06-13 00:50:24 INFO None 4828487: status RUNNING/PENDING
2026-06-13 00:50:24 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:50:25 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:50:25 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:50:25 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:50:25 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:50:25 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:50:25 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:50:25 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:50:27 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:50:27 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:50:27 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:50:27 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:50:27 INFO Jobs still running: ['4828484', '4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:50:42 INFO None 4828484: status RUNNING/PENDING
2026-06-13 00:50:42 INFO None 4828486: status RUNNING/PENDING
2026-06-13 00:50:42 INFO None 4828487: status RUNNING/PENDING
2026-06-13 00:50:42 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:50:42 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:50:42 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:50:42 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:50:42 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:50:42 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:50:42 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:50:42 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:50:42 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:50:42 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:50:42 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:50:42 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:50:42 INFO Jobs still running: ['4828484', '4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:50:57 INFO None 4828484: status RUNNING/PENDING
2026-06-13 00:50:57 INFO None 4828486: status RUNNING/PENDING
2026-06-13 00:50:57 INFO None 4828487: status RUNNING/PENDING
2026-06-13 00:50:57 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:50:57 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:50:57 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:50:57 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:50:57 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:50:57 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:50:59 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:50:59 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:50:59 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:50:59 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:50:59 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:50:59 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:50:59 INFO Jobs still running: ['4828484', '4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:51:14 INFO None 4828484: status RUNNING/PENDING
2026-06-13 00:51:14 INFO None 4828486: status RUNNING/PENDING
2026-06-13 00:51:15 INFO None 4828487: status RUNNING/PENDING
2026-06-13 00:51:15 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:51:15 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:51:15 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:51:15 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:51:15 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:51:15 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:51:15 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:51:15 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:51:15 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:51:15 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:51:15 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:51:15 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:51:15 INFO Jobs still running: ['4828484', '4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:51:30 INFO None 4828484: status RUNNING/PENDING
2026-06-13 00:51:30 INFO None 4828486: status RUNNING/PENDING
2026-06-13 00:51:30 INFO None 4828487: status RUNNING/PENDING
2026-06-13 00:51:30 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:51:30 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:51:30 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:51:32 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:51:32 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:51:32 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:51:32 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:51:32 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:51:32 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:51:32 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:51:32 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:51:32 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:51:32 INFO Jobs still running: ['4828484', '4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:51:47 INFO None 4828484: status RUNNING/PENDING
2026-06-13 00:51:47 INFO None 4828486: status RUNNING/PENDING
2026-06-13 00:51:47 INFO None 4828487: status RUNNING/PENDING
2026-06-13 00:51:47 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:51:47 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:51:47 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:51:47 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:51:47 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:51:47 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:51:47 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:51:47 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:51:47 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:51:47 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:51:47 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:51:47 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:51:47 INFO Jobs still running: ['4828484', '4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:52:02 INFO None 4828484: status RUNNING/PENDING
2026-06-13 00:52:02 INFO None 4828486: status RUNNING/PENDING
2026-06-13 00:52:02 INFO None 4828487: status RUNNING/PENDING
2026-06-13 00:52:02 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:52:02 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:52:03 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:52:03 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:52:03 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:52:03 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:52:03 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:52:03 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:52:03 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:52:03 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:52:03 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:52:03 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:52:03 INFO Jobs still running: ['4828484', '4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:52:18 INFO None 4828484: status RUNNING/PENDING
2026-06-13 00:52:18 INFO None 4828486: status RUNNING/PENDING
2026-06-13 00:52:18 INFO None 4828487: status RUNNING/PENDING
2026-06-13 00:52:18 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:52:18 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:52:18 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:52:18 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:52:18 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:52:18 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:52:18 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:52:18 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:52:18 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:52:18 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:52:18 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:52:18 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:52:18 INFO Jobs still running: ['4828484', '4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:52:34 INFO None 4828484: status RUNNING/PENDING
2026-06-13 00:52:35 INFO None 4828486: status RUNNING/PENDING
2026-06-13 00:52:35 INFO None 4828487: status RUNNING/PENDING
2026-06-13 00:52:35 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:52:35 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:52:35 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:52:35 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:52:35 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:52:35 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:52:35 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:52:35 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:52:35 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:52:35 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:52:35 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:52:35 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:52:35 INFO Jobs still running: ['4828484', '4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:52:50 INFO None 4828484: status RUNNING/PENDING
2026-06-13 00:52:50 INFO None 4828486: status RUNNING/PENDING
2026-06-13 00:52:50 INFO None 4828487: status RUNNING/PENDING
2026-06-13 00:52:50 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:52:50 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:52:50 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:52:50 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:52:50 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:52:50 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:52:50 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:52:50 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:52:50 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:52:50 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:52:50 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:52:52 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:52:52 INFO Jobs still running: ['4828484', '4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:53:07 INFO None 4828484: status RUNNING/PENDING
2026-06-13 00:53:07 INFO None 4828486: status RUNNING/PENDING
2026-06-13 00:53:07 INFO None 4828487: status RUNNING/PENDING
2026-06-13 00:53:07 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:53:07 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:53:07 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:53:07 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:53:07 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:53:07 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:53:07 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:53:07 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:53:07 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:53:07 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:53:07 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:53:07 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:53:07 INFO Jobs still running: ['4828484', '4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:53:22 INFO None 4828484: status RUNNING/PENDING
2026-06-13 00:53:23 INFO None 4828486: status RUNNING/PENDING
2026-06-13 00:53:23 INFO None 4828487: status RUNNING/PENDING
2026-06-13 00:53:23 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:53:23 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:53:23 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:53:23 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:53:23 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:53:23 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:53:23 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:53:23 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:53:23 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:53:25 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:53:25 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:53:25 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:53:25 INFO Jobs still running: ['4828484', '4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:53:40 INFO None 4828484: status RUNNING/PENDING
2026-06-13 00:53:40 INFO None 4828486: status RUNNING/PENDING
2026-06-13 00:53:40 INFO None 4828487: status RUNNING/PENDING
2026-06-13 00:53:40 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:53:40 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:53:40 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:53:40 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:53:40 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:53:40 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:53:40 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:53:40 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:53:40 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:53:40 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:53:40 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:53:40 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:53:40 INFO Jobs still running: ['4828484', '4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:53:55 INFO None 4828484: status RUNNING/PENDING
2026-06-13 00:53:55 INFO None 4828486: status RUNNING/PENDING
2026-06-13 00:53:55 INFO None 4828487: status RUNNING/PENDING
2026-06-13 00:53:55 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:53:55 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:53:55 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:53:55 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:53:55 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:53:55 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:53:55 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:53:55 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:53:55 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:53:55 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:53:55 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:53:55 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:53:55 INFO Jobs still running: ['4828484', '4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:54:10 INFO None 4828484: status RUNNING/PENDING
2026-06-13 00:54:10 INFO None 4828486: status RUNNING/PENDING
2026-06-13 00:54:10 INFO None 4828487: status RUNNING/PENDING
2026-06-13 00:54:10 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:54:10 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:54:10 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:54:10 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:54:10 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:54:10 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:54:11 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:54:11 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:54:11 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:54:11 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:54:11 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:54:11 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:54:11 INFO Jobs still running: ['4828484', '4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:54:26 INFO None 4828484: status RUNNING/PENDING
2026-06-13 00:54:26 INFO None 4828486: status RUNNING/PENDING
2026-06-13 00:54:26 INFO None 4828487: status RUNNING/PENDING
2026-06-13 00:54:26 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:54:26 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:54:26 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:54:26 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:54:26 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:54:26 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:54:26 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:54:26 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:54:26 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:54:26 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:54:26 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:54:26 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:54:26 INFO Jobs still running: ['4828484', '4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:54:42 INFO None 4828484: status RUNNING/PENDING
2026-06-13 00:54:42 INFO None 4828486: status RUNNING/PENDING
2026-06-13 00:54:42 INFO None 4828487: status RUNNING/PENDING
2026-06-13 00:54:42 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:54:42 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:54:42 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:54:42 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:54:42 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:54:42 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:54:42 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:54:42 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:54:42 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:54:42 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:54:42 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:54:42 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:54:42 INFO Jobs still running: ['4828484', '4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:54:57 INFO None 4828484: status RUNNING/PENDING
2026-06-13 00:54:57 INFO None 4828486: status RUNNING/PENDING
2026-06-13 00:54:57 INFO None 4828487: status RUNNING/PENDING
2026-06-13 00:54:57 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:54:57 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:54:57 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:54:57 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:54:57 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:54:57 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:54:57 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:54:57 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:54:57 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:54:57 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:54:57 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:54:57 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:54:57 INFO Jobs still running: ['4828484', '4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:55:12 INFO None 4828484: status RUNNING/PENDING
2026-06-13 00:55:12 INFO None 4828486: status RUNNING/PENDING
2026-06-13 00:55:14 INFO None 4828487: status RUNNING/PENDING
2026-06-13 00:55:14 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:55:14 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:55:14 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:55:14 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:55:14 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:55:14 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:55:14 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:55:14 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:55:14 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:55:14 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:55:14 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:55:14 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:55:14 INFO Jobs still running: ['4828484', '4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:55:30 INFO None 4828484: status RUNNING/PENDING
2026-06-13 00:55:30 INFO None 4828486: status RUNNING/PENDING
2026-06-13 00:55:30 INFO None 4828487: status RUNNING/PENDING
2026-06-13 00:55:30 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:55:30 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:55:30 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:55:30 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:55:30 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:55:30 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:55:30 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:55:30 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:55:30 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:55:32 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:55:32 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:55:32 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:55:32 INFO Jobs still running: ['4828484', '4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:55:47 INFO None 4828484: status RUNNING/PENDING
2026-06-13 00:55:47 INFO None 4828486: status RUNNING/PENDING
2026-06-13 00:55:47 INFO None 4828487: status RUNNING/PENDING
2026-06-13 00:55:47 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:55:47 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:55:47 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:55:47 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:55:47 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:55:47 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:55:47 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:55:47 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:55:47 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:55:47 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:55:47 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:55:47 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:55:47 INFO Jobs still running: ['4828484', '4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:56:02 INFO None 4828484: status RUNNING/PENDING
2026-06-13 00:56:02 INFO None 4828486: status RUNNING/PENDING
2026-06-13 00:56:02 INFO None 4828487: status RUNNING/PENDING
2026-06-13 00:56:02 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:56:02 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:56:02 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:56:02 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:56:02 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:56:02 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:56:02 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:56:02 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:56:02 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:56:02 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:56:02 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:56:04 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:56:04 INFO Jobs still running: ['4828484', '4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:56:19 INFO None 4828484: status RUNNING/PENDING
2026-06-13 00:56:19 INFO None 4828486: status RUNNING/PENDING
2026-06-13 00:56:19 INFO None 4828487: status RUNNING/PENDING
2026-06-13 00:56:19 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:56:19 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:56:20 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:56:20 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:56:20 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:56:20 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:56:20 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:56:20 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:56:20 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:56:20 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:56:20 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:56:20 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:56:20 INFO Jobs still running: ['4828484', '4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:56:35 INFO None 4828484: status RUNNING/PENDING
2026-06-13 00:56:35 INFO None 4828486: status RUNNING/PENDING
2026-06-13 00:56:35 INFO None 4828487: status RUNNING/PENDING
2026-06-13 00:56:35 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:56:35 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:56:35 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:56:35 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:56:35 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:56:35 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:56:35 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:56:35 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:56:35 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:56:35 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:56:35 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:56:35 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:56:35 INFO Jobs still running: ['4828484', '4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:56:52 INFO None 4828484: status RUNNING/PENDING
2026-06-13 00:56:52 INFO None 4828486: status RUNNING/PENDING
2026-06-13 00:56:52 INFO None 4828487: status RUNNING/PENDING
2026-06-13 00:56:52 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:56:52 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:56:52 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:56:52 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:56:52 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:56:52 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:56:52 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:56:52 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:56:52 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:56:52 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:56:52 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:56:52 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:56:52 INFO Jobs still running: ['4828484', '4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:57:07 INFO None 4828484: status FINISHED
2026-06-13 00:57:07 INFO None 4828486: status RUNNING/PENDING
2026-06-13 00:57:07 INFO None 4828487: status RUNNING/PENDING
2026-06-13 00:57:07 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:57:07 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:57:07 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:57:07 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:57:07 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:57:07 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:57:07 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:57:07 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:57:07 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:57:07 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:57:07 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:57:10 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:57:10 INFO Jobs still running: ['4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:57:25 INFO None 4828484: status FINISHED
2026-06-13 00:57:25 INFO None 4828486: status FINISHED
2026-06-13 00:57:25 INFO None 4828487: status FINISHED
2026-06-13 00:57:25 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:57:25 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:57:25 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:57:25 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:57:25 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:57:25 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:57:25 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:57:25 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:57:25 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:57:25 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:57:25 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:57:25 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:57:25 INFO Jobs still running: ['4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:57:40 INFO None 4828484: status FINISHED
2026-06-13 00:57:40 INFO None 4828486: status FINISHED
2026-06-13 00:57:40 INFO None 4828487: status FINISHED
2026-06-13 00:57:40 INFO None 4828488: status RUNNING/PENDING
2026-06-13 00:57:40 INFO None 4828489: status RUNNING/PENDING
2026-06-13 00:57:40 INFO None 4828491: status RUNNING/PENDING
2026-06-13 00:57:40 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:57:40 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:57:40 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:57:40 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:57:42 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:57:42 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:57:42 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:57:42 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:57:42 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:57:42 INFO Jobs still running: ['4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:57:57 INFO None 4828484: status FINISHED
2026-06-13 00:57:57 INFO None 4828486: status FINISHED
2026-06-13 00:57:57 INFO None 4828487: status FINISHED
2026-06-13 00:57:57 INFO None 4828488: status FINISHED
2026-06-13 00:57:57 INFO None 4828489: status FINISHED
2026-06-13 00:57:57 INFO None 4828491: status FINISHED
2026-06-13 00:57:57 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:57:57 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:57:57 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:57:57 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:57:57 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:57:57 INFO None 4828500: status RUNNING/PENDING
2026-06-13 00:57:57 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:57:57 INFO None 4828508: status RUNNING/PENDING
2026-06-13 00:57:57 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:57:57 INFO Jobs still running: ['4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513']. Waiting...
2026-06-13 00:58:13 INFO None 4828484: status FINISHED
2026-06-13 00:58:13 INFO None 4828486: status FINISHED
2026-06-13 00:58:13 INFO None 4828487: status FINISHED
2026-06-13 00:58:13 INFO None 4828488: status FINISHED
2026-06-13 00:58:13 INFO None 4828489: status FINISHED
2026-06-13 00:58:13 INFO None 4828491: status FINISHED
2026-06-13 00:58:13 INFO None 4828492: status RUNNING/PENDING
2026-06-13 00:58:13 INFO None 4828493: status RUNNING/PENDING
2026-06-13 00:58:13 INFO None 4828494: status RUNNING/PENDING
2026-06-13 00:58:15 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:58:15 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:58:15 INFO None 4828500: status FINISHED
2026-06-13 00:58:15 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:58:15 INFO None 4828508: status FINISHED
2026-06-13 00:58:15 INFO None 4828513: status RUNNING/PENDING
2026-06-13 00:58:15 INFO Jobs still running: ['4828492', '4828493', '4828494', '4828495', '4828498', '4828503', '4828513']. Waiting...
2026-06-13 00:58:30 INFO None 4828484: status FINISHED
2026-06-13 00:58:30 INFO None 4828486: status FINISHED
2026-06-13 00:58:30 INFO None 4828487: status FINISHED
2026-06-13 00:58:30 INFO None 4828488: status FINISHED
2026-06-13 00:58:30 INFO None 4828489: status FINISHED
2026-06-13 00:58:30 INFO None 4828491: status FINISHED
2026-06-13 00:58:30 INFO None 4828492: status FINISHED
2026-06-13 00:58:30 INFO None 4828493: status FINISHED
2026-06-13 00:58:30 INFO None 4828494: status FINISHED
2026-06-13 00:58:30 INFO None 4828495: status RUNNING/PENDING
2026-06-13 00:58:30 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:58:30 INFO None 4828500: status FINISHED
2026-06-13 00:58:30 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:58:30 INFO None 4828508: status FINISHED
2026-06-13 00:58:30 INFO None 4828513: status FINISHED
2026-06-13 00:58:30 INFO Jobs still running: ['4828495', '4828498', '4828503']. Waiting...
2026-06-13 00:58:45 INFO None 4828484: status FINISHED
2026-06-13 00:58:45 INFO None 4828486: status FINISHED
2026-06-13 00:58:45 INFO None 4828487: status FINISHED
2026-06-13 00:58:45 INFO None 4828488: status FINISHED
2026-06-13 00:58:45 INFO None 4828489: status FINISHED
2026-06-13 00:58:45 INFO None 4828491: status FINISHED
2026-06-13 00:58:45 INFO None 4828492: status FINISHED
2026-06-13 00:58:45 INFO None 4828493: status FINISHED
2026-06-13 00:58:45 INFO None 4828494: status FINISHED
2026-06-13 00:58:46 INFO None 4828495: status FINISHED
2026-06-13 00:58:46 INFO None 4828498: status RUNNING/PENDING
2026-06-13 00:58:46 INFO None 4828500: status FINISHED
2026-06-13 00:58:46 INFO None 4828503: status RUNNING/PENDING
2026-06-13 00:58:46 INFO None 4828508: status FINISHED
2026-06-13 00:58:46 INFO None 4828513: status FINISHED
2026-06-13 00:58:46 INFO Jobs still running: ['4828498', '4828503']. Waiting...
2026-06-13 00:59:01 INFO None 4828484: status FINISHED
2026-06-13 00:59:01 INFO None 4828486: status FINISHED
2026-06-13 00:59:01 INFO None 4828487: status FINISHED
2026-06-13 00:59:01 INFO None 4828488: status FINISHED
2026-06-13 00:59:01 INFO None 4828489: status FINISHED
2026-06-13 00:59:01 INFO None 4828491: status FINISHED
2026-06-13 00:59:01 INFO None 4828492: status FINISHED
2026-06-13 00:59:01 INFO None 4828493: status FINISHED
2026-06-13 00:59:01 INFO None 4828494: status FINISHED
2026-06-13 00:59:01 INFO None 4828495: status FINISHED
2026-06-13 00:59:01 INFO None 4828498: status FINISHED
2026-06-13 00:59:01 INFO None 4828500: status FINISHED
2026-06-13 00:59:01 INFO None 4828503: status FINISHED
2026-06-13 00:59:01 INFO None 4828508: status FINISHED
2026-06-13 00:59:01 INFO None 4828513: status FINISHED
2026-06-13 00:59:01 INFO Jobs ['4828484', '4828486', '4828487', '4828488', '4828489', '4828491', '4828492', '4828493', '4828494', '4828495', '4828498', '4828500', '4828503', '4828508', '4828513'] have finished
2026-06-13 00:59:01 INFO Checking restart files were created ...
2026-06-13 00:59:01 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS1/end.2020020601_8_ENS1.nc(3005806795 bytes)
2026-06-13 00:59:01 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS2/end.2020020601_8_ENS2.nc(3005806795 bytes)
2026-06-13 00:59:01 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS3/end.2020020601_8_ENS3.nc(3005806795 bytes)
2026-06-13 00:59:01 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS4/end.2020020601_8_ENS4.nc(3005806795 bytes)
2026-06-13 00:59:01 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS5/end.2020020601_8_ENS5.nc(3005806795 bytes)
2026-06-13 00:59:01 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS6/end.2020020601_8_ENS6.nc(3005806795 bytes)
2026-06-13 00:59:01 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS7/end.2020020601_8_ENS7.nc(3005806795 bytes)
2026-06-13 00:59:01 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS8/end.2020020601_8_ENS8.nc(3005806795 bytes)
2026-06-13 00:59:01 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS9/end.2020020601_8_ENS9.nc(3005806795 bytes)
2026-06-13 00:59:01 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS10/end.2020020601_8_ENS10.nc(3005806795 bytes)
2026-06-13 00:59:01 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS11/end.2020020601_8_ENS11.nc(3005806795 bytes)
2026-06-13 00:59:01 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS12/end.2020020601_8_ENS12.nc(3005806795 bytes)
2026-06-13 00:59:01 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS13/end.2020020601_8_ENS13.nc(3005806795 bytes)
2026-06-13 00:59:01 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS14/end.2020020601_8_ENS14.nc(3005806795 bytes)
2026-06-13 00:59:01 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS15/end.2020020601_8_ENS15.nc(3005806795 bytes)
2026-06-13 00:59:01 INFO  Run_model() completed successfully.
2026-06-13 00:59:01 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-13 00:59:01 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 09:00:00 days=153072 seconds=32400
2026-06-13 00:59:01 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-06-13 00:59:01 INFO [TIME] increment current_time 2020-02-06 01:00:00 -> 2020-02-06 09:00:00
2026-06-13 00:59:01 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-13 00:59:01 INFO ---------->>> Running process_satellite_data()
2026-06-13 00:59:01 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12002.nc
2026-06-13 00:59:01 INFO ---------->>> Running run_obs_converter()
2026-06-13 00:59:01 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_33067_153072.out
2026-06-13 00:59:01 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_33067_153072.out
2026-06-13 00:59:01 INFO ---------->>> Running DART
2026-06-13 00:59:01 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-06-13 00:59:01 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS1/chim_ENS1_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS1/chim_ENS1_2020020609_1_out_toDART.nc
2026-06-13 00:59:01 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS2/chim_ENS2_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS2/chim_ENS2_2020020609_1_out_toDART.nc
2026-06-13 00:59:02 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS3/chim_ENS3_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS3/chim_ENS3_2020020609_1_out_toDART.nc
2026-06-13 00:59:02 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS4/chim_ENS4_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS4/chim_ENS4_2020020609_1_out_toDART.nc
2026-06-13 00:59:03 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS5/chim_ENS5_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS5/chim_ENS5_2020020609_1_out_toDART.nc
2026-06-13 00:59:03 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS6/chim_ENS6_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS6/chim_ENS6_2020020609_1_out_toDART.nc
2026-06-13 00:59:03 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS7/chim_ENS7_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS7/chim_ENS7_2020020609_1_out_toDART.nc
2026-06-13 00:59:04 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS8/chim_ENS8_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS8/chim_ENS8_2020020609_1_out_toDART.nc
2026-06-13 00:59:04 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS9/chim_ENS9_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS9/chim_ENS9_2020020609_1_out_toDART.nc
2026-06-13 00:59:04 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS10/chim_ENS10_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS10/chim_ENS10_2020020609_1_out_toDART.nc
2026-06-13 00:59:05 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS11/chim_ENS11_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS11/chim_ENS11_2020020609_1_out_toDART.nc
2026-06-13 00:59:05 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS12/chim_ENS12_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS12/chim_ENS12_2020020609_1_out_toDART.nc
2026-06-13 00:59:05 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS13/chim_ENS13_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS13/chim_ENS13_2020020609_1_out_toDART.nc
2026-06-13 00:59:06 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS14/chim_ENS14_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS14/chim_ENS14_2020020609_1_out_toDART.nc
2026-06-13 00:59:06 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS15/chim_ENS15_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_0629_15m_low_v2/ENS15/chim_ENS15_2020020609_1_out_toDART.nc
2026-06-13 00:59:06 INFO Replacement input_template.nml → input.nml completed successfully.
2026-06-13 00:59:06 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-06-13 00:59:06 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-06-13 00:59:07 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-06-13 00:59:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-06-13 00:59:07 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-06-13 00:59:13 INFO Found: []
2026-06-13 00:59:13 INFO No job id returned by command ./run_filter.bsh
2026-06-13 00:59:13 INFO No monitoring will be performed
2026-06-13 00:59:13 INFO Moving DART output files to analysis and preassim directories for date 2020020609 if present ...
2026-06-13 00:59:13 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/analysis/2020020609'
2026-06-13 00:59:13 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/analysis/2020020609'
2026-06-13 00:59:13 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/analysis/2020020609'
2026-06-13 00:59:13 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/analysis/2020020609'
2026-06-13 00:59:13 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/analysis/2020020609'
2026-06-13 00:59:13 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/preassim/2020020609'
2026-06-13 00:59:13 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/analysis/2020020609'
2026-06-13 00:59:13 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/analysis/2020020609'
2026-06-13 00:59:13 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/analysis/2020020609'
2026-06-13 00:59:13 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/analysis/2020020609'
2026-06-13 00:59:13 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/analysis/2020020609'
2026-06-13 00:59:13 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/preassim/2020020609'
2026-06-13 00:59:13 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/preassim/2020020609'
2026-06-13 00:59:13 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/analysis/2020020609'
2026-06-13 00:59:13 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/preassim/2020020609'
2026-06-13 00:59:13 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/preassim/2020020609'
2026-06-13 00:59:13 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/analysis/2020020609'
2026-06-13 00:59:13 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/preassim/2020020609'
2026-06-13 00:59:13 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/analysis/2020020609'
2026-06-13 00:59:13 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/analysis/2020020609'
2026-06-13 00:59:13 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/preassim/2020020609'
2026-06-13 00:59:13 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/preassim/2020020609'
2026-06-13 00:59:13 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/analysis/2020020609'
2026-06-13 00:59:13 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/analysis/2020020609'
2026-06-13 00:59:13 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/preassim/2020020609'
2026-06-13 00:59:13 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/preassim/2020020609'
2026-06-13 00:59:13 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/preassim/2020020609'
2026-06-13 00:59:13 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/preassim/2020020609'
2026-06-13 00:59:13 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/preassim/2020020609'
2026-06-13 00:59:13 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/preassim/2020020609'
2026-06-13 00:59:13 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/preassim/2020020609'
2026-06-13 00:59:13 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/analysis/2020020609'
2026-06-13 00:59:13 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/preassim/2020020609'
2026-06-13 00:59:13 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_0629_15m_low_v2/preassim/2020020609'
2026-06-13 00:59:13 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-06-13 00:59:16 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-06-13 00:59:16 INFO run_dart() is DONE.
2026-06-13 00:59:16 INFO ---------->>> Running update_pollutant_in_end()
Traceback (most recent call last):
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/site-packages/xarray/core/dataset.py", line 1154, in _construct_dataarray
    variable = self._variables[name]
               ~~~~~~~~~~~~~~~^^^^^^
KeyError: 'EMISA'

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/site-packages/xarray/core/dataset.py", line 1261, in __getitem__
    return self._construct_dataarray(key)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/site-packages/xarray/core/dataset.py", line 1156, in _construct_dataarray
    _, name, variable = _get_virtual_variable(self._variables, name, self.sizes)
                        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/site-packages/xarray/core/dataset_utils.py", line 79, in _get_virtual_variable
    raise KeyError(key)
KeyError: 'EMISA'

The above exception was the direct cause of the following exception:

Traceback (most recent call last):
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/main.py", line 103, in <module>
    pipeline.run_pipeline()
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/pipelines/base_pipeline.py", line 162, in run_pipeline
    self.after_assimilation()
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/pipelines/chimere2023/pipeline.py", line 589, in after_assimilation
    update_pollutant_in_end(dart_file=self.paths.dart_filter_output_list_file(mem, date_ymdH), 
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/orchestrator_utils.py", line 1291, in update_pollutant_in_end
    val_prev = end_ds[pollutant].isel(Time=-2).load().values
               ~~~~~~^^^^^^^^^^^
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/site-packages/xarray/core/dataset.py", line 1274, in __getitem__
    raise KeyError(message) from e
KeyError: "No variable named 'EMISA'. Did you mean one of ('MSA', 'MSIA')?"
+ exit 0
