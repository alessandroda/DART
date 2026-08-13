+ SCRIPT_PID=451316
+ /bin/bash -x /tmp/tmp.utCasKFWnG
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
+ python -u main.py -c config/config_irene_IM_cp2.yaml
2026-07-14 11:15:30 INFO 
███    ███ ██ ███    ███ ███████ ███████ ██
████  ████ ██ ████  ████ ██      ██      ██
██ ████ ██ ██ ██ ████ ██ █████   ███████ ██
██  ██  ██ ██ ██  ██  ██ ██           ██ ██
██      ██ ██ ██      ██ ███████ ███████ ██



2026-07-14 11:15:30 INFO [PIPELINE] =======================================
2026-07-14 11:15:30 INFO [PIPELINE] Starting chimere–DART orchestrator
2026-07-14 11:15:30 INFO [PIPELINE] Config file: config/config_irene_IM_cp2.yaml
2026-07-14 11:15:30 INFO [PIPELINE] Run dir: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart
2026-07-14 11:15:30 INFO [PIPELINE] Log file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/mimesi_orchestrator_logs/chimere_DART_20260714_111530.log
2026-07-14 11:15:30 INFO [PIPELINE] =======================================
2026-07-14 11:15:30 INFO Running assimilation with model_type=ModelType.CHIMERE
2026-07-14 11:15:30 INFO Using scheduler=Scheduler.SLURM, queue=rome
2026-07-14 11:15:30 INFO [STEP] ---- TIME LOOP START ----
2026-07-14 11:15:30 INFO [TIME] step_start current_time=2020-02-06 00:00:00 simulated_time=None dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 11:15:30 INFO [TIME] window start=2020-02-06 00:00:00 end=2020-02-06 01:00:00 run_hours=1 has_assimilation=False
2026-07-14 11:15:30 INFO Creating directories and links for ENS1 to run chimere's parallel part
2026-07-14 11:15:30 INFO Copying EMIS ...
2026-07-14 11:15:31 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens01.nc
2026-07-14 11:15:31 INFO Linking first END ...
2026-07-14 11:15:31 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020500_24_ENS1.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-14 11:15:31 INFO >> Checking links...
2026-07-14 11:15:31 INFO >> All links are good for ENS1  ...
2026-07-14 11:15:31 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 11:15:39 INFO Hourly dataset computed and listing created
2026-07-14 11:15:43 INFO Hourly dataset computed
2026-07-14 11:15:43 INFO Creating directories and links for ENS2 to run chimere's parallel part
2026-07-14 11:15:43 INFO Copying EMIS ...
2026-07-14 11:15:44 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens02.nc
2026-07-14 11:15:44 INFO Linking first END ...
2026-07-14 11:15:44 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020500_24_ENS2.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-14 11:15:44 INFO >> Checking links...
2026-07-14 11:15:44 INFO >> All links are good for ENS2  ...
2026-07-14 11:15:44 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 11:15:45 INFO Hourly dataset computed and listing created
2026-07-14 11:15:45 INFO Hourly dataset computed
2026-07-14 11:15:45 INFO Creating directories and links for ENS3 to run chimere's parallel part
2026-07-14 11:15:45 INFO Copying EMIS ...
2026-07-14 11:15:46 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens03.nc
2026-07-14 11:15:46 INFO Linking first END ...
2026-07-14 11:15:46 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020500_24_ENS3.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-14 11:15:46 INFO >> Checking links...
2026-07-14 11:15:46 INFO >> All links are good for ENS3  ...
2026-07-14 11:15:46 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 11:15:47 INFO Hourly dataset computed and listing created
2026-07-14 11:15:47 INFO Hourly dataset computed
2026-07-14 11:15:47 INFO Creating directories and links for ENS4 to run chimere's parallel part
2026-07-14 11:15:47 INFO Copying EMIS ...
2026-07-14 11:15:48 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens04.nc
2026-07-14 11:15:48 INFO Linking first END ...
2026-07-14 11:15:48 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020500_24_ENS4.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-14 11:15:48 INFO >> Checking links...
2026-07-14 11:15:48 INFO >> All links are good for ENS4  ...
2026-07-14 11:15:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 11:15:49 INFO Hourly dataset computed and listing created
2026-07-14 11:15:49 INFO Hourly dataset computed
2026-07-14 11:15:49 INFO Creating directories and links for ENS5 to run chimere's parallel part
2026-07-14 11:15:49 INFO Copying EMIS ...
2026-07-14 11:15:50 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens05.nc
2026-07-14 11:15:50 INFO Linking first END ...
2026-07-14 11:15:50 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020500_24_ENS5.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-14 11:15:50 INFO >> Checking links...
2026-07-14 11:15:50 INFO >> All links are good for ENS5  ...
2026-07-14 11:15:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 11:15:51 INFO Hourly dataset computed and listing created
2026-07-14 11:15:51 INFO Hourly dataset computed
2026-07-14 11:15:51 INFO Creating directories and links for ENS6 to run chimere's parallel part
2026-07-14 11:15:51 INFO Copying EMIS ...
2026-07-14 11:15:52 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens06.nc
2026-07-14 11:15:52 INFO Linking first END ...
2026-07-14 11:15:52 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020500_24_ENS6.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-14 11:15:52 INFO >> Checking links...
2026-07-14 11:15:52 INFO >> All links are good for ENS6  ...
2026-07-14 11:15:52 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 11:15:53 INFO Hourly dataset computed and listing created
2026-07-14 11:15:53 INFO Hourly dataset computed
2026-07-14 11:15:53 INFO Creating directories and links for ENS7 to run chimere's parallel part
2026-07-14 11:15:53 INFO Copying EMIS ...
2026-07-14 11:15:54 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens07.nc
2026-07-14 11:15:54 INFO Linking first END ...
2026-07-14 11:15:54 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020500_24_ENS7.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-14 11:15:54 INFO >> Checking links...
2026-07-14 11:15:54 INFO >> All links are good for ENS7  ...
2026-07-14 11:15:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 11:15:55 INFO Hourly dataset computed and listing created
2026-07-14 11:15:55 INFO Hourly dataset computed
2026-07-14 11:15:55 INFO Creating directories and links for ENS8 to run chimere's parallel part
2026-07-14 11:15:55 INFO Copying EMIS ...
2026-07-14 11:15:56 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens08.nc
2026-07-14 11:15:56 INFO Linking first END ...
2026-07-14 11:15:56 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020500_24_ENS8.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-14 11:15:56 INFO >> Checking links...
2026-07-14 11:15:56 INFO >> All links are good for ENS8  ...
2026-07-14 11:15:56 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 11:15:57 INFO Hourly dataset computed and listing created
2026-07-14 11:15:57 INFO Hourly dataset computed
2026-07-14 11:15:57 INFO Creating directories and links for ENS9 to run chimere's parallel part
2026-07-14 11:15:57 INFO Copying EMIS ...
2026-07-14 11:15:58 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens09.nc
2026-07-14 11:15:58 INFO Linking first END ...
2026-07-14 11:15:58 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020500_24_ENS9.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-14 11:15:58 INFO >> Checking links...
2026-07-14 11:15:58 INFO >> All links are good for ENS9  ...
2026-07-14 11:15:58 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 11:15:59 INFO Hourly dataset computed and listing created
2026-07-14 11:16:00 INFO Hourly dataset computed
2026-07-14 11:16:00 INFO Creating directories and links for ENS10 to run chimere's parallel part
2026-07-14 11:16:00 INFO Copying EMIS ...
2026-07-14 11:16:00 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens10.nc
2026-07-14 11:16:00 INFO Linking first END ...
2026-07-14 11:16:00 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020500_24_ENS10.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-14 11:16:00 INFO >> Checking links...
2026-07-14 11:16:00 INFO >> All links are good for ENS10  ...
2026-07-14 11:16:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 11:16:01 INFO Hourly dataset computed and listing created
2026-07-14 11:16:02 INFO Hourly dataset computed
2026-07-14 11:16:02 INFO Creating directories and links for ENS11 to run chimere's parallel part
2026-07-14 11:16:02 INFO Copying EMIS ...
2026-07-14 11:18:10 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens11.nc
2026-07-14 11:18:10 INFO Linking first END ...
2026-07-14 11:18:10 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020500_24_ENS11.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-14 11:18:10 INFO >> Checking links...
2026-07-14 11:18:10 INFO >> All links are good for ENS11  ...
2026-07-14 11:18:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 11:18:11 INFO Hourly dataset computed and listing created
2026-07-14 11:18:12 INFO Hourly dataset computed
2026-07-14 11:18:12 INFO Creating directories and links for ENS12 to run chimere's parallel part
2026-07-14 11:18:12 INFO Copying EMIS ...
2026-07-14 11:18:12 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens12.nc
2026-07-14 11:18:12 INFO Linking first END ...
2026-07-14 11:18:12 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020500_24_ENS12.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-14 11:18:12 INFO >> Checking links...
2026-07-14 11:18:12 INFO >> All links are good for ENS12  ...
2026-07-14 11:18:12 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 11:18:13 INFO Hourly dataset computed and listing created
2026-07-14 11:18:14 INFO Hourly dataset computed
2026-07-14 11:18:14 INFO Creating directories and links for ENS13 to run chimere's parallel part
2026-07-14 11:18:14 INFO Copying EMIS ...
2026-07-14 11:18:14 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens13.nc
2026-07-14 11:18:14 INFO Linking first END ...
2026-07-14 11:18:14 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020500_24_ENS13.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-14 11:18:14 INFO >> Checking links...
2026-07-14 11:18:14 INFO >> All links are good for ENS13  ...
2026-07-14 11:18:14 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 11:18:15 INFO Hourly dataset computed and listing created
2026-07-14 11:18:16 INFO Hourly dataset computed
2026-07-14 11:18:16 INFO Creating directories and links for ENS14 to run chimere's parallel part
2026-07-14 11:18:16 INFO Copying EMIS ...
2026-07-14 11:18:16 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens14.nc
2026-07-14 11:18:16 INFO Linking first END ...
2026-07-14 11:18:16 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020500_24_ENS14.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-14 11:18:16 INFO >> Checking links...
2026-07-14 11:18:16 INFO >> All links are good for ENS14  ...
2026-07-14 11:18:16 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 11:18:17 INFO Hourly dataset computed and listing created
2026-07-14 11:18:18 INFO Hourly dataset computed
2026-07-14 11:18:18 INFO Creating directories and links for ENS15 to run chimere's parallel part
2026-07-14 11:18:18 INFO Copying EMIS ...
2026-07-14 11:18:18 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens15.nc
2026-07-14 11:18:18 INFO Linking first END ...
2026-07-14 11:18:18 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020500_24_ENS15.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-14 11:18:18 INFO >> Checking links...
2026-07-14 11:18:18 INFO >> All links are good for ENS15  ...
2026-07-14 11:18:18 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 11:18:19 INFO Hourly dataset computed and listing created
2026-07-14 11:18:20 INFO Hourly dataset computed
2026-07-14 11:18:20 INFO ---------->>> Running CHIMERE model from 2020-02-06 00:00:00 to 2020-02-06 01:00:00
2026-07-14 11:18:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 11:18:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1
2026-07-14 11:18:20 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020500_24_ENS1.nc
2026-07-14 11:18:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-14 11:18:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 11:18:20 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-14 11:18:20 INFO Queuing job for member 1...
2026-07-14 11:18:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 11:18:20 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-14 11:18:21 INFO Found: ['5155171']
2026-07-14 11:18:26 INFO [TGCC-IRENE] Submitted job with ID:['5155171']
2026-07-14 11:18:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 11:18:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2
2026-07-14 11:18:26 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020500_24_ENS2.nc
2026-07-14 11:18:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-14 11:18:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 11:18:26 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-14 11:18:26 INFO Queuing job for member 2...
2026-07-14 11:18:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 11:18:26 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-14 11:18:26 INFO Found: ['5155172']
2026-07-14 11:18:31 INFO [TGCC-IRENE] Submitted job with ID:['5155172']
2026-07-14 11:18:31 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 11:18:31 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3
2026-07-14 11:18:31 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020500_24_ENS3.nc
2026-07-14 11:18:31 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-14 11:18:31 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 11:18:31 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-14 11:18:31 INFO Queuing job for member 3...
2026-07-14 11:18:31 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 11:18:31 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-14 11:18:32 INFO Found: ['5155174']
2026-07-14 11:18:37 INFO [TGCC-IRENE] Submitted job with ID:['5155174']
2026-07-14 11:18:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 11:18:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4
2026-07-14 11:18:37 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020500_24_ENS4.nc
2026-07-14 11:18:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-14 11:18:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 11:18:37 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-14 11:18:37 INFO Queuing job for member 4...
2026-07-14 11:18:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 11:18:37 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-14 11:18:38 INFO Found: ['5155175']
2026-07-14 11:18:43 INFO [TGCC-IRENE] Submitted job with ID:['5155175']
2026-07-14 11:18:43 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 11:18:43 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5
2026-07-14 11:18:43 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020500_24_ENS5.nc
2026-07-14 11:18:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-14 11:18:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 11:18:43 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-14 11:18:43 INFO Queuing job for member 5...
2026-07-14 11:18:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 11:18:43 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-14 11:18:44 INFO Found: ['5155176']
2026-07-14 11:18:49 INFO [TGCC-IRENE] Submitted job with ID:['5155176']
2026-07-14 11:18:49 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 11:18:49 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6
2026-07-14 11:18:49 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020500_24_ENS6.nc
2026-07-14 11:18:49 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-14 11:18:49 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 11:18:49 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-14 11:18:49 INFO Queuing job for member 6...
2026-07-14 11:18:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 11:18:49 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-14 11:18:49 INFO Found: ['5155177']
2026-07-14 11:18:54 INFO [TGCC-IRENE] Submitted job with ID:['5155177']
2026-07-14 11:18:54 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 11:18:54 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7
2026-07-14 11:18:54 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020500_24_ENS7.nc
2026-07-14 11:18:54 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-14 11:18:54 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 11:18:54 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-14 11:18:54 INFO Queuing job for member 7...
2026-07-14 11:18:54 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 11:18:54 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-14 11:18:55 INFO Found: ['5155178']
2026-07-14 11:19:00 INFO [TGCC-IRENE] Submitted job with ID:['5155178']
2026-07-14 11:19:00 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 11:19:00 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8
2026-07-14 11:19:00 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020500_24_ENS8.nc
2026-07-14 11:19:00 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-14 11:19:00 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 11:19:00 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-14 11:19:00 INFO Queuing job for member 8...
2026-07-14 11:19:00 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 11:19:00 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-14 11:19:01 INFO Found: ['5155182']
2026-07-14 11:19:06 INFO [TGCC-IRENE] Submitted job with ID:['5155182']
2026-07-14 11:19:06 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 11:19:06 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9
2026-07-14 11:19:06 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020500_24_ENS9.nc
2026-07-14 11:19:06 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-14 11:19:06 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 11:19:06 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-14 11:19:06 INFO Queuing job for member 9...
2026-07-14 11:19:06 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 11:19:06 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-14 11:19:07 INFO Found: ['5155184']
2026-07-14 11:19:12 INFO [TGCC-IRENE] Submitted job with ID:['5155184']
2026-07-14 11:19:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 11:19:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10
2026-07-14 11:19:12 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020500_24_ENS10.nc
2026-07-14 11:19:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-14 11:19:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 11:19:12 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-14 11:19:12 INFO Queuing job for member 10...
2026-07-14 11:19:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 11:19:12 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-14 11:19:12 INFO Found: ['5155185']
2026-07-14 11:19:17 INFO [TGCC-IRENE] Submitted job with ID:['5155185']
2026-07-14 11:19:17 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 11:19:17 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11
2026-07-14 11:19:17 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020500_24_ENS11.nc
2026-07-14 11:19:17 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-14 11:19:17 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 11:19:17 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-14 11:19:17 INFO Queuing job for member 11...
2026-07-14 11:19:17 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 11:19:17 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-14 11:19:18 INFO Found: ['5155187']
2026-07-14 11:19:23 INFO [TGCC-IRENE] Submitted job with ID:['5155187']
2026-07-14 11:19:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 11:19:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12
2026-07-14 11:19:23 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020500_24_ENS12.nc
2026-07-14 11:19:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-14 11:19:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 11:19:23 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-14 11:19:23 INFO Queuing job for member 12...
2026-07-14 11:19:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 11:19:23 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-14 11:19:24 INFO Found: ['5155188']
2026-07-14 11:19:29 INFO [TGCC-IRENE] Submitted job with ID:['5155188']
2026-07-14 11:19:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 11:19:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13
2026-07-14 11:19:29 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020500_24_ENS13.nc
2026-07-14 11:19:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-14 11:19:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 11:19:29 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-14 11:19:29 INFO Queuing job for member 13...
2026-07-14 11:19:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 11:19:29 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-14 11:19:30 INFO Found: ['5155189']
2026-07-14 11:19:35 INFO [TGCC-IRENE] Submitted job with ID:['5155189']
2026-07-14 11:19:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 11:19:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14
2026-07-14 11:19:35 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020500_24_ENS14.nc
2026-07-14 11:19:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-14 11:19:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 11:19:35 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-14 11:19:35 INFO Queuing job for member 14...
2026-07-14 11:19:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 11:19:35 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-14 11:19:36 INFO Found: ['5155191']
2026-07-14 11:19:41 INFO [TGCC-IRENE] Submitted job with ID:['5155191']
2026-07-14 11:19:41 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 11:19:41 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15
2026-07-14 11:19:41 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020500_24_ENS15.nc
2026-07-14 11:19:41 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-14 11:19:41 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 11:19:41 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-14 11:19:41 INFO Queuing job for member 15...
2026-07-14 11:19:41 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 11:19:41 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-14 11:19:41 INFO Found: ['5155192']
2026-07-14 11:19:46 INFO [TGCC-IRENE] Submitted job with ID:['5155192']
2026-07-14 11:19:46 INFO Checking job status ...
2026-07-14 11:19:46 INFO None 5155171: status RUNNING/PENDING
2026-07-14 11:19:46 INFO None 5155172: status RUNNING/PENDING
2026-07-14 11:19:46 INFO None 5155174: status RUNNING/PENDING
2026-07-14 11:19:46 INFO None 5155175: status RUNNING/PENDING
2026-07-14 11:19:46 INFO None 5155176: status RUNNING/PENDING
2026-07-14 11:19:46 INFO None 5155177: status RUNNING/PENDING
2026-07-14 11:19:46 INFO None 5155178: status RUNNING/PENDING
2026-07-14 11:19:46 INFO None 5155182: status RUNNING/PENDING
2026-07-14 11:19:47 INFO None 5155184: status RUNNING/PENDING
2026-07-14 11:19:47 INFO None 5155185: status RUNNING/PENDING
2026-07-14 11:19:47 INFO None 5155187: status RUNNING/PENDING
2026-07-14 11:19:47 INFO None 5155188: status RUNNING/PENDING
2026-07-14 11:19:47 INFO None 5155189: status RUNNING/PENDING
2026-07-14 11:19:47 INFO None 5155191: status RUNNING/PENDING
2026-07-14 11:19:47 INFO None 5155192: status RUNNING/PENDING
2026-07-14 11:19:47 INFO Jobs still running: ['5155171', '5155172', '5155174', '5155175', '5155176', '5155177', '5155178', '5155182', '5155184', '5155185', '5155187', '5155188', '5155189', '5155191', '5155192']. Waiting...
2026-07-14 11:20:02 INFO None 5155171: status RUNNING/PENDING
2026-07-14 11:20:02 INFO None 5155172: status RUNNING/PENDING
2026-07-14 11:20:02 INFO None 5155174: status RUNNING/PENDING
2026-07-14 11:20:02 INFO None 5155175: status RUNNING/PENDING
2026-07-14 11:20:02 INFO None 5155176: status RUNNING/PENDING
2026-07-14 11:20:02 INFO None 5155177: status RUNNING/PENDING
2026-07-14 11:20:02 INFO None 5155178: status RUNNING/PENDING
2026-07-14 11:20:02 INFO None 5155182: status RUNNING/PENDING
2026-07-14 11:20:02 INFO None 5155184: status RUNNING/PENDING
2026-07-14 11:20:02 INFO None 5155185: status RUNNING/PENDING
2026-07-14 11:20:02 INFO None 5155187: status RUNNING/PENDING
2026-07-14 11:20:02 INFO None 5155188: status RUNNING/PENDING
2026-07-14 11:20:02 INFO None 5155189: status RUNNING/PENDING
2026-07-14 11:20:02 INFO None 5155191: status RUNNING/PENDING
2026-07-14 11:20:02 INFO None 5155192: status RUNNING/PENDING
2026-07-14 11:20:02 INFO Jobs still running: ['5155171', '5155172', '5155174', '5155175', '5155176', '5155177', '5155178', '5155182', '5155184', '5155185', '5155187', '5155188', '5155189', '5155191', '5155192']. Waiting...
2026-07-14 11:20:17 INFO None 5155171: status RUNNING/PENDING
2026-07-14 11:20:17 INFO None 5155172: status RUNNING/PENDING
2026-07-14 11:20:17 INFO None 5155174: status RUNNING/PENDING
2026-07-14 11:20:17 INFO None 5155175: status RUNNING/PENDING
2026-07-14 11:20:17 INFO None 5155176: status RUNNING/PENDING
2026-07-14 11:20:17 INFO None 5155177: status RUNNING/PENDING
2026-07-14 11:20:17 INFO None 5155178: status RUNNING/PENDING
2026-07-14 11:20:17 INFO None 5155182: status RUNNING/PENDING
2026-07-14 11:20:17 INFO None 5155184: status RUNNING/PENDING
2026-07-14 11:20:17 INFO None 5155185: status RUNNING/PENDING
2026-07-14 11:20:17 INFO None 5155187: status RUNNING/PENDING
2026-07-14 11:20:17 INFO None 5155188: status RUNNING/PENDING
2026-07-14 11:20:17 INFO None 5155189: status RUNNING/PENDING
2026-07-14 11:20:17 INFO None 5155191: status RUNNING/PENDING
2026-07-14 11:20:17 INFO None 5155192: status RUNNING/PENDING
2026-07-14 11:20:17 INFO Jobs still running: ['5155171', '5155172', '5155174', '5155175', '5155176', '5155177', '5155178', '5155182', '5155184', '5155185', '5155187', '5155188', '5155189', '5155191', '5155192']. Waiting...
2026-07-14 11:20:32 INFO None 5155171: status RUNNING/PENDING
2026-07-14 11:20:32 INFO None 5155172: status RUNNING/PENDING
2026-07-14 11:20:32 INFO None 5155174: status RUNNING/PENDING
2026-07-14 11:20:32 INFO None 5155175: status RUNNING/PENDING
2026-07-14 11:20:32 INFO None 5155176: status RUNNING/PENDING
2026-07-14 11:20:32 INFO None 5155177: status RUNNING/PENDING
2026-07-14 11:20:32 INFO None 5155178: status RUNNING/PENDING
2026-07-14 11:20:33 INFO None 5155182: status RUNNING/PENDING
2026-07-14 11:20:33 INFO None 5155184: status RUNNING/PENDING
2026-07-14 11:20:33 INFO None 5155185: status RUNNING/PENDING
2026-07-14 11:20:33 INFO None 5155187: status RUNNING/PENDING
2026-07-14 11:20:33 INFO None 5155188: status RUNNING/PENDING
2026-07-14 11:20:33 INFO None 5155189: status RUNNING/PENDING
2026-07-14 11:20:33 INFO None 5155191: status RUNNING/PENDING
2026-07-14 11:20:33 INFO None 5155192: status RUNNING/PENDING
2026-07-14 11:20:33 INFO Jobs still running: ['5155171', '5155172', '5155174', '5155175', '5155176', '5155177', '5155178', '5155182', '5155184', '5155185', '5155187', '5155188', '5155189', '5155191', '5155192']. Waiting...
2026-07-14 11:20:48 INFO None 5155171: status RUNNING/PENDING
2026-07-14 11:20:48 INFO None 5155172: status RUNNING/PENDING
2026-07-14 11:20:48 INFO None 5155174: status RUNNING/PENDING
2026-07-14 11:20:48 INFO None 5155175: status RUNNING/PENDING
2026-07-14 11:20:48 INFO None 5155176: status RUNNING/PENDING
2026-07-14 11:20:48 INFO None 5155177: status RUNNING/PENDING
2026-07-14 11:20:48 INFO None 5155178: status RUNNING/PENDING
2026-07-14 11:20:48 INFO None 5155182: status RUNNING/PENDING
2026-07-14 11:20:48 INFO None 5155184: status RUNNING/PENDING
2026-07-14 11:20:48 INFO None 5155185: status RUNNING/PENDING
2026-07-14 11:20:48 INFO None 5155187: status RUNNING/PENDING
2026-07-14 11:20:48 INFO None 5155188: status RUNNING/PENDING
2026-07-14 11:20:48 INFO None 5155189: status RUNNING/PENDING
2026-07-14 11:20:48 INFO None 5155191: status RUNNING/PENDING
2026-07-14 11:20:48 INFO None 5155192: status RUNNING/PENDING
2026-07-14 11:20:48 INFO Jobs still running: ['5155171', '5155172', '5155174', '5155175', '5155176', '5155177', '5155178', '5155182', '5155184', '5155185', '5155187', '5155188', '5155189', '5155191', '5155192']. Waiting...
2026-07-14 11:21:03 INFO None 5155171: status RUNNING/PENDING
2026-07-14 11:21:03 INFO None 5155172: status RUNNING/PENDING
2026-07-14 11:21:03 INFO None 5155174: status RUNNING/PENDING
2026-07-14 11:21:03 INFO None 5155175: status RUNNING/PENDING
2026-07-14 11:21:03 INFO None 5155176: status RUNNING/PENDING
2026-07-14 11:21:03 INFO None 5155177: status RUNNING/PENDING
2026-07-14 11:21:03 INFO None 5155178: status RUNNING/PENDING
2026-07-14 11:21:03 INFO None 5155182: status RUNNING/PENDING
2026-07-14 11:21:03 INFO None 5155184: status RUNNING/PENDING
2026-07-14 11:21:03 INFO None 5155185: status RUNNING/PENDING
2026-07-14 11:21:03 INFO None 5155187: status RUNNING/PENDING
2026-07-14 11:21:03 INFO None 5155188: status RUNNING/PENDING
2026-07-14 11:21:03 INFO None 5155189: status RUNNING/PENDING
2026-07-14 11:21:03 INFO None 5155191: status RUNNING/PENDING
2026-07-14 11:21:03 INFO None 5155192: status RUNNING/PENDING
2026-07-14 11:21:03 INFO Jobs still running: ['5155171', '5155172', '5155174', '5155175', '5155176', '5155177', '5155178', '5155182', '5155184', '5155185', '5155187', '5155188', '5155189', '5155191', '5155192']. Waiting...
2026-07-14 11:21:18 INFO None 5155171: status RUNNING/PENDING
2026-07-14 11:21:18 INFO None 5155172: status RUNNING/PENDING
2026-07-14 11:21:18 INFO None 5155174: status RUNNING/PENDING
2026-07-14 11:21:18 INFO None 5155175: status RUNNING/PENDING
2026-07-14 11:21:18 INFO None 5155176: status RUNNING/PENDING
2026-07-14 11:21:18 INFO None 5155177: status RUNNING/PENDING
2026-07-14 11:21:18 INFO None 5155178: status RUNNING/PENDING
2026-07-14 11:21:18 INFO None 5155182: status RUNNING/PENDING
2026-07-14 11:21:18 INFO None 5155184: status RUNNING/PENDING
2026-07-14 11:21:18 INFO None 5155185: status RUNNING/PENDING
2026-07-14 11:21:19 INFO None 5155187: status RUNNING/PENDING
2026-07-14 11:21:19 INFO None 5155188: status RUNNING/PENDING
2026-07-14 11:21:19 INFO None 5155189: status RUNNING/PENDING
2026-07-14 11:21:19 INFO None 5155191: status RUNNING/PENDING
2026-07-14 11:21:19 INFO None 5155192: status RUNNING/PENDING
2026-07-14 11:21:19 INFO Jobs still running: ['5155171', '5155172', '5155174', '5155175', '5155176', '5155177', '5155178', '5155182', '5155184', '5155185', '5155187', '5155188', '5155189', '5155191', '5155192']. Waiting...
2026-07-14 11:21:34 INFO None 5155171: status RUNNING/PENDING
2026-07-14 11:21:34 INFO None 5155172: status RUNNING/PENDING
2026-07-14 11:21:34 INFO None 5155174: status RUNNING/PENDING
2026-07-14 11:21:34 INFO None 5155175: status RUNNING/PENDING
2026-07-14 11:21:34 INFO None 5155176: status RUNNING/PENDING
2026-07-14 11:21:34 INFO None 5155177: status RUNNING/PENDING
2026-07-14 11:21:34 INFO None 5155178: status RUNNING/PENDING
2026-07-14 11:21:34 INFO None 5155182: status RUNNING/PENDING
2026-07-14 11:21:34 INFO None 5155184: status RUNNING/PENDING
2026-07-14 11:21:34 INFO None 5155185: status RUNNING/PENDING
2026-07-14 11:21:34 INFO None 5155187: status RUNNING/PENDING
2026-07-14 11:21:34 INFO None 5155188: status RUNNING/PENDING
2026-07-14 11:21:34 INFO None 5155189: status RUNNING/PENDING
2026-07-14 11:21:34 INFO None 5155191: status RUNNING/PENDING
2026-07-14 11:21:34 INFO None 5155192: status RUNNING/PENDING
2026-07-14 11:21:34 INFO Jobs still running: ['5155171', '5155172', '5155174', '5155175', '5155176', '5155177', '5155178', '5155182', '5155184', '5155185', '5155187', '5155188', '5155189', '5155191', '5155192']. Waiting...
2026-07-14 11:21:49 INFO None 5155171: status RUNNING/PENDING
2026-07-14 11:21:49 INFO None 5155172: status RUNNING/PENDING
2026-07-14 11:21:49 INFO None 5155174: status RUNNING/PENDING
2026-07-14 11:21:49 INFO None 5155175: status RUNNING/PENDING
2026-07-14 11:21:49 INFO None 5155176: status RUNNING/PENDING
2026-07-14 11:21:49 INFO None 5155177: status RUNNING/PENDING
2026-07-14 11:21:49 INFO None 5155178: status RUNNING/PENDING
2026-07-14 11:21:49 INFO None 5155182: status RUNNING/PENDING
2026-07-14 11:21:49 INFO None 5155184: status RUNNING/PENDING
2026-07-14 11:21:49 INFO None 5155185: status RUNNING/PENDING
2026-07-14 11:21:49 INFO None 5155187: status RUNNING/PENDING
2026-07-14 11:21:49 INFO None 5155188: status RUNNING/PENDING
2026-07-14 11:21:49 INFO None 5155189: status RUNNING/PENDING
2026-07-14 11:21:49 INFO None 5155191: status RUNNING/PENDING
2026-07-14 11:21:49 INFO None 5155192: status RUNNING/PENDING
2026-07-14 11:21:49 INFO Jobs still running: ['5155171', '5155172', '5155174', '5155175', '5155176', '5155177', '5155178', '5155182', '5155184', '5155185', '5155187', '5155188', '5155189', '5155191', '5155192']. Waiting...
2026-07-14 11:22:04 INFO None 5155171: status RUNNING/PENDING
2026-07-14 11:22:04 INFO None 5155172: status RUNNING/PENDING
2026-07-14 11:22:04 INFO None 5155174: status RUNNING/PENDING
2026-07-14 11:22:04 INFO None 5155175: status RUNNING/PENDING
2026-07-14 11:22:04 INFO None 5155176: status RUNNING/PENDING
2026-07-14 11:22:04 INFO None 5155177: status RUNNING/PENDING
2026-07-14 11:22:04 INFO None 5155178: status RUNNING/PENDING
2026-07-14 11:22:04 INFO None 5155182: status RUNNING/PENDING
2026-07-14 11:22:04 INFO None 5155184: status RUNNING/PENDING
2026-07-14 11:22:04 INFO None 5155185: status RUNNING/PENDING
2026-07-14 11:22:04 INFO None 5155187: status RUNNING/PENDING
2026-07-14 11:22:04 INFO None 5155188: status RUNNING/PENDING
2026-07-14 11:22:04 INFO None 5155189: status RUNNING/PENDING
2026-07-14 11:22:04 INFO None 5155191: status RUNNING/PENDING
2026-07-14 11:22:04 INFO None 5155192: status RUNNING/PENDING
2026-07-14 11:22:04 INFO Jobs still running: ['5155171', '5155172', '5155174', '5155175', '5155176', '5155177', '5155178', '5155182', '5155184', '5155185', '5155187', '5155188', '5155189', '5155191', '5155192']. Waiting...
2026-07-14 11:22:19 INFO None 5155171: status RUNNING/PENDING
2026-07-14 11:23:21 INFO None 5155172: status RUNNING/PENDING
2026-07-14 11:23:21 INFO None 5155174: status RUNNING/PENDING
2026-07-14 11:23:21 INFO None 5155175: status RUNNING/PENDING
2026-07-14 11:23:21 INFO None 5155176: status RUNNING/PENDING
2026-07-14 11:23:21 INFO None 5155177: status RUNNING/PENDING
2026-07-14 11:23:21 INFO None 5155178: status RUNNING/PENDING
2026-07-14 11:23:21 INFO None 5155182: status RUNNING/PENDING
2026-07-14 11:23:21 INFO None 5155184: status RUNNING/PENDING
2026-07-14 11:23:21 INFO None 5155185: status RUNNING/PENDING
2026-07-14 11:23:21 INFO None 5155187: status RUNNING/PENDING
2026-07-14 11:23:21 INFO None 5155188: status RUNNING/PENDING
2026-07-14 11:23:21 INFO None 5155189: status RUNNING/PENDING
2026-07-14 11:23:21 INFO None 5155191: status RUNNING/PENDING
2026-07-14 11:23:21 INFO None 5155192: status RUNNING/PENDING
2026-07-14 11:23:21 INFO Jobs still running: ['5155171', '5155172', '5155174', '5155175', '5155176', '5155177', '5155178', '5155182', '5155184', '5155185', '5155187', '5155188', '5155189', '5155191', '5155192']. Waiting...
2026-07-14 11:23:36 INFO None 5155171: status RUNNING/PENDING
2026-07-14 11:23:36 INFO None 5155172: status RUNNING/PENDING
2026-07-14 11:23:36 INFO None 5155174: status RUNNING/PENDING
2026-07-14 11:23:36 INFO None 5155175: status RUNNING/PENDING
2026-07-14 11:23:36 INFO None 5155176: status RUNNING/PENDING
2026-07-14 11:23:36 INFO None 5155177: status RUNNING/PENDING
2026-07-14 11:23:36 INFO None 5155178: status RUNNING/PENDING
2026-07-14 11:23:36 INFO None 5155182: status RUNNING/PENDING
2026-07-14 11:23:36 INFO None 5155184: status RUNNING/PENDING
2026-07-14 11:23:36 INFO None 5155185: status RUNNING/PENDING
2026-07-14 11:23:36 INFO None 5155187: status RUNNING/PENDING
2026-07-14 11:23:36 INFO None 5155188: status RUNNING/PENDING
2026-07-14 11:23:36 INFO None 5155189: status RUNNING/PENDING
2026-07-14 11:23:36 INFO None 5155191: status RUNNING/PENDING
2026-07-14 11:23:36 INFO None 5155192: status RUNNING/PENDING
2026-07-14 11:23:36 INFO Jobs still running: ['5155171', '5155172', '5155174', '5155175', '5155176', '5155177', '5155178', '5155182', '5155184', '5155185', '5155187', '5155188', '5155189', '5155191', '5155192']. Waiting...
2026-07-14 11:23:51 INFO None 5155171: status RUNNING/PENDING
2026-07-14 11:23:51 INFO None 5155172: status FINISHED
2026-07-14 11:23:51 INFO None 5155174: status FINISHED
2026-07-14 11:23:51 INFO None 5155175: status FINISHED
2026-07-14 11:23:51 INFO None 5155176: status RUNNING/PENDING
2026-07-14 11:23:51 INFO None 5155177: status RUNNING/PENDING
2026-07-14 11:23:51 INFO None 5155178: status RUNNING/PENDING
2026-07-14 11:23:51 INFO None 5155182: status RUNNING/PENDING
2026-07-14 11:23:51 INFO None 5155184: status RUNNING/PENDING
2026-07-14 11:23:51 INFO None 5155185: status RUNNING/PENDING
2026-07-14 11:23:51 INFO None 5155187: status RUNNING/PENDING
2026-07-14 11:23:51 INFO None 5155188: status FINISHED
2026-07-14 11:23:51 INFO None 5155189: status FINISHED
2026-07-14 11:23:52 INFO None 5155191: status FINISHED
2026-07-14 11:23:52 INFO None 5155192: status RUNNING/PENDING
2026-07-14 11:23:52 INFO Jobs still running: ['5155171', '5155176', '5155177', '5155178', '5155182', '5155184', '5155185', '5155187', '5155192']. Waiting...
2026-07-14 11:24:07 INFO None 5155171: status RUNNING/PENDING
2026-07-14 11:24:07 INFO None 5155172: status FINISHED
2026-07-14 11:24:07 INFO None 5155174: status FINISHED
2026-07-14 11:24:07 INFO None 5155175: status FINISHED
2026-07-14 11:24:07 INFO None 5155176: status RUNNING/PENDING
2026-07-14 11:24:07 INFO None 5155177: status RUNNING/PENDING
2026-07-14 11:24:07 INFO None 5155178: status FINISHED
2026-07-14 11:24:07 INFO None 5155182: status RUNNING/PENDING
2026-07-14 11:24:07 INFO None 5155184: status RUNNING/PENDING
2026-07-14 11:24:07 INFO None 5155185: status RUNNING/PENDING
2026-07-14 11:24:07 INFO None 5155187: status FINISHED
2026-07-14 11:24:07 INFO None 5155188: status FINISHED
2026-07-14 11:24:07 INFO None 5155189: status FINISHED
2026-07-14 11:24:07 INFO None 5155191: status FINISHED
2026-07-14 11:24:07 INFO None 5155192: status RUNNING/PENDING
2026-07-14 11:24:07 INFO Jobs still running: ['5155171', '5155176', '5155177', '5155182', '5155184', '5155185', '5155192']. Waiting...
2026-07-14 11:24:22 INFO None 5155171: status RUNNING/PENDING
2026-07-14 11:24:22 INFO None 5155172: status FINISHED
2026-07-14 11:24:22 INFO None 5155174: status FINISHED
2026-07-14 11:24:22 INFO None 5155175: status FINISHED
2026-07-14 11:24:22 INFO None 5155176: status FINISHED
2026-07-14 11:24:22 INFO None 5155177: status FINISHED
2026-07-14 11:24:22 INFO None 5155178: status FINISHED
2026-07-14 11:24:22 INFO None 5155182: status FINISHED
2026-07-14 11:24:22 INFO None 5155184: status FINISHED
2026-07-14 11:24:22 INFO None 5155185: status FINISHED
2026-07-14 11:24:22 INFO None 5155187: status FINISHED
2026-07-14 11:24:22 INFO None 5155188: status FINISHED
2026-07-14 11:24:22 INFO None 5155189: status FINISHED
2026-07-14 11:24:22 INFO None 5155191: status FINISHED
2026-07-14 11:24:22 INFO None 5155192: status RUNNING/PENDING
2026-07-14 11:24:22 INFO Jobs still running: ['5155171', '5155192']. Waiting...
2026-07-14 11:24:37 INFO None 5155171: status RUNNING/PENDING
2026-07-14 11:24:37 INFO None 5155172: status FINISHED
2026-07-14 11:24:37 INFO None 5155174: status FINISHED
2026-07-14 11:24:37 INFO None 5155175: status FINISHED
2026-07-14 11:24:37 INFO None 5155176: status FINISHED
2026-07-14 11:24:37 INFO None 5155177: status FINISHED
2026-07-14 11:24:37 INFO None 5155178: status FINISHED
2026-07-14 11:24:37 INFO None 5155182: status FINISHED
2026-07-14 11:24:37 INFO None 5155184: status FINISHED
2026-07-14 11:24:37 INFO None 5155185: status FINISHED
2026-07-14 11:24:37 INFO None 5155187: status FINISHED
2026-07-14 11:24:37 INFO None 5155188: status FINISHED
2026-07-14 11:24:37 INFO None 5155189: status FINISHED
2026-07-14 11:24:37 INFO None 5155191: status FINISHED
2026-07-14 11:24:37 INFO None 5155192: status FINISHED
2026-07-14 11:24:37 INFO Jobs still running: ['5155171']. Waiting...
2026-07-14 11:24:53 INFO None 5155171: status RUNNING/PENDING
2026-07-14 11:24:53 INFO None 5155172: status FINISHED
2026-07-14 11:24:53 INFO None 5155174: status FINISHED
2026-07-14 11:24:53 INFO None 5155175: status FINISHED
2026-07-14 11:24:53 INFO None 5155176: status FINISHED
2026-07-14 11:24:53 INFO None 5155177: status FINISHED
2026-07-14 11:24:53 INFO None 5155178: status FINISHED
2026-07-14 11:24:53 INFO None 5155182: status FINISHED
2026-07-14 11:24:53 INFO None 5155184: status FINISHED
2026-07-14 11:24:53 INFO None 5155185: status FINISHED
2026-07-14 11:24:53 INFO None 5155187: status FINISHED
2026-07-14 11:24:53 INFO None 5155188: status FINISHED
2026-07-14 11:24:53 INFO None 5155189: status FINISHED
2026-07-14 11:24:53 INFO None 5155191: status FINISHED
2026-07-14 11:24:53 INFO None 5155192: status FINISHED
2026-07-14 11:24:53 INFO Jobs still running: ['5155171']. Waiting...
2026-07-14 11:25:08 INFO None 5155171: status RUNNING/PENDING
2026-07-14 11:25:08 INFO None 5155172: status FINISHED
2026-07-14 11:25:08 INFO None 5155174: status FINISHED
2026-07-14 11:25:08 INFO None 5155175: status FINISHED
2026-07-14 11:25:08 INFO None 5155176: status FINISHED
2026-07-14 11:25:08 INFO None 5155177: status FINISHED
2026-07-14 11:25:08 INFO None 5155178: status FINISHED
2026-07-14 11:25:08 INFO None 5155182: status FINISHED
2026-07-14 11:25:08 INFO None 5155184: status FINISHED
2026-07-14 11:25:08 INFO None 5155185: status FINISHED
2026-07-14 11:25:08 INFO None 5155187: status FINISHED
2026-07-14 11:25:08 INFO None 5155188: status FINISHED
2026-07-14 11:25:08 INFO None 5155189: status FINISHED
2026-07-14 11:25:08 INFO None 5155191: status FINISHED
2026-07-14 11:25:08 INFO None 5155192: status FINISHED
2026-07-14 11:25:08 INFO Jobs still running: ['5155171']. Waiting...
2026-07-14 11:25:25 INFO None 5155171: status FINISHED
2026-07-14 11:25:25 INFO None 5155172: status FINISHED
2026-07-14 11:25:25 INFO None 5155174: status FINISHED
2026-07-14 11:25:25 INFO None 5155175: status FINISHED
2026-07-14 11:25:25 INFO None 5155176: status FINISHED
2026-07-14 11:25:26 INFO None 5155177: status FINISHED
2026-07-14 11:25:26 INFO None 5155178: status FINISHED
2026-07-14 11:25:26 INFO None 5155182: status FINISHED
2026-07-14 11:25:26 INFO None 5155184: status FINISHED
2026-07-14 11:25:26 INFO None 5155185: status FINISHED
2026-07-14 11:25:26 INFO None 5155187: status FINISHED
2026-07-14 11:25:26 INFO None 5155188: status FINISHED
2026-07-14 11:25:26 INFO None 5155189: status FINISHED
2026-07-14 11:25:26 INFO None 5155191: status FINISHED
2026-07-14 11:25:26 INFO None 5155192: status FINISHED
2026-07-14 11:25:26 INFO Jobs ['5155171', '5155172', '5155174', '5155175', '5155176', '5155177', '5155178', '5155182', '5155184', '5155185', '5155187', '5155188', '5155189', '5155191', '5155192'] have finished
2026-07-14 11:25:26 INFO Checking restart files were created ...
2026-07-14 11:25:26 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020600_1_ENS1.nc(668832435 bytes)
2026-07-14 11:25:26 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020600_1_ENS2.nc(668832435 bytes)
2026-07-14 11:25:26 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020600_1_ENS3.nc(668832435 bytes)
2026-07-14 11:25:26 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020600_1_ENS4.nc(668832435 bytes)
2026-07-14 11:25:26 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020600_1_ENS5.nc(668832435 bytes)
2026-07-14 11:25:26 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020600_1_ENS6.nc(668832435 bytes)
2026-07-14 11:25:26 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020600_1_ENS7.nc(668832435 bytes)
2026-07-14 11:25:26 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020600_1_ENS8.nc(668832435 bytes)
2026-07-14 11:25:26 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020600_1_ENS9.nc(668832435 bytes)
2026-07-14 11:25:26 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020600_1_ENS10.nc(668832435 bytes)
2026-07-14 11:25:26 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020600_1_ENS11.nc(668832435 bytes)
2026-07-14 11:25:26 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020600_1_ENS12.nc(668832435 bytes)
2026-07-14 11:25:26 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020600_1_ENS13.nc(668832435 bytes)
2026-07-14 11:25:26 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020600_1_ENS14.nc(668832435 bytes)
2026-07-14 11:25:26 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020600_1_ENS15.nc(668832435 bytes)
2026-07-14 11:25:26 INFO  Run_model() completed successfully.
2026-07-14 11:25:26 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 00:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 11:25:26 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 01:00:00 days=153072 seconds=3600
2026-07-14 11:25:26 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-14 11:25:26 INFO [TIME] increment current_time 2020-02-06 00:00:00 -> 2020-02-06 01:00:00
2026-07-14 11:25:26 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 11:25:26 INFO ---------->>> Running process_satellite_data()
2026-07-14 11:25:26 INFO [DART] No satellite data found, skipping assimilation
2026-07-14 11:25:26 INFO after_assimilation() skipped
2026-07-14 11:25:26 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-14 11:25:26 INFO [TIME] step_end current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 11:25:26 INFO [TIME] step_start current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 11:25:26 INFO [TIME] window start=2020-02-06 01:00:00 end=2020-02-06 09:00:00 run_hours=8 has_assimilation=True
2026-07-14 11:25:26 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 11:25:27 INFO Hourly dataset computed and listing created
2026-07-14 11:25:45 INFO Hourly dataset computed
2026-07-14 11:25:46 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 11:25:47 INFO Hourly dataset computed and listing created
2026-07-14 11:25:49 INFO Hourly dataset computed
2026-07-14 11:25:49 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 11:25:50 INFO Hourly dataset computed and listing created
2026-07-14 11:25:52 INFO Hourly dataset computed
2026-07-14 11:25:52 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 11:25:53 INFO Hourly dataset computed and listing created
2026-07-14 11:25:56 INFO Hourly dataset computed
2026-07-14 11:25:56 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 11:25:57 INFO Hourly dataset computed and listing created
2026-07-14 11:25:59 INFO Hourly dataset computed
2026-07-14 11:26:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 11:26:01 INFO Hourly dataset computed and listing created
2026-07-14 11:28:11 INFO Hourly dataset computed
2026-07-14 11:28:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 11:28:14 INFO Hourly dataset computed and listing created
2026-07-14 11:28:16 INFO Hourly dataset computed
2026-07-14 11:28:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 11:28:17 INFO Hourly dataset computed and listing created
2026-07-14 11:28:19 INFO Hourly dataset computed
2026-07-14 11:28:20 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 11:28:21 INFO Hourly dataset computed and listing created
2026-07-14 11:28:23 INFO Hourly dataset computed
2026-07-14 11:28:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 11:28:24 INFO Hourly dataset computed and listing created
2026-07-14 11:28:26 INFO Hourly dataset computed
2026-07-14 11:28:26 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 11:28:27 INFO Hourly dataset computed and listing created
2026-07-14 11:28:29 INFO Hourly dataset computed
2026-07-14 11:28:29 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 11:28:30 INFO Hourly dataset computed and listing created
2026-07-14 11:28:32 INFO Hourly dataset computed
2026-07-14 11:28:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 11:28:33 INFO Hourly dataset computed and listing created
2026-07-14 11:28:35 INFO Hourly dataset computed
2026-07-14 11:28:35 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 11:28:37 INFO Hourly dataset computed and listing created
2026-07-14 11:28:39 INFO Hourly dataset computed
2026-07-14 11:28:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 11:28:40 INFO Hourly dataset computed and listing created
2026-07-14 11:28:42 INFO Hourly dataset computed
2026-07-14 11:28:42 INFO ---------->>> Running CHIMERE model from 2020-02-06 01:00:00 to 2020-02-06 09:00:00
2026-07-14 11:28:42 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 11:28:42 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1
2026-07-14 11:28:42 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020600_1_ENS1.nc
2026-07-14 11:28:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-14 11:28:42 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 11:28:42 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-14 11:28:42 INFO Queuing job for member 1...
2026-07-14 11:28:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 11:28:42 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-14 11:28:43 INFO Found: ['5155226']
2026-07-14 11:28:48 INFO [TGCC-IRENE] Submitted job with ID:['5155226']
2026-07-14 11:28:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 11:28:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2
2026-07-14 11:28:48 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020600_1_ENS2.nc
2026-07-14 11:28:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-14 11:28:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 11:28:48 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-14 11:28:48 INFO Queuing job for member 2...
2026-07-14 11:28:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 11:28:48 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-14 11:28:48 INFO Found: ['5155227']
2026-07-14 11:28:53 INFO [TGCC-IRENE] Submitted job with ID:['5155227']
2026-07-14 11:28:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 11:28:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3
2026-07-14 11:28:53 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020600_1_ENS3.nc
2026-07-14 11:28:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-14 11:28:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 11:28:53 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-14 11:28:53 INFO Queuing job for member 3...
2026-07-14 11:28:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 11:28:53 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-14 11:28:54 INFO Found: ['5155228']
2026-07-14 11:28:59 INFO [TGCC-IRENE] Submitted job with ID:['5155228']
2026-07-14 11:28:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 11:28:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4
2026-07-14 11:28:59 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020600_1_ENS4.nc
2026-07-14 11:28:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-14 11:28:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 11:28:59 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-14 11:28:59 INFO Queuing job for member 4...
2026-07-14 11:28:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 11:28:59 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-14 11:29:00 INFO Found: ['5155231']
2026-07-14 11:29:05 INFO [TGCC-IRENE] Submitted job with ID:['5155231']
2026-07-14 11:29:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 11:29:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5
2026-07-14 11:29:05 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020600_1_ENS5.nc
2026-07-14 11:29:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-14 11:29:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 11:29:05 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-14 11:29:05 INFO Queuing job for member 5...
2026-07-14 11:29:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 11:29:05 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-14 11:29:06 INFO Found: ['5155239']
2026-07-14 11:29:11 INFO [TGCC-IRENE] Submitted job with ID:['5155239']
2026-07-14 11:29:11 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 11:29:11 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6
2026-07-14 11:29:11 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020600_1_ENS6.nc
2026-07-14 11:29:11 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-14 11:29:11 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 11:29:11 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-14 11:29:11 INFO Queuing job for member 6...
2026-07-14 11:29:11 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 11:29:11 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-14 11:29:12 INFO Found: ['5155240']
2026-07-14 11:29:17 INFO [TGCC-IRENE] Submitted job with ID:['5155240']
2026-07-14 11:29:17 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 11:29:17 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7
2026-07-14 11:29:17 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020600_1_ENS7.nc
2026-07-14 11:29:17 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-14 11:29:17 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 11:29:17 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-14 11:29:17 INFO Queuing job for member 7...
2026-07-14 11:29:17 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 11:29:17 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-14 11:29:17 INFO Found: ['5155241']
2026-07-14 11:29:22 INFO [TGCC-IRENE] Submitted job with ID:['5155241']
2026-07-14 11:29:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 11:29:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8
2026-07-14 11:29:22 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020600_1_ENS8.nc
2026-07-14 11:29:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-14 11:29:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 11:29:22 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-14 11:29:22 INFO Queuing job for member 8...
2026-07-14 11:29:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 11:29:22 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-14 11:29:23 INFO Found: ['5155242']
2026-07-14 11:29:28 INFO [TGCC-IRENE] Submitted job with ID:['5155242']
2026-07-14 11:29:28 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 11:29:28 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9
2026-07-14 11:29:28 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020600_1_ENS9.nc
2026-07-14 11:29:28 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-14 11:29:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 11:29:28 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-14 11:29:28 INFO Queuing job for member 9...
2026-07-14 11:29:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 11:29:28 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-14 11:29:29 INFO Found: ['5155243']
2026-07-14 11:29:34 INFO [TGCC-IRENE] Submitted job with ID:['5155243']
2026-07-14 11:29:34 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 11:29:34 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10
2026-07-14 11:29:34 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020600_1_ENS10.nc
2026-07-14 11:29:34 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-14 11:29:34 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 11:29:34 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-14 11:29:34 INFO Queuing job for member 10...
2026-07-14 11:29:34 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 11:29:34 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-14 11:29:35 INFO Found: ['5155244']
2026-07-14 11:29:40 INFO [TGCC-IRENE] Submitted job with ID:['5155244']
2026-07-14 11:29:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 11:29:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11
2026-07-14 11:29:40 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020600_1_ENS11.nc
2026-07-14 11:29:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-14 11:29:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 11:29:40 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-14 11:29:40 INFO Queuing job for member 11...
2026-07-14 11:29:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 11:29:40 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-14 11:29:40 INFO Found: ['5155245']
2026-07-14 11:29:45 INFO [TGCC-IRENE] Submitted job with ID:['5155245']
2026-07-14 11:29:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 11:29:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12
2026-07-14 11:29:45 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020600_1_ENS12.nc
2026-07-14 11:29:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-14 11:29:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 11:29:45 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-14 11:29:45 INFO Queuing job for member 12...
2026-07-14 11:29:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 11:29:45 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-14 11:29:46 INFO Found: ['5155246']
2026-07-14 11:29:51 INFO [TGCC-IRENE] Submitted job with ID:['5155246']
2026-07-14 11:29:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 11:29:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13
2026-07-14 11:29:51 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020600_1_ENS13.nc
2026-07-14 11:29:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-14 11:29:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 11:29:51 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-14 11:29:51 INFO Queuing job for member 13...
2026-07-14 11:29:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 11:29:51 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-14 11:29:52 INFO Found: ['5155247']
2026-07-14 11:29:57 INFO [TGCC-IRENE] Submitted job with ID:['5155247']
2026-07-14 11:29:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 11:29:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14
2026-07-14 11:29:57 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020600_1_ENS14.nc
2026-07-14 11:29:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-14 11:29:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 11:29:57 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-14 11:29:57 INFO Queuing job for member 14...
2026-07-14 11:29:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 11:29:57 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-14 11:29:58 INFO Found: ['5155249']
2026-07-14 11:30:03 INFO [TGCC-IRENE] Submitted job with ID:['5155249']
2026-07-14 11:30:03 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 11:30:03 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15
2026-07-14 11:30:03 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020600_1_ENS15.nc
2026-07-14 11:30:03 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-14 11:30:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 11:30:03 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-14 11:30:03 INFO Queuing job for member 15...
2026-07-14 11:30:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 11:30:03 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-14 11:30:04 INFO Found: ['5155251']
2026-07-14 11:30:09 INFO [TGCC-IRENE] Submitted job with ID:['5155251']
2026-07-14 11:30:09 INFO Checking job status ...
2026-07-14 11:30:09 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:30:09 INFO None 5155227: status RUNNING/PENDING
2026-07-14 11:30:09 INFO None 5155228: status RUNNING/PENDING
2026-07-14 11:30:09 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:30:09 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:30:09 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:30:09 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:30:09 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:30:09 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:30:09 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:30:09 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:30:09 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:30:09 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:30:09 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:30:09 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:30:09 INFO Jobs still running: ['5155226', '5155227', '5155228', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:30:24 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:30:24 INFO None 5155227: status RUNNING/PENDING
2026-07-14 11:30:24 INFO None 5155228: status RUNNING/PENDING
2026-07-14 11:30:24 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:30:24 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:30:24 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:30:24 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:30:24 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:30:24 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:30:24 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:30:24 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:30:24 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:30:24 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:30:24 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:30:24 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:30:24 INFO Jobs still running: ['5155226', '5155227', '5155228', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:30:39 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:30:39 INFO None 5155227: status RUNNING/PENDING
2026-07-14 11:30:39 INFO None 5155228: status RUNNING/PENDING
2026-07-14 11:30:39 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:30:39 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:30:39 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:30:39 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:30:39 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:30:39 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:30:39 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:30:39 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:30:39 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:30:39 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:30:40 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:30:40 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:30:40 INFO Jobs still running: ['5155226', '5155227', '5155228', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:30:55 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:30:55 INFO None 5155227: status RUNNING/PENDING
2026-07-14 11:30:55 INFO None 5155228: status RUNNING/PENDING
2026-07-14 11:30:55 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:30:55 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:30:55 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:30:55 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:30:55 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:30:55 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:30:55 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:30:55 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:30:55 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:30:55 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:30:55 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:30:55 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:30:55 INFO Jobs still running: ['5155226', '5155227', '5155228', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:31:10 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:31:10 INFO None 5155227: status RUNNING/PENDING
2026-07-14 11:31:10 INFO None 5155228: status RUNNING/PENDING
2026-07-14 11:31:10 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:31:10 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:31:10 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:31:10 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:31:10 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:31:10 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:31:10 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:31:10 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:31:10 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:31:10 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:31:10 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:31:10 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:31:10 INFO Jobs still running: ['5155226', '5155227', '5155228', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:31:25 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:33:16 INFO None 5155227: status RUNNING/PENDING
2026-07-14 11:33:16 INFO None 5155228: status RUNNING/PENDING
2026-07-14 11:33:16 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:33:16 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:33:16 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:33:16 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:33:16 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:33:16 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:33:16 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:33:16 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:33:16 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:33:16 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:33:16 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:33:16 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:33:16 INFO Jobs still running: ['5155226', '5155227', '5155228', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:33:31 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:33:31 INFO None 5155227: status RUNNING/PENDING
2026-07-14 11:33:31 INFO None 5155228: status RUNNING/PENDING
2026-07-14 11:33:31 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:33:31 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:33:31 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:33:32 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:33:32 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:33:32 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:33:32 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:33:32 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:33:32 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:33:32 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:33:32 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:33:32 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:33:32 INFO Jobs still running: ['5155226', '5155227', '5155228', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:33:47 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:33:47 INFO None 5155227: status RUNNING/PENDING
2026-07-14 11:33:47 INFO None 5155228: status RUNNING/PENDING
2026-07-14 11:33:47 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:33:47 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:33:47 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:33:47 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:33:47 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:33:47 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:33:47 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:33:47 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:33:47 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:33:47 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:33:47 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:33:47 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:33:47 INFO Jobs still running: ['5155226', '5155227', '5155228', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:34:02 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:34:02 INFO None 5155227: status RUNNING/PENDING
2026-07-14 11:34:02 INFO None 5155228: status RUNNING/PENDING
2026-07-14 11:34:02 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:34:02 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:34:02 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:34:02 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:34:02 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:34:02 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:34:02 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:34:02 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:34:02 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:34:02 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:34:02 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:34:02 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:34:02 INFO Jobs still running: ['5155226', '5155227', '5155228', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:34:17 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:34:17 INFO None 5155227: status RUNNING/PENDING
2026-07-14 11:34:17 INFO None 5155228: status RUNNING/PENDING
2026-07-14 11:34:17 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:34:17 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:34:17 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:34:17 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:34:17 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:34:18 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:34:18 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:34:18 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:34:18 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:34:18 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:34:18 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:34:18 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:34:18 INFO Jobs still running: ['5155226', '5155227', '5155228', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:34:33 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:34:33 INFO None 5155227: status RUNNING/PENDING
2026-07-14 11:34:33 INFO None 5155228: status RUNNING/PENDING
2026-07-14 11:34:33 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:34:33 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:34:33 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:34:33 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:34:33 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:34:33 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:34:33 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:34:33 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:34:33 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:34:33 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:34:33 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:34:33 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:34:33 INFO Jobs still running: ['5155226', '5155227', '5155228', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:34:48 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:34:48 INFO None 5155227: status RUNNING/PENDING
2026-07-14 11:34:48 INFO None 5155228: status RUNNING/PENDING
2026-07-14 11:34:48 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:34:48 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:34:48 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:34:48 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:34:48 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:34:48 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:34:48 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:34:48 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:34:48 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:34:48 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:34:48 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:34:48 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:34:48 INFO Jobs still running: ['5155226', '5155227', '5155228', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:35:03 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:35:03 INFO None 5155227: status RUNNING/PENDING
2026-07-14 11:35:03 INFO None 5155228: status RUNNING/PENDING
2026-07-14 11:35:03 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:35:03 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:35:03 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:35:03 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:35:03 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:35:03 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:35:03 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:35:04 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:35:04 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:35:04 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:35:04 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:35:04 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:35:04 INFO Jobs still running: ['5155226', '5155227', '5155228', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:35:19 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:35:19 INFO None 5155227: status RUNNING/PENDING
2026-07-14 11:35:19 INFO None 5155228: status RUNNING/PENDING
2026-07-14 11:35:19 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:35:19 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:35:19 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:35:19 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:35:19 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:35:19 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:35:19 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:35:19 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:35:19 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:35:19 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:35:19 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:35:19 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:35:19 INFO Jobs still running: ['5155226', '5155227', '5155228', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:35:34 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:35:34 INFO None 5155227: status RUNNING/PENDING
2026-07-14 11:35:34 INFO None 5155228: status RUNNING/PENDING
2026-07-14 11:35:34 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:35:34 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:35:34 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:35:34 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:35:34 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:35:34 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:35:34 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:35:34 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:35:34 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:35:34 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:35:34 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:35:34 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:35:34 INFO Jobs still running: ['5155226', '5155227', '5155228', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:35:49 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:35:49 INFO None 5155227: status RUNNING/PENDING
2026-07-14 11:35:49 INFO None 5155228: status RUNNING/PENDING
2026-07-14 11:35:49 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:35:49 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:35:49 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:35:49 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:35:49 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:35:49 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:35:49 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:35:49 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:35:49 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:35:49 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:35:49 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:35:50 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:35:50 INFO Jobs still running: ['5155226', '5155227', '5155228', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:36:05 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:36:05 INFO None 5155227: status RUNNING/PENDING
2026-07-14 11:36:05 INFO None 5155228: status RUNNING/PENDING
2026-07-14 11:36:05 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:36:05 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:36:05 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:36:05 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:36:05 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:36:05 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:36:05 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:36:05 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:36:05 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:36:05 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:36:05 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:36:05 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:36:05 INFO Jobs still running: ['5155226', '5155227', '5155228', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:36:20 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:36:20 INFO None 5155227: status RUNNING/PENDING
2026-07-14 11:36:20 INFO None 5155228: status RUNNING/PENDING
2026-07-14 11:36:20 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:36:20 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:36:20 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:36:20 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:36:20 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:36:20 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:36:20 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:36:20 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:36:20 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:36:20 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:36:20 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:36:20 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:36:20 INFO Jobs still running: ['5155226', '5155227', '5155228', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:36:35 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:36:35 INFO None 5155227: status RUNNING/PENDING
2026-07-14 11:36:35 INFO None 5155228: status RUNNING/PENDING
2026-07-14 11:36:35 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:36:35 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:36:35 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:36:35 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:36:35 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:36:35 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:36:35 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:36:35 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:36:35 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:36:35 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:36:35 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:36:35 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:36:35 INFO Jobs still running: ['5155226', '5155227', '5155228', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:36:50 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:36:50 INFO None 5155227: status RUNNING/PENDING
2026-07-14 11:36:50 INFO None 5155228: status RUNNING/PENDING
2026-07-14 11:36:51 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:36:51 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:36:51 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:36:51 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:36:51 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:36:51 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:36:51 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:36:51 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:36:51 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:36:51 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:36:51 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:36:51 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:36:51 INFO Jobs still running: ['5155226', '5155227', '5155228', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:37:06 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:37:06 INFO None 5155227: status RUNNING/PENDING
2026-07-14 11:37:06 INFO None 5155228: status RUNNING/PENDING
2026-07-14 11:37:06 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:37:06 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:37:06 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:37:06 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:37:06 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:37:06 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:37:06 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:37:06 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:37:06 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:37:06 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:37:06 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:37:06 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:37:06 INFO Jobs still running: ['5155226', '5155227', '5155228', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:37:21 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:37:21 INFO None 5155227: status RUNNING/PENDING
2026-07-14 11:37:21 INFO None 5155228: status RUNNING/PENDING
2026-07-14 11:37:21 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:37:21 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:37:21 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:37:21 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:37:21 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:37:21 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:37:21 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:37:21 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:37:21 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:37:21 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:37:21 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:37:21 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:37:21 INFO Jobs still running: ['5155226', '5155227', '5155228', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:37:36 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:38:17 INFO None 5155227: status RUNNING/PENDING
2026-07-14 11:39:42 INFO None 5155228: status RUNNING/PENDING
2026-07-14 11:39:42 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:39:42 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:39:42 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:39:42 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:39:42 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:39:42 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:39:42 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:39:42 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:39:42 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:39:42 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:39:42 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:39:42 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:39:42 INFO Jobs still running: ['5155226', '5155227', '5155228', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:39:57 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:39:57 INFO None 5155227: status RUNNING/PENDING
2026-07-14 11:39:57 INFO None 5155228: status RUNNING/PENDING
2026-07-14 11:39:57 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:39:58 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:39:58 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:39:58 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:39:58 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:39:58 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:39:58 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:39:58 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:39:58 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:39:58 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:39:58 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:39:58 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:39:58 INFO Jobs still running: ['5155226', '5155227', '5155228', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:40:13 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:40:13 INFO None 5155227: status RUNNING/PENDING
2026-07-14 11:40:13 INFO None 5155228: status RUNNING/PENDING
2026-07-14 11:40:13 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:40:13 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:40:13 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:40:13 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:40:13 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:40:13 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:40:13 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:40:13 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:40:13 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:40:13 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:40:13 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:40:13 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:40:13 INFO Jobs still running: ['5155226', '5155227', '5155228', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:40:28 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:40:28 INFO None 5155227: status FINISHED
2026-07-14 11:40:28 INFO None 5155228: status FINISHED
2026-07-14 11:40:28 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:40:28 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:40:28 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:40:28 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:40:28 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:40:28 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:40:28 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:40:28 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:40:28 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:40:28 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:40:28 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:40:28 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:40:28 INFO Jobs still running: ['5155226', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:40:45 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:40:45 INFO None 5155227: status FINISHED
2026-07-14 11:40:45 INFO None 5155228: status FINISHED
2026-07-14 11:40:45 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:40:45 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:40:45 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:40:45 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:40:45 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:40:45 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:40:45 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:40:45 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:40:45 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:40:45 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:40:45 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:40:46 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:40:46 INFO Jobs still running: ['5155226', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:41:01 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:41:01 INFO None 5155227: status FINISHED
2026-07-14 11:41:01 INFO None 5155228: status FINISHED
2026-07-14 11:41:01 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:41:01 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:41:01 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:41:01 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:41:01 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:41:01 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:41:01 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:41:01 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:41:01 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:41:01 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:41:01 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:41:01 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:41:01 INFO Jobs still running: ['5155226', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:41:16 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:43:14 INFO None 5155227: status FINISHED
2026-07-14 11:43:14 INFO None 5155228: status FINISHED
2026-07-14 11:43:14 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:43:14 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:43:14 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:43:14 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:43:14 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:43:14 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:43:14 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:43:14 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:43:14 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:43:14 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:43:14 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:43:14 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:43:14 INFO Jobs still running: ['5155226', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:43:29 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:43:29 INFO None 5155227: status FINISHED
2026-07-14 11:43:29 INFO None 5155228: status FINISHED
2026-07-14 11:43:29 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:43:29 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:43:29 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:43:29 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:43:29 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:43:29 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:43:29 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:43:29 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:43:29 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:43:29 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:43:29 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:43:29 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:43:29 INFO Jobs still running: ['5155226', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:43:44 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:43:44 INFO None 5155227: status FINISHED
2026-07-14 11:43:44 INFO None 5155228: status FINISHED
2026-07-14 11:43:44 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:43:44 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:43:44 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:43:44 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:43:44 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:43:44 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:43:44 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:43:44 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:43:44 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:43:44 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:43:44 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:43:45 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:43:45 INFO Jobs still running: ['5155226', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:44:00 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:44:00 INFO None 5155227: status FINISHED
2026-07-14 11:44:00 INFO None 5155228: status FINISHED
2026-07-14 11:44:00 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:44:00 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:44:00 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:44:00 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:44:00 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:44:00 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:44:00 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:44:00 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:44:00 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:44:00 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:44:00 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:44:00 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:44:00 INFO Jobs still running: ['5155226', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:44:15 INFO None 5155226: status RUNNING/PENDING
2026-07-14 11:44:15 INFO None 5155227: status FINISHED
2026-07-14 11:44:15 INFO None 5155228: status FINISHED
2026-07-14 11:44:15 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:44:15 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:44:15 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:44:15 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:44:15 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:44:15 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:44:15 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:44:15 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:44:15 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:44:15 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:44:15 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:44:15 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:44:15 INFO Jobs still running: ['5155226', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:44:30 INFO None 5155226: status FINISHED
2026-07-14 11:44:30 INFO None 5155227: status FINISHED
2026-07-14 11:44:30 INFO None 5155228: status FINISHED
2026-07-14 11:44:30 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:44:30 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:44:30 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:44:30 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:44:30 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:44:30 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:44:30 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:44:30 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:44:30 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:44:30 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:44:31 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:44:31 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:44:31 INFO Jobs still running: ['5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:44:46 INFO None 5155226: status FINISHED
2026-07-14 11:44:46 INFO None 5155227: status FINISHED
2026-07-14 11:44:46 INFO None 5155228: status FINISHED
2026-07-14 11:44:46 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:44:46 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:44:46 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:44:46 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:44:46 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:44:46 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:44:46 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:44:46 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:44:46 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:44:46 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:44:46 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:44:46 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:44:46 INFO Jobs still running: ['5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:45:01 INFO None 5155226: status FINISHED
2026-07-14 11:45:01 INFO None 5155227: status FINISHED
2026-07-14 11:45:01 INFO None 5155228: status FINISHED
2026-07-14 11:45:01 INFO None 5155231: status RUNNING/PENDING
2026-07-14 11:45:01 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:45:01 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:45:01 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:45:01 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:45:01 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:45:01 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:45:01 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:45:01 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:45:01 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:45:01 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:45:01 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:45:01 INFO Jobs still running: ['5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:45:16 INFO None 5155226: status FINISHED
2026-07-14 11:45:16 INFO None 5155227: status FINISHED
2026-07-14 11:45:16 INFO None 5155228: status FINISHED
2026-07-14 11:45:16 INFO None 5155231: status FINISHED
2026-07-14 11:45:16 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:45:16 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:45:16 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:45:16 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:45:16 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:45:16 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:45:16 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:45:16 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:45:16 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:45:17 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:45:17 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:45:17 INFO Jobs still running: ['5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:45:32 INFO None 5155226: status FINISHED
2026-07-14 11:45:32 INFO None 5155227: status FINISHED
2026-07-14 11:45:32 INFO None 5155228: status FINISHED
2026-07-14 11:45:32 INFO None 5155231: status FINISHED
2026-07-14 11:45:32 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:45:32 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:45:32 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:45:32 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:45:32 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:45:32 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:45:32 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:45:32 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:45:32 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:45:32 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:45:32 INFO None 5155251: status RUNNING/PENDING
2026-07-14 11:45:32 INFO Jobs still running: ['5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251']. Waiting...
2026-07-14 11:45:47 INFO None 5155226: status FINISHED
2026-07-14 11:45:47 INFO None 5155227: status FINISHED
2026-07-14 11:45:47 INFO None 5155228: status FINISHED
2026-07-14 11:45:47 INFO None 5155231: status FINISHED
2026-07-14 11:45:47 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:45:47 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:45:47 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:45:47 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:45:47 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:45:47 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:45:47 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:45:47 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:45:47 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:45:47 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:45:47 INFO None 5155251: status FINISHED
2026-07-14 11:45:47 INFO Jobs still running: ['5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249']. Waiting...
2026-07-14 11:46:02 INFO None 5155226: status FINISHED
2026-07-14 11:46:02 INFO None 5155227: status FINISHED
2026-07-14 11:46:02 INFO None 5155228: status FINISHED
2026-07-14 11:46:02 INFO None 5155231: status FINISHED
2026-07-14 11:46:02 INFO None 5155239: status RUNNING/PENDING
2026-07-14 11:46:02 INFO None 5155240: status RUNNING/PENDING
2026-07-14 11:46:02 INFO None 5155241: status RUNNING/PENDING
2026-07-14 11:46:02 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:46:02 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:46:02 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:46:02 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:46:02 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:46:02 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:46:02 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:46:03 INFO None 5155251: status FINISHED
2026-07-14 11:46:03 INFO Jobs still running: ['5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249']. Waiting...
2026-07-14 11:46:18 INFO None 5155226: status FINISHED
2026-07-14 11:48:16 INFO None 5155227: status FINISHED
2026-07-14 11:48:16 INFO None 5155228: status FINISHED
2026-07-14 11:48:16 INFO None 5155231: status FINISHED
2026-07-14 11:48:16 INFO None 5155239: status FINISHED
2026-07-14 11:48:16 INFO None 5155240: status FINISHED
2026-07-14 11:48:16 INFO None 5155241: status FINISHED
2026-07-14 11:48:16 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:48:16 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:48:16 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:48:16 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:48:16 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:48:16 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:48:16 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:48:16 INFO None 5155251: status FINISHED
2026-07-14 11:48:16 INFO Jobs still running: ['5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249']. Waiting...
2026-07-14 11:48:31 INFO None 5155226: status FINISHED
2026-07-14 11:48:31 INFO None 5155227: status FINISHED
2026-07-14 11:48:31 INFO None 5155228: status FINISHED
2026-07-14 11:48:31 INFO None 5155231: status FINISHED
2026-07-14 11:48:31 INFO None 5155239: status FINISHED
2026-07-14 11:48:31 INFO None 5155240: status FINISHED
2026-07-14 11:48:31 INFO None 5155241: status FINISHED
2026-07-14 11:48:31 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:48:31 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:48:31 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:48:31 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:48:31 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:48:31 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:48:31 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:48:31 INFO None 5155251: status FINISHED
2026-07-14 11:48:31 INFO Jobs still running: ['5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249']. Waiting...
2026-07-14 11:48:46 INFO None 5155226: status FINISHED
2026-07-14 11:48:46 INFO None 5155227: status FINISHED
2026-07-14 11:48:46 INFO None 5155228: status FINISHED
2026-07-14 11:48:46 INFO None 5155231: status FINISHED
2026-07-14 11:48:46 INFO None 5155239: status FINISHED
2026-07-14 11:48:46 INFO None 5155240: status FINISHED
2026-07-14 11:48:46 INFO None 5155241: status FINISHED
2026-07-14 11:48:46 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:48:46 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:48:46 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:48:46 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:48:47 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:48:47 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:48:47 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:48:47 INFO None 5155251: status FINISHED
2026-07-14 11:48:47 INFO Jobs still running: ['5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249']. Waiting...
2026-07-14 11:49:02 INFO None 5155226: status FINISHED
2026-07-14 11:49:02 INFO None 5155227: status FINISHED
2026-07-14 11:49:02 INFO None 5155228: status FINISHED
2026-07-14 11:49:02 INFO None 5155231: status FINISHED
2026-07-14 11:49:02 INFO None 5155239: status FINISHED
2026-07-14 11:49:02 INFO None 5155240: status FINISHED
2026-07-14 11:49:02 INFO None 5155241: status FINISHED
2026-07-14 11:49:02 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:49:02 INFO None 5155243: status RUNNING/PENDING
2026-07-14 11:49:02 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:49:02 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:49:02 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:49:02 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:49:02 INFO None 5155249: status RUNNING/PENDING
2026-07-14 11:49:02 INFO None 5155251: status FINISHED
2026-07-14 11:49:02 INFO Jobs still running: ['5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249']. Waiting...
2026-07-14 11:49:17 INFO None 5155226: status FINISHED
2026-07-14 11:49:17 INFO None 5155227: status FINISHED
2026-07-14 11:49:17 INFO None 5155228: status FINISHED
2026-07-14 11:49:17 INFO None 5155231: status FINISHED
2026-07-14 11:49:17 INFO None 5155239: status FINISHED
2026-07-14 11:49:17 INFO None 5155240: status FINISHED
2026-07-14 11:49:17 INFO None 5155241: status FINISHED
2026-07-14 11:49:17 INFO None 5155242: status RUNNING/PENDING
2026-07-14 11:49:17 INFO None 5155243: status FINISHED
2026-07-14 11:49:17 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:49:17 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:49:17 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:49:17 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:49:17 INFO None 5155249: status FINISHED
2026-07-14 11:49:17 INFO None 5155251: status FINISHED
2026-07-14 11:49:17 INFO Jobs still running: ['5155242', '5155244', '5155245', '5155246', '5155247']. Waiting...
2026-07-14 11:49:32 INFO None 5155226: status FINISHED
2026-07-14 11:49:32 INFO None 5155227: status FINISHED
2026-07-14 11:49:32 INFO None 5155228: status FINISHED
2026-07-14 11:49:32 INFO None 5155231: status FINISHED
2026-07-14 11:49:32 INFO None 5155239: status FINISHED
2026-07-14 11:49:32 INFO None 5155240: status FINISHED
2026-07-14 11:49:32 INFO None 5155241: status FINISHED
2026-07-14 11:49:32 INFO None 5155242: status FINISHED
2026-07-14 11:49:32 INFO None 5155243: status FINISHED
2026-07-14 11:49:33 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:49:33 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:49:33 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:49:33 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:49:33 INFO None 5155249: status FINISHED
2026-07-14 11:49:33 INFO None 5155251: status FINISHED
2026-07-14 11:49:33 INFO Jobs still running: ['5155244', '5155245', '5155246', '5155247']. Waiting...
2026-07-14 11:49:48 INFO None 5155226: status FINISHED
2026-07-14 11:49:48 INFO None 5155227: status FINISHED
2026-07-14 11:49:48 INFO None 5155228: status FINISHED
2026-07-14 11:49:48 INFO None 5155231: status FINISHED
2026-07-14 11:49:48 INFO None 5155239: status FINISHED
2026-07-14 11:49:48 INFO None 5155240: status FINISHED
2026-07-14 11:49:48 INFO None 5155241: status FINISHED
2026-07-14 11:49:48 INFO None 5155242: status FINISHED
2026-07-14 11:49:48 INFO None 5155243: status FINISHED
2026-07-14 11:49:48 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:49:48 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:49:48 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:49:48 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:49:48 INFO None 5155249: status FINISHED
2026-07-14 11:49:48 INFO None 5155251: status FINISHED
2026-07-14 11:49:48 INFO Jobs still running: ['5155244', '5155245', '5155246', '5155247']. Waiting...
2026-07-14 11:50:03 INFO None 5155226: status FINISHED
2026-07-14 11:50:03 INFO None 5155227: status FINISHED (not in squeue)
2026-07-14 11:50:03 INFO None 5155228: status FINISHED (not in squeue)
2026-07-14 11:50:03 INFO None 5155231: status FINISHED
2026-07-14 11:50:03 INFO None 5155239: status FINISHED
2026-07-14 11:50:03 INFO None 5155240: status FINISHED
2026-07-14 11:50:03 INFO None 5155241: status FINISHED
2026-07-14 11:50:03 INFO None 5155242: status FINISHED
2026-07-14 11:50:03 INFO None 5155243: status FINISHED
2026-07-14 11:50:03 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:50:03 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:50:03 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:50:03 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:50:03 INFO None 5155249: status FINISHED
2026-07-14 11:50:03 INFO None 5155251: status FINISHED
2026-07-14 11:50:03 INFO Jobs still running: ['5155244', '5155245', '5155246', '5155247']. Waiting...
2026-07-14 11:50:18 INFO None 5155226: status FINISHED
2026-07-14 11:50:18 INFO None 5155227: status FINISHED (not in squeue)
2026-07-14 11:50:18 INFO None 5155228: status FINISHED (not in squeue)
2026-07-14 11:50:18 INFO None 5155231: status FINISHED
2026-07-14 11:50:18 INFO None 5155239: status FINISHED
2026-07-14 11:50:18 INFO None 5155240: status FINISHED
2026-07-14 11:50:18 INFO None 5155241: status FINISHED
2026-07-14 11:50:18 INFO None 5155242: status FINISHED
2026-07-14 11:50:18 INFO None 5155243: status FINISHED
2026-07-14 11:50:19 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:50:19 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:50:19 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:50:19 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:50:19 INFO None 5155249: status FINISHED
2026-07-14 11:50:19 INFO None 5155251: status FINISHED
2026-07-14 11:50:19 INFO Jobs still running: ['5155244', '5155245', '5155246', '5155247']. Waiting...
2026-07-14 11:50:34 INFO None 5155226: status FINISHED
2026-07-14 11:50:34 INFO None 5155227: status FINISHED (not in squeue)
2026-07-14 11:50:34 INFO None 5155228: status FINISHED (not in squeue)
2026-07-14 11:50:34 INFO None 5155231: status FINISHED
2026-07-14 11:50:34 INFO None 5155239: status FINISHED
2026-07-14 11:50:34 INFO None 5155240: status FINISHED
2026-07-14 11:50:34 INFO None 5155241: status FINISHED
2026-07-14 11:50:34 INFO None 5155242: status FINISHED
2026-07-14 11:50:34 INFO None 5155243: status FINISHED
2026-07-14 11:50:34 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:50:34 INFO None 5155245: status RUNNING/PENDING
2026-07-14 11:50:34 INFO None 5155246: status RUNNING/PENDING
2026-07-14 11:50:34 INFO None 5155247: status RUNNING/PENDING
2026-07-14 11:50:34 INFO None 5155249: status FINISHED
2026-07-14 11:50:34 INFO None 5155251: status FINISHED
2026-07-14 11:50:34 INFO Jobs still running: ['5155244', '5155245', '5155246', '5155247']. Waiting...
2026-07-14 11:50:49 INFO None 5155226: status FINISHED
2026-07-14 11:50:49 INFO None 5155227: status FINISHED (not in squeue)
2026-07-14 11:50:49 INFO None 5155228: status FINISHED (not in squeue)
2026-07-14 11:50:49 INFO None 5155231: status FINISHED
2026-07-14 11:50:49 INFO None 5155239: status FINISHED
2026-07-14 11:50:49 INFO None 5155240: status FINISHED
2026-07-14 11:50:49 INFO None 5155241: status FINISHED
2026-07-14 11:50:49 INFO None 5155242: status FINISHED
2026-07-14 11:50:49 INFO None 5155243: status FINISHED
2026-07-14 11:50:49 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:50:49 INFO None 5155245: status FINISHED
2026-07-14 11:50:49 INFO None 5155246: status FINISHED
2026-07-14 11:50:49 INFO None 5155247: status FINISHED
2026-07-14 11:50:49 INFO None 5155249: status FINISHED
2026-07-14 11:50:49 INFO None 5155251: status FINISHED
2026-07-14 11:50:49 INFO Jobs still running: ['5155244']. Waiting...
2026-07-14 11:51:04 INFO None 5155226: status FINISHED
2026-07-14 11:51:04 INFO None 5155227: status FINISHED (not in squeue)
2026-07-14 11:51:04 INFO None 5155228: status FINISHED (not in squeue)
2026-07-14 11:51:04 INFO None 5155231: status FINISHED
2026-07-14 11:51:04 INFO None 5155239: status FINISHED
2026-07-14 11:51:04 INFO None 5155240: status FINISHED
2026-07-14 11:51:04 INFO None 5155241: status FINISHED
2026-07-14 11:51:04 INFO None 5155242: status FINISHED
2026-07-14 11:51:04 INFO None 5155243: status FINISHED
2026-07-14 11:51:04 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:51:04 INFO None 5155245: status FINISHED
2026-07-14 11:51:05 INFO None 5155246: status FINISHED
2026-07-14 11:51:05 INFO None 5155247: status FINISHED
2026-07-14 11:51:05 INFO None 5155249: status FINISHED
2026-07-14 11:51:05 INFO None 5155251: status FINISHED
2026-07-14 11:51:05 INFO Jobs still running: ['5155244']. Waiting...
2026-07-14 11:51:20 INFO None 5155226: status FINISHED
2026-07-14 11:51:20 INFO None 5155227: status FINISHED (not in squeue)
2026-07-14 11:51:20 INFO None 5155228: status FINISHED (not in squeue)
2026-07-14 11:51:20 INFO None 5155231: status FINISHED
2026-07-14 11:51:20 INFO None 5155239: status FINISHED
2026-07-14 11:51:20 INFO None 5155240: status FINISHED
2026-07-14 11:51:20 INFO None 5155241: status FINISHED
2026-07-14 11:51:20 INFO None 5155242: status FINISHED
2026-07-14 11:51:20 INFO None 5155243: status FINISHED
2026-07-14 11:51:20 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:51:20 INFO None 5155245: status FINISHED
2026-07-14 11:51:20 INFO None 5155246: status FINISHED
2026-07-14 11:51:20 INFO None 5155247: status FINISHED
2026-07-14 11:51:20 INFO None 5155249: status FINISHED
2026-07-14 11:51:20 INFO None 5155251: status FINISHED
2026-07-14 11:51:20 INFO Jobs still running: ['5155244']. Waiting...
2026-07-14 11:51:35 INFO None 5155226: status FINISHED
2026-07-14 11:51:35 INFO None 5155227: status FINISHED (not in squeue)
2026-07-14 11:51:35 INFO None 5155228: status FINISHED (not in squeue)
2026-07-14 11:51:35 INFO None 5155231: status FINISHED
2026-07-14 11:51:35 INFO None 5155239: status FINISHED
2026-07-14 11:51:35 INFO None 5155240: status FINISHED
2026-07-14 11:51:35 INFO None 5155241: status FINISHED
2026-07-14 11:51:35 INFO None 5155242: status FINISHED
2026-07-14 11:51:35 INFO None 5155243: status FINISHED
2026-07-14 11:51:35 INFO None 5155244: status RUNNING/PENDING
2026-07-14 11:51:35 INFO None 5155245: status FINISHED
2026-07-14 11:51:35 INFO None 5155246: status FINISHED
2026-07-14 11:51:35 INFO None 5155247: status FINISHED
2026-07-14 11:51:35 INFO None 5155249: status FINISHED
2026-07-14 11:51:35 INFO None 5155251: status FINISHED
2026-07-14 11:51:35 INFO Jobs still running: ['5155244']. Waiting...
2026-07-14 11:51:50 INFO None 5155226: status FINISHED
2026-07-14 11:53:48 INFO None 5155227: status FINISHED (not in squeue)
2026-07-14 11:53:48 INFO None 5155228: status FINISHED (not in squeue)
2026-07-14 11:53:49 INFO None 5155231: status FINISHED
2026-07-14 11:53:49 INFO None 5155239: status FINISHED
2026-07-14 11:53:49 INFO None 5155240: status FINISHED
2026-07-14 11:53:49 INFO None 5155241: status FINISHED
2026-07-14 11:53:49 INFO None 5155242: status FINISHED
2026-07-14 11:53:49 INFO None 5155243: status FINISHED
2026-07-14 11:53:49 INFO None 5155244: status FINISHED
2026-07-14 11:53:49 INFO None 5155245: status FINISHED
2026-07-14 11:53:49 INFO None 5155246: status FINISHED
2026-07-14 11:53:49 INFO None 5155247: status FINISHED
2026-07-14 11:53:49 INFO None 5155249: status FINISHED
2026-07-14 11:53:49 INFO None 5155251: status FINISHED
2026-07-14 11:53:49 INFO Jobs ['5155226', '5155227', '5155228', '5155231', '5155239', '5155240', '5155241', '5155242', '5155243', '5155244', '5155245', '5155246', '5155247', '5155249', '5155251'] have finished
2026-07-14 11:53:49 INFO Checking restart files were created ...
2026-07-14 11:53:49 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020601_8_ENS1.nc(3005806795 bytes)
2026-07-14 11:53:49 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020601_8_ENS2.nc(3005806795 bytes)
2026-07-14 11:53:49 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020601_8_ENS3.nc(3005806795 bytes)
2026-07-14 11:53:49 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020601_8_ENS4.nc(3005806795 bytes)
2026-07-14 11:53:49 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020601_8_ENS5.nc(3005806795 bytes)
2026-07-14 11:53:49 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020601_8_ENS6.nc(3005806795 bytes)
2026-07-14 11:53:49 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020601_8_ENS7.nc(3005806795 bytes)
2026-07-14 11:53:49 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020601_8_ENS8.nc(3005806795 bytes)
2026-07-14 11:53:49 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020601_8_ENS9.nc(3005806795 bytes)
2026-07-14 11:53:49 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020601_8_ENS10.nc(3005806795 bytes)
2026-07-14 11:53:49 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020601_8_ENS11.nc(3005806795 bytes)
2026-07-14 11:53:49 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020601_8_ENS12.nc(3005806795 bytes)
2026-07-14 11:53:49 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020601_8_ENS13.nc(3005806795 bytes)
2026-07-14 11:53:49 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020601_8_ENS14.nc(3005806795 bytes)
2026-07-14 11:53:49 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020601_8_ENS15.nc(3005806795 bytes)
2026-07-14 11:53:49 INFO  Run_model() completed successfully.
2026-07-14 11:53:49 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 11:53:49 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 09:00:00 days=153072 seconds=32400
2026-07-14 11:53:49 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-14 11:53:49 INFO [TIME] increment current_time 2020-02-06 01:00:00 -> 2020-02-06 09:00:00
2026-07-14 11:53:49 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 11:53:49 INFO ---------->>> Running process_satellite_data()
2026-07-14 11:53:49 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12002.nc
2026-07-14 11:53:49 INFO ---------->>> Running run_obs_converter()
2026-07-14 11:53:49 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_33067_153072.out
2026-07-14 11:53:49 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_33067_153072.out
2026-07-14 11:53:49 INFO ---------->>> Running DART
2026-07-14 11:53:49 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-14 11:53:49 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/chim_ENS1_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/chim_ENS1_2020020609_1_out_toDART.nc
2026-07-14 11:53:49 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/chim_ENS2_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/chim_ENS2_2020020609_1_out_toDART.nc
2026-07-14 11:53:50 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/chim_ENS3_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/chim_ENS3_2020020609_1_out_toDART.nc
2026-07-14 11:53:50 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/chim_ENS4_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/chim_ENS4_2020020609_1_out_toDART.nc
2026-07-14 11:53:50 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/chim_ENS5_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/chim_ENS5_2020020609_1_out_toDART.nc
2026-07-14 11:53:51 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/chim_ENS6_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/chim_ENS6_2020020609_1_out_toDART.nc
2026-07-14 11:53:51 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/chim_ENS7_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/chim_ENS7_2020020609_1_out_toDART.nc
2026-07-14 11:53:52 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/chim_ENS8_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/chim_ENS8_2020020609_1_out_toDART.nc
2026-07-14 11:53:52 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/chim_ENS9_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/chim_ENS9_2020020609_1_out_toDART.nc
2026-07-14 11:53:52 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/chim_ENS10_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/chim_ENS10_2020020609_1_out_toDART.nc
2026-07-14 11:53:53 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/chim_ENS11_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/chim_ENS11_2020020609_1_out_toDART.nc
2026-07-14 11:53:53 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/chim_ENS12_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/chim_ENS12_2020020609_1_out_toDART.nc
2026-07-14 11:53:54 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/chim_ENS13_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/chim_ENS13_2020020609_1_out_toDART.nc
2026-07-14 11:53:54 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/chim_ENS14_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/chim_ENS14_2020020609_1_out_toDART.nc
2026-07-14 11:53:54 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/chim_ENS15_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/chim_ENS15_2020020609_1_out_toDART.nc
2026-07-14 11:53:55 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-14 11:53:55 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-14 11:53:55 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-14 11:53:55 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-14 11:53:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-14 11:53:55 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-14 11:54:04 INFO Found: []
2026-07-14 11:54:04 INFO No job id returned by command ./run_filter.bsh
2026-07-14 11:54:04 INFO No monitoring will be performed
2026-07-14 11:54:04 INFO Moving DART output files to analysis and preassim directories for date 2020020609 if present ...
2026-07-14 11:54:04 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-14 11:54:04 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-14 11:54:04 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-14 11:54:04 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-14 11:54:04 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-14 11:54:04 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-14 11:54:04 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-14 11:54:04 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-14 11:54:04 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-14 11:54:04 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-14 11:54:04 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-14 11:54:04 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-14 11:54:05 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-14 11:54:05 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-14 11:54:05 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-14 11:54:05 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-14 11:54:05 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-14 11:54:05 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-14 11:54:05 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-14 11:54:05 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-14 11:54:05 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-14 11:54:05 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-14 11:54:05 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-14 11:54:05 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-14 11:54:05 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-14 11:54:05 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-14 11:54:05 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-14 11:54:05 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-14 11:54:05 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-14 11:54:05 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-14 11:54:05 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-14 11:54:05 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-14 11:54:05 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-14 11:54:05 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-14 11:54:05 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-14 11:54:08 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-14 11:54:08 INFO run_dart() is DONE.
2026-07-14 11:54:08 INFO ---------->>> Running update_pollutant_in_end()
2026-07-14 11:54:08 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020601_8_ENS1.nc
2026-07-14 11:54:23 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 11:54:23 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020601_8_ENS2.nc
2026-07-14 11:54:38 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 11:54:39 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020601_8_ENS3.nc
2026-07-14 11:54:53 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 11:54:54 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020601_8_ENS4.nc
2026-07-14 11:55:09 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 11:55:09 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020601_8_ENS5.nc
2026-07-14 11:55:23 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 11:55:24 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020601_8_ENS6.nc
2026-07-14 11:55:39 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 11:55:39 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020601_8_ENS7.nc
2026-07-14 11:55:54 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 11:55:54 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020601_8_ENS8.nc
2026-07-14 11:58:30 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 11:58:30 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020601_8_ENS9.nc
2026-07-14 11:58:45 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 11:58:46 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020601_8_ENS10.nc
2026-07-14 11:59:00 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 11:59:01 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020601_8_ENS11.nc
2026-07-14 11:59:16 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 11:59:16 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020601_8_ENS12.nc
2026-07-14 11:59:31 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 11:59:31 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020601_8_ENS13.nc
2026-07-14 11:59:46 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 11:59:46 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020601_8_ENS14.nc
2026-07-14 12:00:01 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 12:00:01 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020601_8_ENS15.nc
2026-07-14 12:00:16 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 12:00:16 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-14 12:00:16 INFO [TIME] step_end current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 12:00:16 INFO [TIME] step_start current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 12:00:16 INFO [TIME] window start=2020-02-06 09:00:00 end=2020-02-06 11:00:00 run_hours=2 has_assimilation=True
2026-07-14 12:00:16 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 12:00:18 INFO Hourly dataset computed and listing created
2026-07-14 12:00:24 INFO Hourly dataset computed
2026-07-14 12:00:24 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 12:00:25 INFO Hourly dataset computed and listing created
2026-07-14 12:00:25 INFO Hourly dataset computed
2026-07-14 12:00:26 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 12:00:26 INFO Hourly dataset computed and listing created
2026-07-14 12:00:27 INFO Hourly dataset computed
2026-07-14 12:00:27 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 12:00:28 INFO Hourly dataset computed and listing created
2026-07-14 12:00:29 INFO Hourly dataset computed
2026-07-14 12:00:29 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 12:00:30 INFO Hourly dataset computed and listing created
2026-07-14 12:00:31 INFO Hourly dataset computed
2026-07-14 12:00:31 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 12:00:32 INFO Hourly dataset computed and listing created
2026-07-14 12:00:33 INFO Hourly dataset computed
2026-07-14 12:00:33 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 12:00:34 INFO Hourly dataset computed and listing created
2026-07-14 12:00:35 INFO Hourly dataset computed
2026-07-14 12:00:35 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 12:00:36 INFO Hourly dataset computed and listing created
2026-07-14 12:00:37 INFO Hourly dataset computed
2026-07-14 12:00:37 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 12:00:38 INFO Hourly dataset computed and listing created
2026-07-14 12:00:39 INFO Hourly dataset computed
2026-07-14 12:00:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 12:00:39 INFO Hourly dataset computed and listing created
2026-07-14 12:00:40 INFO Hourly dataset computed
2026-07-14 12:00:40 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 12:00:41 INFO Hourly dataset computed and listing created
2026-07-14 12:00:42 INFO Hourly dataset computed
2026-07-14 12:00:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 12:00:43 INFO Hourly dataset computed and listing created
2026-07-14 12:00:44 INFO Hourly dataset computed
2026-07-14 12:00:44 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 12:00:45 INFO Hourly dataset computed and listing created
2026-07-14 12:00:46 INFO Hourly dataset computed
2026-07-14 12:00:46 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 12:00:47 INFO Hourly dataset computed and listing created
2026-07-14 12:00:48 INFO Hourly dataset computed
2026-07-14 12:00:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 12:00:49 INFO Hourly dataset computed and listing created
2026-07-14 12:00:50 INFO Hourly dataset computed
2026-07-14 12:00:50 INFO ---------->>> Running CHIMERE model from 2020-02-06 09:00:00 to 2020-02-06 11:00:00
2026-07-14 12:00:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 12:00:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1
2026-07-14 12:00:50 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020601_8_ENS1.nc
2026-07-14 12:00:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-14 12:00:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 12:00:50 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-14 12:00:50 INFO Queuing job for member 1...
2026-07-14 12:00:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 12:00:50 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-14 12:00:51 INFO Found: ['5155319']
2026-07-14 12:00:56 INFO [TGCC-IRENE] Submitted job with ID:['5155319']
2026-07-14 12:00:56 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 12:00:56 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2
2026-07-14 12:00:56 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020601_8_ENS2.nc
2026-07-14 12:00:56 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-14 12:00:56 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 12:00:56 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-14 12:00:56 INFO Queuing job for member 2...
2026-07-14 12:00:56 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 12:00:56 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-14 12:00:56 INFO Found: ['5155320']
2026-07-14 12:01:01 INFO [TGCC-IRENE] Submitted job with ID:['5155320']
2026-07-14 12:01:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 12:01:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3
2026-07-14 12:01:01 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020601_8_ENS3.nc
2026-07-14 12:03:09 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-14 12:03:09 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 12:03:09 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-14 12:03:09 INFO Queuing job for member 3...
2026-07-14 12:03:09 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 12:03:09 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-14 12:03:09 INFO Found: ['5155328']
2026-07-14 12:03:14 INFO [TGCC-IRENE] Submitted job with ID:['5155328']
2026-07-14 12:03:14 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 12:03:14 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4
2026-07-14 12:03:14 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020601_8_ENS4.nc
2026-07-14 12:03:14 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-14 12:03:14 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 12:03:14 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-14 12:03:14 INFO Queuing job for member 4...
2026-07-14 12:03:14 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 12:03:14 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-14 12:03:15 INFO Found: ['5155329']
2026-07-14 12:03:20 INFO [TGCC-IRENE] Submitted job with ID:['5155329']
2026-07-14 12:03:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 12:03:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5
2026-07-14 12:03:20 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020601_8_ENS5.nc
2026-07-14 12:03:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-14 12:03:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 12:03:20 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-14 12:03:20 INFO Queuing job for member 5...
2026-07-14 12:03:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 12:03:20 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-14 12:03:21 INFO Found: ['5155330']
2026-07-14 12:03:26 INFO [TGCC-IRENE] Submitted job with ID:['5155330']
2026-07-14 12:03:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 12:03:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6
2026-07-14 12:03:26 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020601_8_ENS6.nc
2026-07-14 12:03:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-14 12:03:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 12:03:26 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-14 12:03:26 INFO Queuing job for member 6...
2026-07-14 12:03:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 12:03:26 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-14 12:03:27 INFO Found: ['5155331']
2026-07-14 12:03:32 INFO [TGCC-IRENE] Submitted job with ID:['5155331']
2026-07-14 12:03:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 12:03:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7
2026-07-14 12:03:32 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020601_8_ENS7.nc
2026-07-14 12:03:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-14 12:03:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 12:03:32 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-14 12:03:32 INFO Queuing job for member 7...
2026-07-14 12:03:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 12:03:32 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-14 12:03:32 INFO Found: ['5155332']
2026-07-14 12:03:37 INFO [TGCC-IRENE] Submitted job with ID:['5155332']
2026-07-14 12:03:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 12:03:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8
2026-07-14 12:03:37 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020601_8_ENS8.nc
2026-07-14 12:03:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-14 12:03:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 12:03:37 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-14 12:03:37 INFO Queuing job for member 8...
2026-07-14 12:03:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 12:03:37 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-14 12:03:38 INFO Found: ['5155333']
2026-07-14 12:03:43 INFO [TGCC-IRENE] Submitted job with ID:['5155333']
2026-07-14 12:03:43 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 12:03:43 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9
2026-07-14 12:03:43 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020601_8_ENS9.nc
2026-07-14 12:03:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-14 12:03:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 12:03:43 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-14 12:03:43 INFO Queuing job for member 9...
2026-07-14 12:03:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 12:03:43 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-14 12:03:44 INFO Found: ['5155334']
2026-07-14 12:03:49 INFO [TGCC-IRENE] Submitted job with ID:['5155334']
2026-07-14 12:03:49 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 12:03:49 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10
2026-07-14 12:03:49 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020601_8_ENS10.nc
2026-07-14 12:03:49 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-14 12:03:49 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 12:03:49 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-14 12:03:49 INFO Queuing job for member 10...
2026-07-14 12:03:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 12:03:49 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-14 12:03:49 INFO Found: ['5155335']
2026-07-14 12:03:54 INFO [TGCC-IRENE] Submitted job with ID:['5155335']
2026-07-14 12:03:54 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 12:03:54 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11
2026-07-14 12:03:54 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020601_8_ENS11.nc
2026-07-14 12:03:54 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-14 12:03:54 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 12:03:54 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-14 12:03:54 INFO Queuing job for member 11...
2026-07-14 12:03:54 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 12:03:54 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-14 12:03:55 INFO Found: ['5155337']
2026-07-14 12:04:00 INFO [TGCC-IRENE] Submitted job with ID:['5155337']
2026-07-14 12:04:00 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 12:04:00 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12
2026-07-14 12:04:00 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020601_8_ENS12.nc
2026-07-14 12:04:00 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-14 12:04:00 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 12:04:00 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-14 12:04:00 INFO Queuing job for member 12...
2026-07-14 12:04:00 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 12:04:00 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-14 12:04:01 INFO Found: ['5155338']
2026-07-14 12:04:06 INFO [TGCC-IRENE] Submitted job with ID:['5155338']
2026-07-14 12:04:06 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 12:04:06 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13
2026-07-14 12:04:06 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020601_8_ENS13.nc
2026-07-14 12:04:06 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-14 12:04:06 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 12:04:06 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-14 12:04:06 INFO Queuing job for member 13...
2026-07-14 12:04:06 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 12:04:06 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-14 12:04:07 INFO Found: ['5155340']
2026-07-14 12:04:12 INFO [TGCC-IRENE] Submitted job with ID:['5155340']
2026-07-14 12:04:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 12:04:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14
2026-07-14 12:04:12 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020601_8_ENS14.nc
2026-07-14 12:04:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-14 12:04:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 12:04:12 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-14 12:04:12 INFO Queuing job for member 14...
2026-07-14 12:04:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 12:04:12 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-14 12:04:13 INFO Found: ['5155341']
2026-07-14 12:04:18 INFO [TGCC-IRENE] Submitted job with ID:['5155341']
2026-07-14 12:04:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 12:04:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15
2026-07-14 12:04:18 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020601_8_ENS15.nc
2026-07-14 12:04:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-14 12:04:18 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 12:04:18 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-14 12:04:18 INFO Queuing job for member 15...
2026-07-14 12:04:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 12:04:18 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-14 12:04:18 INFO Found: ['5155342']
2026-07-14 12:04:23 INFO [TGCC-IRENE] Submitted job with ID:['5155342']
2026-07-14 12:04:23 INFO Checking job status ...
2026-07-14 12:04:23 INFO None 5155319: status RUNNING/PENDING
2026-07-14 12:04:23 INFO None 5155320: status RUNNING/PENDING
2026-07-14 12:04:23 INFO None 5155328: status RUNNING/PENDING
2026-07-14 12:04:23 INFO None 5155329: status RUNNING/PENDING
2026-07-14 12:04:23 INFO None 5155330: status RUNNING/PENDING
2026-07-14 12:04:24 INFO None 5155331: status RUNNING/PENDING
2026-07-14 12:04:24 INFO None 5155332: status RUNNING/PENDING
2026-07-14 12:04:24 INFO None 5155333: status RUNNING/PENDING
2026-07-14 12:04:24 INFO None 5155334: status RUNNING/PENDING
2026-07-14 12:04:24 INFO None 5155335: status RUNNING/PENDING
2026-07-14 12:04:24 INFO None 5155337: status RUNNING/PENDING
2026-07-14 12:04:24 INFO None 5155338: status RUNNING/PENDING
2026-07-14 12:04:24 INFO None 5155340: status RUNNING/PENDING
2026-07-14 12:04:24 INFO None 5155341: status RUNNING/PENDING
2026-07-14 12:04:24 INFO None 5155342: status RUNNING/PENDING
2026-07-14 12:04:24 INFO Jobs still running: ['5155319', '5155320', '5155328', '5155329', '5155330', '5155331', '5155332', '5155333', '5155334', '5155335', '5155337', '5155338', '5155340', '5155341', '5155342']. Waiting...
2026-07-14 12:04:39 INFO None 5155319: status RUNNING/PENDING
2026-07-14 12:04:39 INFO None 5155320: status RUNNING/PENDING
2026-07-14 12:04:39 INFO None 5155328: status RUNNING/PENDING
2026-07-14 12:04:39 INFO None 5155329: status RUNNING/PENDING
2026-07-14 12:04:39 INFO None 5155330: status RUNNING/PENDING
2026-07-14 12:04:39 INFO None 5155331: status RUNNING/PENDING
2026-07-14 12:04:39 INFO None 5155332: status RUNNING/PENDING
2026-07-14 12:04:39 INFO None 5155333: status RUNNING/PENDING
2026-07-14 12:04:39 INFO None 5155334: status RUNNING/PENDING
2026-07-14 12:04:39 INFO None 5155335: status RUNNING/PENDING
2026-07-14 12:04:39 INFO None 5155337: status RUNNING/PENDING
2026-07-14 12:04:39 INFO None 5155338: status RUNNING/PENDING
2026-07-14 12:04:39 INFO None 5155340: status RUNNING/PENDING
2026-07-14 12:04:39 INFO None 5155341: status RUNNING/PENDING
2026-07-14 12:04:39 INFO None 5155342: status RUNNING/PENDING
2026-07-14 12:04:39 INFO Jobs still running: ['5155319', '5155320', '5155328', '5155329', '5155330', '5155331', '5155332', '5155333', '5155334', '5155335', '5155337', '5155338', '5155340', '5155341', '5155342']. Waiting...
2026-07-14 12:04:54 INFO None 5155319: status RUNNING/PENDING
2026-07-14 12:04:54 INFO None 5155320: status RUNNING/PENDING
2026-07-14 12:04:54 INFO None 5155328: status RUNNING/PENDING
2026-07-14 12:04:54 INFO None 5155329: status RUNNING/PENDING
2026-07-14 12:04:54 INFO None 5155330: status RUNNING/PENDING
2026-07-14 12:04:54 INFO None 5155331: status RUNNING/PENDING
2026-07-14 12:04:54 INFO None 5155332: status RUNNING/PENDING
2026-07-14 12:04:54 INFO None 5155333: status RUNNING/PENDING
2026-07-14 12:04:54 INFO None 5155334: status RUNNING/PENDING
2026-07-14 12:04:54 INFO None 5155335: status RUNNING/PENDING
2026-07-14 12:04:54 INFO None 5155337: status RUNNING/PENDING
2026-07-14 12:04:54 INFO None 5155338: status RUNNING/PENDING
2026-07-14 12:04:54 INFO None 5155340: status RUNNING/PENDING
2026-07-14 12:04:54 INFO None 5155341: status RUNNING/PENDING
2026-07-14 12:04:54 INFO None 5155342: status RUNNING/PENDING
2026-07-14 12:04:54 INFO Jobs still running: ['5155319', '5155320', '5155328', '5155329', '5155330', '5155331', '5155332', '5155333', '5155334', '5155335', '5155337', '5155338', '5155340', '5155341', '5155342']. Waiting...
2026-07-14 12:05:09 INFO None 5155319: status RUNNING/PENDING
2026-07-14 12:05:09 INFO None 5155320: status RUNNING/PENDING
2026-07-14 12:05:09 INFO None 5155328: status RUNNING/PENDING
2026-07-14 12:05:09 INFO None 5155329: status RUNNING/PENDING
2026-07-14 12:05:09 INFO None 5155330: status RUNNING/PENDING
2026-07-14 12:05:09 INFO None 5155331: status RUNNING/PENDING
2026-07-14 12:05:10 INFO None 5155332: status RUNNING/PENDING
2026-07-14 12:05:10 INFO None 5155333: status RUNNING/PENDING
2026-07-14 12:05:10 INFO None 5155334: status RUNNING/PENDING
2026-07-14 12:05:10 INFO None 5155335: status RUNNING/PENDING
2026-07-14 12:05:10 INFO None 5155337: status RUNNING/PENDING
2026-07-14 12:05:10 INFO None 5155338: status RUNNING/PENDING
2026-07-14 12:05:10 INFO None 5155340: status RUNNING/PENDING
2026-07-14 12:05:10 INFO None 5155341: status RUNNING/PENDING
2026-07-14 12:05:10 INFO None 5155342: status RUNNING/PENDING
2026-07-14 12:05:10 INFO Jobs still running: ['5155319', '5155320', '5155328', '5155329', '5155330', '5155331', '5155332', '5155333', '5155334', '5155335', '5155337', '5155338', '5155340', '5155341', '5155342']. Waiting...
2026-07-14 12:05:25 INFO None 5155319: status RUNNING/PENDING
2026-07-14 12:05:25 INFO None 5155320: status RUNNING/PENDING
2026-07-14 12:05:25 INFO None 5155328: status RUNNING/PENDING
2026-07-14 12:05:25 INFO None 5155329: status RUNNING/PENDING
2026-07-14 12:05:25 INFO None 5155330: status RUNNING/PENDING
2026-07-14 12:05:25 INFO None 5155331: status RUNNING/PENDING
2026-07-14 12:05:25 INFO None 5155332: status RUNNING/PENDING
2026-07-14 12:05:25 INFO None 5155333: status RUNNING/PENDING
2026-07-14 12:05:25 INFO None 5155334: status RUNNING/PENDING
2026-07-14 12:05:25 INFO None 5155335: status RUNNING/PENDING
2026-07-14 12:05:25 INFO None 5155337: status RUNNING/PENDING
2026-07-14 12:05:25 INFO None 5155338: status RUNNING/PENDING
2026-07-14 12:05:25 INFO None 5155340: status RUNNING/PENDING
2026-07-14 12:05:25 INFO None 5155341: status RUNNING/PENDING
2026-07-14 12:05:25 INFO None 5155342: status RUNNING/PENDING
2026-07-14 12:05:25 INFO Jobs still running: ['5155319', '5155320', '5155328', '5155329', '5155330', '5155331', '5155332', '5155333', '5155334', '5155335', '5155337', '5155338', '5155340', '5155341', '5155342']. Waiting...
2026-07-14 12:05:40 INFO None 5155319: status RUNNING/PENDING
2026-07-14 12:05:40 INFO None 5155320: status RUNNING/PENDING
2026-07-14 12:05:40 INFO None 5155328: status RUNNING/PENDING
2026-07-14 12:05:40 INFO None 5155329: status RUNNING/PENDING
2026-07-14 12:05:40 INFO None 5155330: status RUNNING/PENDING
2026-07-14 12:05:40 INFO None 5155331: status RUNNING/PENDING
2026-07-14 12:05:40 INFO None 5155332: status RUNNING/PENDING
2026-07-14 12:05:40 INFO None 5155333: status RUNNING/PENDING
2026-07-14 12:05:40 INFO None 5155334: status RUNNING/PENDING
2026-07-14 12:05:40 INFO None 5155335: status RUNNING/PENDING
2026-07-14 12:05:40 INFO None 5155337: status RUNNING/PENDING
2026-07-14 12:05:40 INFO None 5155338: status RUNNING/PENDING
2026-07-14 12:05:40 INFO None 5155340: status RUNNING/PENDING
2026-07-14 12:05:40 INFO None 5155341: status RUNNING/PENDING
2026-07-14 12:05:40 INFO None 5155342: status RUNNING/PENDING
2026-07-14 12:05:40 INFO Jobs still running: ['5155319', '5155320', '5155328', '5155329', '5155330', '5155331', '5155332', '5155333', '5155334', '5155335', '5155337', '5155338', '5155340', '5155341', '5155342']. Waiting...
2026-07-14 12:05:55 INFO None 5155319: status RUNNING/PENDING
2026-07-14 12:05:55 INFO None 5155320: status RUNNING/PENDING
2026-07-14 12:05:55 INFO None 5155328: status RUNNING/PENDING
2026-07-14 12:05:55 INFO None 5155329: status RUNNING/PENDING
2026-07-14 12:05:55 INFO None 5155330: status RUNNING/PENDING
2026-07-14 12:05:55 INFO None 5155331: status RUNNING/PENDING
2026-07-14 12:05:55 INFO None 5155332: status RUNNING/PENDING
2026-07-14 12:05:55 INFO None 5155333: status RUNNING/PENDING
2026-07-14 12:05:56 INFO None 5155334: status RUNNING/PENDING
2026-07-14 12:05:56 INFO None 5155335: status RUNNING/PENDING
2026-07-14 12:05:56 INFO None 5155337: status RUNNING/PENDING
2026-07-14 12:05:56 INFO None 5155338: status RUNNING/PENDING
2026-07-14 12:05:56 INFO None 5155340: status RUNNING/PENDING
2026-07-14 12:05:56 INFO None 5155341: status RUNNING/PENDING
2026-07-14 12:05:56 INFO None 5155342: status RUNNING/PENDING
2026-07-14 12:05:56 INFO Jobs still running: ['5155319', '5155320', '5155328', '5155329', '5155330', '5155331', '5155332', '5155333', '5155334', '5155335', '5155337', '5155338', '5155340', '5155341', '5155342']. Waiting...
2026-07-14 12:06:11 INFO None 5155319: status RUNNING/PENDING
2026-07-14 12:06:11 INFO None 5155320: status RUNNING/PENDING
2026-07-14 12:06:11 INFO None 5155328: status RUNNING/PENDING
2026-07-14 12:06:11 INFO None 5155329: status RUNNING/PENDING
2026-07-14 12:06:11 INFO None 5155330: status RUNNING/PENDING
2026-07-14 12:06:11 INFO None 5155331: status RUNNING/PENDING
2026-07-14 12:06:11 INFO None 5155332: status RUNNING/PENDING
2026-07-14 12:06:11 INFO None 5155333: status RUNNING/PENDING
2026-07-14 12:06:11 INFO None 5155334: status RUNNING/PENDING
2026-07-14 12:06:11 INFO None 5155335: status RUNNING/PENDING
2026-07-14 12:06:11 INFO None 5155337: status RUNNING/PENDING
2026-07-14 12:06:11 INFO None 5155338: status RUNNING/PENDING
2026-07-14 12:06:11 INFO None 5155340: status RUNNING/PENDING
2026-07-14 12:06:11 INFO None 5155341: status RUNNING/PENDING
2026-07-14 12:06:11 INFO None 5155342: status RUNNING/PENDING
2026-07-14 12:06:11 INFO Jobs still running: ['5155319', '5155320', '5155328', '5155329', '5155330', '5155331', '5155332', '5155333', '5155334', '5155335', '5155337', '5155338', '5155340', '5155341', '5155342']. Waiting...
2026-07-14 12:06:26 INFO None 5155319: status RUNNING/PENDING
2026-07-14 12:06:26 INFO None 5155320: status RUNNING/PENDING
2026-07-14 12:06:26 INFO None 5155328: status RUNNING/PENDING
2026-07-14 12:06:26 INFO None 5155329: status RUNNING/PENDING
2026-07-14 12:06:26 INFO None 5155330: status RUNNING/PENDING
2026-07-14 12:06:26 INFO None 5155331: status RUNNING/PENDING
2026-07-14 12:06:26 INFO None 5155332: status RUNNING/PENDING
2026-07-14 12:06:26 INFO None 5155333: status RUNNING/PENDING
2026-07-14 12:06:26 INFO None 5155334: status RUNNING/PENDING
2026-07-14 12:06:26 INFO None 5155335: status FINISHED
2026-07-14 12:06:26 INFO None 5155337: status RUNNING/PENDING
2026-07-14 12:06:26 INFO None 5155338: status RUNNING/PENDING
2026-07-14 12:06:26 INFO None 5155340: status RUNNING/PENDING
2026-07-14 12:06:26 INFO None 5155341: status RUNNING/PENDING
2026-07-14 12:06:26 INFO None 5155342: status RUNNING/PENDING
2026-07-14 12:06:26 INFO Jobs still running: ['5155319', '5155320', '5155328', '5155329', '5155330', '5155331', '5155332', '5155333', '5155334', '5155337', '5155338', '5155340', '5155341', '5155342']. Waiting...
2026-07-14 12:06:41 INFO None 5155319: status RUNNING/PENDING
2026-07-14 12:06:41 INFO None 5155320: status RUNNING/PENDING
2026-07-14 12:06:41 INFO None 5155328: status RUNNING/PENDING
2026-07-14 12:06:41 INFO None 5155329: status RUNNING/PENDING
2026-07-14 12:06:41 INFO None 5155330: status RUNNING/PENDING
2026-07-14 12:06:41 INFO None 5155331: status RUNNING/PENDING
2026-07-14 12:06:41 INFO None 5155332: status RUNNING/PENDING
2026-07-14 12:06:41 INFO None 5155333: status RUNNING/PENDING
2026-07-14 12:06:41 INFO None 5155334: status RUNNING/PENDING
2026-07-14 12:06:41 INFO None 5155335: status FINISHED
2026-07-14 12:06:41 INFO None 5155337: status RUNNING/PENDING
2026-07-14 12:06:42 INFO None 5155338: status RUNNING/PENDING
2026-07-14 12:06:42 INFO None 5155340: status RUNNING/PENDING
2026-07-14 12:06:42 INFO None 5155341: status RUNNING/PENDING
2026-07-14 12:06:42 INFO None 5155342: status RUNNING/PENDING
2026-07-14 12:06:42 INFO Jobs still running: ['5155319', '5155320', '5155328', '5155329', '5155330', '5155331', '5155332', '5155333', '5155334', '5155337', '5155338', '5155340', '5155341', '5155342']. Waiting...
2026-07-14 12:06:57 INFO None 5155319: status RUNNING/PENDING
2026-07-14 12:06:57 INFO None 5155320: status RUNNING/PENDING
2026-07-14 12:06:57 INFO None 5155328: status RUNNING/PENDING
2026-07-14 12:06:57 INFO None 5155329: status RUNNING/PENDING
2026-07-14 12:06:57 INFO None 5155330: status RUNNING/PENDING
2026-07-14 12:06:57 INFO None 5155331: status RUNNING/PENDING
2026-07-14 12:06:57 INFO None 5155332: status RUNNING/PENDING
2026-07-14 12:06:57 INFO None 5155333: status RUNNING/PENDING
2026-07-14 12:06:57 INFO None 5155334: status RUNNING/PENDING
2026-07-14 12:06:57 INFO None 5155335: status FINISHED
2026-07-14 12:06:57 INFO None 5155337: status RUNNING/PENDING
2026-07-14 12:06:57 INFO None 5155338: status RUNNING/PENDING
2026-07-14 12:06:57 INFO None 5155340: status RUNNING/PENDING
2026-07-14 12:06:57 INFO None 5155341: status RUNNING/PENDING
2026-07-14 12:06:57 INFO None 5155342: status RUNNING/PENDING
2026-07-14 12:06:57 INFO Jobs still running: ['5155319', '5155320', '5155328', '5155329', '5155330', '5155331', '5155332', '5155333', '5155334', '5155337', '5155338', '5155340', '5155341', '5155342']. Waiting...
2026-07-14 12:07:12 INFO None 5155319: status RUNNING/PENDING
2026-07-14 12:07:12 INFO None 5155320: status RUNNING/PENDING
2026-07-14 12:07:12 INFO None 5155328: status RUNNING/PENDING
2026-07-14 12:07:12 INFO None 5155329: status FINISHED
2026-07-14 12:07:12 INFO None 5155330: status FINISHED
2026-07-14 12:07:12 INFO None 5155331: status FINISHED
2026-07-14 12:07:12 INFO None 5155332: status FINISHED
2026-07-14 12:07:12 INFO None 5155333: status FINISHED
2026-07-14 12:07:12 INFO None 5155334: status FINISHED
2026-07-14 12:07:12 INFO None 5155335: status FINISHED
2026-07-14 12:07:12 INFO None 5155337: status RUNNING/PENDING
2026-07-14 12:07:12 INFO None 5155338: status RUNNING/PENDING
2026-07-14 12:07:12 INFO None 5155340: status RUNNING/PENDING
2026-07-14 12:07:12 INFO None 5155341: status RUNNING/PENDING
2026-07-14 12:07:12 INFO None 5155342: status FINISHED
2026-07-14 12:07:12 INFO Jobs still running: ['5155319', '5155320', '5155328', '5155337', '5155338', '5155340', '5155341']. Waiting...
2026-07-14 12:07:27 INFO None 5155319: status RUNNING/PENDING
2026-07-14 12:09:25 INFO None 5155320: status FINISHED
2026-07-14 12:09:25 INFO None 5155328: status FINISHED
2026-07-14 12:09:25 INFO None 5155329: status FINISHED
2026-07-14 12:09:26 INFO None 5155330: status FINISHED
2026-07-14 12:09:26 INFO None 5155331: status FINISHED
2026-07-14 12:09:26 INFO None 5155332: status FINISHED
2026-07-14 12:09:26 INFO None 5155333: status FINISHED
2026-07-14 12:09:26 INFO None 5155334: status FINISHED
2026-07-14 12:09:26 INFO None 5155335: status FINISHED
2026-07-14 12:09:26 INFO None 5155337: status RUNNING/PENDING
2026-07-14 12:09:26 INFO None 5155338: status RUNNING/PENDING
2026-07-14 12:09:26 INFO None 5155340: status RUNNING/PENDING
2026-07-14 12:09:26 INFO None 5155341: status RUNNING/PENDING
2026-07-14 12:09:26 INFO None 5155342: status FINISHED
2026-07-14 12:09:26 INFO Jobs still running: ['5155319', '5155337', '5155338', '5155340', '5155341']. Waiting...
2026-07-14 12:09:41 INFO None 5155319: status RUNNING/PENDING
2026-07-14 12:09:41 INFO None 5155320: status FINISHED
2026-07-14 12:09:41 INFO None 5155328: status FINISHED
2026-07-14 12:09:41 INFO None 5155329: status FINISHED
2026-07-14 12:09:41 INFO None 5155330: status FINISHED
2026-07-14 12:09:41 INFO None 5155331: status FINISHED
2026-07-14 12:09:41 INFO None 5155332: status FINISHED
2026-07-14 12:09:41 INFO None 5155333: status FINISHED
2026-07-14 12:09:41 INFO None 5155334: status FINISHED
2026-07-14 12:09:41 INFO None 5155335: status FINISHED
2026-07-14 12:09:41 INFO None 5155337: status RUNNING/PENDING
2026-07-14 12:09:41 INFO None 5155338: status RUNNING/PENDING
2026-07-14 12:09:41 INFO None 5155340: status FINISHED
2026-07-14 12:09:41 INFO None 5155341: status RUNNING/PENDING
2026-07-14 12:09:41 INFO None 5155342: status FINISHED
2026-07-14 12:09:41 INFO Jobs still running: ['5155319', '5155337', '5155338', '5155341']. Waiting...
2026-07-14 12:09:56 INFO None 5155319: status RUNNING/PENDING
2026-07-14 12:09:56 INFO None 5155320: status FINISHED
2026-07-14 12:09:56 INFO None 5155328: status FINISHED
2026-07-14 12:09:56 INFO None 5155329: status FINISHED
2026-07-14 12:09:56 INFO None 5155330: status FINISHED
2026-07-14 12:09:56 INFO None 5155331: status FINISHED
2026-07-14 12:09:56 INFO None 5155332: status FINISHED
2026-07-14 12:09:56 INFO None 5155333: status FINISHED
2026-07-14 12:09:56 INFO None 5155334: status FINISHED
2026-07-14 12:09:56 INFO None 5155335: status FINISHED
2026-07-14 12:09:56 INFO None 5155337: status RUNNING/PENDING
2026-07-14 12:09:56 INFO None 5155338: status RUNNING/PENDING
2026-07-14 12:09:56 INFO None 5155340: status FINISHED
2026-07-14 12:09:56 INFO None 5155341: status RUNNING/PENDING
2026-07-14 12:09:56 INFO None 5155342: status FINISHED
2026-07-14 12:09:56 INFO Jobs still running: ['5155319', '5155337', '5155338', '5155341']. Waiting...
2026-07-14 12:10:11 INFO None 5155319: status RUNNING/PENDING
2026-07-14 12:10:11 INFO None 5155320: status FINISHED
2026-07-14 12:10:11 INFO None 5155328: status FINISHED
2026-07-14 12:10:11 INFO None 5155329: status FINISHED
2026-07-14 12:10:11 INFO None 5155330: status FINISHED
2026-07-14 12:10:11 INFO None 5155331: status FINISHED
2026-07-14 12:10:11 INFO None 5155332: status FINISHED
2026-07-14 12:10:11 INFO None 5155333: status FINISHED
2026-07-14 12:10:12 INFO None 5155334: status FINISHED
2026-07-14 12:10:12 INFO None 5155335: status FINISHED
2026-07-14 12:10:12 INFO None 5155337: status FINISHED
2026-07-14 12:10:12 INFO None 5155338: status FINISHED
2026-07-14 12:10:12 INFO None 5155340: status FINISHED
2026-07-14 12:10:12 INFO None 5155341: status FINISHED
2026-07-14 12:10:12 INFO None 5155342: status FINISHED
2026-07-14 12:10:12 INFO Jobs still running: ['5155319']. Waiting...
2026-07-14 12:10:27 INFO None 5155319: status FINISHED
2026-07-14 12:10:27 INFO None 5155320: status FINISHED
2026-07-14 12:10:27 INFO None 5155328: status FINISHED
2026-07-14 12:10:27 INFO None 5155329: status FINISHED
2026-07-14 12:10:27 INFO None 5155330: status FINISHED
2026-07-14 12:10:27 INFO None 5155331: status FINISHED
2026-07-14 12:10:27 INFO None 5155332: status FINISHED
2026-07-14 12:10:27 INFO None 5155333: status FINISHED
2026-07-14 12:10:27 INFO None 5155334: status FINISHED
2026-07-14 12:10:27 INFO None 5155335: status FINISHED
2026-07-14 12:10:27 INFO None 5155337: status FINISHED
2026-07-14 12:10:27 INFO None 5155338: status FINISHED
2026-07-14 12:10:27 INFO None 5155340: status FINISHED
2026-07-14 12:10:27 INFO None 5155341: status FINISHED
2026-07-14 12:10:27 INFO None 5155342: status FINISHED
2026-07-14 12:10:27 INFO Jobs ['5155319', '5155320', '5155328', '5155329', '5155330', '5155331', '5155332', '5155333', '5155334', '5155335', '5155337', '5155338', '5155340', '5155341', '5155342'] have finished
2026-07-14 12:10:27 INFO Checking restart files were created ...
2026-07-14 12:10:27 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc(1002685915 bytes)
2026-07-14 12:10:27 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc(1002685915 bytes)
2026-07-14 12:10:27 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc(1002685915 bytes)
2026-07-14 12:10:27 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc(1002685915 bytes)
2026-07-14 12:10:27 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc(1002685915 bytes)
2026-07-14 12:10:27 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc(1002685915 bytes)
2026-07-14 12:10:27 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc(1002685915 bytes)
2026-07-14 12:10:27 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc(1002685915 bytes)
2026-07-14 12:10:27 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc(1002685915 bytes)
2026-07-14 12:10:27 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc(1002685915 bytes)
2026-07-14 12:10:27 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc(1002685915 bytes)
2026-07-14 12:10:27 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc(1002685915 bytes)
2026-07-14 12:10:27 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc(1002685915 bytes)
2026-07-14 12:10:27 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc(1002685915 bytes)
2026-07-14 12:10:27 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc(1002685915 bytes)
2026-07-14 12:10:27 INFO  Run_model() completed successfully.
2026-07-14 12:10:27 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 12:10:27 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 11:00:00 days=153072 seconds=39600
2026-07-14 12:10:27 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-14 12:10:27 INFO [TIME] increment current_time 2020-02-06 09:00:00 -> 2020-02-06 11:00:00
2026-07-14 12:10:27 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 12:10:27 INFO ---------->>> Running process_satellite_data()
2026-07-14 12:10:27 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12003.nc
2026-07-14 12:10:27 INFO ---------->>> Running run_obs_converter()
2026-07-14 12:10:27 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_39045_153072.out
2026-07-14 12:10:27 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_39045_153072.out
2026-07-14 12:10:27 INFO ---------->>> Running DART
2026-07-14 12:10:27 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-14 12:10:27 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/chim_ENS1_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/chim_ENS1_2020020611_1_out_toDART.nc
2026-07-14 12:10:27 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/chim_ENS2_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/chim_ENS2_2020020611_1_out_toDART.nc
2026-07-14 12:10:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/chim_ENS3_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/chim_ENS3_2020020611_1_out_toDART.nc
2026-07-14 12:10:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/chim_ENS4_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/chim_ENS4_2020020611_1_out_toDART.nc
2026-07-14 12:10:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/chim_ENS5_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/chim_ENS5_2020020611_1_out_toDART.nc
2026-07-14 12:10:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/chim_ENS6_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/chim_ENS6_2020020611_1_out_toDART.nc
2026-07-14 12:10:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/chim_ENS7_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/chim_ENS7_2020020611_1_out_toDART.nc
2026-07-14 12:10:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/chim_ENS8_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/chim_ENS8_2020020611_1_out_toDART.nc
2026-07-14 12:10:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/chim_ENS9_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/chim_ENS9_2020020611_1_out_toDART.nc
2026-07-14 12:10:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/chim_ENS10_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/chim_ENS10_2020020611_1_out_toDART.nc
2026-07-14 12:10:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/chim_ENS11_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/chim_ENS11_2020020611_1_out_toDART.nc
2026-07-14 12:10:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/chim_ENS12_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/chim_ENS12_2020020611_1_out_toDART.nc
2026-07-14 12:10:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/chim_ENS13_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/chim_ENS13_2020020611_1_out_toDART.nc
2026-07-14 12:10:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/chim_ENS14_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/chim_ENS14_2020020611_1_out_toDART.nc
2026-07-14 12:10:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/chim_ENS15_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/chim_ENS15_2020020611_1_out_toDART.nc
2026-07-14 12:10:32 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-14 12:10:32 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-14 12:10:32 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-14 12:10:32 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-14 12:10:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-14 12:10:32 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-14 12:10:46 INFO Found: []
2026-07-14 12:10:46 INFO No job id returned by command ./run_filter.bsh
2026-07-14 12:10:46 INFO No monitoring will be performed
2026-07-14 12:10:46 INFO Moving DART output files to analysis and preassim directories for date 2020020611 if present ...
2026-07-14 12:10:46 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 12:10:46 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 12:10:46 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 12:10:46 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 12:10:46 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 12:10:46 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 12:10:46 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 12:10:46 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 12:10:47 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 12:10:47 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 12:10:47 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 12:10:47 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 12:10:47 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 12:10:47 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 12:10:47 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 12:10:47 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 12:10:47 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 12:10:47 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 12:10:47 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 12:10:47 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 12:10:47 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 12:10:47 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 12:10:47 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 12:10:47 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 12:10:47 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 12:10:47 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 12:10:47 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 12:10:47 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 12:10:47 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 12:10:47 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 12:10:47 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 12:10:47 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 12:10:47 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 12:10:47 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 12:10:47 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-14 12:10:47 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-14 12:10:47 INFO run_dart() is DONE.
2026-07-14 12:10:47 INFO ---------->>> Running update_pollutant_in_end()
2026-07-14 12:10:47 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc
2026-07-14 12:10:53 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 12:10:53 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc
2026-07-14 12:10:58 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 12:10:58 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc
--- Logging error ---
Traceback (most recent call last):
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/logging/__init__.py", line 1164, in emit
    self.flush()
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/logging/__init__.py", line 1144, in flush
    self.stream.flush()
BrokenPipeError: [Errno 108] Cannot send after transport endpoint shutdown
Call stack:
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/main.py", line 103, in <module>
    pipeline.run_pipeline()
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/pipelines/base_pipeline.py", line 176, in run_pipeline
    self.after_assimilation()
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/pipelines/chimere2023/pipeline.py", line 607, in after_assimilation
    logger.info(f"Computing differences between posterior vs. original CHIMERE outputs ...")
Message: 'Computing differences between posterior vs. original CHIMERE outputs ...'
Arguments: ()
2026-07-14 12:13:11 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 12:13:50 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc
2026-07-14 12:13:55 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 12:13:55 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc
2026-07-14 12:14:00 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 12:14:01 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc
2026-07-14 12:14:06 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 12:14:06 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc
2026-07-14 12:14:11 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 12:14:11 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc
2026-07-14 12:14:16 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 12:14:17 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc
2026-07-14 12:14:22 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 12:14:22 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc
2026-07-14 12:14:27 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 12:14:27 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc
2026-07-14 12:14:32 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 12:14:33 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc
2026-07-14 12:14:38 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 12:14:38 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc
2026-07-14 12:14:43 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 12:14:43 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc
2026-07-14 12:14:48 INFO Computing differences between posterior vs. original CHIMERE outputs ...
/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/site-packages/xarray/backends/api.py:696: RuntimeWarning: 'netcdf4' fails while guessing
  engine = plugins.guess_engine(filename_or_obj)
Traceback (most recent call last):
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/site-packages/xarray/backends/file_manager.py", line 211, in _acquire_with_cache_info
    file = self._cache[self._key]
           ~~~~~~~~~~~^^^^^^^^^^^
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/site-packages/xarray/backends/lru_cache.py", line 56, in __getitem__
    value = self._cache[key]
            ~~~~~~~~~~~^^^^^
KeyError: [<function _open_scipy_netcdf at 0x2adebaaa5300>, ('/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/posteriors/2020020611/chim_ENS15_2020020611_1_out_fromDART.nc',), 'r', (('mmap', None), ('version', 2)), '446756b3-c4ce-4014-9ff1-17e49f09b56d']

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/main.py", line 103, in <module>
    pipeline.run_pipeline()
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/pipelines/base_pipeline.py", line 176, in run_pipeline
    self.after_assimilation()
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/pipelines/chimere2023/pipeline.py", line 596, in after_assimilation
    update_pollutant_in_end(dart_file=self.paths.dart_filter_output_list_file(mem, date_ymdH), 
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/orchestrator_utils.py", line 1302, in update_pollutant_in_end
    with xr.open_dataset(dart_file, decode_timedelta=True) as dart_ds, xr.open_dataset(end_file, decode_timedelta=True) as end_ds, xr.open_dataset(out_file, decode_timedelta=True) as out_ds, xr.open_dataset(emis_file, decode_timedelta=True) as emis_ds,  xr.open_dataset(dart_in_file, decode_timedelta=True) as dart_in_ds:
         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/site-packages/xarray/backends/api.py", line 715, in open_dataset
    backend_ds = backend.open_dataset(
                 ^^^^^^^^^^^^^^^^^^^^^
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/site-packages/xarray/backends/scipy_.py", line 333, in open_dataset
    ds = store_entrypoint.open_dataset(
         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/site-packages/xarray/backends/store.py", line 46, in open_dataset
    vars, attrs = filename_or_obj.load()
                  ^^^^^^^^^^^^^^^^^^^^^^
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/site-packages/xarray/backends/common.py", line 312, in load
    (_decode_variable_name(k), v) for k, v in self.get_variables().items()
                                              ^^^^^^^^^^^^^^^^^^^^
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/site-packages/xarray/backends/scipy_.py", line 199, in get_variables
    (k, self.open_store_variable(k, v)) for k, v in self.ds.variables.items()
                                                    ^^^^^^^
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/site-packages/xarray/backends/scipy_.py", line 188, in ds
    return self._manager.acquire()
           ^^^^^^^^^^^^^^^^^^^^^^^
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/site-packages/xarray/backends/file_manager.py", line 193, in acquire
    file, _ = self._acquire_with_cache_info(needs_lock)
              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/site-packages/xarray/backends/file_manager.py", line 217, in _acquire_with_cache_info
    file = self._opener(*self._args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/site-packages/xarray/backends/scipy_.py", line 127, in _open_scipy_netcdf
    return scipy.io.netcdf_file(filename, mode=mode, mmap=mmap, version=version)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/site-packages/scipy/io/_netcdf.py", line 283, in __init__
    self._read()
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/site-packages/scipy/io/_netcdf.py", line 615, in _read
    self._read_var_array()
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/site-packages/scipy/io/_netcdf.py", line 699, in _read_var_array
    data.shape = shape
    ^^^^^^^^^^
ValueError: cannot reshape array of size 90608 into shape (1,20,85,101)
+ exit 0
