+ SCRIPT_PID=2311094
+ /bin/bash -x /tmp/tmp.6kTp0fT3iK
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
2026-07-15 14:15:51 INFO 
███    ███ ██ ███    ███ ███████ ███████ ██
████  ████ ██ ████  ████ ██      ██      ██
██ ████ ██ ██ ██ ████ ██ █████   ███████ ██
██  ██  ██ ██ ██  ██  ██ ██           ██ ██
██      ██ ██ ██      ██ ███████ ███████ ██



2026-07-15 14:15:51 INFO [PIPELINE] =======================================
2026-07-15 14:15:51 INFO [PIPELINE] Starting chimere–DART orchestrator
2026-07-15 14:15:51 INFO [PIPELINE] Config file: config/config_irene_IM_cp2.yaml
2026-07-15 14:15:51 INFO [PIPELINE] Run dir: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart
2026-07-15 14:15:51 INFO [PIPELINE] Log file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/mimesi_orchestrator_logs/chimere_DART_20260715_141551.log
2026-07-15 14:15:51 INFO [PIPELINE] =======================================
2026-07-15 14:15:51 INFO Running assimilation with model_type=ModelType.CHIMERE
2026-07-15 14:15:51 INFO Using scheduler=Scheduler.SLURM, queue=rome
2026-07-15 14:15:51 INFO [STEP] ---- TIME LOOP START ----
2026-07-15 14:15:51 INFO [TIME] step_start current_time=2020-02-06 00:00:00 simulated_time=None dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 14:15:51 INFO [TIME] window start=2020-02-06 00:00:00 end=2020-02-06 01:00:00 run_hours=1 has_assimilation=False
2026-07-15 14:15:51 INFO Creating directories and links for ENS1 to run chimere's parallel part
2026-07-15 14:15:51 INFO Copying EMIS ...
2026-07-15 14:15:52 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens01.nc
2026-07-15 14:15:52 INFO Linking first END ...
2026-07-15 14:15:52 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020500_24_ENS1.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-15 14:15:52 INFO >> Checking links...
2026-07-15 14:15:52 INFO >> All links are good for ENS1  ...
2026-07-15 14:15:52 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:15:59 INFO Hourly dataset computed and listing created
2026-07-15 14:16:03 INFO Hourly dataset computed
2026-07-15 14:16:03 INFO Creating directories and links for ENS2 to run chimere's parallel part
2026-07-15 14:16:03 INFO Copying EMIS ...
2026-07-15 14:16:04 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens02.nc
2026-07-15 14:16:04 INFO Linking first END ...
2026-07-15 14:16:04 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020500_24_ENS2.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-15 14:16:04 INFO >> Checking links...
2026-07-15 14:16:04 INFO >> All links are good for ENS2  ...
2026-07-15 14:16:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:16:05 INFO Hourly dataset computed and listing created
2026-07-15 14:16:06 INFO Hourly dataset computed
2026-07-15 14:16:06 INFO Creating directories and links for ENS3 to run chimere's parallel part
2026-07-15 14:16:06 INFO Copying EMIS ...
2026-07-15 14:16:07 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens03.nc
2026-07-15 14:16:07 INFO Linking first END ...
2026-07-15 14:16:07 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020500_24_ENS3.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-15 14:16:07 INFO >> Checking links...
2026-07-15 14:16:07 INFO >> All links are good for ENS3  ...
2026-07-15 14:16:07 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:16:08 INFO Hourly dataset computed and listing created
2026-07-15 14:16:09 INFO Hourly dataset computed
2026-07-15 14:16:10 INFO Creating directories and links for ENS4 to run chimere's parallel part
2026-07-15 14:16:10 INFO Copying EMIS ...
2026-07-15 14:16:10 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens04.nc
2026-07-15 14:16:10 INFO Linking first END ...
2026-07-15 14:16:10 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020500_24_ENS4.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-15 14:16:10 INFO >> Checking links...
2026-07-15 14:16:10 INFO >> All links are good for ENS4  ...
2026-07-15 14:16:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:16:11 INFO Hourly dataset computed and listing created
2026-07-15 14:16:13 INFO Hourly dataset computed
2026-07-15 14:16:13 INFO Creating directories and links for ENS5 to run chimere's parallel part
2026-07-15 14:16:13 INFO Copying EMIS ...
2026-07-15 14:16:13 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens05.nc
2026-07-15 14:16:13 INFO Linking first END ...
2026-07-15 14:16:13 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020500_24_ENS5.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-15 14:16:13 INFO >> Checking links...
2026-07-15 14:16:13 INFO >> All links are good for ENS5  ...
2026-07-15 14:16:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:16:14 INFO Hourly dataset computed and listing created
2026-07-15 14:16:16 INFO Hourly dataset computed
2026-07-15 14:16:16 INFO Creating directories and links for ENS6 to run chimere's parallel part
2026-07-15 14:16:16 INFO Copying EMIS ...
2026-07-15 14:16:17 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens06.nc
2026-07-15 14:16:17 INFO Linking first END ...
2026-07-15 14:16:17 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020500_24_ENS6.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-15 14:16:17 INFO >> Checking links...
2026-07-15 14:16:17 INFO >> All links are good for ENS6  ...
2026-07-15 14:16:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:16:18 INFO Hourly dataset computed and listing created
2026-07-15 14:16:19 INFO Hourly dataset computed
2026-07-15 14:16:19 INFO Creating directories and links for ENS7 to run chimere's parallel part
2026-07-15 14:16:19 INFO Copying EMIS ...
2026-07-15 14:16:20 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens07.nc
2026-07-15 14:16:20 INFO Linking first END ...
2026-07-15 14:16:20 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020500_24_ENS7.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-15 14:16:20 INFO >> Checking links...
2026-07-15 14:16:20 INFO >> All links are good for ENS7  ...
2026-07-15 14:16:20 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:16:21 INFO Hourly dataset computed and listing created
2026-07-15 14:16:22 INFO Hourly dataset computed
2026-07-15 14:16:22 INFO Creating directories and links for ENS8 to run chimere's parallel part
2026-07-15 14:16:22 INFO Copying EMIS ...
2026-07-15 14:16:23 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens08.nc
2026-07-15 14:16:23 INFO Linking first END ...
2026-07-15 14:16:23 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020500_24_ENS8.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-15 14:16:23 INFO >> Checking links...
2026-07-15 14:16:23 INFO >> All links are good for ENS8  ...
2026-07-15 14:16:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:16:24 INFO Hourly dataset computed and listing created
2026-07-15 14:16:26 INFO Hourly dataset computed
2026-07-15 14:16:26 INFO Creating directories and links for ENS9 to run chimere's parallel part
2026-07-15 14:16:26 INFO Copying EMIS ...
2026-07-15 14:16:26 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens09.nc
2026-07-15 14:16:26 INFO Linking first END ...
2026-07-15 14:16:26 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020500_24_ENS9.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-15 14:16:26 INFO >> Checking links...
2026-07-15 14:16:26 INFO >> All links are good for ENS9  ...
2026-07-15 14:16:26 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:16:27 INFO Hourly dataset computed and listing created
2026-07-15 14:16:29 INFO Hourly dataset computed
2026-07-15 14:16:29 INFO Creating directories and links for ENS10 to run chimere's parallel part
2026-07-15 14:16:29 INFO Copying EMIS ...
2026-07-15 14:16:29 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens10.nc
2026-07-15 14:16:29 INFO Linking first END ...
2026-07-15 14:16:29 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020500_24_ENS10.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-15 14:16:29 INFO >> Checking links...
2026-07-15 14:16:30 INFO >> All links are good for ENS10  ...
2026-07-15 14:16:30 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:16:31 INFO Hourly dataset computed and listing created
2026-07-15 14:16:32 INFO Hourly dataset computed
2026-07-15 14:16:32 INFO Creating directories and links for ENS11 to run chimere's parallel part
2026-07-15 14:16:32 INFO Copying EMIS ...
2026-07-15 14:16:33 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens11.nc
2026-07-15 14:16:33 INFO Linking first END ...
2026-07-15 14:16:33 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020500_24_ENS11.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-15 14:16:33 INFO >> Checking links...
2026-07-15 14:16:33 INFO >> All links are good for ENS11  ...
2026-07-15 14:16:33 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:16:34 INFO Hourly dataset computed and listing created
2026-07-15 14:16:35 INFO Hourly dataset computed
2026-07-15 14:16:35 INFO Creating directories and links for ENS12 to run chimere's parallel part
2026-07-15 14:16:35 INFO Copying EMIS ...
2026-07-15 14:16:36 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens12.nc
2026-07-15 14:16:36 INFO Linking first END ...
2026-07-15 14:16:36 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020500_24_ENS12.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-15 14:16:36 INFO >> Checking links...
2026-07-15 14:16:36 INFO >> All links are good for ENS12  ...
2026-07-15 14:16:36 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:16:37 INFO Hourly dataset computed and listing created
2026-07-15 14:16:39 INFO Hourly dataset computed
2026-07-15 14:16:39 INFO Creating directories and links for ENS13 to run chimere's parallel part
2026-07-15 14:16:39 INFO Copying EMIS ...
2026-07-15 14:16:39 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens13.nc
2026-07-15 14:16:39 INFO Linking first END ...
2026-07-15 14:16:39 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020500_24_ENS13.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-15 14:16:39 INFO >> Checking links...
2026-07-15 14:16:39 INFO >> All links are good for ENS13  ...
2026-07-15 14:16:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:16:41 INFO Hourly dataset computed and listing created
2026-07-15 14:16:43 INFO Hourly dataset computed
2026-07-15 14:16:43 INFO Creating directories and links for ENS14 to run chimere's parallel part
2026-07-15 14:16:43 INFO Copying EMIS ...
2026-07-15 14:16:43 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens14.nc
2026-07-15 14:16:43 INFO Linking first END ...
2026-07-15 14:16:43 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020500_24_ENS14.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-15 14:16:43 INFO >> Checking links...
2026-07-15 14:16:43 INFO >> All links are good for ENS14  ...
2026-07-15 14:16:43 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:16:44 INFO Hourly dataset computed and listing created
2026-07-15 14:16:45 INFO Hourly dataset computed
2026-07-15 14:16:45 INFO Creating directories and links for ENS15 to run chimere's parallel part
2026-07-15 14:16:45 INFO Copying EMIS ...
2026-07-15 14:16:46 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens15.nc
2026-07-15 14:16:46 INFO Linking first END ...
2026-07-15 14:16:46 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020500_24_ENS15.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-15 14:16:46 INFO >> Checking links...
2026-07-15 14:16:46 INFO >> All links are good for ENS15  ...
2026-07-15 14:16:46 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:16:47 INFO Hourly dataset computed and listing created
2026-07-15 14:16:47 INFO Hourly dataset computed
2026-07-15 14:16:47 INFO ---------->>> Running CHIMERE model from 2020-02-06 00:00:00 to 2020-02-06 01:00:00
2026-07-15 14:16:47 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:16:47 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1
2026-07-15 14:16:47 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020500_24_ENS1.nc
2026-07-15 14:16:47 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-15 14:16:47 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:16:47 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-15 14:16:47 INFO Queuing job for member 1...
2026-07-15 14:16:47 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:16:47 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-15 14:16:48 INFO Found: ['5163918']
2026-07-15 14:16:53 INFO [TGCC-IRENE] Submitted job with ID:['5163918']
2026-07-15 14:16:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:16:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2
2026-07-15 14:16:53 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020500_24_ENS2.nc
2026-07-15 14:16:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-15 14:16:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:16:53 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-15 14:16:53 INFO Queuing job for member 2...
2026-07-15 14:16:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:16:53 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-15 14:16:54 INFO Found: ['5163919']
2026-07-15 14:16:59 INFO [TGCC-IRENE] Submitted job with ID:['5163919']
2026-07-15 14:16:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:16:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3
2026-07-15 14:16:59 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020500_24_ENS3.nc
2026-07-15 14:16:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-15 14:16:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:16:59 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-15 14:16:59 INFO Queuing job for member 3...
2026-07-15 14:16:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:16:59 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-15 14:17:00 INFO Found: ['5163920']
2026-07-15 14:17:05 INFO [TGCC-IRENE] Submitted job with ID:['5163920']
2026-07-15 14:17:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:17:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4
2026-07-15 14:17:05 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020500_24_ENS4.nc
2026-07-15 14:17:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-15 14:17:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:17:05 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-15 14:17:05 INFO Queuing job for member 4...
2026-07-15 14:17:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:17:05 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-15 14:17:05 INFO Found: ['5163922']
2026-07-15 14:17:10 INFO [TGCC-IRENE] Submitted job with ID:['5163922']
2026-07-15 14:17:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:17:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5
2026-07-15 14:17:10 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020500_24_ENS5.nc
2026-07-15 14:17:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-15 14:17:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:17:10 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-15 14:17:10 INFO Queuing job for member 5...
2026-07-15 14:17:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:17:10 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-15 14:17:11 INFO Found: ['5163923']
2026-07-15 14:17:16 INFO [TGCC-IRENE] Submitted job with ID:['5163923']
2026-07-15 14:17:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:17:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6
2026-07-15 14:17:16 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020500_24_ENS6.nc
2026-07-15 14:17:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-15 14:17:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:17:16 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-15 14:17:16 INFO Queuing job for member 6...
2026-07-15 14:17:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:17:16 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-15 14:17:17 INFO Found: ['5163924']
2026-07-15 14:17:22 INFO [TGCC-IRENE] Submitted job with ID:['5163924']
2026-07-15 14:17:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:17:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7
2026-07-15 14:17:22 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020500_24_ENS7.nc
2026-07-15 14:17:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-15 14:17:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:17:22 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-15 14:17:22 INFO Queuing job for member 7...
2026-07-15 14:17:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:17:22 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-15 14:17:24 INFO Found: ['5163925']
2026-07-15 14:17:29 INFO [TGCC-IRENE] Submitted job with ID:['5163925']
2026-07-15 14:17:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:17:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8
2026-07-15 14:17:29 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020500_24_ENS8.nc
2026-07-15 14:17:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-15 14:17:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:17:29 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-15 14:17:29 INFO Queuing job for member 8...
2026-07-15 14:17:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:17:29 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-15 14:17:31 INFO Found: ['5163926']
2026-07-15 14:17:36 INFO [TGCC-IRENE] Submitted job with ID:['5163926']
2026-07-15 14:17:36 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:17:36 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9
2026-07-15 14:17:36 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020500_24_ENS9.nc
2026-07-15 14:17:36 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-15 14:17:36 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:17:36 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-15 14:17:36 INFO Queuing job for member 9...
2026-07-15 14:17:36 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:17:36 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-15 14:17:37 INFO Found: ['5163927']
2026-07-15 14:17:42 INFO [TGCC-IRENE] Submitted job with ID:['5163927']
2026-07-15 14:17:42 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:17:42 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10
2026-07-15 14:17:42 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020500_24_ENS10.nc
2026-07-15 14:17:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-15 14:17:42 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:17:42 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-15 14:17:42 INFO Queuing job for member 10...
2026-07-15 14:17:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:17:42 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-15 14:17:43 INFO Found: ['5163929']
2026-07-15 14:17:48 INFO [TGCC-IRENE] Submitted job with ID:['5163929']
2026-07-15 14:17:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:17:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11
2026-07-15 14:17:48 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020500_24_ENS11.nc
2026-07-15 14:17:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-15 14:17:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:17:48 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-15 14:17:48 INFO Queuing job for member 11...
2026-07-15 14:17:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:17:48 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-15 14:17:49 INFO Found: ['5163931']
2026-07-15 14:17:54 INFO [TGCC-IRENE] Submitted job with ID:['5163931']
2026-07-15 14:17:54 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:17:54 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12
2026-07-15 14:17:54 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020500_24_ENS12.nc
2026-07-15 14:17:54 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-15 14:17:54 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:17:54 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-15 14:17:54 INFO Queuing job for member 12...
2026-07-15 14:17:54 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:17:54 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-15 14:17:54 INFO Found: ['5163932']
2026-07-15 14:17:59 INFO [TGCC-IRENE] Submitted job with ID:['5163932']
2026-07-15 14:17:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:17:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13
2026-07-15 14:17:59 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020500_24_ENS13.nc
2026-07-15 14:17:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-15 14:17:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:17:59 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-15 14:17:59 INFO Queuing job for member 13...
2026-07-15 14:17:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:17:59 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-15 14:18:00 INFO Found: ['5163933']
2026-07-15 14:18:05 INFO [TGCC-IRENE] Submitted job with ID:['5163933']
2026-07-15 14:18:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:18:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14
2026-07-15 14:18:05 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020500_24_ENS14.nc
2026-07-15 14:18:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-15 14:18:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:18:05 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-15 14:18:05 INFO Queuing job for member 14...
2026-07-15 14:18:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:18:05 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-15 14:18:06 INFO Found: ['5163935']
2026-07-15 14:18:11 INFO [TGCC-IRENE] Submitted job with ID:['5163935']
2026-07-15 14:18:11 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:18:11 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15
2026-07-15 14:18:11 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020500_24_ENS15.nc
2026-07-15 14:18:11 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-15 14:18:11 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:18:11 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-15 14:18:11 INFO Queuing job for member 15...
2026-07-15 14:18:11 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:18:11 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-15 14:18:12 INFO Found: ['5163937']
2026-07-15 14:18:17 INFO [TGCC-IRENE] Submitted job with ID:['5163937']
2026-07-15 14:18:17 INFO Checking job status ...
2026-07-15 14:18:18 INFO None 5163918: status RUNNING/PENDING
2026-07-15 14:18:18 INFO None 5163919: status RUNNING/PENDING
2026-07-15 14:18:18 INFO None 5163920: status RUNNING/PENDING
2026-07-15 14:18:18 INFO None 5163922: status RUNNING/PENDING
2026-07-15 14:18:18 INFO None 5163923: status RUNNING/PENDING
2026-07-15 14:18:18 INFO None 5163924: status RUNNING/PENDING
2026-07-15 14:18:18 INFO None 5163925: status RUNNING/PENDING
2026-07-15 14:18:18 INFO None 5163926: status RUNNING/PENDING
2026-07-15 14:18:18 INFO None 5163927: status RUNNING/PENDING
2026-07-15 14:18:18 INFO None 5163929: status RUNNING/PENDING
2026-07-15 14:18:18 INFO None 5163931: status RUNNING/PENDING
2026-07-15 14:18:18 INFO None 5163932: status RUNNING/PENDING
2026-07-15 14:18:18 INFO None 5163933: status RUNNING/PENDING
2026-07-15 14:18:18 INFO None 5163935: status RUNNING/PENDING
2026-07-15 14:18:18 INFO None 5163937: status RUNNING/PENDING
2026-07-15 14:18:18 INFO Jobs still running: ['5163918', '5163919', '5163920', '5163922', '5163923', '5163924', '5163925', '5163926', '5163927', '5163929', '5163931', '5163932', '5163933', '5163935', '5163937']. Waiting...
2026-07-15 14:18:33 INFO None 5163918: status RUNNING/PENDING
2026-07-15 14:18:33 INFO None 5163919: status RUNNING/PENDING
2026-07-15 14:18:34 INFO None 5163920: status RUNNING/PENDING
2026-07-15 14:18:34 INFO None 5163922: status RUNNING/PENDING
2026-07-15 14:18:34 INFO None 5163923: status RUNNING/PENDING
2026-07-15 14:18:34 INFO None 5163924: status RUNNING/PENDING
2026-07-15 14:18:34 INFO None 5163925: status RUNNING/PENDING
2026-07-15 14:18:34 INFO None 5163926: status RUNNING/PENDING
2026-07-15 14:18:34 INFO None 5163927: status RUNNING/PENDING
2026-07-15 14:18:34 INFO None 5163929: status RUNNING/PENDING
2026-07-15 14:18:34 INFO None 5163931: status RUNNING/PENDING
2026-07-15 14:18:34 INFO None 5163932: status RUNNING/PENDING
2026-07-15 14:18:34 INFO None 5163933: status RUNNING/PENDING
2026-07-15 14:18:34 INFO None 5163935: status RUNNING/PENDING
2026-07-15 14:18:34 INFO None 5163937: status RUNNING/PENDING
2026-07-15 14:18:34 INFO Jobs still running: ['5163918', '5163919', '5163920', '5163922', '5163923', '5163924', '5163925', '5163926', '5163927', '5163929', '5163931', '5163932', '5163933', '5163935', '5163937']. Waiting...
2026-07-15 14:18:49 INFO None 5163918: status RUNNING/PENDING
2026-07-15 14:18:49 INFO None 5163919: status RUNNING/PENDING
2026-07-15 14:18:49 INFO None 5163920: status RUNNING/PENDING
2026-07-15 14:18:49 INFO None 5163922: status RUNNING/PENDING
2026-07-15 14:18:49 INFO None 5163923: status RUNNING/PENDING
2026-07-15 14:18:49 INFO None 5163924: status RUNNING/PENDING
2026-07-15 14:18:49 INFO None 5163925: status RUNNING/PENDING
2026-07-15 14:18:49 INFO None 5163926: status RUNNING/PENDING
2026-07-15 14:18:49 INFO None 5163927: status RUNNING/PENDING
2026-07-15 14:18:49 INFO None 5163929: status RUNNING/PENDING
2026-07-15 14:18:49 INFO None 5163931: status RUNNING/PENDING
2026-07-15 14:18:49 INFO None 5163932: status RUNNING/PENDING
2026-07-15 14:18:49 INFO None 5163933: status RUNNING/PENDING
2026-07-15 14:18:49 INFO None 5163935: status RUNNING/PENDING
2026-07-15 14:18:49 INFO None 5163937: status RUNNING/PENDING
2026-07-15 14:18:49 INFO Jobs still running: ['5163918', '5163919', '5163920', '5163922', '5163923', '5163924', '5163925', '5163926', '5163927', '5163929', '5163931', '5163932', '5163933', '5163935', '5163937']. Waiting...
2026-07-15 14:19:04 INFO None 5163918: status RUNNING/PENDING
2026-07-15 14:19:04 INFO None 5163919: status RUNNING/PENDING
2026-07-15 14:19:04 INFO None 5163920: status RUNNING/PENDING
2026-07-15 14:19:04 INFO None 5163922: status RUNNING/PENDING
2026-07-15 14:19:04 INFO None 5163923: status RUNNING/PENDING
2026-07-15 14:19:04 INFO None 5163924: status RUNNING/PENDING
2026-07-15 14:19:04 INFO None 5163925: status RUNNING/PENDING
2026-07-15 14:19:04 INFO None 5163926: status RUNNING/PENDING
2026-07-15 14:19:04 INFO None 5163927: status RUNNING/PENDING
2026-07-15 14:19:04 INFO None 5163929: status RUNNING/PENDING
2026-07-15 14:19:04 INFO None 5163931: status RUNNING/PENDING
2026-07-15 14:19:04 INFO None 5163932: status RUNNING/PENDING
2026-07-15 14:19:04 INFO None 5163933: status RUNNING/PENDING
2026-07-15 14:19:04 INFO None 5163935: status RUNNING/PENDING
2026-07-15 14:19:04 INFO None 5163937: status RUNNING/PENDING
2026-07-15 14:19:04 INFO Jobs still running: ['5163918', '5163919', '5163920', '5163922', '5163923', '5163924', '5163925', '5163926', '5163927', '5163929', '5163931', '5163932', '5163933', '5163935', '5163937']. Waiting...
2026-07-15 14:19:20 INFO None 5163918: status RUNNING/PENDING
2026-07-15 14:19:20 INFO None 5163919: status RUNNING/PENDING
2026-07-15 14:19:20 INFO None 5163920: status RUNNING/PENDING
2026-07-15 14:19:20 INFO None 5163922: status RUNNING/PENDING
2026-07-15 14:19:20 INFO None 5163923: status RUNNING/PENDING
2026-07-15 14:19:20 INFO None 5163924: status RUNNING/PENDING
2026-07-15 14:19:20 INFO None 5163925: status RUNNING/PENDING
2026-07-15 14:19:20 INFO None 5163926: status RUNNING/PENDING
2026-07-15 14:19:20 INFO None 5163927: status RUNNING/PENDING
2026-07-15 14:19:20 INFO None 5163929: status RUNNING/PENDING
2026-07-15 14:19:20 INFO None 5163931: status RUNNING/PENDING
2026-07-15 14:19:20 INFO None 5163932: status RUNNING/PENDING
2026-07-15 14:19:20 INFO None 5163933: status RUNNING/PENDING
2026-07-15 14:19:20 INFO None 5163935: status RUNNING/PENDING
2026-07-15 14:19:20 INFO None 5163937: status RUNNING/PENDING
2026-07-15 14:19:20 INFO Jobs still running: ['5163918', '5163919', '5163920', '5163922', '5163923', '5163924', '5163925', '5163926', '5163927', '5163929', '5163931', '5163932', '5163933', '5163935', '5163937']. Waiting...
2026-07-15 14:19:35 INFO None 5163918: status RUNNING/PENDING
2026-07-15 14:19:35 INFO None 5163919: status RUNNING/PENDING
2026-07-15 14:19:35 INFO None 5163920: status RUNNING/PENDING
2026-07-15 14:19:35 INFO None 5163922: status RUNNING/PENDING
2026-07-15 14:19:35 INFO None 5163923: status FINISHED
2026-07-15 14:19:35 INFO None 5163924: status FINISHED
2026-07-15 14:19:36 INFO None 5163925: status RUNNING/PENDING
2026-07-15 14:19:36 INFO None 5163926: status RUNNING/PENDING
2026-07-15 14:19:36 INFO None 5163927: status RUNNING/PENDING
2026-07-15 14:19:36 INFO None 5163929: status RUNNING/PENDING
2026-07-15 14:19:36 INFO None 5163931: status RUNNING/PENDING
2026-07-15 14:19:36 INFO None 5163932: status RUNNING/PENDING
2026-07-15 14:19:36 INFO None 5163933: status RUNNING/PENDING
2026-07-15 14:19:36 INFO None 5163935: status RUNNING/PENDING
2026-07-15 14:19:36 INFO None 5163937: status RUNNING/PENDING
2026-07-15 14:19:36 INFO Jobs still running: ['5163918', '5163919', '5163920', '5163922', '5163925', '5163926', '5163927', '5163929', '5163931', '5163932', '5163933', '5163935', '5163937']. Waiting...
2026-07-15 14:19:51 INFO None 5163918: status RUNNING/PENDING
2026-07-15 14:19:51 INFO None 5163919: status RUNNING/PENDING
2026-07-15 14:19:51 INFO None 5163920: status RUNNING/PENDING
2026-07-15 14:19:51 INFO None 5163922: status RUNNING/PENDING
2026-07-15 14:19:51 INFO None 5163923: status FINISHED
2026-07-15 14:19:51 INFO None 5163924: status FINISHED
2026-07-15 14:19:51 INFO None 5163925: status RUNNING/PENDING
2026-07-15 14:19:51 INFO None 5163926: status RUNNING/PENDING
2026-07-15 14:19:51 INFO None 5163927: status RUNNING/PENDING
2026-07-15 14:19:51 INFO None 5163929: status RUNNING/PENDING
2026-07-15 14:19:51 INFO None 5163931: status RUNNING/PENDING
2026-07-15 14:19:51 INFO None 5163932: status RUNNING/PENDING
2026-07-15 14:19:51 INFO None 5163933: status RUNNING/PENDING
2026-07-15 14:19:51 INFO None 5163935: status RUNNING/PENDING
2026-07-15 14:19:51 INFO None 5163937: status RUNNING/PENDING
2026-07-15 14:19:51 INFO Jobs still running: ['5163918', '5163919', '5163920', '5163922', '5163925', '5163926', '5163927', '5163929', '5163931', '5163932', '5163933', '5163935', '5163937']. Waiting...
2026-07-15 14:20:06 INFO None 5163918: status RUNNING/PENDING
2026-07-15 14:20:06 INFO None 5163919: status RUNNING/PENDING
2026-07-15 14:20:06 INFO None 5163920: status RUNNING/PENDING
2026-07-15 14:20:06 INFO None 5163922: status RUNNING/PENDING
2026-07-15 14:20:06 INFO None 5163923: status FINISHED
2026-07-15 14:20:06 INFO None 5163924: status FINISHED
2026-07-15 14:20:06 INFO None 5163925: status RUNNING/PENDING
2026-07-15 14:20:06 INFO None 5163926: status RUNNING/PENDING
2026-07-15 14:20:06 INFO None 5163927: status RUNNING/PENDING
2026-07-15 14:20:06 INFO None 5163929: status RUNNING/PENDING
2026-07-15 14:20:06 INFO None 5163931: status RUNNING/PENDING
2026-07-15 14:20:06 INFO None 5163932: status RUNNING/PENDING
2026-07-15 14:20:06 INFO None 5163933: status RUNNING/PENDING
2026-07-15 14:20:06 INFO None 5163935: status RUNNING/PENDING
2026-07-15 14:20:06 INFO None 5163937: status RUNNING/PENDING
2026-07-15 14:20:06 INFO Jobs still running: ['5163918', '5163919', '5163920', '5163922', '5163925', '5163926', '5163927', '5163929', '5163931', '5163932', '5163933', '5163935', '5163937']. Waiting...
2026-07-15 14:20:23 INFO None 5163918: status RUNNING/PENDING
2026-07-15 14:20:23 INFO None 5163919: status RUNNING/PENDING
2026-07-15 14:20:23 INFO None 5163920: status RUNNING/PENDING
2026-07-15 14:20:23 INFO None 5163922: status RUNNING/PENDING
2026-07-15 14:20:23 INFO None 5163923: status FINISHED
2026-07-15 14:20:23 INFO None 5163924: status FINISHED
2026-07-15 14:20:23 INFO None 5163925: status RUNNING/PENDING
2026-07-15 14:20:23 INFO None 5163926: status RUNNING/PENDING
2026-07-15 14:20:23 INFO None 5163927: status RUNNING/PENDING
2026-07-15 14:20:23 INFO None 5163929: status RUNNING/PENDING
2026-07-15 14:20:23 INFO None 5163931: status RUNNING/PENDING
2026-07-15 14:20:23 INFO None 5163932: status RUNNING/PENDING
2026-07-15 14:20:23 INFO None 5163933: status RUNNING/PENDING
2026-07-15 14:20:23 INFO None 5163935: status RUNNING/PENDING
2026-07-15 14:20:23 INFO None 5163937: status RUNNING/PENDING
2026-07-15 14:20:23 INFO Jobs still running: ['5163918', '5163919', '5163920', '5163922', '5163925', '5163926', '5163927', '5163929', '5163931', '5163932', '5163933', '5163935', '5163937']. Waiting...
2026-07-15 14:20:38 INFO None 5163918: status RUNNING/PENDING
2026-07-15 14:20:38 INFO None 5163919: status RUNNING/PENDING
2026-07-15 14:20:38 INFO None 5163920: status RUNNING/PENDING
2026-07-15 14:20:38 INFO None 5163922: status RUNNING/PENDING
2026-07-15 14:20:38 INFO None 5163923: status FINISHED
2026-07-15 14:20:38 INFO None 5163924: status FINISHED
2026-07-15 14:20:38 INFO None 5163925: status RUNNING/PENDING
2026-07-15 14:20:38 INFO None 5163926: status RUNNING/PENDING
2026-07-15 14:20:38 INFO None 5163927: status RUNNING/PENDING
2026-07-15 14:20:38 INFO None 5163929: status RUNNING/PENDING
2026-07-15 14:20:38 INFO None 5163931: status RUNNING/PENDING
2026-07-15 14:20:38 INFO None 5163932: status RUNNING/PENDING
2026-07-15 14:20:38 INFO None 5163933: status RUNNING/PENDING
2026-07-15 14:20:38 INFO None 5163935: status RUNNING/PENDING
2026-07-15 14:20:38 INFO None 5163937: status RUNNING/PENDING
2026-07-15 14:20:38 INFO Jobs still running: ['5163918', '5163919', '5163920', '5163922', '5163925', '5163926', '5163927', '5163929', '5163931', '5163932', '5163933', '5163935', '5163937']. Waiting...
2026-07-15 14:20:53 INFO None 5163918: status RUNNING/PENDING
2026-07-15 14:20:53 INFO None 5163919: status RUNNING/PENDING
2026-07-15 14:20:53 INFO None 5163920: status RUNNING/PENDING
2026-07-15 14:20:53 INFO None 5163922: status FINISHED
2026-07-15 14:20:53 INFO None 5163923: status FINISHED
2026-07-15 14:20:53 INFO None 5163924: status FINISHED
2026-07-15 14:20:53 INFO None 5163925: status RUNNING/PENDING
2026-07-15 14:20:53 INFO None 5163926: status RUNNING/PENDING
2026-07-15 14:20:53 INFO None 5163927: status RUNNING/PENDING
2026-07-15 14:20:53 INFO None 5163929: status RUNNING/PENDING
2026-07-15 14:20:53 INFO None 5163931: status RUNNING/PENDING
2026-07-15 14:20:54 INFO None 5163932: status RUNNING/PENDING
2026-07-15 14:20:54 INFO None 5163933: status RUNNING/PENDING
2026-07-15 14:20:54 INFO None 5163935: status RUNNING/PENDING
2026-07-15 14:20:54 INFO None 5163937: status RUNNING/PENDING
2026-07-15 14:20:54 INFO Jobs still running: ['5163918', '5163919', '5163920', '5163925', '5163926', '5163927', '5163929', '5163931', '5163932', '5163933', '5163935', '5163937']. Waiting...
2026-07-15 14:21:09 INFO None 5163918: status RUNNING/PENDING
2026-07-15 14:21:09 INFO None 5163919: status RUNNING/PENDING
2026-07-15 14:21:09 INFO None 5163920: status RUNNING/PENDING
2026-07-15 14:21:09 INFO None 5163922: status FINISHED
2026-07-15 14:21:09 INFO None 5163923: status FINISHED
2026-07-15 14:21:09 INFO None 5163924: status FINISHED
2026-07-15 14:21:09 INFO None 5163925: status RUNNING/PENDING
2026-07-15 14:21:09 INFO None 5163926: status RUNNING/PENDING
2026-07-15 14:21:09 INFO None 5163927: status FINISHED
2026-07-15 14:21:09 INFO None 5163929: status RUNNING/PENDING
2026-07-15 14:21:09 INFO None 5163931: status RUNNING/PENDING
2026-07-15 14:21:09 INFO None 5163932: status RUNNING/PENDING
2026-07-15 14:21:09 INFO None 5163933: status RUNNING/PENDING
2026-07-15 14:21:09 INFO None 5163935: status RUNNING/PENDING
2026-07-15 14:21:09 INFO None 5163937: status RUNNING/PENDING
2026-07-15 14:21:09 INFO Jobs still running: ['5163918', '5163919', '5163920', '5163925', '5163926', '5163929', '5163931', '5163932', '5163933', '5163935', '5163937']. Waiting...
2026-07-15 14:21:26 INFO None 5163918: status RUNNING/PENDING
2026-07-15 14:21:26 INFO None 5163919: status FINISHED
2026-07-15 14:21:26 INFO None 5163920: status RUNNING/PENDING
2026-07-15 14:21:26 INFO None 5163922: status FINISHED
2026-07-15 14:21:26 INFO None 5163923: status FINISHED
2026-07-15 14:21:26 INFO None 5163924: status FINISHED
2026-07-15 14:21:26 INFO None 5163925: status FINISHED
2026-07-15 14:21:26 INFO None 5163926: status RUNNING/PENDING
2026-07-15 14:21:26 INFO None 5163927: status FINISHED
2026-07-15 14:21:26 INFO None 5163929: status RUNNING/PENDING
2026-07-15 14:21:26 INFO None 5163931: status RUNNING/PENDING
2026-07-15 14:21:26 INFO None 5163932: status RUNNING/PENDING
2026-07-15 14:21:26 INFO None 5163933: status RUNNING/PENDING
2026-07-15 14:21:26 INFO None 5163935: status RUNNING/PENDING
2026-07-15 14:21:26 INFO None 5163937: status RUNNING/PENDING
2026-07-15 14:21:26 INFO Jobs still running: ['5163918', '5163920', '5163926', '5163929', '5163931', '5163932', '5163933', '5163935', '5163937']. Waiting...
2026-07-15 14:21:41 INFO None 5163918: status RUNNING/PENDING
2026-07-15 14:21:41 INFO None 5163919: status FINISHED
2026-07-15 14:21:41 INFO None 5163920: status RUNNING/PENDING
2026-07-15 14:21:41 INFO None 5163922: status FINISHED
2026-07-15 14:21:41 INFO None 5163923: status FINISHED
2026-07-15 14:21:41 INFO None 5163924: status FINISHED
2026-07-15 14:21:41 INFO None 5163925: status FINISHED
2026-07-15 14:21:41 INFO None 5163926: status RUNNING/PENDING
2026-07-15 14:21:41 INFO None 5163927: status FINISHED
2026-07-15 14:21:41 INFO None 5163929: status RUNNING/PENDING
2026-07-15 14:21:41 INFO None 5163931: status RUNNING/PENDING
2026-07-15 14:21:41 INFO None 5163932: status RUNNING/PENDING
2026-07-15 14:21:41 INFO None 5163933: status FINISHED
2026-07-15 14:21:41 INFO None 5163935: status FINISHED
2026-07-15 14:21:41 INFO None 5163937: status RUNNING/PENDING
2026-07-15 14:21:41 INFO Jobs still running: ['5163918', '5163920', '5163926', '5163929', '5163931', '5163932', '5163937']. Waiting...
2026-07-15 14:21:56 INFO None 5163918: status RUNNING/PENDING
2026-07-15 14:21:56 INFO None 5163919: status FINISHED
2026-07-15 14:21:56 INFO None 5163920: status RUNNING/PENDING
2026-07-15 14:21:57 INFO None 5163922: status FINISHED
2026-07-15 14:21:57 INFO None 5163923: status FINISHED
2026-07-15 14:21:57 INFO None 5163924: status FINISHED
2026-07-15 14:21:57 INFO None 5163925: status FINISHED
2026-07-15 14:21:57 INFO None 5163926: status FINISHED
2026-07-15 14:21:57 INFO None 5163927: status FINISHED
2026-07-15 14:21:57 INFO None 5163929: status RUNNING/PENDING
2026-07-15 14:21:57 INFO None 5163931: status RUNNING/PENDING
2026-07-15 14:21:57 INFO None 5163932: status RUNNING/PENDING
2026-07-15 14:21:57 INFO None 5163933: status FINISHED
2026-07-15 14:21:57 INFO None 5163935: status FINISHED
2026-07-15 14:21:57 INFO None 5163937: status RUNNING/PENDING
2026-07-15 14:21:57 INFO Jobs still running: ['5163918', '5163920', '5163929', '5163931', '5163932', '5163937']. Waiting...
2026-07-15 14:22:12 INFO None 5163918: status FINISHED
2026-07-15 14:22:12 INFO None 5163919: status FINISHED
2026-07-15 14:22:12 INFO None 5163920: status RUNNING/PENDING
2026-07-15 14:22:12 INFO None 5163922: status FINISHED
2026-07-15 14:22:12 INFO None 5163923: status FINISHED
2026-07-15 14:22:12 INFO None 5163924: status FINISHED
2026-07-15 14:22:12 INFO None 5163925: status FINISHED
2026-07-15 14:22:12 INFO None 5163926: status FINISHED
2026-07-15 14:22:12 INFO None 5163927: status FINISHED
2026-07-15 14:22:12 INFO None 5163929: status FINISHED
2026-07-15 14:22:12 INFO None 5163931: status FINISHED
2026-07-15 14:22:12 INFO None 5163932: status FINISHED
2026-07-15 14:22:12 INFO None 5163933: status FINISHED
2026-07-15 14:22:12 INFO None 5163935: status FINISHED
2026-07-15 14:22:12 INFO None 5163937: status RUNNING/PENDING
2026-07-15 14:22:12 INFO Jobs still running: ['5163920', '5163937']. Waiting...
2026-07-15 14:22:28 INFO None 5163918: status FINISHED
2026-07-15 14:22:28 INFO None 5163919: status FINISHED
2026-07-15 14:22:28 INFO None 5163920: status RUNNING/PENDING
2026-07-15 14:22:28 INFO None 5163922: status FINISHED
2026-07-15 14:22:29 INFO None 5163923: status FINISHED
2026-07-15 14:22:29 INFO None 5163924: status FINISHED
2026-07-15 14:22:29 INFO None 5163925: status FINISHED
2026-07-15 14:22:29 INFO None 5163926: status FINISHED
2026-07-15 14:22:29 INFO None 5163927: status FINISHED
2026-07-15 14:22:29 INFO None 5163929: status FINISHED
2026-07-15 14:22:29 INFO None 5163931: status FINISHED
2026-07-15 14:22:29 INFO None 5163932: status FINISHED
2026-07-15 14:22:29 INFO None 5163933: status FINISHED
2026-07-15 14:22:29 INFO None 5163935: status FINISHED
2026-07-15 14:22:29 INFO None 5163937: status RUNNING/PENDING
2026-07-15 14:22:29 INFO Jobs still running: ['5163920', '5163937']. Waiting...
2026-07-15 14:22:44 INFO None 5163918: status FINISHED
2026-07-15 14:22:44 INFO None 5163919: status FINISHED
2026-07-15 14:22:44 INFO None 5163920: status RUNNING/PENDING
2026-07-15 14:22:44 INFO None 5163922: status FINISHED
2026-07-15 14:22:44 INFO None 5163923: status FINISHED
2026-07-15 14:22:44 INFO None 5163924: status FINISHED
2026-07-15 14:22:44 INFO None 5163925: status FINISHED
2026-07-15 14:22:44 INFO None 5163926: status FINISHED
2026-07-15 14:22:44 INFO None 5163927: status FINISHED
2026-07-15 14:22:44 INFO None 5163929: status FINISHED
2026-07-15 14:22:44 INFO None 5163931: status FINISHED
2026-07-15 14:22:44 INFO None 5163932: status FINISHED
2026-07-15 14:22:44 INFO None 5163933: status FINISHED
2026-07-15 14:22:44 INFO None 5163935: status FINISHED
2026-07-15 14:22:44 INFO None 5163937: status RUNNING/PENDING
2026-07-15 14:22:44 INFO Jobs still running: ['5163920', '5163937']. Waiting...
2026-07-15 14:22:59 INFO None 5163918: status FINISHED
2026-07-15 14:22:59 INFO None 5163919: status FINISHED
2026-07-15 14:22:59 INFO None 5163920: status RUNNING/PENDING
2026-07-15 14:22:59 INFO None 5163922: status FINISHED
2026-07-15 14:22:59 INFO None 5163923: status FINISHED
2026-07-15 14:22:59 INFO None 5163924: status FINISHED
2026-07-15 14:22:59 INFO None 5163925: status FINISHED
2026-07-15 14:22:59 INFO None 5163926: status FINISHED
2026-07-15 14:22:59 INFO None 5163927: status FINISHED
2026-07-15 14:22:59 INFO None 5163929: status FINISHED
2026-07-15 14:22:59 INFO None 5163931: status FINISHED
2026-07-15 14:22:59 INFO None 5163932: status FINISHED
2026-07-15 14:22:59 INFO None 5163933: status FINISHED
2026-07-15 14:22:59 INFO None 5163935: status FINISHED
2026-07-15 14:22:59 INFO None 5163937: status FINISHED
2026-07-15 14:22:59 INFO Jobs still running: ['5163920']. Waiting...
2026-07-15 14:23:15 INFO None 5163918: status FINISHED
2026-07-15 14:23:16 INFO None 5163919: status FINISHED
2026-07-15 14:23:16 INFO None 5163920: status RUNNING/PENDING
2026-07-15 14:23:16 INFO None 5163922: status FINISHED
2026-07-15 14:23:16 INFO None 5163923: status FINISHED
2026-07-15 14:23:16 INFO None 5163924: status FINISHED
2026-07-15 14:23:16 INFO None 5163925: status FINISHED
2026-07-15 14:23:16 INFO None 5163926: status FINISHED
2026-07-15 14:23:16 INFO None 5163927: status FINISHED
2026-07-15 14:23:16 INFO None 5163929: status FINISHED
2026-07-15 14:23:16 INFO None 5163931: status FINISHED
2026-07-15 14:23:16 INFO None 5163932: status FINISHED
2026-07-15 14:23:16 INFO None 5163933: status FINISHED
2026-07-15 14:23:16 INFO None 5163935: status FINISHED
2026-07-15 14:23:16 INFO None 5163937: status FINISHED
2026-07-15 14:23:16 INFO Jobs still running: ['5163920']. Waiting...
2026-07-15 14:23:31 INFO None 5163918: status FINISHED
2026-07-15 14:23:31 INFO None 5163919: status FINISHED
2026-07-15 14:23:31 INFO None 5163920: status RUNNING/PENDING
2026-07-15 14:23:31 INFO None 5163922: status FINISHED
2026-07-15 14:23:31 INFO None 5163923: status FINISHED
2026-07-15 14:23:31 INFO None 5163924: status FINISHED
2026-07-15 14:23:31 INFO None 5163925: status FINISHED
2026-07-15 14:23:31 INFO None 5163926: status FINISHED
2026-07-15 14:23:31 INFO None 5163927: status FINISHED
2026-07-15 14:23:31 INFO None 5163929: status FINISHED
2026-07-15 14:23:31 INFO None 5163931: status FINISHED
2026-07-15 14:23:31 INFO None 5163932: status FINISHED
2026-07-15 14:23:31 INFO None 5163933: status FINISHED
2026-07-15 14:23:31 INFO None 5163935: status FINISHED
2026-07-15 14:23:31 INFO None 5163937: status FINISHED
2026-07-15 14:23:31 INFO Jobs still running: ['5163920']. Waiting...
2026-07-15 14:23:46 INFO None 5163918: status FINISHED
2026-07-15 14:23:46 INFO None 5163919: status FINISHED
2026-07-15 14:23:46 INFO None 5163920: status RUNNING/PENDING
2026-07-15 14:23:46 INFO None 5163922: status FINISHED
2026-07-15 14:23:46 INFO None 5163923: status FINISHED
2026-07-15 14:23:46 INFO None 5163924: status FINISHED
2026-07-15 14:23:46 INFO None 5163925: status FINISHED
2026-07-15 14:23:46 INFO None 5163926: status FINISHED
2026-07-15 14:23:46 INFO None 5163927: status FINISHED
2026-07-15 14:23:46 INFO None 5163929: status FINISHED
2026-07-15 14:23:46 INFO None 5163931: status FINISHED
2026-07-15 14:23:46 INFO None 5163932: status FINISHED
2026-07-15 14:23:46 INFO None 5163933: status FINISHED
2026-07-15 14:23:46 INFO None 5163935: status FINISHED
2026-07-15 14:23:46 INFO None 5163937: status FINISHED
2026-07-15 14:23:46 INFO Jobs still running: ['5163920']. Waiting...
2026-07-15 14:24:01 INFO None 5163918: status FINISHED
2026-07-15 14:24:01 INFO None 5163919: status FINISHED
2026-07-15 14:24:01 INFO None 5163920: status RUNNING/PENDING
2026-07-15 14:24:01 INFO None 5163922: status FINISHED
2026-07-15 14:24:01 INFO None 5163923: status FINISHED
2026-07-15 14:24:01 INFO None 5163924: status FINISHED
2026-07-15 14:24:01 INFO None 5163925: status FINISHED
2026-07-15 14:24:01 INFO None 5163926: status FINISHED
2026-07-15 14:24:02 INFO None 5163927: status FINISHED
2026-07-15 14:24:02 INFO None 5163929: status FINISHED
2026-07-15 14:24:02 INFO None 5163931: status FINISHED
2026-07-15 14:24:02 INFO None 5163932: status FINISHED
2026-07-15 14:24:02 INFO None 5163933: status FINISHED
2026-07-15 14:24:02 INFO None 5163935: status FINISHED
2026-07-15 14:24:02 INFO None 5163937: status FINISHED
2026-07-15 14:24:02 INFO Jobs still running: ['5163920']. Waiting...
2026-07-15 14:24:18 INFO None 5163918: status FINISHED
2026-07-15 14:24:18 INFO None 5163919: status FINISHED
2026-07-15 14:24:19 INFO None 5163920: status RUNNING/PENDING
2026-07-15 14:24:19 INFO None 5163922: status FINISHED
2026-07-15 14:24:19 INFO None 5163923: status FINISHED
2026-07-15 14:24:19 INFO None 5163924: status FINISHED
2026-07-15 14:24:19 INFO None 5163925: status FINISHED
2026-07-15 14:24:19 INFO None 5163926: status FINISHED
2026-07-15 14:24:19 INFO None 5163927: status FINISHED
2026-07-15 14:24:19 INFO None 5163929: status FINISHED
2026-07-15 14:24:19 INFO None 5163931: status FINISHED
2026-07-15 14:24:19 INFO None 5163932: status FINISHED
2026-07-15 14:24:19 INFO None 5163933: status FINISHED
2026-07-15 14:24:19 INFO None 5163935: status FINISHED
2026-07-15 14:24:19 INFO None 5163937: status FINISHED
2026-07-15 14:24:19 INFO Jobs still running: ['5163920']. Waiting...
2026-07-15 14:24:34 INFO None 5163918: status FINISHED
2026-07-15 14:24:34 INFO None 5163919: status FINISHED
2026-07-15 14:24:34 INFO None 5163920: status RUNNING/PENDING
2026-07-15 14:24:34 INFO None 5163922: status FINISHED
2026-07-15 14:24:34 INFO None 5163923: status FINISHED
2026-07-15 14:24:34 INFO None 5163924: status FINISHED
2026-07-15 14:24:34 INFO None 5163925: status FINISHED
2026-07-15 14:24:34 INFO None 5163926: status FINISHED
2026-07-15 14:24:34 INFO None 5163927: status FINISHED
2026-07-15 14:24:34 INFO None 5163929: status FINISHED
2026-07-15 14:24:34 INFO None 5163931: status FINISHED
2026-07-15 14:24:34 INFO None 5163932: status FINISHED
2026-07-15 14:24:34 INFO None 5163933: status FINISHED
2026-07-15 14:24:34 INFO None 5163935: status FINISHED
2026-07-15 14:24:34 INFO None 5163937: status FINISHED
2026-07-15 14:24:34 INFO Jobs still running: ['5163920']. Waiting...
2026-07-15 14:24:49 INFO None 5163918: status FINISHED
2026-07-15 14:24:49 INFO None 5163919: status FINISHED
2026-07-15 14:24:49 INFO None 5163920: status FINISHED
2026-07-15 14:24:49 INFO None 5163922: status FINISHED
2026-07-15 14:24:49 INFO None 5163923: status FINISHED
2026-07-15 14:24:49 INFO None 5163924: status FINISHED
2026-07-15 14:24:49 INFO None 5163925: status FINISHED
2026-07-15 14:24:49 INFO None 5163926: status FINISHED
2026-07-15 14:24:49 INFO None 5163927: status FINISHED
2026-07-15 14:24:49 INFO None 5163929: status FINISHED
2026-07-15 14:24:49 INFO None 5163931: status FINISHED
2026-07-15 14:24:49 INFO None 5163932: status FINISHED
2026-07-15 14:24:49 INFO None 5163933: status FINISHED
2026-07-15 14:24:49 INFO None 5163935: status FINISHED
2026-07-15 14:24:49 INFO None 5163937: status FINISHED
2026-07-15 14:24:49 INFO Jobs ['5163918', '5163919', '5163920', '5163922', '5163923', '5163924', '5163925', '5163926', '5163927', '5163929', '5163931', '5163932', '5163933', '5163935', '5163937'] have finished
2026-07-15 14:24:49 INFO Checking restart files were created ...
2026-07-15 14:24:49 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020600_1_ENS1.nc(668832435 bytes)
2026-07-15 14:24:49 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020600_1_ENS2.nc(668832435 bytes)
2026-07-15 14:24:49 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020600_1_ENS3.nc(668832435 bytes)
2026-07-15 14:24:49 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020600_1_ENS4.nc(668832435 bytes)
2026-07-15 14:24:49 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020600_1_ENS5.nc(668832435 bytes)
2026-07-15 14:24:49 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020600_1_ENS6.nc(668832435 bytes)
2026-07-15 14:24:49 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020600_1_ENS7.nc(668832435 bytes)
2026-07-15 14:24:49 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020600_1_ENS8.nc(668832435 bytes)
2026-07-15 14:24:49 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020600_1_ENS9.nc(668832435 bytes)
2026-07-15 14:24:49 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020600_1_ENS10.nc(668832435 bytes)
2026-07-15 14:24:49 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020600_1_ENS11.nc(668832435 bytes)
2026-07-15 14:24:49 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020600_1_ENS12.nc(668832435 bytes)
2026-07-15 14:24:49 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020600_1_ENS13.nc(668832435 bytes)
2026-07-15 14:24:49 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020600_1_ENS14.nc(668832435 bytes)
2026-07-15 14:24:49 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020600_1_ENS15.nc(668832435 bytes)
2026-07-15 14:24:49 INFO  Run_model() completed successfully.
2026-07-15 14:24:49 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 00:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 14:24:49 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 01:00:00 days=153072 seconds=3600
2026-07-15 14:24:49 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-15 14:24:49 INFO [TIME] increment current_time 2020-02-06 00:00:00 -> 2020-02-06 01:00:00
2026-07-15 14:24:49 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 14:24:49 INFO ---------->>> Running process_satellite_data()
2026-07-15 14:24:49 INFO [DART] No satellite data found, skipping assimilation
2026-07-15 14:24:49 INFO after_assimilation() skipped
2026-07-15 14:24:49 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-15 14:24:49 INFO [TIME] step_end current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 14:24:49 INFO [TIME] step_start current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 14:24:49 INFO [TIME] window start=2020-02-06 01:00:00 end=2020-02-06 09:00:00 run_hours=8 has_assimilation=True
2026-07-15 14:24:49 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:24:51 INFO Hourly dataset computed and listing created
2026-07-15 14:25:07 INFO Hourly dataset computed
2026-07-15 14:25:07 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:25:09 INFO Hourly dataset computed and listing created
2026-07-15 14:25:21 INFO Hourly dataset computed
2026-07-15 14:25:21 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:25:23 INFO Hourly dataset computed and listing created
2026-07-15 14:25:35 INFO Hourly dataset computed
2026-07-15 14:25:35 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:25:37 INFO Hourly dataset computed and listing created
2026-07-15 14:25:50 INFO Hourly dataset computed
2026-07-15 14:25:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:25:52 INFO Hourly dataset computed and listing created
2026-07-15 14:26:05 INFO Hourly dataset computed
2026-07-15 14:26:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:26:07 INFO Hourly dataset computed and listing created
2026-07-15 14:26:20 INFO Hourly dataset computed
2026-07-15 14:26:20 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:26:21 INFO Hourly dataset computed and listing created
2026-07-15 14:26:37 INFO Hourly dataset computed
2026-07-15 14:26:37 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:26:38 INFO Hourly dataset computed and listing created
2026-07-15 14:26:52 INFO Hourly dataset computed
2026-07-15 14:26:52 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:26:54 INFO Hourly dataset computed and listing created
2026-07-15 14:27:10 INFO Hourly dataset computed
2026-07-15 14:27:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:27:12 INFO Hourly dataset computed and listing created
2026-07-15 14:27:26 INFO Hourly dataset computed
2026-07-15 14:27:26 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:27:27 INFO Hourly dataset computed and listing created
2026-07-15 14:27:39 INFO Hourly dataset computed
2026-07-15 14:27:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:27:41 INFO Hourly dataset computed and listing created
2026-07-15 14:27:54 INFO Hourly dataset computed
2026-07-15 14:27:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:27:56 INFO Hourly dataset computed and listing created
2026-07-15 14:28:09 INFO Hourly dataset computed
2026-07-15 14:28:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:28:11 INFO Hourly dataset computed and listing created
2026-07-15 14:28:24 INFO Hourly dataset computed
2026-07-15 14:28:24 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:28:26 INFO Hourly dataset computed and listing created
2026-07-15 14:28:39 INFO Hourly dataset computed
2026-07-15 14:28:39 INFO ---------->>> Running CHIMERE model from 2020-02-06 01:00:00 to 2020-02-06 09:00:00
2026-07-15 14:28:39 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:28:39 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1
2026-07-15 14:28:39 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020600_1_ENS1.nc
2026-07-15 14:28:39 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-15 14:28:39 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:28:39 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-15 14:28:39 INFO Queuing job for member 1...
2026-07-15 14:28:39 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:28:39 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-15 14:28:40 INFO Found: ['5163985']
2026-07-15 14:28:45 INFO [TGCC-IRENE] Submitted job with ID:['5163985']
2026-07-15 14:28:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:28:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2
2026-07-15 14:28:45 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020600_1_ENS2.nc
2026-07-15 14:28:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-15 14:28:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:28:45 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-15 14:28:45 INFO Queuing job for member 2...
2026-07-15 14:28:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:28:45 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-15 14:28:46 INFO Found: ['5163986']
2026-07-15 14:28:51 INFO [TGCC-IRENE] Submitted job with ID:['5163986']
2026-07-15 14:28:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:28:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3
2026-07-15 14:28:51 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020600_1_ENS3.nc
2026-07-15 14:28:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-15 14:28:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:28:51 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-15 14:28:51 INFO Queuing job for member 3...
2026-07-15 14:28:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:28:51 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-15 14:28:51 INFO Found: ['5163988']
2026-07-15 14:28:56 INFO [TGCC-IRENE] Submitted job with ID:['5163988']
2026-07-15 14:28:56 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:28:56 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4
2026-07-15 14:28:56 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020600_1_ENS4.nc
2026-07-15 14:28:56 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-15 14:28:56 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:28:56 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-15 14:28:56 INFO Queuing job for member 4...
2026-07-15 14:28:56 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:28:56 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-15 14:28:57 INFO Found: ['5163990']
2026-07-15 14:29:02 INFO [TGCC-IRENE] Submitted job with ID:['5163990']
2026-07-15 14:29:02 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:29:02 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5
2026-07-15 14:29:02 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020600_1_ENS5.nc
2026-07-15 14:29:02 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-15 14:29:02 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:29:02 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-15 14:29:02 INFO Queuing job for member 5...
2026-07-15 14:29:02 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:29:02 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-15 14:29:03 INFO Found: ['5163991']
2026-07-15 14:29:08 INFO [TGCC-IRENE] Submitted job with ID:['5163991']
2026-07-15 14:29:08 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:29:08 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6
2026-07-15 14:29:08 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020600_1_ENS6.nc
2026-07-15 14:29:08 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-15 14:29:08 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:29:08 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-15 14:29:08 INFO Queuing job for member 6...
2026-07-15 14:29:08 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:29:08 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-15 14:29:11 INFO Found: ['5163994']
2026-07-15 14:29:16 INFO [TGCC-IRENE] Submitted job with ID:['5163994']
2026-07-15 14:29:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:29:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7
2026-07-15 14:29:16 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020600_1_ENS7.nc
2026-07-15 14:29:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-15 14:29:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:29:16 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-15 14:29:16 INFO Queuing job for member 7...
2026-07-15 14:29:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:29:16 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-15 14:29:16 INFO Found: ['5163996']
2026-07-15 14:29:21 INFO [TGCC-IRENE] Submitted job with ID:['5163996']
2026-07-15 14:29:21 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:29:21 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8
2026-07-15 14:29:21 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020600_1_ENS8.nc
2026-07-15 14:29:21 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-15 14:29:21 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:29:21 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-15 14:29:21 INFO Queuing job for member 8...
2026-07-15 14:29:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:29:21 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-15 14:29:22 INFO Found: ['5163997']
2026-07-15 14:29:27 INFO [TGCC-IRENE] Submitted job with ID:['5163997']
2026-07-15 14:29:27 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:29:27 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9
2026-07-15 14:29:27 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020600_1_ENS9.nc
2026-07-15 14:29:27 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-15 14:29:27 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:29:27 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-15 14:29:27 INFO Queuing job for member 9...
2026-07-15 14:29:27 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:29:27 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-15 14:29:28 INFO Found: ['5163998']
2026-07-15 14:29:33 INFO [TGCC-IRENE] Submitted job with ID:['5163998']
2026-07-15 14:29:33 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:29:33 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10
2026-07-15 14:29:33 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020600_1_ENS10.nc
2026-07-15 14:29:33 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-15 14:29:33 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:29:33 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-15 14:29:33 INFO Queuing job for member 10...
2026-07-15 14:29:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:29:33 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-15 14:29:34 INFO Found: ['5163999']
2026-07-15 14:29:39 INFO [TGCC-IRENE] Submitted job with ID:['5163999']
2026-07-15 14:29:39 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:29:39 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11
2026-07-15 14:29:39 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020600_1_ENS11.nc
2026-07-15 14:29:39 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-15 14:29:39 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:29:39 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-15 14:29:39 INFO Queuing job for member 11...
2026-07-15 14:29:39 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:29:39 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-15 14:29:40 INFO Found: ['5164000']
2026-07-15 14:29:45 INFO [TGCC-IRENE] Submitted job with ID:['5164000']
2026-07-15 14:29:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:29:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12
2026-07-15 14:29:45 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020600_1_ENS12.nc
2026-07-15 14:29:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-15 14:29:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:29:45 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-15 14:29:45 INFO Queuing job for member 12...
2026-07-15 14:29:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:29:45 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-15 14:29:45 INFO Found: ['5164001']
2026-07-15 14:29:50 INFO [TGCC-IRENE] Submitted job with ID:['5164001']
2026-07-15 14:29:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:29:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13
2026-07-15 14:29:50 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020600_1_ENS13.nc
2026-07-15 14:29:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-15 14:29:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:29:50 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-15 14:29:50 INFO Queuing job for member 13...
2026-07-15 14:29:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:29:50 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-15 14:29:51 INFO Found: ['5164002']
2026-07-15 14:29:56 INFO [TGCC-IRENE] Submitted job with ID:['5164002']
2026-07-15 14:29:56 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:29:56 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14
2026-07-15 14:29:56 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020600_1_ENS14.nc
2026-07-15 14:29:56 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-15 14:29:56 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:29:56 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-15 14:29:56 INFO Queuing job for member 14...
2026-07-15 14:29:56 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:29:56 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-15 14:29:57 INFO Found: ['5164003']
2026-07-15 14:30:02 INFO [TGCC-IRENE] Submitted job with ID:['5164003']
2026-07-15 14:30:02 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:30:02 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15
2026-07-15 14:30:02 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020600_1_ENS15.nc
2026-07-15 14:30:02 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-15 14:30:02 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:30:02 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-15 14:30:02 INFO Queuing job for member 15...
2026-07-15 14:30:02 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:30:02 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-15 14:30:03 INFO Found: ['5164005']
2026-07-15 14:30:08 INFO [TGCC-IRENE] Submitted job with ID:['5164005']
2026-07-15 14:30:08 INFO Checking job status ...
2026-07-15 14:30:08 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:30:08 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:30:08 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:30:08 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:30:08 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:30:08 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:30:08 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:30:09 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:30:09 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:30:09 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:30:09 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:30:09 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:30:09 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:30:09 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:30:09 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:30:09 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:30:24 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:30:24 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:30:24 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:30:24 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:30:24 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:30:24 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:30:24 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:30:24 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:30:24 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:30:24 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:30:24 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:30:24 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:30:24 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:30:24 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:30:24 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:30:24 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:30:39 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:30:39 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:30:39 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:30:39 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:30:39 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:30:39 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:30:39 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:30:39 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:30:39 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:30:39 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:30:39 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:30:39 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:30:39 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:30:39 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:30:39 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:30:39 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:30:54 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:30:54 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:30:54 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:30:54 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:30:54 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:30:54 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:30:54 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:30:54 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:30:54 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:30:54 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:30:54 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:30:54 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:30:54 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:30:54 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:30:54 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:30:54 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:31:10 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:31:10 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:31:10 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:31:10 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:31:10 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:31:10 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:31:10 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:31:10 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:31:10 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:31:10 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:31:10 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:31:10 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:31:10 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:31:10 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:31:10 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:31:10 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:31:25 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:31:25 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:31:25 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:31:25 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:31:25 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:31:25 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:31:25 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:31:25 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:31:25 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:31:26 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:31:26 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:31:26 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:31:26 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:31:26 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:31:26 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:31:26 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:31:41 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:31:41 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:31:41 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:31:41 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:31:41 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:31:41 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:31:41 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:31:41 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:31:41 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:31:41 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:31:41 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:31:41 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:31:41 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:31:41 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:31:41 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:31:41 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:31:57 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:31:57 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:31:57 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:31:57 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:31:57 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:31:57 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:31:57 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:31:57 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:31:57 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:31:57 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:31:57 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:31:57 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:31:57 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:31:57 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:31:57 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:31:57 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:32:13 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:32:13 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:32:13 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:32:13 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:32:13 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:32:13 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:32:13 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:32:13 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:32:13 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:32:13 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:32:13 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:32:13 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:32:13 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:32:13 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:32:13 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:32:13 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:32:28 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:32:28 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:32:28 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:32:28 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:32:28 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:32:28 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:32:28 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:32:28 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:32:28 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:32:28 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:32:28 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:32:28 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:32:28 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:32:28 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:32:28 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:32:28 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:32:43 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:32:43 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:32:43 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:32:43 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:32:43 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:32:43 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:32:43 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:32:43 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:32:43 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:32:43 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:32:43 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:32:43 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:32:43 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:32:43 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:32:43 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:32:43 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:32:59 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:32:59 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:32:59 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:32:59 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:32:59 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:32:59 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:32:59 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:32:59 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:32:59 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:32:59 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:32:59 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:32:59 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:32:59 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:32:59 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:32:59 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:32:59 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:33:14 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:33:14 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:33:14 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:33:14 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:33:14 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:33:14 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:33:14 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:33:14 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:33:14 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:33:14 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:33:14 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:33:14 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:33:14 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:33:14 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:33:14 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:33:14 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:33:29 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:33:29 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:33:29 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:33:29 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:33:29 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:33:29 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:33:29 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:33:29 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:33:29 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:33:29 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:33:29 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:33:29 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:33:30 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:33:30 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:33:30 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:33:30 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:33:45 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:33:45 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:33:45 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:33:45 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:33:45 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:33:45 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:33:45 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:33:45 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:33:45 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:33:45 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:33:45 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:33:45 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:33:45 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:33:45 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:33:45 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:33:45 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:34:00 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:34:00 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:34:00 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:34:02 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:34:02 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:34:02 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:34:02 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:34:02 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:34:02 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:34:02 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:34:02 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:34:02 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:34:02 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:34:02 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:34:02 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:34:02 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:34:17 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:34:17 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:34:17 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:34:17 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:34:17 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:34:17 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:34:17 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:34:17 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:34:17 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:34:17 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:34:17 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:34:17 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:34:17 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:34:17 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:34:17 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:34:17 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:34:32 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:34:32 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:34:33 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:34:33 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:34:33 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:34:33 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:34:33 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:34:33 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:34:33 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:34:33 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:34:33 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:34:33 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:34:33 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:34:33 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:34:33 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:34:33 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:34:48 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:34:48 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:34:48 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:34:48 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:34:48 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:34:48 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:34:48 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:34:48 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:34:48 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:34:48 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:34:48 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:34:48 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:34:48 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:34:48 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:34:48 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:34:48 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:35:03 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:35:03 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:35:03 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:35:03 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:35:03 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:35:03 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:35:03 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:35:03 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:35:03 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:35:03 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:35:03 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:35:03 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:35:03 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:35:03 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:35:03 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:35:03 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:35:19 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:35:19 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:35:19 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:35:19 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:35:19 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:35:19 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:35:19 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:35:19 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:35:19 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:35:19 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:35:19 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:35:19 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:35:19 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:35:19 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:35:19 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:35:19 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:35:34 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:35:34 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:35:34 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:35:34 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:35:34 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:35:34 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:35:34 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:35:34 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:35:34 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:35:34 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:35:34 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:35:34 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:35:34 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:35:34 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:35:34 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:35:34 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:35:49 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:35:49 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:35:49 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:35:49 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:35:49 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:35:49 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:35:49 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:35:49 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:35:49 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:35:49 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:35:49 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:35:49 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:35:49 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:35:49 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:35:49 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:35:49 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:36:06 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:36:06 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:36:06 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:36:06 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:36:06 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:36:06 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:36:06 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:36:07 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:36:07 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:36:07 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:36:07 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:36:07 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:36:07 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:36:07 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:36:07 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:36:07 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:36:22 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:36:22 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:36:22 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:36:22 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:36:22 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:36:22 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:36:22 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:36:22 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:36:22 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:36:22 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:36:22 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:36:22 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:36:22 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:36:22 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:36:22 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:36:22 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:36:37 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:36:37 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:36:37 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:36:37 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:36:37 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:36:37 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:36:37 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:36:37 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:36:37 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:36:37 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:36:37 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:36:37 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:36:37 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:36:37 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:36:37 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:36:37 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:36:52 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:36:52 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:36:52 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:36:52 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:36:52 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:36:52 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:36:52 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:36:52 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:36:52 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:36:52 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:36:52 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:36:52 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:36:52 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:36:53 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:36:53 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:36:53 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:37:08 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:37:08 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:37:08 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:37:10 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:37:10 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:37:10 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:37:10 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:37:10 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:37:10 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:37:10 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:37:10 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:37:10 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:37:10 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:37:10 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:37:10 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:37:10 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:37:25 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:37:25 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:37:25 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:37:25 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:37:25 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:37:25 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:37:25 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:37:25 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:37:25 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:37:25 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:37:25 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:37:25 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:37:25 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:37:25 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:37:25 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:37:25 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:37:40 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:37:40 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:37:40 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:37:40 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:37:40 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:37:40 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:37:40 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:37:40 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:37:40 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:37:40 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:37:40 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:37:40 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:37:40 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:37:40 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:37:40 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:37:40 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:37:55 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:37:55 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:37:55 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:37:55 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:37:55 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:37:56 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:37:56 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:37:56 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:37:56 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:37:56 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:37:56 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:37:56 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:37:56 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:37:56 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:37:56 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:37:56 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:38:11 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:38:11 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:38:11 INFO None 5163988: status RUNNING/PENDING
2026-07-15 14:38:11 INFO None 5163990: status RUNNING/PENDING
2026-07-15 14:38:11 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:38:11 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:38:11 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:38:11 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:38:11 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:38:11 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:38:11 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:38:11 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:38:11 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:38:11 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:38:11 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:38:11 INFO Jobs still running: ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:38:26 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:38:26 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:38:26 INFO None 5163988: status FINISHED
2026-07-15 14:38:26 INFO None 5163990: status FINISHED
2026-07-15 14:38:26 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:38:26 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:38:26 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:38:26 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:38:26 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:38:26 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:38:26 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:38:26 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:38:26 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:38:26 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:38:26 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:38:26 INFO Jobs still running: ['5163985', '5163986', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:38:41 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:38:41 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:38:41 INFO None 5163988: status FINISHED
2026-07-15 14:38:41 INFO None 5163990: status FINISHED
2026-07-15 14:38:41 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:38:41 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:38:41 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:38:41 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:38:41 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:38:41 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:38:41 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:38:41 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:38:42 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:38:42 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:38:42 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:38:42 INFO Jobs still running: ['5163985', '5163986', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:38:57 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:38:57 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:38:57 INFO None 5163988: status FINISHED
2026-07-15 14:38:57 INFO None 5163990: status FINISHED
2026-07-15 14:38:57 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:38:57 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:38:57 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:38:57 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:38:57 INFO None 5163998: status RUNNING/PENDING
2026-07-15 14:38:57 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:38:57 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:38:57 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:38:57 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:38:57 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:38:57 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:38:57 INFO Jobs still running: ['5163985', '5163986', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:39:12 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:39:12 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:39:12 INFO None 5163988: status FINISHED
2026-07-15 14:39:12 INFO None 5163990: status FINISHED
2026-07-15 14:39:12 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:39:12 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:39:12 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:39:12 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:39:12 INFO None 5163998: status FINISHED
2026-07-15 14:39:12 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:39:12 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:39:12 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:39:12 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:39:12 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:39:12 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:39:12 INFO Jobs still running: ['5163985', '5163986', '5163991', '5163994', '5163996', '5163997', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:39:27 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:39:27 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:39:27 INFO None 5163988: status FINISHED
2026-07-15 14:39:27 INFO None 5163990: status FINISHED
2026-07-15 14:39:27 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:39:27 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:39:27 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:39:27 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:39:27 INFO None 5163998: status FINISHED
2026-07-15 14:39:27 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:39:27 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:39:27 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:39:27 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:39:28 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:39:28 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:39:28 INFO Jobs still running: ['5163985', '5163986', '5163991', '5163994', '5163996', '5163997', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:39:43 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:39:43 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:39:43 INFO None 5163988: status FINISHED
2026-07-15 14:39:43 INFO None 5163990: status FINISHED
2026-07-15 14:39:43 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:39:43 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:39:43 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:39:43 INFO None 5163997: status RUNNING/PENDING
2026-07-15 14:39:43 INFO None 5163998: status FINISHED
2026-07-15 14:39:43 INFO None 5163999: status RUNNING/PENDING
2026-07-15 14:39:43 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:39:43 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:39:43 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:39:43 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:39:43 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:39:43 INFO Jobs still running: ['5163985', '5163986', '5163991', '5163994', '5163996', '5163997', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:39:59 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:39:59 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:39:59 INFO None 5163988: status FINISHED
2026-07-15 14:39:59 INFO None 5163990: status FINISHED
2026-07-15 14:39:59 INFO None 5163991: status RUNNING/PENDING
2026-07-15 14:39:59 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:39:59 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:39:59 INFO None 5163997: status FINISHED
2026-07-15 14:39:59 INFO None 5163998: status FINISHED
2026-07-15 14:39:59 INFO None 5163999: status FINISHED
2026-07-15 14:39:59 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:39:59 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:40:00 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:40:00 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:40:00 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:40:00 INFO Jobs still running: ['5163985', '5163986', '5163991', '5163994', '5163996', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:40:15 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:40:15 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:40:15 INFO None 5163988: status FINISHED
2026-07-15 14:40:15 INFO None 5163990: status FINISHED
2026-07-15 14:40:15 INFO None 5163991: status FINISHED
2026-07-15 14:40:15 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:40:15 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:40:15 INFO None 5163997: status FINISHED
2026-07-15 14:40:15 INFO None 5163998: status FINISHED
2026-07-15 14:40:15 INFO None 5163999: status FINISHED
2026-07-15 14:40:15 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:40:15 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:40:15 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:40:15 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:40:15 INFO None 5164005: status RUNNING/PENDING
2026-07-15 14:40:15 INFO Jobs still running: ['5163985', '5163986', '5163994', '5163996', '5164000', '5164001', '5164002', '5164003', '5164005']. Waiting...
2026-07-15 14:40:30 INFO None 5163985: status RUNNING/PENDING
2026-07-15 14:40:30 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:40:30 INFO None 5163988: status FINISHED
2026-07-15 14:40:30 INFO None 5163990: status FINISHED
2026-07-15 14:40:30 INFO None 5163991: status FINISHED
2026-07-15 14:40:30 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:40:30 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:40:30 INFO None 5163997: status FINISHED
2026-07-15 14:40:30 INFO None 5163998: status FINISHED
2026-07-15 14:40:30 INFO None 5163999: status FINISHED
2026-07-15 14:40:30 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:40:30 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:40:30 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:40:30 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:40:30 INFO None 5164005: status FINISHED
2026-07-15 14:40:30 INFO Jobs still running: ['5163985', '5163986', '5163994', '5163996', '5164000', '5164001', '5164002', '5164003']. Waiting...
2026-07-15 14:40:45 INFO None 5163985: status FINISHED
2026-07-15 14:40:45 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:40:45 INFO None 5163988: status FINISHED
2026-07-15 14:40:45 INFO None 5163990: status FINISHED
2026-07-15 14:40:45 INFO None 5163991: status FINISHED
2026-07-15 14:40:45 INFO None 5163994: status RUNNING/PENDING
2026-07-15 14:40:45 INFO None 5163996: status RUNNING/PENDING
2026-07-15 14:40:45 INFO None 5163997: status FINISHED
2026-07-15 14:40:45 INFO None 5163998: status FINISHED
2026-07-15 14:40:45 INFO None 5163999: status FINISHED
2026-07-15 14:40:45 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:40:45 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:40:45 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:40:45 INFO None 5164003: status RUNNING/PENDING
2026-07-15 14:40:45 INFO None 5164005: status FINISHED
2026-07-15 14:40:45 INFO Jobs still running: ['5163986', '5163994', '5163996', '5164000', '5164001', '5164002', '5164003']. Waiting...
2026-07-15 14:41:01 INFO None 5163985: status FINISHED
2026-07-15 14:41:01 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:41:01 INFO None 5163988: status FINISHED
2026-07-15 14:41:01 INFO None 5163990: status FINISHED
2026-07-15 14:41:01 INFO None 5163991: status FINISHED
2026-07-15 14:41:01 INFO None 5163994: status FINISHED
2026-07-15 14:41:01 INFO None 5163996: status FINISHED
2026-07-15 14:41:01 INFO None 5163997: status FINISHED
2026-07-15 14:41:01 INFO None 5163998: status FINISHED
2026-07-15 14:41:01 INFO None 5163999: status FINISHED
2026-07-15 14:41:01 INFO None 5164000: status RUNNING/PENDING
2026-07-15 14:41:01 INFO None 5164001: status RUNNING/PENDING
2026-07-15 14:41:01 INFO None 5164002: status RUNNING/PENDING
2026-07-15 14:41:01 INFO None 5164003: status FINISHED
2026-07-15 14:41:01 INFO None 5164005: status FINISHED
2026-07-15 14:41:01 INFO Jobs still running: ['5163986', '5164000', '5164001', '5164002']. Waiting...
2026-07-15 14:41:16 INFO None 5163985: status FINISHED
2026-07-15 14:41:16 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:41:16 INFO None 5163988: status FINISHED
2026-07-15 14:41:16 INFO None 5163990: status FINISHED
2026-07-15 14:41:16 INFO None 5163991: status FINISHED
2026-07-15 14:41:16 INFO None 5163994: status FINISHED
2026-07-15 14:41:16 INFO None 5163996: status FINISHED
2026-07-15 14:41:16 INFO None 5163997: status FINISHED
2026-07-15 14:41:16 INFO None 5163998: status FINISHED
2026-07-15 14:41:16 INFO None 5163999: status FINISHED
2026-07-15 14:41:16 INFO None 5164000: status FINISHED
2026-07-15 14:41:16 INFO None 5164001: status FINISHED
2026-07-15 14:41:16 INFO None 5164002: status FINISHED
2026-07-15 14:41:16 INFO None 5164003: status FINISHED
2026-07-15 14:41:16 INFO None 5164005: status FINISHED
2026-07-15 14:41:16 INFO Jobs still running: ['5163986']. Waiting...
2026-07-15 14:41:31 INFO None 5163985: status FINISHED
2026-07-15 14:41:32 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:41:32 INFO None 5163988: status FINISHED
2026-07-15 14:41:32 INFO None 5163990: status FINISHED
2026-07-15 14:41:32 INFO None 5163991: status FINISHED
2026-07-15 14:41:32 INFO None 5163994: status FINISHED
2026-07-15 14:41:32 INFO None 5163996: status FINISHED
2026-07-15 14:41:32 INFO None 5163997: status FINISHED
2026-07-15 14:41:32 INFO None 5163998: status FINISHED
2026-07-15 14:41:32 INFO None 5163999: status FINISHED
2026-07-15 14:41:32 INFO None 5164000: status FINISHED
2026-07-15 14:41:32 INFO None 5164001: status FINISHED
2026-07-15 14:41:32 INFO None 5164002: status FINISHED
2026-07-15 14:41:32 INFO None 5164003: status FINISHED
2026-07-15 14:41:32 INFO None 5164005: status FINISHED
2026-07-15 14:41:32 INFO Jobs still running: ['5163986']. Waiting...
2026-07-15 14:41:47 INFO None 5163985: status FINISHED
2026-07-15 14:41:47 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:41:47 INFO None 5163988: status FINISHED
2026-07-15 14:41:47 INFO None 5163990: status FINISHED
2026-07-15 14:41:47 INFO None 5163991: status FINISHED
2026-07-15 14:41:47 INFO None 5163994: status FINISHED
2026-07-15 14:41:47 INFO None 5163996: status FINISHED
2026-07-15 14:41:47 INFO None 5163997: status FINISHED
2026-07-15 14:41:47 INFO None 5163998: status FINISHED
2026-07-15 14:41:47 INFO None 5163999: status FINISHED
2026-07-15 14:41:47 INFO None 5164000: status FINISHED
2026-07-15 14:41:47 INFO None 5164001: status FINISHED
2026-07-15 14:41:47 INFO None 5164002: status FINISHED
2026-07-15 14:41:47 INFO None 5164003: status FINISHED
2026-07-15 14:41:47 INFO None 5164005: status FINISHED
2026-07-15 14:41:47 INFO Jobs still running: ['5163986']. Waiting...
2026-07-15 14:42:02 INFO None 5163985: status FINISHED
2026-07-15 14:42:02 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:42:02 INFO None 5163988: status FINISHED
2026-07-15 14:42:02 INFO None 5163990: status FINISHED
2026-07-15 14:42:02 INFO None 5163991: status FINISHED
2026-07-15 14:42:02 INFO None 5163994: status FINISHED
2026-07-15 14:42:02 INFO None 5163996: status FINISHED
2026-07-15 14:42:02 INFO None 5163997: status FINISHED
2026-07-15 14:42:02 INFO None 5163998: status FINISHED
2026-07-15 14:42:04 INFO None 5163999: status FINISHED
2026-07-15 14:42:04 INFO None 5164000: status FINISHED
2026-07-15 14:42:04 INFO None 5164001: status FINISHED
2026-07-15 14:42:04 INFO None 5164002: status FINISHED
2026-07-15 14:42:04 INFO None 5164003: status FINISHED
2026-07-15 14:42:04 INFO None 5164005: status FINISHED
2026-07-15 14:42:04 INFO Jobs still running: ['5163986']. Waiting...
2026-07-15 14:42:19 INFO None 5163985: status FINISHED
2026-07-15 14:42:19 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:42:19 INFO None 5163988: status FINISHED
2026-07-15 14:42:19 INFO None 5163990: status FINISHED
2026-07-15 14:42:19 INFO None 5163991: status FINISHED
2026-07-15 14:42:19 INFO None 5163994: status FINISHED
2026-07-15 14:42:19 INFO None 5163996: status FINISHED
2026-07-15 14:42:19 INFO None 5163997: status FINISHED
2026-07-15 14:42:20 INFO None 5163998: status FINISHED
2026-07-15 14:42:20 INFO None 5163999: status FINISHED
2026-07-15 14:42:20 INFO None 5164000: status FINISHED
2026-07-15 14:42:20 INFO None 5164001: status FINISHED
2026-07-15 14:42:20 INFO None 5164002: status FINISHED
2026-07-15 14:42:20 INFO None 5164003: status FINISHED
2026-07-15 14:42:20 INFO None 5164005: status FINISHED
2026-07-15 14:42:20 INFO Jobs still running: ['5163986']. Waiting...
2026-07-15 14:42:35 INFO None 5163985: status FINISHED
2026-07-15 14:42:35 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:42:35 INFO None 5163988: status FINISHED
2026-07-15 14:42:35 INFO None 5163990: status FINISHED
2026-07-15 14:42:35 INFO None 5163991: status FINISHED
2026-07-15 14:42:35 INFO None 5163994: status FINISHED
2026-07-15 14:42:35 INFO None 5163996: status FINISHED
2026-07-15 14:42:35 INFO None 5163997: status FINISHED
2026-07-15 14:42:35 INFO None 5163998: status FINISHED
2026-07-15 14:42:35 INFO None 5163999: status FINISHED
2026-07-15 14:42:35 INFO None 5164000: status FINISHED
2026-07-15 14:42:35 INFO None 5164001: status FINISHED
2026-07-15 14:42:35 INFO None 5164002: status FINISHED
2026-07-15 14:42:35 INFO None 5164003: status FINISHED
2026-07-15 14:42:35 INFO None 5164005: status FINISHED
2026-07-15 14:42:35 INFO Jobs still running: ['5163986']. Waiting...
2026-07-15 14:42:50 INFO None 5163985: status FINISHED
2026-07-15 14:42:50 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:42:52 INFO None 5163988: status FINISHED
2026-07-15 14:42:52 INFO None 5163990: status FINISHED
2026-07-15 14:42:52 INFO None 5163991: status FINISHED
2026-07-15 14:42:52 INFO None 5163994: status FINISHED
2026-07-15 14:42:52 INFO None 5163996: status FINISHED
2026-07-15 14:42:52 INFO None 5163997: status FINISHED
2026-07-15 14:42:52 INFO None 5163998: status FINISHED
2026-07-15 14:42:52 INFO None 5163999: status FINISHED
2026-07-15 14:42:52 INFO None 5164000: status FINISHED
2026-07-15 14:42:52 INFO None 5164001: status FINISHED
2026-07-15 14:42:52 INFO None 5164002: status FINISHED
2026-07-15 14:42:52 INFO None 5164003: status FINISHED
2026-07-15 14:42:52 INFO None 5164005: status FINISHED
2026-07-15 14:42:52 INFO Jobs still running: ['5163986']. Waiting...
2026-07-15 14:43:07 INFO None 5163985: status FINISHED
2026-07-15 14:43:07 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:43:07 INFO None 5163988: status FINISHED
2026-07-15 14:43:07 INFO None 5163990: status FINISHED
2026-07-15 14:43:07 INFO None 5163991: status FINISHED
2026-07-15 14:43:07 INFO None 5163994: status FINISHED
2026-07-15 14:43:07 INFO None 5163996: status FINISHED
2026-07-15 14:43:07 INFO None 5163997: status FINISHED
2026-07-15 14:43:07 INFO None 5163998: status FINISHED
2026-07-15 14:43:07 INFO None 5163999: status FINISHED
2026-07-15 14:43:07 INFO None 5164000: status FINISHED
2026-07-15 14:43:07 INFO None 5164001: status FINISHED
2026-07-15 14:43:07 INFO None 5164002: status FINISHED
2026-07-15 14:43:07 INFO None 5164003: status FINISHED
2026-07-15 14:43:07 INFO None 5164005: status FINISHED
2026-07-15 14:43:07 INFO Jobs still running: ['5163986']. Waiting...
2026-07-15 14:43:22 INFO None 5163985: status FINISHED
2026-07-15 14:43:23 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:43:23 INFO None 5163988: status FINISHED
2026-07-15 14:43:23 INFO None 5163990: status FINISHED
2026-07-15 14:43:23 INFO None 5163991: status FINISHED
2026-07-15 14:43:23 INFO None 5163994: status FINISHED
2026-07-15 14:43:23 INFO None 5163996: status FINISHED
2026-07-15 14:43:23 INFO None 5163997: status FINISHED
2026-07-15 14:43:23 INFO None 5163998: status FINISHED
2026-07-15 14:43:23 INFO None 5163999: status FINISHED
2026-07-15 14:43:23 INFO None 5164000: status FINISHED
2026-07-15 14:43:23 INFO None 5164001: status FINISHED
2026-07-15 14:43:23 INFO None 5164002: status FINISHED
2026-07-15 14:43:23 INFO None 5164003: status FINISHED
2026-07-15 14:43:23 INFO None 5164005: status FINISHED
2026-07-15 14:43:23 INFO Jobs still running: ['5163986']. Waiting...
2026-07-15 14:43:38 INFO None 5163985: status FINISHED
2026-07-15 14:43:38 INFO None 5163986: status RUNNING/PENDING
2026-07-15 14:43:38 INFO None 5163988: status FINISHED
2026-07-15 14:43:38 INFO None 5163990: status FINISHED
2026-07-15 14:43:38 INFO None 5163991: status FINISHED
2026-07-15 14:43:38 INFO None 5163994: status FINISHED
2026-07-15 14:43:38 INFO None 5163996: status FINISHED
2026-07-15 14:43:38 INFO None 5163997: status FINISHED
2026-07-15 14:43:38 INFO None 5163998: status FINISHED
2026-07-15 14:43:38 INFO None 5163999: status FINISHED
2026-07-15 14:43:38 INFO None 5164000: status FINISHED
2026-07-15 14:43:38 INFO None 5164001: status FINISHED
2026-07-15 14:43:38 INFO None 5164002: status FINISHED
2026-07-15 14:43:38 INFO None 5164003: status FINISHED
2026-07-15 14:43:38 INFO None 5164005: status FINISHED
2026-07-15 14:43:38 INFO Jobs still running: ['5163986']. Waiting...
2026-07-15 14:43:54 INFO None 5163985: status FINISHED
2026-07-15 14:43:54 INFO None 5163986: status FINISHED
2026-07-15 14:43:54 INFO None 5163988: status FINISHED
2026-07-15 14:43:54 INFO None 5163990: status FINISHED
2026-07-15 14:43:54 INFO None 5163991: status FINISHED
2026-07-15 14:43:54 INFO None 5163994: status FINISHED
2026-07-15 14:43:54 INFO None 5163996: status FINISHED
2026-07-15 14:43:54 INFO None 5163997: status FINISHED
2026-07-15 14:43:54 INFO None 5163998: status FINISHED
2026-07-15 14:43:54 INFO None 5163999: status FINISHED
2026-07-15 14:43:54 INFO None 5164000: status FINISHED
2026-07-15 14:43:54 INFO None 5164001: status FINISHED
2026-07-15 14:43:54 INFO None 5164002: status FINISHED
2026-07-15 14:43:54 INFO None 5164003: status FINISHED
2026-07-15 14:43:54 INFO None 5164005: status FINISHED
2026-07-15 14:43:54 INFO Jobs ['5163985', '5163986', '5163988', '5163990', '5163991', '5163994', '5163996', '5163997', '5163998', '5163999', '5164000', '5164001', '5164002', '5164003', '5164005'] have finished
2026-07-15 14:43:54 INFO Checking restart files were created ...
2026-07-15 14:43:54 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020601_8_ENS1.nc(3005806795 bytes)
2026-07-15 14:43:54 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020601_8_ENS2.nc(3005806795 bytes)
2026-07-15 14:43:54 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020601_8_ENS3.nc(3005806795 bytes)
2026-07-15 14:43:54 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020601_8_ENS4.nc(3005806795 bytes)
2026-07-15 14:43:54 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020601_8_ENS5.nc(3005806795 bytes)
2026-07-15 14:43:54 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020601_8_ENS6.nc(3005806795 bytes)
2026-07-15 14:43:54 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020601_8_ENS7.nc(3005806795 bytes)
2026-07-15 14:43:54 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020601_8_ENS8.nc(3005806795 bytes)
2026-07-15 14:43:54 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020601_8_ENS9.nc(3005806795 bytes)
2026-07-15 14:43:54 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020601_8_ENS10.nc(3005806795 bytes)
2026-07-15 14:43:54 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020601_8_ENS11.nc(3005806795 bytes)
2026-07-15 14:43:54 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020601_8_ENS12.nc(3005806795 bytes)
2026-07-15 14:43:54 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020601_8_ENS13.nc(3005806795 bytes)
2026-07-15 14:43:54 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020601_8_ENS14.nc(3005806795 bytes)
2026-07-15 14:43:54 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020601_8_ENS15.nc(3005806795 bytes)
2026-07-15 14:43:54 INFO  Run_model() completed successfully.
2026-07-15 14:43:54 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 14:43:54 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 09:00:00 days=153072 seconds=32400
2026-07-15 14:43:54 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-15 14:43:54 INFO [TIME] increment current_time 2020-02-06 01:00:00 -> 2020-02-06 09:00:00
2026-07-15 14:43:54 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 14:43:54 INFO ---------->>> Running process_satellite_data()
2026-07-15 14:43:54 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12002.nc
2026-07-15 14:43:54 INFO ---------->>> Running run_obs_converter()
2026-07-15 14:43:54 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_33067_153072.out
2026-07-15 14:43:54 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_33067_153072.out
2026-07-15 14:43:54 INFO ---------->>> Running DART
2026-07-15 14:43:54 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-15 14:43:54 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/chim_ENS1_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/chim_ENS1_2020020609_1_out_toDART.nc
2026-07-15 14:43:54 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/chim_ENS2_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/chim_ENS2_2020020609_1_out_toDART.nc
2026-07-15 14:43:55 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/chim_ENS3_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/chim_ENS3_2020020609_1_out_toDART.nc
2026-07-15 14:43:55 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/chim_ENS4_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/chim_ENS4_2020020609_1_out_toDART.nc
2026-07-15 14:43:55 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/chim_ENS5_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/chim_ENS5_2020020609_1_out_toDART.nc
2026-07-15 14:43:56 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/chim_ENS6_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/chim_ENS6_2020020609_1_out_toDART.nc
2026-07-15 14:43:56 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/chim_ENS7_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/chim_ENS7_2020020609_1_out_toDART.nc
2026-07-15 14:43:57 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/chim_ENS8_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/chim_ENS8_2020020609_1_out_toDART.nc
2026-07-15 14:43:57 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/chim_ENS9_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/chim_ENS9_2020020609_1_out_toDART.nc
2026-07-15 14:43:57 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/chim_ENS10_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/chim_ENS10_2020020609_1_out_toDART.nc
2026-07-15 14:43:58 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/chim_ENS11_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/chim_ENS11_2020020609_1_out_toDART.nc
2026-07-15 14:43:58 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/chim_ENS12_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/chim_ENS12_2020020609_1_out_toDART.nc
2026-07-15 14:43:58 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/chim_ENS13_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/chim_ENS13_2020020609_1_out_toDART.nc
2026-07-15 14:43:59 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/chim_ENS14_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/chim_ENS14_2020020609_1_out_toDART.nc
2026-07-15 14:43:59 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/chim_ENS15_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/chim_ENS15_2020020609_1_out_toDART.nc
2026-07-15 14:44:00 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-15 14:44:00 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-15 14:44:00 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-15 14:44:00 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-15 14:44:00 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-15 14:44:00 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-15 14:44:07 INFO Found: []
2026-07-15 14:44:07 INFO No job id returned by command ./run_filter.bsh
2026-07-15 14:44:07 INFO No monitoring will be performed
2026-07-15 14:44:07 INFO Moving DART output files to analysis and preassim directories for date 2020020609 if present ...
2026-07-15 14:44:07 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-15 14:44:07 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-15 14:44:07 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-15 14:44:07 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-15 14:44:07 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-15 14:44:07 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-15 14:44:07 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-15 14:44:07 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-15 14:44:07 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-15 14:44:07 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-15 14:44:07 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-15 14:44:07 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-15 14:44:07 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-15 14:44:07 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-15 14:44:07 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-15 14:44:07 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-15 14:44:07 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-15 14:44:07 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-15 14:44:07 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-15 14:44:07 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-15 14:44:07 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-15 14:44:07 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-15 14:44:07 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-15 14:44:07 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-15 14:44:07 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-15 14:44:07 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-15 14:44:07 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-15 14:44:07 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-15 14:44:07 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-15 14:44:07 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-15 14:44:07 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-15 14:44:07 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020609'
2026-07-15 14:44:07 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-15 14:44:07 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020609'
2026-07-15 14:44:07 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-15 14:44:10 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-15 14:44:10 INFO run_dart() is DONE.
2026-07-15 14:44:10 INFO ---------->>> Running update_pollutant_in_end()
2026-07-15 14:44:10 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020601_8_ENS1.nc
2026-07-15 14:44:24 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 14:44:25 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020601_8_ENS2.nc
2026-07-15 14:44:39 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 14:44:39 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020601_8_ENS3.nc
2026-07-15 14:44:53 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 14:44:54 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020601_8_ENS4.nc
2026-07-15 14:45:09 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 14:45:09 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020601_8_ENS5.nc
2026-07-15 14:45:24 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 14:45:24 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020601_8_ENS6.nc
2026-07-15 14:45:38 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 14:45:39 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020601_8_ENS7.nc
2026-07-15 14:45:53 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 14:45:53 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020601_8_ENS8.nc
2026-07-15 14:46:07 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 14:46:08 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020601_8_ENS9.nc
2026-07-15 14:46:22 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 14:46:22 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020601_8_ENS10.nc
2026-07-15 14:46:36 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 14:46:37 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020601_8_ENS11.nc
2026-07-15 14:46:50 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 14:46:50 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020601_8_ENS12.nc
2026-07-15 14:47:04 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 14:47:05 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020601_8_ENS13.nc
2026-07-15 14:47:18 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 14:47:19 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020601_8_ENS14.nc
2026-07-15 14:47:33 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 14:47:34 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020601_8_ENS15.nc
2026-07-15 14:47:48 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 14:47:48 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-15 14:47:48 INFO [TIME] step_end current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 14:47:48 INFO [TIME] step_start current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 14:47:48 INFO [TIME] window start=2020-02-06 09:00:00 end=2020-02-06 11:00:00 run_hours=2 has_assimilation=True
2026-07-15 14:47:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:47:49 INFO Hourly dataset computed and listing created
2026-07-15 14:47:55 INFO Hourly dataset computed
2026-07-15 14:47:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:47:56 INFO Hourly dataset computed and listing created
2026-07-15 14:48:00 INFO Hourly dataset computed
2026-07-15 14:48:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:48:01 INFO Hourly dataset computed and listing created
2026-07-15 14:48:04 INFO Hourly dataset computed
2026-07-15 14:48:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:48:05 INFO Hourly dataset computed and listing created
2026-07-15 14:48:08 INFO Hourly dataset computed
2026-07-15 14:48:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:48:09 INFO Hourly dataset computed and listing created
2026-07-15 14:48:13 INFO Hourly dataset computed
2026-07-15 14:48:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:48:14 INFO Hourly dataset computed and listing created
2026-07-15 14:48:17 INFO Hourly dataset computed
2026-07-15 14:48:18 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:48:19 INFO Hourly dataset computed and listing created
2026-07-15 14:48:21 INFO Hourly dataset computed
2026-07-15 14:48:21 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:48:22 INFO Hourly dataset computed and listing created
2026-07-15 14:48:24 INFO Hourly dataset computed
2026-07-15 14:48:24 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:48:26 INFO Hourly dataset computed and listing created
2026-07-15 14:48:28 INFO Hourly dataset computed
2026-07-15 14:48:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:48:29 INFO Hourly dataset computed and listing created
2026-07-15 14:48:32 INFO Hourly dataset computed
2026-07-15 14:48:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:48:34 INFO Hourly dataset computed and listing created
2026-07-15 14:48:36 INFO Hourly dataset computed
2026-07-15 14:48:36 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:48:37 INFO Hourly dataset computed and listing created
2026-07-15 14:48:39 INFO Hourly dataset computed
2026-07-15 14:48:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:48:41 INFO Hourly dataset computed and listing created
2026-07-15 14:48:43 INFO Hourly dataset computed
2026-07-15 14:48:43 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:48:44 INFO Hourly dataset computed and listing created
2026-07-15 14:48:47 INFO Hourly dataset computed
2026-07-15 14:48:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:48:48 INFO Hourly dataset computed and listing created
2026-07-15 14:48:50 INFO Hourly dataset computed
2026-07-15 14:48:50 INFO ---------->>> Running CHIMERE model from 2020-02-06 09:00:00 to 2020-02-06 11:00:00
2026-07-15 14:48:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:48:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1
2026-07-15 14:48:50 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020601_8_ENS1.nc
2026-07-15 14:48:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-15 14:48:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:48:50 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-15 14:48:51 INFO Queuing job for member 1...
2026-07-15 14:48:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:48:51 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-15 14:48:52 INFO Found: ['5164196']
2026-07-15 14:48:57 INFO [TGCC-IRENE] Submitted job with ID:['5164196']
2026-07-15 14:48:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:48:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2
2026-07-15 14:48:57 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020601_8_ENS2.nc
2026-07-15 14:48:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-15 14:48:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:48:57 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-15 14:48:57 INFO Queuing job for member 2...
2026-07-15 14:48:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:48:57 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-15 14:48:59 INFO Found: ['5164197']
2026-07-15 14:49:04 INFO [TGCC-IRENE] Submitted job with ID:['5164197']
2026-07-15 14:49:04 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:49:04 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3
2026-07-15 14:49:04 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020601_8_ENS3.nc
2026-07-15 14:49:04 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-15 14:49:04 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:49:04 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-15 14:49:04 INFO Queuing job for member 3...
2026-07-15 14:49:04 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:49:04 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-15 14:49:06 INFO Found: ['5164199']
2026-07-15 14:49:11 INFO [TGCC-IRENE] Submitted job with ID:['5164199']
2026-07-15 14:49:11 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:49:11 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4
2026-07-15 14:49:11 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020601_8_ENS4.nc
2026-07-15 14:49:11 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-15 14:49:11 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:49:11 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-15 14:49:11 INFO Queuing job for member 4...
2026-07-15 14:49:11 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:49:11 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-15 14:49:11 INFO Found: ['5164200']
2026-07-15 14:49:16 INFO [TGCC-IRENE] Submitted job with ID:['5164200']
2026-07-15 14:49:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:49:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5
2026-07-15 14:49:16 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020601_8_ENS5.nc
2026-07-15 14:49:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-15 14:49:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:49:16 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-15 14:49:16 INFO Queuing job for member 5...
2026-07-15 14:49:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:49:16 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-15 14:49:17 INFO Found: ['5164211']
2026-07-15 14:49:22 INFO [TGCC-IRENE] Submitted job with ID:['5164211']
2026-07-15 14:49:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:49:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6
2026-07-15 14:49:22 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020601_8_ENS6.nc
2026-07-15 14:49:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-15 14:49:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:49:22 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-15 14:49:22 INFO Queuing job for member 6...
2026-07-15 14:49:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:49:22 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-15 14:49:23 INFO Found: ['5164212']
2026-07-15 14:49:28 INFO [TGCC-IRENE] Submitted job with ID:['5164212']
2026-07-15 14:49:28 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:49:28 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7
2026-07-15 14:49:28 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020601_8_ENS7.nc
2026-07-15 14:49:28 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-15 14:49:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:49:28 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-15 14:49:28 INFO Queuing job for member 7...
2026-07-15 14:49:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:49:28 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-15 14:49:29 INFO Found: ['5164213']
2026-07-15 14:49:34 INFO [TGCC-IRENE] Submitted job with ID:['5164213']
2026-07-15 14:49:34 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:49:34 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8
2026-07-15 14:49:34 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020601_8_ENS8.nc
2026-07-15 14:49:34 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-15 14:49:34 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:49:34 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-15 14:49:34 INFO Queuing job for member 8...
2026-07-15 14:49:34 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:49:34 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-15 14:49:34 INFO Found: ['5164215']
2026-07-15 14:49:39 INFO [TGCC-IRENE] Submitted job with ID:['5164215']
2026-07-15 14:49:39 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:49:39 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9
2026-07-15 14:49:39 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020601_8_ENS9.nc
2026-07-15 14:49:39 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-15 14:49:39 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:49:39 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-15 14:49:39 INFO Queuing job for member 9...
2026-07-15 14:49:39 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:49:39 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-15 14:49:40 INFO Found: ['5164217']
2026-07-15 14:49:45 INFO [TGCC-IRENE] Submitted job with ID:['5164217']
2026-07-15 14:49:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:49:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10
2026-07-15 14:49:45 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020601_8_ENS10.nc
2026-07-15 14:49:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-15 14:49:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:49:45 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-15 14:49:45 INFO Queuing job for member 10...
2026-07-15 14:49:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:49:45 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-15 14:49:46 INFO Found: ['5164219']
2026-07-15 14:49:51 INFO [TGCC-IRENE] Submitted job with ID:['5164219']
2026-07-15 14:49:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:49:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11
2026-07-15 14:49:51 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020601_8_ENS11.nc
2026-07-15 14:49:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-15 14:49:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:49:51 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-15 14:49:51 INFO Queuing job for member 11...
2026-07-15 14:49:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:49:51 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-15 14:49:53 INFO Found: ['5164221']
2026-07-15 14:49:58 INFO [TGCC-IRENE] Submitted job with ID:['5164221']
2026-07-15 14:49:58 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:49:58 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12
2026-07-15 14:49:58 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020601_8_ENS12.nc
2026-07-15 14:49:58 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-15 14:49:58 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:49:58 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-15 14:49:58 INFO Queuing job for member 12...
2026-07-15 14:49:58 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:49:58 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-15 14:50:00 INFO Found: ['5164223']
2026-07-15 14:50:05 INFO [TGCC-IRENE] Submitted job with ID:['5164223']
2026-07-15 14:50:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:50:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13
2026-07-15 14:50:05 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020601_8_ENS13.nc
2026-07-15 14:50:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-15 14:50:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:50:05 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-15 14:50:05 INFO Queuing job for member 13...
2026-07-15 14:50:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:50:05 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-15 14:50:08 INFO Found: ['5164226']
2026-07-15 14:50:13 INFO [TGCC-IRENE] Submitted job with ID:['5164226']
2026-07-15 14:50:13 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:50:13 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14
2026-07-15 14:50:13 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020601_8_ENS14.nc
2026-07-15 14:50:13 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-15 14:50:13 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:50:13 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-15 14:50:13 INFO Queuing job for member 14...
2026-07-15 14:50:13 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:50:13 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-15 14:50:13 INFO Found: ['5164228']
2026-07-15 14:50:18 INFO [TGCC-IRENE] Submitted job with ID:['5164228']
2026-07-15 14:50:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:50:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15
2026-07-15 14:50:18 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020601_8_ENS15.nc
2026-07-15 14:50:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-15 14:50:18 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:50:18 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-15 14:50:18 INFO Queuing job for member 15...
2026-07-15 14:50:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:50:18 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-15 14:50:19 INFO Found: ['5164229']
2026-07-15 14:50:24 INFO [TGCC-IRENE] Submitted job with ID:['5164229']
2026-07-15 14:50:24 INFO Checking job status ...
2026-07-15 14:50:24 INFO None 5164196: status RUNNING/PENDING
2026-07-15 14:50:24 INFO None 5164197: status RUNNING/PENDING
2026-07-15 14:50:24 INFO None 5164199: status RUNNING/PENDING
2026-07-15 14:50:24 INFO None 5164200: status RUNNING/PENDING
2026-07-15 14:50:24 INFO None 5164211: status RUNNING/PENDING
2026-07-15 14:50:24 INFO None 5164212: status RUNNING/PENDING
2026-07-15 14:50:24 INFO None 5164213: status RUNNING/PENDING
2026-07-15 14:50:24 INFO None 5164215: status RUNNING/PENDING
2026-07-15 14:50:24 INFO None 5164217: status RUNNING/PENDING
2026-07-15 14:50:24 INFO None 5164219: status RUNNING/PENDING
2026-07-15 14:50:24 INFO None 5164221: status RUNNING/PENDING
2026-07-15 14:50:24 INFO None 5164223: status RUNNING/PENDING
2026-07-15 14:50:24 INFO None 5164226: status RUNNING/PENDING
2026-07-15 14:50:24 INFO None 5164228: status RUNNING/PENDING
2026-07-15 14:50:24 INFO None 5164229: status RUNNING/PENDING
2026-07-15 14:50:24 INFO Jobs still running: ['5164196', '5164197', '5164199', '5164200', '5164211', '5164212', '5164213', '5164215', '5164217', '5164219', '5164221', '5164223', '5164226', '5164228', '5164229']. Waiting...
2026-07-15 14:50:40 INFO None 5164196: status RUNNING/PENDING
2026-07-15 14:50:40 INFO None 5164197: status RUNNING/PENDING
2026-07-15 14:50:40 INFO None 5164199: status RUNNING/PENDING
2026-07-15 14:50:40 INFO None 5164200: status RUNNING/PENDING
2026-07-15 14:50:40 INFO None 5164211: status RUNNING/PENDING
2026-07-15 14:50:40 INFO None 5164212: status RUNNING/PENDING
2026-07-15 14:50:40 INFO None 5164213: status RUNNING/PENDING
2026-07-15 14:50:40 INFO None 5164215: status RUNNING/PENDING
2026-07-15 14:50:40 INFO None 5164217: status RUNNING/PENDING
2026-07-15 14:50:40 INFO None 5164219: status RUNNING/PENDING
2026-07-15 14:50:40 INFO None 5164221: status RUNNING/PENDING
2026-07-15 14:50:40 INFO None 5164223: status RUNNING/PENDING
2026-07-15 14:50:40 INFO None 5164226: status RUNNING/PENDING
2026-07-15 14:50:40 INFO None 5164228: status RUNNING/PENDING
2026-07-15 14:50:40 INFO None 5164229: status RUNNING/PENDING
2026-07-15 14:50:40 INFO Jobs still running: ['5164196', '5164197', '5164199', '5164200', '5164211', '5164212', '5164213', '5164215', '5164217', '5164219', '5164221', '5164223', '5164226', '5164228', '5164229']. Waiting...
2026-07-15 14:50:55 INFO None 5164196: status RUNNING/PENDING
2026-07-15 14:50:55 INFO None 5164197: status RUNNING/PENDING
2026-07-15 14:50:55 INFO None 5164199: status RUNNING/PENDING
2026-07-15 14:50:55 INFO None 5164200: status RUNNING/PENDING
2026-07-15 14:50:55 INFO None 5164211: status RUNNING/PENDING
2026-07-15 14:50:55 INFO None 5164212: status RUNNING/PENDING
2026-07-15 14:50:55 INFO None 5164213: status RUNNING/PENDING
2026-07-15 14:50:55 INFO None 5164215: status RUNNING/PENDING
2026-07-15 14:50:55 INFO None 5164217: status RUNNING/PENDING
2026-07-15 14:50:55 INFO None 5164219: status RUNNING/PENDING
2026-07-15 14:50:55 INFO None 5164221: status RUNNING/PENDING
2026-07-15 14:50:55 INFO None 5164223: status RUNNING/PENDING
2026-07-15 14:50:55 INFO None 5164226: status RUNNING/PENDING
2026-07-15 14:50:55 INFO None 5164228: status RUNNING/PENDING
2026-07-15 14:50:55 INFO None 5164229: status RUNNING/PENDING
2026-07-15 14:50:55 INFO Jobs still running: ['5164196', '5164197', '5164199', '5164200', '5164211', '5164212', '5164213', '5164215', '5164217', '5164219', '5164221', '5164223', '5164226', '5164228', '5164229']. Waiting...
2026-07-15 14:51:11 INFO None 5164196: status RUNNING/PENDING
2026-07-15 14:51:12 INFO None 5164197: status RUNNING/PENDING
2026-07-15 14:51:12 INFO None 5164199: status RUNNING/PENDING
2026-07-15 14:51:12 INFO None 5164200: status RUNNING/PENDING
2026-07-15 14:51:12 INFO None 5164211: status RUNNING/PENDING
2026-07-15 14:51:12 INFO None 5164212: status RUNNING/PENDING
2026-07-15 14:51:12 INFO None 5164213: status RUNNING/PENDING
2026-07-15 14:51:12 INFO None 5164215: status RUNNING/PENDING
2026-07-15 14:51:12 INFO None 5164217: status RUNNING/PENDING
2026-07-15 14:51:12 INFO None 5164219: status RUNNING/PENDING
2026-07-15 14:51:12 INFO None 5164221: status RUNNING/PENDING
2026-07-15 14:51:12 INFO None 5164223: status RUNNING/PENDING
2026-07-15 14:51:12 INFO None 5164226: status RUNNING/PENDING
2026-07-15 14:51:12 INFO None 5164228: status RUNNING/PENDING
2026-07-15 14:51:12 INFO None 5164229: status RUNNING/PENDING
2026-07-15 14:51:12 INFO Jobs still running: ['5164196', '5164197', '5164199', '5164200', '5164211', '5164212', '5164213', '5164215', '5164217', '5164219', '5164221', '5164223', '5164226', '5164228', '5164229']. Waiting...
2026-07-15 14:51:27 INFO None 5164196: status RUNNING/PENDING
2026-07-15 14:51:27 INFO None 5164197: status RUNNING/PENDING
2026-07-15 14:51:27 INFO None 5164199: status RUNNING/PENDING
2026-07-15 14:51:27 INFO None 5164200: status RUNNING/PENDING
2026-07-15 14:51:27 INFO None 5164211: status RUNNING/PENDING
2026-07-15 14:51:27 INFO None 5164212: status RUNNING/PENDING
2026-07-15 14:51:27 INFO None 5164213: status RUNNING/PENDING
2026-07-15 14:51:27 INFO None 5164215: status RUNNING/PENDING
2026-07-15 14:51:27 INFO None 5164217: status RUNNING/PENDING
2026-07-15 14:51:27 INFO None 5164219: status RUNNING/PENDING
2026-07-15 14:51:27 INFO None 5164221: status RUNNING/PENDING
2026-07-15 14:51:27 INFO None 5164223: status RUNNING/PENDING
2026-07-15 14:51:27 INFO None 5164226: status RUNNING/PENDING
2026-07-15 14:51:27 INFO None 5164228: status RUNNING/PENDING
2026-07-15 14:51:27 INFO None 5164229: status RUNNING/PENDING
2026-07-15 14:51:27 INFO Jobs still running: ['5164196', '5164197', '5164199', '5164200', '5164211', '5164212', '5164213', '5164215', '5164217', '5164219', '5164221', '5164223', '5164226', '5164228', '5164229']. Waiting...
2026-07-15 14:51:42 INFO None 5164196: status RUNNING/PENDING
2026-07-15 14:51:42 INFO None 5164197: status RUNNING/PENDING
2026-07-15 14:51:42 INFO None 5164199: status RUNNING/PENDING
2026-07-15 14:51:42 INFO None 5164200: status RUNNING/PENDING
2026-07-15 14:51:42 INFO None 5164211: status RUNNING/PENDING
2026-07-15 14:51:42 INFO None 5164212: status RUNNING/PENDING
2026-07-15 14:51:42 INFO None 5164213: status RUNNING/PENDING
2026-07-15 14:51:42 INFO None 5164215: status RUNNING/PENDING
2026-07-15 14:51:42 INFO None 5164217: status RUNNING/PENDING
2026-07-15 14:51:42 INFO None 5164219: status RUNNING/PENDING
2026-07-15 14:51:42 INFO None 5164221: status RUNNING/PENDING
2026-07-15 14:51:42 INFO None 5164223: status RUNNING/PENDING
2026-07-15 14:51:42 INFO None 5164226: status RUNNING/PENDING
2026-07-15 14:51:42 INFO None 5164228: status RUNNING/PENDING
2026-07-15 14:51:42 INFO None 5164229: status RUNNING/PENDING
2026-07-15 14:51:42 INFO Jobs still running: ['5164196', '5164197', '5164199', '5164200', '5164211', '5164212', '5164213', '5164215', '5164217', '5164219', '5164221', '5164223', '5164226', '5164228', '5164229']. Waiting...
2026-07-15 14:51:57 INFO None 5164196: status RUNNING/PENDING
2026-07-15 14:51:57 INFO None 5164197: status RUNNING/PENDING
2026-07-15 14:51:57 INFO None 5164199: status RUNNING/PENDING
2026-07-15 14:51:57 INFO None 5164200: status RUNNING/PENDING
2026-07-15 14:51:57 INFO None 5164211: status RUNNING/PENDING
2026-07-15 14:51:57 INFO None 5164212: status RUNNING/PENDING
2026-07-15 14:51:57 INFO None 5164213: status RUNNING/PENDING
2026-07-15 14:51:58 INFO None 5164215: status RUNNING/PENDING
2026-07-15 14:51:58 INFO None 5164217: status RUNNING/PENDING
2026-07-15 14:51:58 INFO None 5164219: status RUNNING/PENDING
2026-07-15 14:51:58 INFO None 5164221: status RUNNING/PENDING
2026-07-15 14:51:58 INFO None 5164223: status RUNNING/PENDING
2026-07-15 14:51:58 INFO None 5164226: status RUNNING/PENDING
2026-07-15 14:51:58 INFO None 5164228: status RUNNING/PENDING
2026-07-15 14:51:58 INFO None 5164229: status RUNNING/PENDING
2026-07-15 14:51:58 INFO Jobs still running: ['5164196', '5164197', '5164199', '5164200', '5164211', '5164212', '5164213', '5164215', '5164217', '5164219', '5164221', '5164223', '5164226', '5164228', '5164229']. Waiting...
2026-07-15 14:52:13 INFO None 5164196: status RUNNING/PENDING
2026-07-15 14:52:13 INFO None 5164197: status RUNNING/PENDING
2026-07-15 14:52:13 INFO None 5164199: status RUNNING/PENDING
2026-07-15 14:52:13 INFO None 5164200: status RUNNING/PENDING
2026-07-15 14:52:14 INFO None 5164211: status RUNNING/PENDING
2026-07-15 14:52:14 INFO None 5164212: status RUNNING/PENDING
2026-07-15 14:52:14 INFO None 5164213: status RUNNING/PENDING
2026-07-15 14:52:14 INFO None 5164215: status RUNNING/PENDING
2026-07-15 14:52:14 INFO None 5164217: status RUNNING/PENDING
2026-07-15 14:52:14 INFO None 5164219: status RUNNING/PENDING
2026-07-15 14:52:14 INFO None 5164221: status RUNNING/PENDING
2026-07-15 14:52:14 INFO None 5164223: status RUNNING/PENDING
2026-07-15 14:52:14 INFO None 5164226: status RUNNING/PENDING
2026-07-15 14:52:14 INFO None 5164228: status RUNNING/PENDING
2026-07-15 14:52:14 INFO None 5164229: status RUNNING/PENDING
2026-07-15 14:52:14 INFO Jobs still running: ['5164196', '5164197', '5164199', '5164200', '5164211', '5164212', '5164213', '5164215', '5164217', '5164219', '5164221', '5164223', '5164226', '5164228', '5164229']. Waiting...
2026-07-15 14:52:29 INFO None 5164196: status RUNNING/PENDING
2026-07-15 14:52:29 INFO None 5164197: status RUNNING/PENDING
2026-07-15 14:52:29 INFO None 5164199: status RUNNING/PENDING
2026-07-15 14:52:29 INFO None 5164200: status RUNNING/PENDING
2026-07-15 14:52:29 INFO None 5164211: status RUNNING/PENDING
2026-07-15 14:52:29 INFO None 5164212: status RUNNING/PENDING
2026-07-15 14:52:29 INFO None 5164213: status RUNNING/PENDING
2026-07-15 14:52:29 INFO None 5164215: status RUNNING/PENDING
2026-07-15 14:52:29 INFO None 5164217: status RUNNING/PENDING
2026-07-15 14:52:29 INFO None 5164219: status RUNNING/PENDING
2026-07-15 14:52:29 INFO None 5164221: status RUNNING/PENDING
2026-07-15 14:52:29 INFO None 5164223: status RUNNING/PENDING
2026-07-15 14:52:29 INFO None 5164226: status RUNNING/PENDING
2026-07-15 14:52:29 INFO None 5164228: status RUNNING/PENDING
2026-07-15 14:52:29 INFO None 5164229: status RUNNING/PENDING
2026-07-15 14:52:29 INFO Jobs still running: ['5164196', '5164197', '5164199', '5164200', '5164211', '5164212', '5164213', '5164215', '5164217', '5164219', '5164221', '5164223', '5164226', '5164228', '5164229']. Waiting...
2026-07-15 14:52:44 INFO None 5164196: status RUNNING/PENDING
2026-07-15 14:52:44 INFO None 5164197: status RUNNING/PENDING
2026-07-15 14:52:44 INFO None 5164199: status RUNNING/PENDING
2026-07-15 14:52:44 INFO None 5164200: status RUNNING/PENDING
2026-07-15 14:52:44 INFO None 5164211: status RUNNING/PENDING
2026-07-15 14:52:44 INFO None 5164212: status RUNNING/PENDING
2026-07-15 14:52:44 INFO None 5164213: status RUNNING/PENDING
2026-07-15 14:52:44 INFO None 5164215: status RUNNING/PENDING
2026-07-15 14:52:44 INFO None 5164217: status RUNNING/PENDING
2026-07-15 14:52:44 INFO None 5164219: status RUNNING/PENDING
2026-07-15 14:52:44 INFO None 5164221: status RUNNING/PENDING
2026-07-15 14:52:44 INFO None 5164223: status RUNNING/PENDING
2026-07-15 14:52:44 INFO None 5164226: status RUNNING/PENDING
2026-07-15 14:52:44 INFO None 5164228: status RUNNING/PENDING
2026-07-15 14:52:44 INFO None 5164229: status RUNNING/PENDING
2026-07-15 14:52:44 INFO Jobs still running: ['5164196', '5164197', '5164199', '5164200', '5164211', '5164212', '5164213', '5164215', '5164217', '5164219', '5164221', '5164223', '5164226', '5164228', '5164229']. Waiting...
2026-07-15 14:52:59 INFO None 5164196: status RUNNING/PENDING
2026-07-15 14:52:59 INFO None 5164197: status RUNNING/PENDING
2026-07-15 14:52:59 INFO None 5164199: status RUNNING/PENDING
2026-07-15 14:52:59 INFO None 5164200: status RUNNING/PENDING
2026-07-15 14:52:59 INFO None 5164211: status RUNNING/PENDING
2026-07-15 14:52:59 INFO None 5164212: status RUNNING/PENDING
2026-07-15 14:52:59 INFO None 5164213: status RUNNING/PENDING
2026-07-15 14:52:59 INFO None 5164215: status RUNNING/PENDING
2026-07-15 14:52:59 INFO None 5164217: status FINISHED
2026-07-15 14:52:59 INFO None 5164219: status FINISHED
2026-07-15 14:53:00 INFO None 5164221: status RUNNING/PENDING
2026-07-15 14:53:00 INFO None 5164223: status RUNNING/PENDING
2026-07-15 14:53:00 INFO None 5164226: status RUNNING/PENDING
2026-07-15 14:53:00 INFO None 5164228: status RUNNING/PENDING
2026-07-15 14:53:00 INFO None 5164229: status RUNNING/PENDING
2026-07-15 14:53:00 INFO Jobs still running: ['5164196', '5164197', '5164199', '5164200', '5164211', '5164212', '5164213', '5164215', '5164221', '5164223', '5164226', '5164228', '5164229']. Waiting...
2026-07-15 14:53:15 INFO None 5164196: status RUNNING/PENDING
2026-07-15 14:53:15 INFO None 5164197: status RUNNING/PENDING
2026-07-15 14:53:15 INFO None 5164199: status RUNNING/PENDING
2026-07-15 14:53:15 INFO None 5164200: status RUNNING/PENDING
2026-07-15 14:53:15 INFO None 5164211: status RUNNING/PENDING
2026-07-15 14:53:15 INFO None 5164212: status RUNNING/PENDING
2026-07-15 14:53:15 INFO None 5164213: status RUNNING/PENDING
2026-07-15 14:53:15 INFO None 5164215: status RUNNING/PENDING
2026-07-15 14:53:15 INFO None 5164217: status FINISHED
2026-07-15 14:53:15 INFO None 5164219: status FINISHED
2026-07-15 14:53:15 INFO None 5164221: status RUNNING/PENDING
2026-07-15 14:53:15 INFO None 5164223: status RUNNING/PENDING
2026-07-15 14:53:15 INFO None 5164226: status RUNNING/PENDING
2026-07-15 14:53:15 INFO None 5164228: status RUNNING/PENDING
2026-07-15 14:53:15 INFO None 5164229: status RUNNING/PENDING
2026-07-15 14:53:15 INFO Jobs still running: ['5164196', '5164197', '5164199', '5164200', '5164211', '5164212', '5164213', '5164215', '5164221', '5164223', '5164226', '5164228', '5164229']. Waiting...
2026-07-15 14:53:30 INFO None 5164196: status RUNNING/PENDING
2026-07-15 14:53:30 INFO None 5164197: status RUNNING/PENDING
2026-07-15 14:53:30 INFO None 5164199: status RUNNING/PENDING
2026-07-15 14:53:30 INFO None 5164200: status RUNNING/PENDING
2026-07-15 14:53:31 INFO None 5164211: status RUNNING/PENDING
2026-07-15 14:53:31 INFO None 5164212: status RUNNING/PENDING
2026-07-15 14:53:31 INFO None 5164213: status RUNNING/PENDING
2026-07-15 14:53:31 INFO None 5164215: status RUNNING/PENDING
2026-07-15 14:53:31 INFO None 5164217: status FINISHED
2026-07-15 14:53:31 INFO None 5164219: status FINISHED
2026-07-15 14:53:33 INFO None 5164221: status RUNNING/PENDING
2026-07-15 14:53:33 INFO None 5164223: status RUNNING/PENDING
2026-07-15 14:53:33 INFO None 5164226: status RUNNING/PENDING
2026-07-15 14:53:33 INFO None 5164228: status RUNNING/PENDING
2026-07-15 14:53:33 INFO None 5164229: status RUNNING/PENDING
2026-07-15 14:53:33 INFO Jobs still running: ['5164196', '5164197', '5164199', '5164200', '5164211', '5164212', '5164213', '5164215', '5164221', '5164223', '5164226', '5164228', '5164229']. Waiting...
2026-07-15 14:53:48 INFO None 5164196: status RUNNING/PENDING
2026-07-15 14:53:48 INFO None 5164197: status RUNNING/PENDING
2026-07-15 14:53:48 INFO None 5164199: status FINISHED
2026-07-15 14:53:48 INFO None 5164200: status FINISHED
2026-07-15 14:53:48 INFO None 5164211: status FINISHED
2026-07-15 14:53:48 INFO None 5164212: status FINISHED
2026-07-15 14:53:48 INFO None 5164213: status FINISHED
2026-07-15 14:53:48 INFO None 5164215: status RUNNING/PENDING
2026-07-15 14:53:48 INFO None 5164217: status FINISHED
2026-07-15 14:53:48 INFO None 5164219: status FINISHED
2026-07-15 14:53:48 INFO None 5164221: status RUNNING/PENDING
2026-07-15 14:53:48 INFO None 5164223: status RUNNING/PENDING
2026-07-15 14:53:48 INFO None 5164226: status RUNNING/PENDING
2026-07-15 14:53:48 INFO None 5164228: status RUNNING/PENDING
2026-07-15 14:53:48 INFO None 5164229: status RUNNING/PENDING
2026-07-15 14:53:48 INFO Jobs still running: ['5164196', '5164197', '5164215', '5164221', '5164223', '5164226', '5164228', '5164229']. Waiting...
2026-07-15 14:54:03 INFO None 5164196: status RUNNING/PENDING
2026-07-15 14:54:03 INFO None 5164197: status RUNNING/PENDING
2026-07-15 14:54:03 INFO None 5164199: status FINISHED
2026-07-15 14:54:03 INFO None 5164200: status FINISHED
2026-07-15 14:54:03 INFO None 5164211: status FINISHED
2026-07-15 14:54:03 INFO None 5164212: status FINISHED
2026-07-15 14:54:03 INFO None 5164213: status FINISHED
2026-07-15 14:54:03 INFO None 5164215: status RUNNING/PENDING
2026-07-15 14:54:03 INFO None 5164217: status FINISHED
2026-07-15 14:54:03 INFO None 5164219: status FINISHED
2026-07-15 14:54:03 INFO None 5164221: status RUNNING/PENDING
2026-07-15 14:54:03 INFO None 5164223: status RUNNING/PENDING
2026-07-15 14:54:03 INFO None 5164226: status RUNNING/PENDING
2026-07-15 14:54:03 INFO None 5164228: status RUNNING/PENDING
2026-07-15 14:54:03 INFO None 5164229: status RUNNING/PENDING
2026-07-15 14:54:03 INFO Jobs still running: ['5164196', '5164197', '5164215', '5164221', '5164223', '5164226', '5164228', '5164229']. Waiting...
2026-07-15 14:54:20 INFO None 5164196: status RUNNING/PENDING
2026-07-15 14:54:20 INFO None 5164197: status FINISHED
2026-07-15 14:54:20 INFO None 5164199: status FINISHED
2026-07-15 14:54:20 INFO None 5164200: status FINISHED
2026-07-15 14:54:20 INFO None 5164211: status FINISHED
2026-07-15 14:54:20 INFO None 5164212: status FINISHED
2026-07-15 14:54:20 INFO None 5164213: status FINISHED
2026-07-15 14:54:20 INFO None 5164215: status RUNNING/PENDING
2026-07-15 14:54:20 INFO None 5164217: status FINISHED
2026-07-15 14:54:20 INFO None 5164219: status FINISHED
2026-07-15 14:54:20 INFO None 5164221: status RUNNING/PENDING
2026-07-15 14:54:20 INFO None 5164223: status RUNNING/PENDING
2026-07-15 14:54:20 INFO None 5164226: status RUNNING/PENDING
2026-07-15 14:54:20 INFO None 5164228: status RUNNING/PENDING
2026-07-15 14:54:20 INFO None 5164229: status RUNNING/PENDING
2026-07-15 14:54:20 INFO Jobs still running: ['5164196', '5164215', '5164221', '5164223', '5164226', '5164228', '5164229']. Waiting...
2026-07-15 14:54:35 INFO None 5164196: status RUNNING/PENDING
2026-07-15 14:54:35 INFO None 5164197: status FINISHED
2026-07-15 14:54:35 INFO None 5164199: status FINISHED
2026-07-15 14:54:35 INFO None 5164200: status FINISHED
2026-07-15 14:54:35 INFO None 5164211: status FINISHED
2026-07-15 14:54:35 INFO None 5164212: status FINISHED
2026-07-15 14:54:35 INFO None 5164213: status FINISHED
2026-07-15 14:54:35 INFO None 5164215: status RUNNING/PENDING
2026-07-15 14:54:35 INFO None 5164217: status FINISHED
2026-07-15 14:54:35 INFO None 5164219: status FINISHED
2026-07-15 14:54:35 INFO None 5164221: status RUNNING/PENDING
2026-07-15 14:54:35 INFO None 5164223: status RUNNING/PENDING
2026-07-15 14:54:35 INFO None 5164226: status RUNNING/PENDING
2026-07-15 14:54:37 INFO None 5164228: status RUNNING/PENDING
2026-07-15 14:54:37 INFO None 5164229: status RUNNING/PENDING
2026-07-15 14:54:37 INFO Jobs still running: ['5164196', '5164215', '5164221', '5164223', '5164226', '5164228', '5164229']. Waiting...
2026-07-15 14:54:52 INFO None 5164196: status RUNNING/PENDING
2026-07-15 14:54:52 INFO None 5164197: status FINISHED
2026-07-15 14:54:52 INFO None 5164199: status FINISHED
2026-07-15 14:54:52 INFO None 5164200: status FINISHED
2026-07-15 14:54:52 INFO None 5164211: status FINISHED
2026-07-15 14:54:52 INFO None 5164212: status FINISHED
2026-07-15 14:54:52 INFO None 5164213: status FINISHED
2026-07-15 14:54:52 INFO None 5164215: status RUNNING/PENDING
2026-07-15 14:54:52 INFO None 5164217: status FINISHED
2026-07-15 14:54:52 INFO None 5164219: status FINISHED
2026-07-15 14:54:52 INFO None 5164221: status RUNNING/PENDING
2026-07-15 14:54:52 INFO None 5164223: status RUNNING/PENDING
2026-07-15 14:54:52 INFO None 5164226: status RUNNING/PENDING
2026-07-15 14:54:52 INFO None 5164228: status RUNNING/PENDING
2026-07-15 14:54:52 INFO None 5164229: status FINISHED
2026-07-15 14:54:52 INFO Jobs still running: ['5164196', '5164215', '5164221', '5164223', '5164226', '5164228']. Waiting...
2026-07-15 14:55:07 INFO None 5164196: status RUNNING/PENDING
2026-07-15 14:55:07 INFO None 5164197: status FINISHED
2026-07-15 14:55:07 INFO None 5164199: status FINISHED
2026-07-15 14:55:07 INFO None 5164200: status FINISHED
2026-07-15 14:55:08 INFO None 5164211: status FINISHED
2026-07-15 14:55:08 INFO None 5164212: status FINISHED
2026-07-15 14:55:08 INFO None 5164213: status FINISHED
2026-07-15 14:55:08 INFO None 5164215: status RUNNING/PENDING
2026-07-15 14:55:08 INFO None 5164217: status FINISHED
2026-07-15 14:55:08 INFO None 5164219: status FINISHED
2026-07-15 14:55:08 INFO None 5164221: status FINISHED
2026-07-15 14:55:08 INFO None 5164223: status FINISHED
2026-07-15 14:55:08 INFO None 5164226: status RUNNING/PENDING
2026-07-15 14:55:08 INFO None 5164228: status RUNNING/PENDING
2026-07-15 14:55:08 INFO None 5164229: status FINISHED
2026-07-15 14:55:08 INFO Jobs still running: ['5164196', '5164215', '5164226', '5164228']. Waiting...
2026-07-15 14:55:23 INFO None 5164196: status RUNNING/PENDING
2026-07-15 14:55:23 INFO None 5164197: status FINISHED
2026-07-15 14:55:23 INFO None 5164199: status FINISHED
2026-07-15 14:55:23 INFO None 5164200: status FINISHED
2026-07-15 14:55:23 INFO None 5164211: status FINISHED
2026-07-15 14:55:23 INFO None 5164212: status FINISHED
2026-07-15 14:55:23 INFO None 5164213: status FINISHED
2026-07-15 14:55:23 INFO None 5164215: status RUNNING/PENDING
2026-07-15 14:55:23 INFO None 5164217: status FINISHED
2026-07-15 14:55:23 INFO None 5164219: status FINISHED
2026-07-15 14:55:23 INFO None 5164221: status FINISHED
2026-07-15 14:55:23 INFO None 5164223: status FINISHED
2026-07-15 14:55:23 INFO None 5164226: status FINISHED
2026-07-15 14:55:23 INFO None 5164228: status RUNNING/PENDING
2026-07-15 14:55:23 INFO None 5164229: status FINISHED
2026-07-15 14:55:23 INFO Jobs still running: ['5164196', '5164215', '5164228']. Waiting...
2026-07-15 14:55:39 INFO None 5164196: status RUNNING/PENDING
2026-07-15 14:55:39 INFO None 5164197: status FINISHED
2026-07-15 14:55:39 INFO None 5164199: status FINISHED
2026-07-15 14:55:39 INFO None 5164200: status FINISHED
2026-07-15 14:55:39 INFO None 5164211: status FINISHED
2026-07-15 14:55:39 INFO None 5164212: status FINISHED
2026-07-15 14:55:39 INFO None 5164213: status FINISHED
2026-07-15 14:55:39 INFO None 5164215: status FINISHED
2026-07-15 14:55:39 INFO None 5164217: status FINISHED
2026-07-15 14:55:39 INFO None 5164219: status FINISHED
2026-07-15 14:55:39 INFO None 5164221: status FINISHED
2026-07-15 14:55:39 INFO None 5164223: status FINISHED
2026-07-15 14:55:39 INFO None 5164226: status FINISHED
2026-07-15 14:55:39 INFO None 5164228: status RUNNING/PENDING
2026-07-15 14:55:39 INFO None 5164229: status FINISHED
2026-07-15 14:55:39 INFO Jobs still running: ['5164196', '5164228']. Waiting...
2026-07-15 14:55:54 INFO None 5164196: status FINISHED
2026-07-15 14:55:54 INFO None 5164197: status FINISHED
2026-07-15 14:55:54 INFO None 5164199: status FINISHED
2026-07-15 14:55:54 INFO None 5164200: status FINISHED
2026-07-15 14:55:54 INFO None 5164211: status FINISHED
2026-07-15 14:55:54 INFO None 5164212: status FINISHED
2026-07-15 14:55:54 INFO None 5164213: status FINISHED
2026-07-15 14:55:54 INFO None 5164215: status FINISHED
2026-07-15 14:55:54 INFO None 5164217: status FINISHED
2026-07-15 14:55:54 INFO None 5164219: status FINISHED
2026-07-15 14:55:54 INFO None 5164221: status FINISHED
2026-07-15 14:55:55 INFO None 5164223: status FINISHED
2026-07-15 14:55:55 INFO None 5164226: status FINISHED
2026-07-15 14:55:55 INFO None 5164228: status RUNNING/PENDING
2026-07-15 14:55:55 INFO None 5164229: status FINISHED
2026-07-15 14:55:55 INFO Jobs still running: ['5164228']. Waiting...
2026-07-15 14:56:10 INFO None 5164196: status FINISHED
2026-07-15 14:56:10 INFO None 5164197: status FINISHED
2026-07-15 14:56:10 INFO None 5164199: status FINISHED
2026-07-15 14:56:10 INFO None 5164200: status FINISHED
2026-07-15 14:56:10 INFO None 5164211: status FINISHED
2026-07-15 14:56:10 INFO None 5164212: status FINISHED
2026-07-15 14:56:10 INFO None 5164213: status FINISHED
2026-07-15 14:56:10 INFO None 5164215: status FINISHED
2026-07-15 14:56:10 INFO None 5164217: status FINISHED
2026-07-15 14:56:10 INFO None 5164219: status FINISHED
2026-07-15 14:56:10 INFO None 5164221: status FINISHED
2026-07-15 14:56:10 INFO None 5164223: status FINISHED
2026-07-15 14:56:10 INFO None 5164226: status FINISHED
2026-07-15 14:56:10 INFO None 5164228: status RUNNING/PENDING
2026-07-15 14:56:10 INFO None 5164229: status FINISHED
2026-07-15 14:56:10 INFO Jobs still running: ['5164228']. Waiting...
2026-07-15 14:56:25 INFO None 5164196: status FINISHED
2026-07-15 14:56:25 INFO None 5164197: status FINISHED
2026-07-15 14:56:25 INFO None 5164199: status FINISHED
2026-07-15 14:56:25 INFO None 5164200: status FINISHED
2026-07-15 14:56:25 INFO None 5164211: status FINISHED
2026-07-15 14:56:25 INFO None 5164212: status FINISHED
2026-07-15 14:56:25 INFO None 5164213: status FINISHED
2026-07-15 14:56:25 INFO None 5164215: status FINISHED
2026-07-15 14:56:25 INFO None 5164217: status FINISHED
2026-07-15 14:56:25 INFO None 5164219: status FINISHED
2026-07-15 14:56:25 INFO None 5164221: status FINISHED
2026-07-15 14:56:25 INFO None 5164223: status FINISHED
2026-07-15 14:56:25 INFO None 5164226: status FINISHED
2026-07-15 14:56:25 INFO None 5164228: status FINISHED
2026-07-15 14:56:25 INFO None 5164229: status FINISHED
2026-07-15 14:56:25 INFO Jobs ['5164196', '5164197', '5164199', '5164200', '5164211', '5164212', '5164213', '5164215', '5164217', '5164219', '5164221', '5164223', '5164226', '5164228', '5164229'] have finished
2026-07-15 14:56:25 INFO Checking restart files were created ...
2026-07-15 14:56:25 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc(1002685915 bytes)
2026-07-15 14:56:25 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc(1002685915 bytes)
2026-07-15 14:56:25 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc(1002685915 bytes)
2026-07-15 14:56:25 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc(1002685915 bytes)
2026-07-15 14:56:25 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc(1002685915 bytes)
2026-07-15 14:56:25 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc(1002685915 bytes)
2026-07-15 14:56:25 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc(1002685915 bytes)
2026-07-15 14:56:25 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc(1002685915 bytes)
2026-07-15 14:56:25 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc(1002685915 bytes)
2026-07-15 14:56:25 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc(1002685915 bytes)
2026-07-15 14:56:25 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc(1002685915 bytes)
2026-07-15 14:56:25 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc(1002685915 bytes)
2026-07-15 14:56:25 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc(1002685915 bytes)
2026-07-15 14:56:25 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc(1002685915 bytes)
2026-07-15 14:56:25 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc(1002685915 bytes)
2026-07-15 14:56:25 INFO  Run_model() completed successfully.
2026-07-15 14:56:25 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 14:56:25 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 11:00:00 days=153072 seconds=39600
2026-07-15 14:56:25 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-15 14:56:25 INFO [TIME] increment current_time 2020-02-06 09:00:00 -> 2020-02-06 11:00:00
2026-07-15 14:56:25 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 14:56:25 INFO ---------->>> Running process_satellite_data()
2026-07-15 14:56:25 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12003.nc
2026-07-15 14:56:25 INFO ---------->>> Running run_obs_converter()
2026-07-15 14:56:25 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_39045_153072.out
2026-07-15 14:56:25 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_39045_153072.out
2026-07-15 14:56:25 INFO ---------->>> Running DART
2026-07-15 14:56:25 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-15 14:56:25 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/chim_ENS1_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/chim_ENS1_2020020611_1_out_toDART.nc
2026-07-15 14:56:26 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/chim_ENS2_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/chim_ENS2_2020020611_1_out_toDART.nc
2026-07-15 14:56:26 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/chim_ENS3_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/chim_ENS3_2020020611_1_out_toDART.nc
2026-07-15 14:56:26 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/chim_ENS4_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/chim_ENS4_2020020611_1_out_toDART.nc
2026-07-15 14:56:27 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/chim_ENS5_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/chim_ENS5_2020020611_1_out_toDART.nc
2026-07-15 14:56:27 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/chim_ENS6_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/chim_ENS6_2020020611_1_out_toDART.nc
2026-07-15 14:56:27 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/chim_ENS7_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/chim_ENS7_2020020611_1_out_toDART.nc
2026-07-15 14:56:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/chim_ENS8_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/chim_ENS8_2020020611_1_out_toDART.nc
2026-07-15 14:56:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/chim_ENS9_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/chim_ENS9_2020020611_1_out_toDART.nc
2026-07-15 14:56:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/chim_ENS10_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/chim_ENS10_2020020611_1_out_toDART.nc
2026-07-15 14:56:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/chim_ENS11_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/chim_ENS11_2020020611_1_out_toDART.nc
2026-07-15 14:56:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/chim_ENS12_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/chim_ENS12_2020020611_1_out_toDART.nc
2026-07-15 14:56:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/chim_ENS13_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/chim_ENS13_2020020611_1_out_toDART.nc
2026-07-15 14:56:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/chim_ENS14_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/chim_ENS14_2020020611_1_out_toDART.nc
2026-07-15 14:56:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/chim_ENS15_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/chim_ENS15_2020020611_1_out_toDART.nc
2026-07-15 14:56:30 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-15 14:56:30 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-15 14:56:30 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-15 14:56:30 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-15 14:56:30 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-15 14:56:30 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-15 14:56:43 INFO Found: []
2026-07-15 14:56:43 INFO No job id returned by command ./run_filter.bsh
2026-07-15 14:56:43 INFO No monitoring will be performed
2026-07-15 14:56:43 INFO Moving DART output files to analysis and preassim directories for date 2020020611 if present ...
2026-07-15 14:56:43 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-15 14:56:43 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-15 14:56:43 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-15 14:56:43 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-15 14:56:43 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-15 14:56:43 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-15 14:56:43 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-15 14:56:43 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-15 14:56:43 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-15 14:56:43 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-15 14:56:43 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-15 14:56:43 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-15 14:56:43 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-15 14:56:43 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-15 14:56:43 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-15 14:56:43 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-15 14:56:43 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-15 14:56:43 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-15 14:56:43 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-15 14:56:43 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-15 14:56:43 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-15 14:56:43 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-15 14:56:43 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-15 14:56:43 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-15 14:56:43 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-15 14:56:43 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-15 14:56:43 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-15 14:56:43 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-15 14:56:43 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-15 14:56:44 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-15 14:56:44 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-15 14:56:44 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-15 14:56:44 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-15 14:56:44 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-15 14:56:44 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-15 14:56:44 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-15 14:56:44 INFO run_dart() is DONE.
2026-07-15 14:56:44 INFO ---------->>> Running update_pollutant_in_end()
2026-07-15 14:56:44 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc
2026-07-15 14:56:49 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 14:56:49 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc
2026-07-15 14:56:55 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 14:56:55 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc
2026-07-15 14:57:00 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 14:57:01 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc
2026-07-15 14:57:06 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 14:57:06 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc
2026-07-15 14:57:12 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 14:57:12 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc
2026-07-15 14:57:18 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 14:57:18 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc
2026-07-15 14:57:23 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 14:57:24 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc
2026-07-15 14:57:29 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 14:57:29 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc
2026-07-15 14:57:35 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 14:57:35 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc
2026-07-15 14:57:40 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 14:57:41 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc
2026-07-15 14:57:47 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 14:57:47 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc
2026-07-15 14:57:54 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 14:57:54 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc
2026-07-15 14:58:00 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 14:58:00 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc
2026-07-15 14:58:05 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 14:58:06 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc
2026-07-15 14:58:11 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 14:58:11 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-15 14:58:11 INFO [TIME] step_end current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 14:58:11 INFO [TIME] step_start current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 14:58:11 INFO [TIME] window start=2020-02-06 11:00:00 end=2020-02-06 13:00:00 run_hours=2 has_assimilation=True
2026-07-15 14:58:11 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:58:13 INFO Hourly dataset computed and listing created
2026-07-15 14:58:19 INFO Hourly dataset computed
2026-07-15 14:58:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:58:20 INFO Hourly dataset computed and listing created
2026-07-15 14:58:23 INFO Hourly dataset computed
2026-07-15 14:58:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:58:24 INFO Hourly dataset computed and listing created
2026-07-15 14:58:26 INFO Hourly dataset computed
2026-07-15 14:58:26 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:58:28 INFO Hourly dataset computed and listing created
2026-07-15 14:58:31 INFO Hourly dataset computed
2026-07-15 14:58:31 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:58:32 INFO Hourly dataset computed and listing created
2026-07-15 14:58:35 INFO Hourly dataset computed
2026-07-15 14:58:35 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:58:36 INFO Hourly dataset computed and listing created
2026-07-15 14:58:39 INFO Hourly dataset computed
2026-07-15 14:58:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:58:40 INFO Hourly dataset computed and listing created
2026-07-15 14:58:43 INFO Hourly dataset computed
2026-07-15 14:58:43 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:58:44 INFO Hourly dataset computed and listing created
2026-07-15 14:58:47 INFO Hourly dataset computed
2026-07-15 14:58:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:58:48 INFO Hourly dataset computed and listing created
2026-07-15 14:58:51 INFO Hourly dataset computed
2026-07-15 14:58:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:58:52 INFO Hourly dataset computed and listing created
2026-07-15 14:58:54 INFO Hourly dataset computed
2026-07-15 14:58:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:58:55 INFO Hourly dataset computed and listing created
2026-07-15 14:58:59 INFO Hourly dataset computed
2026-07-15 14:58:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:59:00 INFO Hourly dataset computed and listing created
2026-07-15 14:59:02 INFO Hourly dataset computed
2026-07-15 14:59:02 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:59:03 INFO Hourly dataset computed and listing created
2026-07-15 14:59:06 INFO Hourly dataset computed
2026-07-15 14:59:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:59:07 INFO Hourly dataset computed and listing created
2026-07-15 14:59:09 INFO Hourly dataset computed
2026-07-15 14:59:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 14:59:10 INFO Hourly dataset computed and listing created
2026-07-15 14:59:12 INFO Hourly dataset computed
2026-07-15 14:59:12 INFO ---------->>> Running CHIMERE model from 2020-02-06 11:00:00 to 2020-02-06 13:00:00
2026-07-15 14:59:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:59:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1
2026-07-15 14:59:12 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc
2026-07-15 14:59:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-15 14:59:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:59:12 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-15 14:59:12 INFO Queuing job for member 1...
2026-07-15 14:59:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:59:12 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-15 14:59:14 INFO Found: ['5164278']
2026-07-15 14:59:19 INFO [TGCC-IRENE] Submitted job with ID:['5164278']
2026-07-15 14:59:19 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:59:19 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2
2026-07-15 14:59:19 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc
2026-07-15 14:59:19 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-15 14:59:19 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:59:19 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-15 14:59:19 INFO Queuing job for member 2...
2026-07-15 14:59:19 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:59:19 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-15 14:59:20 INFO Found: ['5164279']
2026-07-15 14:59:25 INFO [TGCC-IRENE] Submitted job with ID:['5164279']
2026-07-15 14:59:25 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:59:25 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3
2026-07-15 14:59:25 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc
2026-07-15 14:59:25 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-15 14:59:25 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:59:25 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-15 14:59:25 INFO Queuing job for member 3...
2026-07-15 14:59:25 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:59:25 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-15 14:59:26 INFO Found: ['5164280']
2026-07-15 14:59:31 INFO [TGCC-IRENE] Submitted job with ID:['5164280']
2026-07-15 14:59:31 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:59:31 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4
2026-07-15 14:59:31 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc
2026-07-15 14:59:31 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-15 14:59:31 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:59:31 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-15 14:59:31 INFO Queuing job for member 4...
2026-07-15 14:59:31 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:59:31 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-15 14:59:31 INFO Found: ['5164281']
2026-07-15 14:59:36 INFO [TGCC-IRENE] Submitted job with ID:['5164281']
2026-07-15 14:59:36 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:59:36 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5
2026-07-15 14:59:36 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc
2026-07-15 14:59:36 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-15 14:59:36 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:59:36 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-15 14:59:36 INFO Queuing job for member 5...
2026-07-15 14:59:36 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:59:36 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-15 14:59:37 INFO Found: ['5164282']
2026-07-15 14:59:42 INFO [TGCC-IRENE] Submitted job with ID:['5164282']
2026-07-15 14:59:42 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:59:42 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6
2026-07-15 14:59:42 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc
2026-07-15 14:59:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-15 14:59:42 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:59:42 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-15 14:59:42 INFO Queuing job for member 6...
2026-07-15 14:59:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:59:42 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-15 14:59:43 INFO Found: ['5164283']
2026-07-15 14:59:48 INFO [TGCC-IRENE] Submitted job with ID:['5164283']
2026-07-15 14:59:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:59:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7
2026-07-15 14:59:48 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc
2026-07-15 14:59:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-15 14:59:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:59:48 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-15 14:59:48 INFO Queuing job for member 7...
2026-07-15 14:59:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:59:48 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-15 14:59:49 INFO Found: ['5164284']
2026-07-15 14:59:54 INFO [TGCC-IRENE] Submitted job with ID:['5164284']
2026-07-15 14:59:54 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:59:54 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8
2026-07-15 14:59:54 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc
2026-07-15 14:59:54 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-15 14:59:54 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 14:59:54 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-15 14:59:54 INFO Queuing job for member 8...
2026-07-15 14:59:54 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 14:59:54 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-15 14:59:54 INFO Found: ['5164285']
2026-07-15 14:59:59 INFO [TGCC-IRENE] Submitted job with ID:['5164285']
2026-07-15 14:59:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 14:59:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9
2026-07-15 14:59:59 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc
2026-07-15 14:59:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-15 15:00:00 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 15:00:00 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-15 15:00:00 INFO Queuing job for member 9...
2026-07-15 15:00:00 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 15:00:00 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-15 15:00:01 INFO Found: ['5164286']
2026-07-15 15:00:06 INFO [TGCC-IRENE] Submitted job with ID:['5164286']
2026-07-15 15:00:06 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 15:00:06 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10
2026-07-15 15:00:06 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc
2026-07-15 15:00:06 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-15 15:00:06 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 15:00:06 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-15 15:00:06 INFO Queuing job for member 10...
2026-07-15 15:00:06 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 15:00:06 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-15 15:00:09 INFO Found: ['5164288']
2026-07-15 15:00:14 INFO [TGCC-IRENE] Submitted job with ID:['5164288']
2026-07-15 15:00:14 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 15:00:14 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11
2026-07-15 15:00:14 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc
2026-07-15 15:00:14 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-15 15:00:14 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 15:00:14 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-15 15:00:14 INFO Queuing job for member 11...
2026-07-15 15:00:14 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 15:00:14 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-15 15:00:16 INFO Found: ['5164289']
2026-07-15 15:00:21 INFO [TGCC-IRENE] Submitted job with ID:['5164289']
2026-07-15 15:00:21 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 15:00:21 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12
2026-07-15 15:00:21 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc
2026-07-15 15:00:21 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-15 15:00:21 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 15:00:21 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-15 15:00:21 INFO Queuing job for member 12...
2026-07-15 15:00:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 15:00:21 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-15 15:01:52 INFO Found: []
Traceback (most recent call last):
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/main.py", line 103, in <module>
    pipeline.run_pipeline()
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/pipelines/base_pipeline.py", line 138, in run_pipeline
    self.run_model()
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/pipelines/chimere2023/pipeline.py", line 317, in run_model
    job_ids.append(submit_irene(CommandSpec(
                   ^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/orchestrator_utils.py", line 268, in submit_irene
    raise SchedulerError(
pipeline_errors.SchedulerError: Submission command failed: ccc_msub ./submit_p_12.sh (return code 1)
+ exit 0
