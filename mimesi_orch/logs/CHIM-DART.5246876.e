+ SCRIPT_PID=3354747
+ /bin/bash -x /tmp/tmp.YTOF07HRmL
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
2026-07-22 12:17:41 INFO 
███    ███ ██ ███    ███ ███████ ███████ ██
████  ████ ██ ████  ████ ██      ██      ██
██ ████ ██ ██ ██ ████ ██ █████   ███████ ██
██  ██  ██ ██ ██  ██  ██ ██           ██ ██
██      ██ ██ ██      ██ ███████ ███████ ██



2026-07-22 12:17:41 INFO [PIPELINE] =======================================
2026-07-22 12:17:41 INFO [PIPELINE] Starting chimere–DART orchestrator
2026-07-22 12:17:41 INFO [PIPELINE] Config file: config/config_irene_IM.yaml
2026-07-22 12:17:41 INFO [PIPELINE] Run dir: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart
2026-07-22 12:17:41 INFO [PIPELINE] Log file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/mimesi_orchestrator_logs/chimere_DART_20260722_121741.log
2026-07-22 12:17:41 INFO [PIPELINE] =======================================
2026-07-22 12:17:41 INFO Running assimilation with model_type=ModelType.CHIMERE
2026-07-22 12:17:41 INFO Using scheduler=Scheduler.SLURM, queue=rome
2026-07-22 12:17:41 INFO [STEP] ---- TIME LOOP START ----
2026-07-22 12:17:41 INFO [TIME] step_start current_time=2020-02-06 00:00:00 simulated_time=None dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 12:17:41 INFO [TIME] window start=2020-02-06 00:00:00 end=2020-02-06 01:00:00 run_hours=1 has_assimilation=False
2026-07-22 12:17:42 INFO Asked to restart from control run ...
2026-07-22 12:17:42 INFO Copying EMIS ...
2026-07-22 12:17:42 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens01.nc
2026-07-22 12:17:42 INFO Linking first END ...
2026-07-22 12:17:42 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020500_24_ENS1.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-22 12:17:42 INFO >> Checking links...
2026-07-22 12:17:42 INFO >> All links are good for ENS1  ...
2026-07-22 12:17:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 12:17:49 INFO Hourly dataset computed and listing created
2026-07-22 12:17:54 INFO Hourly dataset computed
2026-07-22 12:17:54 INFO Asked to restart from control run ...
2026-07-22 12:17:54 INFO Copying EMIS ...
2026-07-22 12:17:54 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens02.nc
2026-07-22 12:17:54 INFO Linking first END ...
2026-07-22 12:17:54 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020500_24_ENS2.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-22 12:17:54 INFO >> Checking links...
2026-07-22 12:17:54 INFO >> All links are good for ENS2  ...
2026-07-22 12:17:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 12:17:56 INFO Hourly dataset computed and listing created
2026-07-22 12:17:58 INFO Hourly dataset computed
2026-07-22 12:17:58 INFO Asked to restart from control run ...
2026-07-22 12:17:58 INFO Copying EMIS ...
2026-07-22 12:17:59 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens03.nc
2026-07-22 12:17:59 INFO Linking first END ...
2026-07-22 12:17:59 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020500_24_ENS3.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-22 12:17:59 INFO >> Checking links...
2026-07-22 12:17:59 INFO >> All links are good for ENS3  ...
2026-07-22 12:17:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 12:18:00 INFO Hourly dataset computed and listing created
2026-07-22 12:18:00 INFO Hourly dataset computed
2026-07-22 12:18:00 INFO Asked to restart from control run ...
2026-07-22 12:18:00 INFO Copying EMIS ...
2026-07-22 12:18:01 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens04.nc
2026-07-22 12:18:01 INFO Linking first END ...
2026-07-22 12:18:01 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020500_24_ENS4.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-22 12:18:01 INFO >> Checking links...
2026-07-22 12:18:01 INFO >> All links are good for ENS4  ...
2026-07-22 12:18:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 12:18:01 INFO Hourly dataset computed and listing created
2026-07-22 12:18:02 INFO Hourly dataset computed
2026-07-22 12:18:02 INFO Asked to restart from control run ...
2026-07-22 12:18:02 INFO Copying EMIS ...
2026-07-22 12:18:03 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens05.nc
2026-07-22 12:18:03 INFO Linking first END ...
2026-07-22 12:18:03 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020500_24_ENS5.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-22 12:18:03 INFO >> Checking links...
2026-07-22 12:18:03 INFO >> All links are good for ENS5  ...
2026-07-22 12:18:03 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 12:18:04 INFO Hourly dataset computed and listing created
2026-07-22 12:18:04 INFO Hourly dataset computed
2026-07-22 12:18:04 INFO Asked to restart from control run ...
2026-07-22 12:18:04 INFO Copying EMIS ...
2026-07-22 12:18:05 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens06.nc
2026-07-22 12:18:05 INFO Linking first END ...
2026-07-22 12:18:05 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020500_24_ENS6.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-22 12:18:05 INFO >> Checking links...
2026-07-22 12:18:05 INFO >> All links are good for ENS6  ...
2026-07-22 12:18:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 12:18:06 INFO Hourly dataset computed and listing created
2026-07-22 12:18:06 INFO Hourly dataset computed
2026-07-22 12:18:06 INFO Asked to restart from control run ...
2026-07-22 12:18:06 INFO Copying EMIS ...
2026-07-22 12:18:07 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens07.nc
2026-07-22 12:18:07 INFO Linking first END ...
2026-07-22 12:18:07 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020500_24_ENS7.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-22 12:18:07 INFO >> Checking links...
2026-07-22 12:18:07 INFO >> All links are good for ENS7  ...
2026-07-22 12:18:07 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 12:18:08 INFO Hourly dataset computed and listing created
2026-07-22 12:18:08 INFO Hourly dataset computed
2026-07-22 12:18:08 INFO Asked to restart from control run ...
2026-07-22 12:18:08 INFO Copying EMIS ...
2026-07-22 12:18:09 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens08.nc
2026-07-22 12:18:09 INFO Linking first END ...
2026-07-22 12:18:09 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020500_24_ENS8.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-22 12:18:09 INFO >> Checking links...
2026-07-22 12:18:09 INFO >> All links are good for ENS8  ...
2026-07-22 12:18:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 12:18:10 INFO Hourly dataset computed and listing created
2026-07-22 12:18:10 INFO Hourly dataset computed
2026-07-22 12:18:10 INFO Asked to restart from control run ...
2026-07-22 12:18:10 INFO Copying EMIS ...
2026-07-22 12:18:11 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens09.nc
2026-07-22 12:18:11 INFO Linking first END ...
2026-07-22 12:18:11 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020500_24_ENS9.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-22 12:18:11 INFO >> Checking links...
2026-07-22 12:18:11 INFO >> All links are good for ENS9  ...
2026-07-22 12:18:11 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 12:18:12 INFO Hourly dataset computed and listing created
2026-07-22 12:18:12 INFO Hourly dataset computed
2026-07-22 12:18:12 INFO Asked to restart from control run ...
2026-07-22 12:18:12 INFO Copying EMIS ...
2026-07-22 12:18:13 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens10.nc
2026-07-22 12:18:13 INFO Linking first END ...
2026-07-22 12:18:13 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020500_24_ENS10.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-22 12:18:13 INFO >> Checking links...
2026-07-22 12:18:13 INFO >> All links are good for ENS10  ...
2026-07-22 12:18:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 12:18:14 INFO Hourly dataset computed and listing created
2026-07-22 12:18:14 INFO Hourly dataset computed
2026-07-22 12:18:15 INFO Asked to restart from control run ...
2026-07-22 12:18:15 INFO Copying EMIS ...
2026-07-22 12:18:15 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens11.nc
2026-07-22 12:18:15 INFO Linking first END ...
2026-07-22 12:18:15 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020500_24_ENS11.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-22 12:18:15 INFO >> Checking links...
2026-07-22 12:18:15 INFO >> All links are good for ENS11  ...
2026-07-22 12:18:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 12:18:16 INFO Hourly dataset computed and listing created
2026-07-22 12:18:16 INFO Hourly dataset computed
2026-07-22 12:18:17 INFO Asked to restart from control run ...
2026-07-22 12:18:17 INFO Copying EMIS ...
2026-07-22 12:18:17 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens12.nc
2026-07-22 12:18:17 INFO Linking first END ...
2026-07-22 12:18:17 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020500_24_ENS12.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-22 12:18:17 INFO >> Checking links...
2026-07-22 12:18:17 INFO >> All links are good for ENS12  ...
2026-07-22 12:18:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 12:18:18 INFO Hourly dataset computed and listing created
2026-07-22 12:18:19 INFO Hourly dataset computed
2026-07-22 12:18:19 INFO Asked to restart from control run ...
2026-07-22 12:18:19 INFO Copying EMIS ...
2026-07-22 12:18:19 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens13.nc
2026-07-22 12:18:19 INFO Linking first END ...
2026-07-22 12:18:19 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020500_24_ENS13.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-22 12:18:19 INFO >> Checking links...
2026-07-22 12:18:19 INFO >> All links are good for ENS13  ...
2026-07-22 12:18:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 12:18:20 INFO Hourly dataset computed and listing created
2026-07-22 12:18:21 INFO Hourly dataset computed
2026-07-22 12:18:21 INFO Asked to restart from control run ...
2026-07-22 12:18:21 INFO Copying EMIS ...
2026-07-22 12:18:21 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens14.nc
2026-07-22 12:18:21 INFO Linking first END ...
2026-07-22 12:18:21 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020500_24_ENS14.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-22 12:18:21 INFO >> Checking links...
2026-07-22 12:18:21 INFO >> All links are good for ENS14  ...
2026-07-22 12:18:21 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 12:18:22 INFO Hourly dataset computed and listing created
2026-07-22 12:18:22 INFO Hourly dataset computed
2026-07-22 12:18:23 INFO Asked to restart from control run ...
2026-07-22 12:18:23 INFO Copying EMIS ...
2026-07-22 12:18:23 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens15.nc
2026-07-22 12:18:23 INFO Linking first END ...
2026-07-22 12:18:23 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020500_24_ENS15.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-22 12:18:23 INFO >> Checking links...
2026-07-22 12:18:23 INFO >> All links are good for ENS15  ...
2026-07-22 12:18:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 12:18:24 INFO Hourly dataset computed and listing created
2026-07-22 12:18:25 INFO Hourly dataset computed
2026-07-22 12:18:25 INFO ---------->>> Running CHIMERE model from 2020-02-06 00:00:00 to 2020-02-06 01:00:00
2026-07-22 12:18:25 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 12:18:25 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1
2026-07-22 12:18:25 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020500_24_ENS1.nc
2026-07-22 12:18:25 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-22 12:18:25 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 12:18:25 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-22 12:18:25 INFO Queuing job for member 1...
2026-07-22 12:18:25 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 12:18:25 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-22 12:18:27 INFO Found: ['5246944']
2026-07-22 12:18:32 INFO [TGCC-IRENE] Submitted job with ID:['5246944']
2026-07-22 12:18:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 12:18:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2
2026-07-22 12:18:32 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020500_24_ENS2.nc
2026-07-22 12:18:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-22 12:18:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 12:18:32 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-22 12:18:32 INFO Queuing job for member 2...
2026-07-22 12:18:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 12:18:32 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-22 12:18:32 INFO Found: ['5246945']
2026-07-22 12:18:37 INFO [TGCC-IRENE] Submitted job with ID:['5246945']
2026-07-22 12:18:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 12:18:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3
2026-07-22 12:18:37 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020500_24_ENS3.nc
2026-07-22 12:18:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-22 12:18:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 12:18:37 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-22 12:18:37 INFO Queuing job for member 3...
2026-07-22 12:18:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 12:18:37 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-22 12:18:38 INFO Found: ['5246947']
2026-07-22 12:18:43 INFO [TGCC-IRENE] Submitted job with ID:['5246947']
2026-07-22 12:18:43 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 12:18:43 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4
2026-07-22 12:18:43 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020500_24_ENS4.nc
2026-07-22 12:18:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-22 12:18:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 12:18:43 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-22 12:18:43 INFO Queuing job for member 4...
2026-07-22 12:18:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 12:18:43 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-22 12:18:44 INFO Found: ['5246948']
2026-07-22 12:18:49 INFO [TGCC-IRENE] Submitted job with ID:['5246948']
2026-07-22 12:18:49 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 12:18:49 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5
2026-07-22 12:18:49 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020500_24_ENS5.nc
2026-07-22 12:18:49 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-22 12:18:49 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 12:18:49 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-22 12:18:49 INFO Queuing job for member 5...
2026-07-22 12:18:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 12:18:49 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-22 12:18:50 INFO Found: ['5246949']
2026-07-22 12:18:55 INFO [TGCC-IRENE] Submitted job with ID:['5246949']
2026-07-22 12:18:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 12:18:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6
2026-07-22 12:18:55 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020500_24_ENS6.nc
2026-07-22 12:18:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-22 12:18:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 12:18:55 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-22 12:18:55 INFO Queuing job for member 6...
2026-07-22 12:18:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 12:18:55 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-22 12:18:55 INFO Found: ['5246950']
2026-07-22 12:19:00 INFO [TGCC-IRENE] Submitted job with ID:['5246950']
2026-07-22 12:19:00 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 12:19:00 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7
2026-07-22 12:19:00 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020500_24_ENS7.nc
2026-07-22 12:19:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-22 12:19:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 12:19:01 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-22 12:19:01 INFO Queuing job for member 7...
2026-07-22 12:19:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 12:19:01 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-22 12:19:01 INFO Found: ['5246952']
2026-07-22 12:19:06 INFO [TGCC-IRENE] Submitted job with ID:['5246952']
2026-07-22 12:19:06 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 12:19:06 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8
2026-07-22 12:19:06 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020500_24_ENS8.nc
2026-07-22 12:19:06 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-22 12:19:06 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 12:19:06 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-22 12:19:06 INFO Queuing job for member 8...
2026-07-22 12:19:06 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 12:19:06 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-22 12:19:07 INFO Found: ['5246953']
2026-07-22 12:19:12 INFO [TGCC-IRENE] Submitted job with ID:['5246953']
2026-07-22 12:19:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 12:19:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9
2026-07-22 12:19:12 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020500_24_ENS9.nc
2026-07-22 12:19:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-22 12:19:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 12:19:12 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-22 12:19:12 INFO Queuing job for member 9...
2026-07-22 12:19:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 12:19:12 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-22 12:19:13 INFO Found: ['5246955']
2026-07-22 12:19:18 INFO [TGCC-IRENE] Submitted job with ID:['5246955']
2026-07-22 12:19:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 12:19:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10
2026-07-22 12:19:18 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020500_24_ENS10.nc
2026-07-22 12:19:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-22 12:19:18 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 12:19:18 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-22 12:19:18 INFO Queuing job for member 10...
2026-07-22 12:19:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 12:19:18 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-22 12:19:19 INFO Found: ['5246959']
2026-07-22 12:19:24 INFO [TGCC-IRENE] Submitted job with ID:['5246959']
2026-07-22 12:19:24 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 12:19:24 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11
2026-07-22 12:19:24 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020500_24_ENS11.nc
2026-07-22 12:19:24 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-22 12:19:24 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 12:19:24 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-22 12:19:24 INFO Queuing job for member 11...
2026-07-22 12:19:24 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 12:19:24 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-22 12:19:27 INFO Found: ['5246963']
2026-07-22 12:19:32 INFO [TGCC-IRENE] Submitted job with ID:['5246963']
2026-07-22 12:19:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 12:19:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12
2026-07-22 12:19:32 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020500_24_ENS12.nc
2026-07-22 12:19:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-22 12:19:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 12:19:32 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-22 12:19:32 INFO Queuing job for member 12...
2026-07-22 12:19:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 12:19:32 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-22 12:19:32 INFO Found: ['5246968']
2026-07-22 12:19:37 INFO [TGCC-IRENE] Submitted job with ID:['5246968']
2026-07-22 12:19:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 12:19:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13
2026-07-22 12:19:37 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020500_24_ENS13.nc
2026-07-22 12:19:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-22 12:19:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 12:19:37 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-22 12:19:37 INFO Queuing job for member 13...
2026-07-22 12:19:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 12:19:37 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-22 12:19:38 INFO Found: ['5246969']
2026-07-22 12:19:43 INFO [TGCC-IRENE] Submitted job with ID:['5246969']
2026-07-22 12:19:43 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 12:19:43 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14
2026-07-22 12:19:43 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020500_24_ENS14.nc
2026-07-22 12:19:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-22 12:19:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 12:19:43 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-22 12:19:43 INFO Queuing job for member 14...
2026-07-22 12:19:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 12:19:43 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-22 12:19:44 INFO Found: ['5246975']
2026-07-22 12:19:49 INFO [TGCC-IRENE] Submitted job with ID:['5246975']
2026-07-22 12:19:49 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 12:19:49 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15
2026-07-22 12:19:49 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020500_24_ENS15.nc
2026-07-22 12:19:49 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-22 12:19:49 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 12:19:49 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-22 12:19:49 INFO Queuing job for member 15...
2026-07-22 12:19:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 12:19:49 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-22 12:19:50 INFO Found: ['5246985']
2026-07-22 12:19:55 INFO [TGCC-IRENE] Submitted job with ID:['5246985']
2026-07-22 12:19:55 INFO Checking job status ...
2026-07-22 12:19:55 INFO None 5246944: status RUNNING/PENDING
2026-07-22 12:19:55 INFO None 5246945: status RUNNING/PENDING
2026-07-22 12:19:55 INFO None 5246947: status RUNNING/PENDING
2026-07-22 12:19:55 INFO None 5246948: status RUNNING/PENDING
2026-07-22 12:19:55 INFO None 5246949: status RUNNING/PENDING
2026-07-22 12:19:55 INFO None 5246950: status RUNNING/PENDING
2026-07-22 12:19:55 INFO None 5246952: status RUNNING/PENDING
2026-07-22 12:19:55 INFO None 5246953: status RUNNING/PENDING
2026-07-22 12:19:55 INFO None 5246955: status RUNNING/PENDING
2026-07-22 12:19:55 INFO None 5246959: status RUNNING/PENDING
2026-07-22 12:19:55 INFO None 5246963: status RUNNING/PENDING
2026-07-22 12:19:55 INFO None 5246968: status RUNNING/PENDING
2026-07-22 12:19:55 INFO None 5246969: status RUNNING/PENDING
2026-07-22 12:19:55 INFO None 5246975: status RUNNING/PENDING
2026-07-22 12:19:55 INFO None 5246985: status RUNNING/PENDING
2026-07-22 12:19:55 INFO Jobs still running: ['5246944', '5246945', '5246947', '5246948', '5246949', '5246950', '5246952', '5246953', '5246955', '5246959', '5246963', '5246968', '5246969', '5246975', '5246985']. Waiting...
2026-07-22 12:20:10 INFO None 5246944: status RUNNING/PENDING
2026-07-22 12:20:10 INFO None 5246945: status RUNNING/PENDING
2026-07-22 12:20:10 INFO None 5246947: status RUNNING/PENDING
2026-07-22 12:20:10 INFO None 5246948: status RUNNING/PENDING
2026-07-22 12:20:10 INFO None 5246949: status RUNNING/PENDING
2026-07-22 12:20:10 INFO None 5246950: status RUNNING/PENDING
2026-07-22 12:20:10 INFO None 5246952: status RUNNING/PENDING
2026-07-22 12:20:10 INFO None 5246953: status RUNNING/PENDING
2026-07-22 12:20:10 INFO None 5246955: status RUNNING/PENDING
2026-07-22 12:20:10 INFO None 5246959: status RUNNING/PENDING
2026-07-22 12:20:10 INFO None 5246963: status RUNNING/PENDING
2026-07-22 12:20:10 INFO None 5246968: status RUNNING/PENDING
2026-07-22 12:20:10 INFO None 5246969: status RUNNING/PENDING
2026-07-22 12:20:10 INFO None 5246975: status RUNNING/PENDING
2026-07-22 12:20:10 INFO None 5246985: status RUNNING/PENDING
2026-07-22 12:20:10 INFO Jobs still running: ['5246944', '5246945', '5246947', '5246948', '5246949', '5246950', '5246952', '5246953', '5246955', '5246959', '5246963', '5246968', '5246969', '5246975', '5246985']. Waiting...
2026-07-22 12:20:27 INFO None 5246944: status RUNNING/PENDING
2026-07-22 12:20:27 INFO None 5246945: status RUNNING/PENDING
2026-07-22 12:20:27 INFO None 5246947: status RUNNING/PENDING
2026-07-22 12:20:27 INFO None 5246948: status RUNNING/PENDING
2026-07-22 12:20:27 INFO None 5246949: status RUNNING/PENDING
2026-07-22 12:20:27 INFO None 5246950: status RUNNING/PENDING
2026-07-22 12:20:27 INFO None 5246952: status RUNNING/PENDING
2026-07-22 12:20:27 INFO None 5246953: status RUNNING/PENDING
2026-07-22 12:20:27 INFO None 5246955: status RUNNING/PENDING
2026-07-22 12:20:27 INFO None 5246959: status RUNNING/PENDING
2026-07-22 12:20:27 INFO None 5246963: status RUNNING/PENDING
2026-07-22 12:20:27 INFO None 5246968: status RUNNING/PENDING
2026-07-22 12:20:27 INFO None 5246969: status RUNNING/PENDING
2026-07-22 12:20:27 INFO None 5246975: status RUNNING/PENDING
2026-07-22 12:20:27 INFO None 5246985: status RUNNING/PENDING
2026-07-22 12:20:27 INFO Jobs still running: ['5246944', '5246945', '5246947', '5246948', '5246949', '5246950', '5246952', '5246953', '5246955', '5246959', '5246963', '5246968', '5246969', '5246975', '5246985']. Waiting...
2026-07-22 12:20:42 INFO None 5246944: status RUNNING/PENDING
2026-07-22 12:20:42 INFO None 5246945: status RUNNING/PENDING
2026-07-22 12:20:42 INFO None 5246947: status RUNNING/PENDING
2026-07-22 12:20:42 INFO None 5246948: status RUNNING/PENDING
2026-07-22 12:20:42 INFO None 5246949: status RUNNING/PENDING
2026-07-22 12:20:42 INFO None 5246950: status RUNNING/PENDING
2026-07-22 12:20:42 INFO None 5246952: status RUNNING/PENDING
2026-07-22 12:20:42 INFO None 5246953: status RUNNING/PENDING
2026-07-22 12:20:42 INFO None 5246955: status RUNNING/PENDING
2026-07-22 12:20:42 INFO None 5246959: status RUNNING/PENDING
2026-07-22 12:20:42 INFO None 5246963: status RUNNING/PENDING
2026-07-22 12:20:42 INFO None 5246968: status RUNNING/PENDING
2026-07-22 12:20:42 INFO None 5246969: status RUNNING/PENDING
2026-07-22 12:20:42 INFO None 5246975: status RUNNING/PENDING
2026-07-22 12:20:42 INFO None 5246985: status RUNNING/PENDING
2026-07-22 12:20:42 INFO Jobs still running: ['5246944', '5246945', '5246947', '5246948', '5246949', '5246950', '5246952', '5246953', '5246955', '5246959', '5246963', '5246968', '5246969', '5246975', '5246985']. Waiting...
2026-07-22 12:20:57 INFO None 5246944: status RUNNING/PENDING
2026-07-22 12:20:57 INFO None 5246945: status RUNNING/PENDING
2026-07-22 12:20:57 INFO None 5246947: status RUNNING/PENDING
2026-07-22 12:20:57 INFO None 5246948: status RUNNING/PENDING
2026-07-22 12:20:58 INFO None 5246949: status RUNNING/PENDING
2026-07-22 12:20:58 INFO None 5246950: status RUNNING/PENDING
2026-07-22 12:20:58 INFO None 5246952: status RUNNING/PENDING
2026-07-22 12:20:58 INFO None 5246953: status RUNNING/PENDING
2026-07-22 12:20:58 INFO None 5246955: status RUNNING/PENDING
2026-07-22 12:20:58 INFO None 5246959: status RUNNING/PENDING
2026-07-22 12:20:58 INFO None 5246963: status RUNNING/PENDING
2026-07-22 12:20:58 INFO None 5246968: status RUNNING/PENDING
2026-07-22 12:20:58 INFO None 5246969: status RUNNING/PENDING
2026-07-22 12:20:58 INFO None 5246975: status RUNNING/PENDING
2026-07-22 12:20:58 INFO None 5246985: status RUNNING/PENDING
2026-07-22 12:20:58 INFO Jobs still running: ['5246944', '5246945', '5246947', '5246948', '5246949', '5246950', '5246952', '5246953', '5246955', '5246959', '5246963', '5246968', '5246969', '5246975', '5246985']. Waiting...
2026-07-22 12:21:13 INFO None 5246944: status RUNNING/PENDING
2026-07-22 12:21:13 INFO None 5246945: status RUNNING/PENDING
2026-07-22 12:21:13 INFO None 5246947: status RUNNING/PENDING
2026-07-22 12:21:13 INFO None 5246948: status RUNNING/PENDING
2026-07-22 12:21:13 INFO None 5246949: status RUNNING/PENDING
2026-07-22 12:21:13 INFO None 5246950: status RUNNING/PENDING
2026-07-22 12:21:13 INFO None 5246952: status RUNNING/PENDING
2026-07-22 12:21:13 INFO None 5246953: status RUNNING/PENDING
2026-07-22 12:21:13 INFO None 5246955: status RUNNING/PENDING
2026-07-22 12:21:13 INFO None 5246959: status RUNNING/PENDING
2026-07-22 12:21:13 INFO None 5246963: status RUNNING/PENDING
2026-07-22 12:21:13 INFO None 5246968: status RUNNING/PENDING
2026-07-22 12:21:13 INFO None 5246969: status RUNNING/PENDING
2026-07-22 12:21:13 INFO None 5246975: status RUNNING/PENDING
2026-07-22 12:21:13 INFO None 5246985: status RUNNING/PENDING
2026-07-22 12:21:13 INFO Jobs still running: ['5246944', '5246945', '5246947', '5246948', '5246949', '5246950', '5246952', '5246953', '5246955', '5246959', '5246963', '5246968', '5246969', '5246975', '5246985']. Waiting...
2026-07-22 12:21:30 INFO None 5246944: status FINISHED
2026-07-22 12:21:30 INFO None 5246945: status RUNNING/PENDING
2026-07-22 12:21:30 INFO None 5246947: status RUNNING/PENDING
2026-07-22 12:21:30 INFO None 5246948: status RUNNING/PENDING
2026-07-22 12:21:30 INFO None 5246949: status RUNNING/PENDING
2026-07-22 12:21:30 INFO None 5246950: status RUNNING/PENDING
2026-07-22 12:21:30 INFO None 5246952: status RUNNING/PENDING
2026-07-22 12:21:30 INFO None 5246953: status FINISHED
2026-07-22 12:21:30 INFO None 5246955: status FINISHED
2026-07-22 12:21:30 INFO None 5246959: status RUNNING/PENDING
2026-07-22 12:21:30 INFO None 5246963: status RUNNING/PENDING
2026-07-22 12:21:30 INFO None 5246968: status RUNNING/PENDING
2026-07-22 12:21:30 INFO None 5246969: status RUNNING/PENDING
2026-07-22 12:21:30 INFO None 5246975: status RUNNING/PENDING
2026-07-22 12:21:30 INFO None 5246985: status RUNNING/PENDING
2026-07-22 12:21:30 INFO Jobs still running: ['5246945', '5246947', '5246948', '5246949', '5246950', '5246952', '5246959', '5246963', '5246968', '5246969', '5246975', '5246985']. Waiting...
2026-07-22 12:21:45 INFO None 5246944: status FINISHED
2026-07-22 12:21:45 INFO None 5246945: status FINISHED
2026-07-22 12:21:45 INFO None 5246947: status RUNNING/PENDING
2026-07-22 12:21:45 INFO None 5246948: status RUNNING/PENDING
2026-07-22 12:21:45 INFO None 5246949: status RUNNING/PENDING
2026-07-22 12:21:45 INFO None 5246950: status RUNNING/PENDING
2026-07-22 12:21:45 INFO None 5246952: status RUNNING/PENDING
2026-07-22 12:21:45 INFO None 5246953: status FINISHED
2026-07-22 12:21:45 INFO None 5246955: status FINISHED
2026-07-22 12:21:45 INFO None 5246959: status RUNNING/PENDING
2026-07-22 12:21:45 INFO None 5246963: status RUNNING/PENDING
2026-07-22 12:21:45 INFO None 5246968: status RUNNING/PENDING
2026-07-22 12:21:45 INFO None 5246969: status RUNNING/PENDING
2026-07-22 12:21:45 INFO None 5246975: status RUNNING/PENDING
2026-07-22 12:21:45 INFO None 5246985: status RUNNING/PENDING
2026-07-22 12:21:45 INFO Jobs still running: ['5246947', '5246948', '5246949', '5246950', '5246952', '5246959', '5246963', '5246968', '5246969', '5246975', '5246985']. Waiting...
2026-07-22 12:22:00 INFO None 5246944: status FINISHED
2026-07-22 12:22:00 INFO None 5246945: status FINISHED
2026-07-22 12:22:00 INFO None 5246947: status RUNNING/PENDING
2026-07-22 12:22:00 INFO None 5246948: status RUNNING/PENDING
2026-07-22 12:22:00 INFO None 5246949: status RUNNING/PENDING
2026-07-22 12:22:00 INFO None 5246950: status RUNNING/PENDING
2026-07-22 12:22:00 INFO None 5246952: status RUNNING/PENDING
2026-07-22 12:22:00 INFO None 5246953: status FINISHED
2026-07-22 12:22:00 INFO None 5246955: status FINISHED
2026-07-22 12:22:00 INFO None 5246959: status RUNNING/PENDING
2026-07-22 12:22:00 INFO None 5246963: status RUNNING/PENDING
2026-07-22 12:22:00 INFO None 5246968: status RUNNING/PENDING
2026-07-22 12:22:00 INFO None 5246969: status RUNNING/PENDING
2026-07-22 12:22:00 INFO None 5246975: status RUNNING/PENDING
2026-07-22 12:22:00 INFO None 5246985: status RUNNING/PENDING
2026-07-22 12:22:00 INFO Jobs still running: ['5246947', '5246948', '5246949', '5246950', '5246952', '5246959', '5246963', '5246968', '5246969', '5246975', '5246985']. Waiting...
2026-07-22 12:22:15 INFO None 5246944: status FINISHED
2026-07-22 12:22:15 INFO None 5246945: status FINISHED
2026-07-22 12:22:16 INFO None 5246947: status RUNNING/PENDING
2026-07-22 12:22:16 INFO None 5246948: status RUNNING/PENDING
2026-07-22 12:22:16 INFO None 5246949: status RUNNING/PENDING
2026-07-22 12:22:16 INFO None 5246950: status RUNNING/PENDING
2026-07-22 12:22:16 INFO None 5246952: status RUNNING/PENDING
2026-07-22 12:22:16 INFO None 5246953: status FINISHED
2026-07-22 12:22:16 INFO None 5246955: status FINISHED
2026-07-22 12:22:16 INFO None 5246959: status RUNNING/PENDING
2026-07-22 12:22:16 INFO None 5246963: status RUNNING/PENDING
2026-07-22 12:22:16 INFO None 5246968: status RUNNING/PENDING
2026-07-22 12:22:16 INFO None 5246969: status RUNNING/PENDING
2026-07-22 12:22:16 INFO None 5246975: status RUNNING/PENDING
2026-07-22 12:22:16 INFO None 5246985: status RUNNING/PENDING
2026-07-22 12:22:16 INFO Jobs still running: ['5246947', '5246948', '5246949', '5246950', '5246952', '5246959', '5246963', '5246968', '5246969', '5246975', '5246985']. Waiting...
2026-07-22 12:22:31 INFO None 5246944: status FINISHED
2026-07-22 12:22:32 INFO None 5246945: status FINISHED
2026-07-22 12:22:32 INFO None 5246947: status RUNNING/PENDING
2026-07-22 12:22:32 INFO None 5246948: status RUNNING/PENDING
2026-07-22 12:22:32 INFO None 5246949: status RUNNING/PENDING
2026-07-22 12:22:32 INFO None 5246950: status RUNNING/PENDING
2026-07-22 12:22:32 INFO None 5246952: status RUNNING/PENDING
2026-07-22 12:22:32 INFO None 5246953: status FINISHED
2026-07-22 12:22:32 INFO None 5246955: status FINISHED
2026-07-22 12:22:32 INFO None 5246959: status RUNNING/PENDING
2026-07-22 12:22:32 INFO None 5246963: status RUNNING/PENDING
2026-07-22 12:22:32 INFO None 5246968: status RUNNING/PENDING
2026-07-22 12:22:32 INFO None 5246969: status RUNNING/PENDING
2026-07-22 12:22:32 INFO None 5246975: status RUNNING/PENDING
2026-07-22 12:22:32 INFO None 5246985: status RUNNING/PENDING
2026-07-22 12:22:32 INFO Jobs still running: ['5246947', '5246948', '5246949', '5246950', '5246952', '5246959', '5246963', '5246968', '5246969', '5246975', '5246985']. Waiting...
2026-07-22 12:22:47 INFO None 5246944: status FINISHED
2026-07-22 12:22:47 INFO None 5246945: status FINISHED
2026-07-22 12:22:47 INFO None 5246947: status RUNNING/PENDING
2026-07-22 12:22:47 INFO None 5246948: status RUNNING/PENDING
2026-07-22 12:22:47 INFO None 5246949: status RUNNING/PENDING
2026-07-22 12:22:47 INFO None 5246950: status RUNNING/PENDING
2026-07-22 12:22:47 INFO None 5246952: status RUNNING/PENDING
2026-07-22 12:22:47 INFO None 5246953: status FINISHED
2026-07-22 12:22:47 INFO None 5246955: status FINISHED
2026-07-22 12:22:47 INFO None 5246959: status RUNNING/PENDING
2026-07-22 12:22:47 INFO None 5246963: status RUNNING/PENDING
2026-07-22 12:22:47 INFO None 5246968: status RUNNING/PENDING
2026-07-22 12:22:47 INFO None 5246969: status RUNNING/PENDING
2026-07-22 12:22:47 INFO None 5246975: status RUNNING/PENDING
2026-07-22 12:22:47 INFO None 5246985: status RUNNING/PENDING
2026-07-22 12:22:47 INFO Jobs still running: ['5246947', '5246948', '5246949', '5246950', '5246952', '5246959', '5246963', '5246968', '5246969', '5246975', '5246985']. Waiting...
2026-07-22 12:23:02 INFO None 5246944: status FINISHED
2026-07-22 12:23:02 INFO None 5246945: status FINISHED
2026-07-22 12:23:02 INFO None 5246947: status RUNNING/PENDING
2026-07-22 12:23:02 INFO None 5246948: status FINISHED
2026-07-22 12:23:02 INFO None 5246949: status FINISHED
2026-07-22 12:23:02 INFO None 5246950: status RUNNING/PENDING
2026-07-22 12:23:02 INFO None 5246952: status RUNNING/PENDING
2026-07-22 12:23:02 INFO None 5246953: status FINISHED
2026-07-22 12:23:02 INFO None 5246955: status FINISHED
2026-07-22 12:23:02 INFO None 5246959: status RUNNING/PENDING
2026-07-22 12:23:02 INFO None 5246963: status RUNNING/PENDING
2026-07-22 12:23:02 INFO None 5246968: status FINISHED
2026-07-22 12:23:02 INFO None 5246969: status RUNNING/PENDING
2026-07-22 12:23:02 INFO None 5246975: status RUNNING/PENDING
2026-07-22 12:23:02 INFO None 5246985: status RUNNING/PENDING
2026-07-22 12:23:02 INFO Jobs still running: ['5246947', '5246950', '5246952', '5246959', '5246963', '5246969', '5246975', '5246985']. Waiting...
2026-07-22 12:23:17 INFO None 5246944: status FINISHED
2026-07-22 12:23:17 INFO None 5246945: status FINISHED
2026-07-22 12:23:17 INFO None 5246947: status FINISHED
2026-07-22 12:23:18 INFO None 5246948: status FINISHED
2026-07-22 12:23:18 INFO None 5246949: status FINISHED
2026-07-22 12:23:18 INFO None 5246950: status FINISHED
2026-07-22 12:23:18 INFO None 5246952: status FINISHED
2026-07-22 12:23:18 INFO None 5246953: status FINISHED
2026-07-22 12:23:18 INFO None 5246955: status FINISHED
2026-07-22 12:23:18 INFO None 5246959: status RUNNING/PENDING
2026-07-22 12:23:18 INFO None 5246963: status RUNNING/PENDING
2026-07-22 12:23:18 INFO None 5246968: status FINISHED
2026-07-22 12:23:18 INFO None 5246969: status RUNNING/PENDING
2026-07-22 12:23:18 INFO None 5246975: status RUNNING/PENDING
2026-07-22 12:23:18 INFO None 5246985: status RUNNING/PENDING
2026-07-22 12:23:18 INFO Jobs still running: ['5246959', '5246963', '5246969', '5246975', '5246985']. Waiting...
2026-07-22 12:23:33 INFO None 5246944: status FINISHED
2026-07-22 12:23:33 INFO None 5246945: status FINISHED
2026-07-22 12:23:33 INFO None 5246947: status FINISHED
2026-07-22 12:23:33 INFO None 5246948: status FINISHED
2026-07-22 12:23:33 INFO None 5246949: status FINISHED
2026-07-22 12:23:34 INFO None 5246950: status FINISHED
2026-07-22 12:23:34 INFO None 5246952: status FINISHED
2026-07-22 12:23:34 INFO None 5246953: status FINISHED
2026-07-22 12:23:34 INFO None 5246955: status FINISHED
2026-07-22 12:23:34 INFO None 5246959: status RUNNING/PENDING
2026-07-22 12:23:34 INFO None 5246963: status RUNNING/PENDING
2026-07-22 12:23:34 INFO None 5246968: status FINISHED
2026-07-22 12:23:34 INFO None 5246969: status RUNNING/PENDING
2026-07-22 12:23:34 INFO None 5246975: status RUNNING/PENDING
2026-07-22 12:23:34 INFO None 5246985: status RUNNING/PENDING
2026-07-22 12:23:34 INFO Jobs still running: ['5246959', '5246963', '5246969', '5246975', '5246985']. Waiting...
2026-07-22 12:23:49 INFO None 5246944: status FINISHED
2026-07-22 12:23:49 INFO None 5246945: status FINISHED
2026-07-22 12:23:49 INFO None 5246947: status FINISHED
2026-07-22 12:23:49 INFO None 5246948: status FINISHED
2026-07-22 12:23:49 INFO None 5246949: status FINISHED
2026-07-22 12:23:49 INFO None 5246950: status FINISHED
2026-07-22 12:23:49 INFO None 5246952: status FINISHED
2026-07-22 12:23:49 INFO None 5246953: status FINISHED
2026-07-22 12:23:49 INFO None 5246955: status FINISHED
2026-07-22 12:23:49 INFO None 5246959: status FINISHED
2026-07-22 12:23:49 INFO None 5246963: status RUNNING/PENDING
2026-07-22 12:23:49 INFO None 5246968: status FINISHED
2026-07-22 12:23:49 INFO None 5246969: status FINISHED
2026-07-22 12:23:49 INFO None 5246975: status RUNNING/PENDING
2026-07-22 12:23:49 INFO None 5246985: status RUNNING/PENDING
2026-07-22 12:23:49 INFO Jobs still running: ['5246963', '5246975', '5246985']. Waiting...
2026-07-22 12:24:04 INFO None 5246944: status FINISHED
2026-07-22 12:24:04 INFO None 5246945: status FINISHED
2026-07-22 12:24:04 INFO None 5246947: status FINISHED
2026-07-22 12:24:04 INFO None 5246948: status FINISHED
2026-07-22 12:24:04 INFO None 5246949: status FINISHED
2026-07-22 12:24:04 INFO None 5246950: status FINISHED
2026-07-22 12:24:04 INFO None 5246952: status FINISHED
2026-07-22 12:24:04 INFO None 5246953: status FINISHED
2026-07-22 12:24:04 INFO None 5246955: status FINISHED
2026-07-22 12:24:04 INFO None 5246959: status FINISHED
2026-07-22 12:24:04 INFO None 5246963: status FINISHED
2026-07-22 12:24:04 INFO None 5246968: status FINISHED
2026-07-22 12:24:04 INFO None 5246969: status FINISHED
2026-07-22 12:24:04 INFO None 5246975: status FINISHED
2026-07-22 12:24:04 INFO None 5246985: status FINISHED
2026-07-22 12:24:04 INFO Jobs ['5246944', '5246945', '5246947', '5246948', '5246949', '5246950', '5246952', '5246953', '5246955', '5246959', '5246963', '5246968', '5246969', '5246975', '5246985'] have finished
2026-07-22 12:24:04 INFO Checking restart files were created ...
2026-07-22 12:24:04 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020600_1_ENS1.nc(668832435 bytes)
2026-07-22 12:24:04 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020600_1_ENS2.nc(668832435 bytes)
2026-07-22 12:24:04 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020600_1_ENS3.nc(668832435 bytes)
2026-07-22 12:24:04 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020600_1_ENS4.nc(668832435 bytes)
2026-07-22 12:24:04 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020600_1_ENS5.nc(668832435 bytes)
2026-07-22 12:24:04 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020600_1_ENS6.nc(668832435 bytes)
2026-07-22 12:24:04 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020600_1_ENS7.nc(668832435 bytes)
2026-07-22 12:24:04 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020600_1_ENS8.nc(668832435 bytes)
2026-07-22 12:24:04 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020600_1_ENS9.nc(668832435 bytes)
2026-07-22 12:24:04 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020600_1_ENS10.nc(668832435 bytes)
2026-07-22 12:24:04 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020600_1_ENS11.nc(668832435 bytes)
2026-07-22 12:24:04 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020600_1_ENS12.nc(668832435 bytes)
2026-07-22 12:24:04 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020600_1_ENS13.nc(668832435 bytes)
2026-07-22 12:24:04 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020600_1_ENS14.nc(668832435 bytes)
2026-07-22 12:24:04 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020600_1_ENS15.nc(668832435 bytes)
2026-07-22 12:24:04 INFO  Run_model() completed successfully.
2026-07-22 12:24:04 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 00:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 12:24:04 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 01:00:00 days=153072 seconds=3600
2026-07-22 12:24:04 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-22 12:24:04 INFO [TIME] increment current_time 2020-02-06 00:00:00 -> 2020-02-06 01:00:00
2026-07-22 12:24:04 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 12:24:04 INFO ---------->>> Running process_satellite_data()
2026-07-22 12:24:04 INFO [DART] No satellite data found, skipping assimilation
2026-07-22 12:24:04 INFO after_assimilation() skipped
2026-07-22 12:24:04 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-22 12:24:04 INFO [TIME] step_end current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 12:24:04 INFO [TIME] step_start current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 12:24:04 INFO [TIME] window start=2020-02-06 01:00:00 end=2020-02-06 09:00:00 run_hours=8 has_assimilation=True
2026-07-22 12:24:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 12:24:06 INFO Hourly dataset computed and listing created
2026-07-22 12:24:18 INFO Hourly dataset computed
2026-07-22 12:24:18 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 12:24:20 INFO Hourly dataset computed and listing created
2026-07-22 12:24:33 INFO Hourly dataset computed
2026-07-22 12:24:33 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 12:24:35 INFO Hourly dataset computed and listing created
2026-07-22 12:24:54 INFO Hourly dataset computed
2026-07-22 12:24:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 12:24:55 INFO Hourly dataset computed and listing created
2026-07-22 12:25:12 INFO Hourly dataset computed
2026-07-22 12:25:12 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 12:25:14 INFO Hourly dataset computed and listing created
2026-07-22 12:25:35 INFO Hourly dataset computed
2026-07-22 12:25:35 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 12:25:37 INFO Hourly dataset computed and listing created
2026-07-22 12:25:55 INFO Hourly dataset computed
2026-07-22 12:25:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 12:25:56 INFO Hourly dataset computed and listing created
2026-07-22 12:26:15 INFO Hourly dataset computed
2026-07-22 12:26:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 12:26:16 INFO Hourly dataset computed and listing created
2026-07-22 12:26:32 INFO Hourly dataset computed
2026-07-22 12:26:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 12:26:34 INFO Hourly dataset computed and listing created
2026-07-22 12:26:50 INFO Hourly dataset computed
2026-07-22 12:26:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 12:26:52 INFO Hourly dataset computed and listing created
2026-07-22 12:27:10 INFO Hourly dataset computed
2026-07-22 12:27:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 12:27:12 INFO Hourly dataset computed and listing created
2026-07-22 12:27:28 INFO Hourly dataset computed
2026-07-22 12:27:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 12:27:29 INFO Hourly dataset computed and listing created
2026-07-22 12:27:47 INFO Hourly dataset computed
2026-07-22 12:27:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 12:27:49 INFO Hourly dataset computed and listing created
2026-07-22 12:28:04 INFO Hourly dataset computed
2026-07-22 12:28:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 12:28:06 INFO Hourly dataset computed and listing created
2026-07-22 12:28:22 INFO Hourly dataset computed
2026-07-22 12:28:22 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 12:28:23 INFO Hourly dataset computed and listing created
2026-07-22 12:28:34 INFO Hourly dataset computed
2026-07-22 12:28:34 INFO ---------->>> Running CHIMERE model from 2020-02-06 01:00:00 to 2020-02-06 09:00:00
2026-07-22 12:28:34 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 12:28:34 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1
2026-07-22 12:28:34 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020600_1_ENS1.nc
2026-07-22 12:28:34 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-22 12:28:34 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 12:28:34 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-22 12:28:34 INFO Queuing job for member 1...
2026-07-22 12:28:34 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 12:28:34 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-22 12:28:35 INFO Found: ['5247397']
2026-07-22 12:28:40 INFO [TGCC-IRENE] Submitted job with ID:['5247397']
2026-07-22 12:28:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 12:28:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2
2026-07-22 12:28:40 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020600_1_ENS2.nc
2026-07-22 12:28:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-22 12:28:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 12:28:40 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-22 12:28:40 INFO Queuing job for member 2...
2026-07-22 12:28:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 12:28:40 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-22 12:28:40 INFO Found: ['5247403']
2026-07-22 12:28:45 INFO [TGCC-IRENE] Submitted job with ID:['5247403']
2026-07-22 12:28:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 12:28:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3
2026-07-22 12:28:45 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020600_1_ENS3.nc
2026-07-22 12:28:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-22 12:28:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 12:28:45 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-22 12:28:45 INFO Queuing job for member 3...
2026-07-22 12:28:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 12:28:45 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-22 12:28:46 INFO Found: ['5247410']
2026-07-22 12:28:51 INFO [TGCC-IRENE] Submitted job with ID:['5247410']
2026-07-22 12:28:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 12:28:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4
2026-07-22 12:28:51 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020600_1_ENS4.nc
2026-07-22 12:28:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-22 12:28:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 12:28:51 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-22 12:28:51 INFO Queuing job for member 4...
2026-07-22 12:28:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 12:28:51 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-22 12:28:54 INFO Found: ['5247413']
2026-07-22 12:28:59 INFO [TGCC-IRENE] Submitted job with ID:['5247413']
2026-07-22 12:28:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 12:28:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5
2026-07-22 12:28:59 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020600_1_ENS5.nc
2026-07-22 12:28:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-22 12:28:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 12:28:59 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-22 12:28:59 INFO Queuing job for member 5...
2026-07-22 12:28:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 12:28:59 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-22 12:29:01 INFO Found: ['5247416']
2026-07-22 12:29:06 INFO [TGCC-IRENE] Submitted job with ID:['5247416']
2026-07-22 12:29:06 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 12:29:06 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6
2026-07-22 12:29:06 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020600_1_ENS6.nc
2026-07-22 12:29:06 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-22 12:29:06 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 12:29:06 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-22 12:29:06 INFO Queuing job for member 6...
2026-07-22 12:29:06 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 12:29:06 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-22 12:29:07 INFO Found: ['5247418']
2026-07-22 12:29:12 INFO [TGCC-IRENE] Submitted job with ID:['5247418']
2026-07-22 12:29:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 12:29:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7
2026-07-22 12:29:12 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020600_1_ENS7.nc
2026-07-22 12:29:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-22 12:29:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 12:29:12 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-22 12:29:12 INFO Queuing job for member 7...
2026-07-22 12:29:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 12:29:12 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-22 12:29:12 INFO Found: ['5247419']
2026-07-22 12:29:17 INFO [TGCC-IRENE] Submitted job with ID:['5247419']
2026-07-22 12:29:17 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 12:29:17 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8
2026-07-22 12:29:17 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020600_1_ENS8.nc
2026-07-22 12:29:17 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-22 12:29:17 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 12:29:18 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-22 12:29:18 INFO Queuing job for member 8...
2026-07-22 12:29:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 12:29:18 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-22 12:29:18 INFO Found: ['5247421']
2026-07-22 12:29:23 INFO [TGCC-IRENE] Submitted job with ID:['5247421']
2026-07-22 12:29:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 12:29:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9
2026-07-22 12:29:23 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020600_1_ENS9.nc
2026-07-22 12:29:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-22 12:29:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 12:29:23 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-22 12:29:23 INFO Queuing job for member 9...
2026-07-22 12:29:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 12:29:23 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-22 12:29:24 INFO Found: ['5247422']
2026-07-22 12:29:29 INFO [TGCC-IRENE] Submitted job with ID:['5247422']
2026-07-22 12:29:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 12:29:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10
2026-07-22 12:29:29 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020600_1_ENS10.nc
2026-07-22 12:29:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-22 12:29:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 12:29:29 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-22 12:29:29 INFO Queuing job for member 10...
2026-07-22 12:29:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 12:29:29 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-22 12:29:30 INFO Found: ['5247423']
2026-07-22 12:29:35 INFO [TGCC-IRENE] Submitted job with ID:['5247423']
2026-07-22 12:29:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 12:29:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11
2026-07-22 12:29:35 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020600_1_ENS11.nc
2026-07-22 12:29:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-22 12:29:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 12:29:35 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-22 12:29:35 INFO Queuing job for member 11...
2026-07-22 12:29:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 12:29:35 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-22 12:29:35 INFO Found: ['5247424']
2026-07-22 12:29:40 INFO [TGCC-IRENE] Submitted job with ID:['5247424']
2026-07-22 12:29:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 12:29:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12
2026-07-22 12:29:40 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020600_1_ENS12.nc
2026-07-22 12:29:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-22 12:29:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 12:29:41 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-22 12:29:41 INFO Queuing job for member 12...
2026-07-22 12:29:41 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 12:29:41 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-22 12:29:41 INFO Found: ['5247425']
2026-07-22 12:29:46 INFO [TGCC-IRENE] Submitted job with ID:['5247425']
2026-07-22 12:29:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 12:29:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13
2026-07-22 12:29:46 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020600_1_ENS13.nc
2026-07-22 12:29:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-22 12:29:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 12:29:46 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-22 12:29:46 INFO Queuing job for member 13...
2026-07-22 12:29:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 12:29:46 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-22 12:29:47 INFO Found: ['5247426']
2026-07-22 12:29:52 INFO [TGCC-IRENE] Submitted job with ID:['5247426']
2026-07-22 12:29:52 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 12:29:52 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14
2026-07-22 12:29:52 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020600_1_ENS14.nc
2026-07-22 12:29:52 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-22 12:29:52 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 12:29:52 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-22 12:29:52 INFO Queuing job for member 14...
2026-07-22 12:29:52 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 12:29:52 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-22 12:29:53 INFO Found: ['5247430']
2026-07-22 12:29:58 INFO [TGCC-IRENE] Submitted job with ID:['5247430']
2026-07-22 12:29:58 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 12:29:58 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15
2026-07-22 12:29:58 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020600_1_ENS15.nc
2026-07-22 12:29:58 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-22 12:29:58 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 12:29:58 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-22 12:29:58 INFO Queuing job for member 15...
2026-07-22 12:29:58 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 12:29:58 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-22 12:30:00 INFO Found: ['5247433']
2026-07-22 12:30:05 INFO [TGCC-IRENE] Submitted job with ID:['5247433']
2026-07-22 12:30:05 INFO Checking job status ...
2026-07-22 12:30:05 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:30:05 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:30:05 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:30:05 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:30:05 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:30:05 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:30:05 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:30:05 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:30:05 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:30:05 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:30:05 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:30:05 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:30:05 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:30:05 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:30:05 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:30:05 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:30:20 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:30:20 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:30:20 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:30:20 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:30:20 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:30:20 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:30:20 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:30:20 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:30:20 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:30:20 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:30:20 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:30:20 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:30:20 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:30:20 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:30:20 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:30:20 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:30:35 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:30:35 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:30:36 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:30:36 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:30:36 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:30:36 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:30:36 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:30:36 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:30:36 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:30:36 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:30:36 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:30:36 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:30:36 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:30:36 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:30:36 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:30:36 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:30:51 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:30:51 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:30:51 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:30:51 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:30:51 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:30:51 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:30:51 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:30:51 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:30:51 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:30:51 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:30:51 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:30:51 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:30:51 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:30:51 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:30:51 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:30:51 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:31:08 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:31:08 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:31:08 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:31:08 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:31:08 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:31:08 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:31:08 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:31:08 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:31:08 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:31:08 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:31:08 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:31:08 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:31:08 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:31:08 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:31:08 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:31:08 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:31:23 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:31:23 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:31:23 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:31:23 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:31:23 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:31:23 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:31:23 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:31:23 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:31:23 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:31:23 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:31:23 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:31:23 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:31:23 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:31:23 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:31:23 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:31:23 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:31:38 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:31:38 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:31:38 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:31:38 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:31:38 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:31:38 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:31:38 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:31:38 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:31:38 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:31:38 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:31:38 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:31:38 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:31:39 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:31:39 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:31:39 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:31:39 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:31:54 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:31:54 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:31:54 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:31:54 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:31:54 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:31:54 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:31:54 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:31:54 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:31:54 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:31:54 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:31:54 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:31:54 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:31:54 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:31:54 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:31:54 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:31:54 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:32:09 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:32:09 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:32:09 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:32:09 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:32:09 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:32:09 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:32:09 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:32:09 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:32:09 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:32:09 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:32:09 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:32:09 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:32:09 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:32:09 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:32:09 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:32:09 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:32:24 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:32:24 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:32:24 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:32:24 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:32:24 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:32:24 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:32:25 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:32:27 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:32:27 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:32:27 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:32:27 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:32:27 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:32:27 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:32:27 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:32:27 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:32:27 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:32:42 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:32:42 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:32:42 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:32:42 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:32:42 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:32:42 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:32:42 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:32:42 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:32:42 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:32:42 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:32:42 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:32:42 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:32:42 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:32:42 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:32:42 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:32:42 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:32:57 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:32:57 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:32:57 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:32:57 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:32:57 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:32:57 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:32:57 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:32:57 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:32:57 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:32:57 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:32:57 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:32:57 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:32:57 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:32:57 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:32:57 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:32:57 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:33:14 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:33:14 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:33:14 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:33:14 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:33:15 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:33:15 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:33:15 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:33:17 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:33:17 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:33:17 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:33:17 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:33:17 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:33:17 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:33:17 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:33:17 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:33:17 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:33:32 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:33:32 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:33:32 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:33:32 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:33:32 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:33:32 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:33:32 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:33:32 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:33:32 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:33:32 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:33:32 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:33:32 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:33:32 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:33:32 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:33:32 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:33:32 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:33:47 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:33:47 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:33:47 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:33:47 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:33:47 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:33:47 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:33:47 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:33:47 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:33:47 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:33:47 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:33:47 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:33:48 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:33:48 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:33:48 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:33:48 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:33:48 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:34:03 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:34:03 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:34:03 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:34:03 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:34:03 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:34:03 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:34:03 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:34:03 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:34:03 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:34:03 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:34:03 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:34:03 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:34:03 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:34:03 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:34:03 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:34:03 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:34:18 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:34:18 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:34:18 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:34:18 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:34:18 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:34:18 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:34:18 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:34:18 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:34:18 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:34:18 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:34:18 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:34:18 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:34:18 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:34:18 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:34:18 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:34:18 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:34:33 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:34:33 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:34:33 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:34:33 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:34:33 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:34:33 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:34:33 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:34:33 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:34:35 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:34:35 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:34:35 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:34:35 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:34:35 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:34:35 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:34:35 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:34:35 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:34:51 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:34:51 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:34:51 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:34:51 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:34:51 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:34:51 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:34:51 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:34:51 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:34:51 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:34:51 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:34:51 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:34:51 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:34:51 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:34:51 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:34:51 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:34:51 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:35:06 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:35:06 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:35:06 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:35:06 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:35:06 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:35:06 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:35:06 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:35:06 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:35:06 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:35:06 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:35:06 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:35:06 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:35:06 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:35:06 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:35:06 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:35:06 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:35:21 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:35:21 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:35:21 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:35:21 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:35:21 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:35:21 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:35:21 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:35:21 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:35:21 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:35:21 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:35:21 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:35:21 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:35:21 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:35:21 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:35:21 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:35:21 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:35:37 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:35:37 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:35:37 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:35:37 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:35:37 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:35:37 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:35:37 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:35:37 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:35:37 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:35:37 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:35:37 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:35:37 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:35:37 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:35:39 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:35:39 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:35:39 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:35:54 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:35:54 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:35:54 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:35:54 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:35:54 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:35:54 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:35:54 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:35:54 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:35:54 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:35:54 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:35:54 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:35:54 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:35:54 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:35:54 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:35:54 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:35:54 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:36:09 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:36:09 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:36:09 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:36:09 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:36:09 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:36:09 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:36:09 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:36:09 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:36:09 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:36:09 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:36:09 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:36:09 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:36:09 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:36:09 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:36:09 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:36:09 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:36:24 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:36:25 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:36:25 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:36:25 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:36:25 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:36:25 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:36:25 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:36:25 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:36:25 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:36:25 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:36:25 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:36:25 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:36:25 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:36:25 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:36:25 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:36:25 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:36:40 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:36:40 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:36:40 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:36:42 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:36:42 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:36:42 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:36:42 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:36:42 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:36:42 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:36:42 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:36:42 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:36:42 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:36:42 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:36:42 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:36:42 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:36:42 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:36:57 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:36:57 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:36:57 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:36:57 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:36:57 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:36:57 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:36:57 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:36:57 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:36:57 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:36:57 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:36:57 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:36:57 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:36:57 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:36:59 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:36:59 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:36:59 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:37:14 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:37:14 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:37:14 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:37:14 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:37:14 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:37:14 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:37:14 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:37:14 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:37:14 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:37:14 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:37:14 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:37:14 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:37:14 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:37:14 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:37:14 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:37:14 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:37:29 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:37:29 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:37:29 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:37:29 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:37:30 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:37:30 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:37:30 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:37:30 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:37:30 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:37:30 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:37:30 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:37:30 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:37:30 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:37:30 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:37:30 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:37:30 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:37:46 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:37:46 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:37:46 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:37:46 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:37:46 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:37:46 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:37:46 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:37:46 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:37:46 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:37:46 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:37:46 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:37:46 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:37:46 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:37:46 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:37:46 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:37:46 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:38:01 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:38:01 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:38:01 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:38:01 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:38:02 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:38:02 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:38:02 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:38:02 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:38:02 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:38:02 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:38:02 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:38:02 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:38:02 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:38:02 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:38:04 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:38:04 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:38:19 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:38:19 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:38:19 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:38:19 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:38:19 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:38:19 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:38:19 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:38:19 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:38:19 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:38:19 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:38:19 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:38:19 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:38:19 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:38:19 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:38:19 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:38:19 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:38:34 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:38:34 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:38:34 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:38:34 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:38:34 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:38:34 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:38:34 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:38:34 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:38:34 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:38:34 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:38:34 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:38:34 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:38:34 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:38:34 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:38:34 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:38:34 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:38:49 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:38:50 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:38:50 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:38:50 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:38:50 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:38:50 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:38:50 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:38:50 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:38:50 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:38:50 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:38:50 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:38:50 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:38:50 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:38:50 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:38:50 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:38:50 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:39:07 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:39:07 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:39:07 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:39:07 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:39:07 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:39:07 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:39:07 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:39:07 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:39:07 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:39:07 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:39:07 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:39:07 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:39:07 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:39:07 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:39:07 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:39:07 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:39:22 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:39:22 INFO None 5247403: status RUNNING/PENDING
2026-07-22 12:39:22 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:39:22 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:39:22 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:39:22 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:39:22 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:39:22 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:39:22 INFO None 5247422: status RUNNING/PENDING
2026-07-22 12:39:22 INFO None 5247423: status RUNNING/PENDING
2026-07-22 12:39:22 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:39:22 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:39:22 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:39:22 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:39:22 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:39:22 INFO Jobs still running: ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:39:37 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:39:37 INFO None 5247403: status FINISHED
2026-07-22 12:39:37 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:39:37 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:39:37 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:39:37 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:39:37 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:39:37 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:39:38 INFO None 5247422: status FINISHED
2026-07-22 12:39:38 INFO None 5247423: status FINISHED
2026-07-22 12:39:38 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:39:38 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:39:38 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:39:38 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:39:38 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:39:38 INFO Jobs still running: ['5247397', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:39:53 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:39:53 INFO None 5247403: status FINISHED
2026-07-22 12:39:53 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:39:53 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:39:53 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:39:53 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:39:53 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:39:53 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:39:53 INFO None 5247422: status FINISHED
2026-07-22 12:39:53 INFO None 5247423: status FINISHED
2026-07-22 12:39:53 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:39:53 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:39:53 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:39:53 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:39:53 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:39:53 INFO Jobs still running: ['5247397', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:40:08 INFO None 5247397: status RUNNING/PENDING
2026-07-22 12:40:08 INFO None 5247403: status FINISHED
2026-07-22 12:40:08 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:40:08 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:40:08 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:40:08 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:40:08 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:40:08 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:40:08 INFO None 5247422: status FINISHED
2026-07-22 12:40:08 INFO None 5247423: status FINISHED
2026-07-22 12:40:08 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:40:08 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:40:08 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:40:08 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:40:08 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:40:08 INFO Jobs still running: ['5247397', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:40:23 INFO None 5247397: status FINISHED
2026-07-22 12:40:23 INFO None 5247403: status FINISHED
2026-07-22 12:40:23 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:40:23 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:40:23 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:40:24 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:40:24 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:40:24 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:40:24 INFO None 5247422: status FINISHED
2026-07-22 12:40:24 INFO None 5247423: status FINISHED
2026-07-22 12:40:24 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:40:24 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:40:24 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:40:24 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:40:26 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:40:26 INFO Jobs still running: ['5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:40:41 INFO None 5247397: status FINISHED
2026-07-22 12:40:41 INFO None 5247403: status FINISHED
2026-07-22 12:40:41 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:40:41 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:40:41 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:40:41 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:40:41 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:40:41 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:40:41 INFO None 5247422: status FINISHED
2026-07-22 12:40:41 INFO None 5247423: status FINISHED
2026-07-22 12:40:41 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:40:41 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:40:41 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:40:41 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:40:41 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:40:41 INFO Jobs still running: ['5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:40:56 INFO None 5247397: status FINISHED
2026-07-22 12:40:56 INFO None 5247403: status FINISHED
2026-07-22 12:40:56 INFO None 5247410: status RUNNING/PENDING
2026-07-22 12:40:56 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:40:56 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:40:56 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:40:56 INFO None 5247419: status RUNNING/PENDING
2026-07-22 12:40:56 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:40:56 INFO None 5247422: status FINISHED
2026-07-22 12:40:56 INFO None 5247423: status FINISHED
2026-07-22 12:40:56 INFO None 5247424: status RUNNING/PENDING
2026-07-22 12:40:56 INFO None 5247425: status RUNNING/PENDING
2026-07-22 12:40:56 INFO None 5247426: status RUNNING/PENDING
2026-07-22 12:40:56 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:40:56 INFO None 5247433: status RUNNING/PENDING
2026-07-22 12:40:56 INFO Jobs still running: ['5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247424', '5247425', '5247426', '5247430', '5247433']. Waiting...
2026-07-22 12:41:11 INFO None 5247397: status FINISHED
2026-07-22 12:41:11 INFO None 5247403: status FINISHED
2026-07-22 12:41:11 INFO None 5247410: status FINISHED
2026-07-22 12:41:11 INFO None 5247413: status RUNNING/PENDING
2026-07-22 12:41:11 INFO None 5247416: status RUNNING/PENDING
2026-07-22 12:41:12 INFO None 5247418: status RUNNING/PENDING
2026-07-22 12:41:12 INFO None 5247419: status FINISHED
2026-07-22 12:41:12 INFO None 5247421: status RUNNING/PENDING
2026-07-22 12:41:12 INFO None 5247422: status FINISHED
2026-07-22 12:41:12 INFO None 5247423: status FINISHED
2026-07-22 12:41:12 INFO None 5247424: status FINISHED
2026-07-22 12:41:12 INFO None 5247425: status FINISHED
2026-07-22 12:41:12 INFO None 5247426: status FINISHED
2026-07-22 12:41:12 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:41:12 INFO None 5247433: status FINISHED
2026-07-22 12:41:12 INFO Jobs still running: ['5247413', '5247416', '5247418', '5247421', '5247430']. Waiting...
2026-07-22 12:41:27 INFO None 5247397: status FINISHED
2026-07-22 12:41:27 INFO None 5247403: status FINISHED
2026-07-22 12:41:27 INFO None 5247410: status FINISHED
2026-07-22 12:41:27 INFO None 5247413: status FINISHED
2026-07-22 12:41:27 INFO None 5247416: status FINISHED
2026-07-22 12:41:27 INFO None 5247418: status FINISHED
2026-07-22 12:41:27 INFO None 5247419: status FINISHED
2026-07-22 12:41:27 INFO None 5247421: status FINISHED
2026-07-22 12:41:27 INFO None 5247422: status FINISHED
2026-07-22 12:41:27 INFO None 5247423: status FINISHED
2026-07-22 12:41:27 INFO None 5247424: status FINISHED
2026-07-22 12:41:27 INFO None 5247425: status FINISHED
2026-07-22 12:41:27 INFO None 5247426: status FINISHED
2026-07-22 12:41:27 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:41:27 INFO None 5247433: status FINISHED
2026-07-22 12:41:27 INFO Jobs still running: ['5247430']. Waiting...
2026-07-22 12:41:42 INFO None 5247397: status FINISHED
2026-07-22 12:41:42 INFO None 5247403: status FINISHED
2026-07-22 12:41:42 INFO None 5247410: status FINISHED
2026-07-22 12:41:42 INFO None 5247413: status FINISHED
2026-07-22 12:41:42 INFO None 5247416: status FINISHED
2026-07-22 12:41:44 INFO None 5247418: status FINISHED
2026-07-22 12:41:44 INFO None 5247419: status FINISHED
2026-07-22 12:41:44 INFO None 5247421: status FINISHED
2026-07-22 12:41:44 INFO None 5247422: status FINISHED
2026-07-22 12:41:44 INFO None 5247423: status FINISHED
2026-07-22 12:41:44 INFO None 5247424: status FINISHED
2026-07-22 12:41:44 INFO None 5247425: status FINISHED
2026-07-22 12:41:44 INFO None 5247426: status FINISHED
2026-07-22 12:41:44 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:41:44 INFO None 5247433: status FINISHED
2026-07-22 12:41:44 INFO Jobs still running: ['5247430']. Waiting...
2026-07-22 12:41:59 INFO None 5247397: status FINISHED
2026-07-22 12:41:59 INFO None 5247403: status FINISHED
2026-07-22 12:41:59 INFO None 5247410: status FINISHED
2026-07-22 12:42:00 INFO None 5247413: status FINISHED
2026-07-22 12:42:00 INFO None 5247416: status FINISHED
2026-07-22 12:42:00 INFO None 5247418: status FINISHED
2026-07-22 12:42:00 INFO None 5247419: status FINISHED
2026-07-22 12:42:00 INFO None 5247421: status FINISHED
2026-07-22 12:42:00 INFO None 5247422: status FINISHED
2026-07-22 12:42:00 INFO None 5247423: status FINISHED
2026-07-22 12:42:00 INFO None 5247424: status FINISHED
2026-07-22 12:42:00 INFO None 5247425: status FINISHED
2026-07-22 12:42:00 INFO None 5247426: status FINISHED
2026-07-22 12:42:00 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:42:00 INFO None 5247433: status FINISHED
2026-07-22 12:42:00 INFO Jobs still running: ['5247430']. Waiting...
2026-07-22 12:42:15 INFO None 5247397: status FINISHED
2026-07-22 12:42:15 INFO None 5247403: status FINISHED
2026-07-22 12:42:15 INFO None 5247410: status FINISHED
2026-07-22 12:42:15 INFO None 5247413: status FINISHED
2026-07-22 12:42:15 INFO None 5247416: status FINISHED
2026-07-22 12:42:15 INFO None 5247418: status FINISHED
2026-07-22 12:42:15 INFO None 5247419: status FINISHED
2026-07-22 12:42:15 INFO None 5247421: status FINISHED
2026-07-22 12:42:15 INFO None 5247422: status FINISHED
2026-07-22 12:42:15 INFO None 5247423: status FINISHED
2026-07-22 12:42:15 INFO None 5247424: status FINISHED
2026-07-22 12:42:15 INFO None 5247425: status FINISHED
2026-07-22 12:42:15 INFO None 5247426: status FINISHED
2026-07-22 12:42:15 INFO None 5247430: status RUNNING/PENDING
2026-07-22 12:42:15 INFO None 5247433: status FINISHED
2026-07-22 12:42:15 INFO Jobs still running: ['5247430']. Waiting...
2026-07-22 12:42:30 INFO None 5247397: status FINISHED
2026-07-22 12:42:30 INFO None 5247403: status FINISHED
2026-07-22 12:42:30 INFO None 5247410: status FINISHED
2026-07-22 12:42:30 INFO None 5247413: status FINISHED
2026-07-22 12:42:30 INFO None 5247416: status FINISHED
2026-07-22 12:42:30 INFO None 5247418: status FINISHED
2026-07-22 12:42:30 INFO None 5247419: status FINISHED
2026-07-22 12:42:30 INFO None 5247421: status FINISHED
2026-07-22 12:42:30 INFO None 5247422: status FINISHED
2026-07-22 12:42:30 INFO None 5247423: status FINISHED
2026-07-22 12:42:30 INFO None 5247424: status FINISHED
2026-07-22 12:42:30 INFO None 5247425: status FINISHED
2026-07-22 12:42:30 INFO None 5247426: status FINISHED
2026-07-22 12:42:30 INFO None 5247430: status FINISHED
2026-07-22 12:42:30 INFO None 5247433: status FINISHED
2026-07-22 12:42:30 INFO Jobs ['5247397', '5247403', '5247410', '5247413', '5247416', '5247418', '5247419', '5247421', '5247422', '5247423', '5247424', '5247425', '5247426', '5247430', '5247433'] have finished
2026-07-22 12:42:30 INFO Checking restart files were created ...
2026-07-22 12:42:30 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020601_8_ENS1.nc(3005806795 bytes)
2026-07-22 12:42:30 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020601_8_ENS2.nc(3005806795 bytes)
2026-07-22 12:42:30 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020601_8_ENS3.nc(3005806795 bytes)
2026-07-22 12:42:30 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020601_8_ENS4.nc(3005806795 bytes)
2026-07-22 12:42:30 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020601_8_ENS5.nc(3005806795 bytes)
2026-07-22 12:42:30 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020601_8_ENS6.nc(3005806795 bytes)
2026-07-22 12:42:30 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020601_8_ENS7.nc(3005806795 bytes)
2026-07-22 12:42:30 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020601_8_ENS8.nc(3005806795 bytes)
2026-07-22 12:42:30 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020601_8_ENS9.nc(3005806795 bytes)
2026-07-22 12:42:30 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020601_8_ENS10.nc(3005806795 bytes)
2026-07-22 12:42:30 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020601_8_ENS11.nc(3005806795 bytes)
2026-07-22 12:42:30 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020601_8_ENS12.nc(3005806795 bytes)
2026-07-22 12:42:30 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020601_8_ENS13.nc(3005806795 bytes)
2026-07-22 12:42:30 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020601_8_ENS14.nc(3005806795 bytes)
2026-07-22 12:42:30 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020601_8_ENS15.nc(3005806795 bytes)
2026-07-22 12:42:30 INFO  Run_model() completed successfully.
2026-07-22 12:42:30 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 12:42:30 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 09:00:00 days=153072 seconds=32400
2026-07-22 12:42:30 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-22 12:42:30 INFO [TIME] increment current_time 2020-02-06 01:00:00 -> 2020-02-06 09:00:00
2026-07-22 12:42:30 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 12:42:30 INFO ---------->>> Running process_satellite_data()
2026-07-22 12:42:30 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12002.nc
2026-07-22 12:42:30 INFO ---------->>> Running run_obs_converter()
2026-07-22 12:42:30 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_33067_153072.out
2026-07-22 12:42:30 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_33067_153072.out
2026-07-22 12:42:30 INFO ---------->>> Running DART
2026-07-22 12:42:30 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-22 12:42:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/chim_ENS1_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/chim_ENS1_2020020609_1_out_toDART.nc
2026-07-22 12:42:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/chim_ENS2_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/chim_ENS2_2020020609_1_out_toDART.nc
2026-07-22 12:42:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/chim_ENS3_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/chim_ENS3_2020020609_1_out_toDART.nc
2026-07-22 12:42:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/chim_ENS4_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/chim_ENS4_2020020609_1_out_toDART.nc
2026-07-22 12:42:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/chim_ENS5_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/chim_ENS5_2020020609_1_out_toDART.nc
2026-07-22 12:42:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/chim_ENS6_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/chim_ENS6_2020020609_1_out_toDART.nc
2026-07-22 12:42:33 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/chim_ENS7_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/chim_ENS7_2020020609_1_out_toDART.nc
2026-07-22 12:42:33 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/chim_ENS8_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/chim_ENS8_2020020609_1_out_toDART.nc
2026-07-22 12:42:33 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/chim_ENS9_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/chim_ENS9_2020020609_1_out_toDART.nc
2026-07-22 12:42:34 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/chim_ENS10_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/chim_ENS10_2020020609_1_out_toDART.nc
2026-07-22 12:42:34 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/chim_ENS11_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/chim_ENS11_2020020609_1_out_toDART.nc
2026-07-22 12:42:34 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/chim_ENS12_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/chim_ENS12_2020020609_1_out_toDART.nc
2026-07-22 12:42:35 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/chim_ENS13_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/chim_ENS13_2020020609_1_out_toDART.nc
2026-07-22 12:42:35 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/chim_ENS14_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/chim_ENS14_2020020609_1_out_toDART.nc
2026-07-22 12:42:35 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/chim_ENS15_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/chim_ENS15_2020020609_1_out_toDART.nc
2026-07-22 12:42:36 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-22 12:42:36 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-22 12:42:36 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-22 12:42:36 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-22 12:42:36 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-22 12:42:36 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-22 12:42:43 INFO Found: []
2026-07-22 12:42:43 INFO No job id returned by command ./run_filter.bsh
2026-07-22 12:42:43 INFO No monitoring will be performed
2026-07-22 12:42:43 INFO Moving DART output files to analysis and preassim directories for date 2020020609 if present ...
2026-07-22 12:42:43 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 12:42:43 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 12:42:43 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 12:42:43 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 12:42:43 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 12:42:43 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 12:42:43 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 12:42:43 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 12:42:43 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 12:42:43 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 12:42:43 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 12:42:43 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 12:42:43 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 12:42:43 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 12:42:43 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 12:42:43 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 12:42:43 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 12:42:43 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 12:42:43 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 12:42:43 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 12:42:43 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 12:42:43 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 12:42:43 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 12:42:43 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 12:42:43 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 12:42:43 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 12:42:43 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 12:42:43 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 12:42:43 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 12:42:43 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 12:42:43 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 12:42:43 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 12:42:43 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 12:42:43 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 12:42:43 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-22 12:42:46 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-22 12:42:46 INFO run_dart() is DONE.
2026-07-22 12:42:46 INFO ---------->>> Running update_pollutant_in_end()
2026-07-22 12:42:46 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 12:42:46 INFO No previous orbit memory found.
2026-07-22 12:42:47 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc
Traceback (most recent call last):
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/site-packages/xarray/backends/file_manager.py", line 211, in _acquire_with_cache_info
    file = self._cache[self._key]
           ~~~~~~~~~~~^^^^^^^^^^^
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/site-packages/xarray/backends/lru_cache.py", line 56, in __getitem__
    value = self._cache[key]
            ~~~~~~~~~~~^^^^^
KeyError: [<class 'netCDF4._netCDF4.Dataset'>, ('/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc',), 'r', (('clobber', True), ('diskless', False), ('format', 'NETCDF4'), ('persist', False)), '2fd2ad14-ea0f-40f1-bc5e-0ea02f7de47a']

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/main.py", line 103, in <module>
    pipeline.run_pipeline()
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/pipelines/base_pipeline.py", line 176, in run_pipeline
    self.after_assimilation()
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/pipelines/chimere2023/pipeline.py", line 604, in after_assimilation
    update_pollutant_in_end(dart_file=self.paths.dart_filter_output_list_file(mem, date_ymdH), 
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/orchestrator_utils.py", line 1459, in update_pollutant_in_end
    with xr.open_dataset(next_emis_file, decode_timedelta=True) as next_emis_ds:
         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/site-packages/xarray/backends/api.py", line 715, in open_dataset
    backend_ds = backend.open_dataset(
                 ^^^^^^^^^^^^^^^^^^^^^
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/site-packages/xarray/backends/netCDF4_.py", line 671, in open_dataset
    store = NetCDF4DataStore.open(
            ^^^^^^^^^^^^^^^^^^^^^^
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/site-packages/xarray/backends/netCDF4_.py", line 457, in open
    return cls(manager, group=group, mode=mode, lock=lock, autoclose=autoclose)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/site-packages/xarray/backends/netCDF4_.py", line 398, in __init__
    self.format = self.ds.data_model
                  ^^^^^^^
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/site-packages/xarray/backends/netCDF4_.py", line 466, in ds
    return self._acquire()
           ^^^^^^^^^^^^^^^
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/site-packages/xarray/backends/netCDF4_.py", line 460, in _acquire
    with self._manager.acquire_context(needs_lock) as root:
         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/contextlib.py", line 137, in __enter__
    return next(self.gen)
           ^^^^^^^^^^^^^^
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/site-packages/xarray/backends/file_manager.py", line 199, in acquire_context
    file, cached = self._acquire_with_cache_info(needs_lock)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/ccc/products2/python3-3.12/Rhel_8__x86_64/system/default/lib/python3.12/site-packages/xarray/backends/file_manager.py", line 217, in _acquire_with_cache_info
    file = self._opener(*self._args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "src/netCDF4/_netCDF4.pyx", line 2521, in netCDF4._netCDF4.Dataset.__init__
  File "src/netCDF4/_netCDF4.pyx", line 2158, in netCDF4._netCDF4._ensure_nc_success
FileNotFoundError: [Errno 2] No such file or directory: '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc'
+ exit 0
