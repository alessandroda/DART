+ SCRIPT_PID=3401118
+ /bin/bash -x /tmp/tmp.ZkC9aKeR1j
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
2026-07-22 14:21:51 INFO 
███    ███ ██ ███    ███ ███████ ███████ ██
████  ████ ██ ████  ████ ██      ██      ██
██ ████ ██ ██ ██ ████ ██ █████   ███████ ██
██  ██  ██ ██ ██  ██  ██ ██           ██ ██
██      ██ ██ ██      ██ ███████ ███████ ██



2026-07-22 14:21:51 INFO [PIPELINE] =======================================
2026-07-22 14:21:51 INFO [PIPELINE] Starting chimere–DART orchestrator
2026-07-22 14:21:51 INFO [PIPELINE] Config file: config/config_irene_IM.yaml
2026-07-22 14:21:51 INFO [PIPELINE] Run dir: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart
2026-07-22 14:21:51 INFO [PIPELINE] Log file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/mimesi_orchestrator_logs/chimere_DART_20260722_142151.log
2026-07-22 14:21:51 INFO [PIPELINE] =======================================
2026-07-22 14:21:51 INFO Running assimilation with model_type=ModelType.CHIMERE
2026-07-22 14:21:51 INFO Using scheduler=Scheduler.SLURM, queue=rome
2026-07-22 14:21:51 INFO [STEP] ---- TIME LOOP START ----
2026-07-22 14:21:51 INFO [TIME] step_start current_time=2020-02-06 00:00:00 simulated_time=None dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 14:21:51 INFO [TIME] window start=2020-02-06 00:00:00 end=2020-02-06 01:00:00 run_hours=1 has_assimilation=False
2026-07-22 14:21:51 INFO Asked to restart from control run ...
2026-07-22 14:21:51 INFO Copying EMIS ...
2026-07-22 14:21:52 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens01.nc
2026-07-22 14:21:52 INFO Linking first END ...
2026-07-22 14:21:52 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020500_24_ENS1.nc
2026-07-22 14:21:52 INFO >> Checking links...
2026-07-22 14:21:52 INFO >> All links are good for ENS1  ...
2026-07-22 14:21:52 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:22:00 INFO Hourly dataset computed and listing created
2026-07-22 14:22:03 INFO Hourly dataset computed
2026-07-22 14:22:03 INFO Asked to restart from control run ...
2026-07-22 14:22:03 INFO Copying EMIS ...
2026-07-22 14:22:03 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens02.nc
2026-07-22 14:22:03 INFO Linking first END ...
2026-07-22 14:22:03 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020500_24_ENS2.nc
2026-07-22 14:22:03 INFO >> Checking links...
2026-07-22 14:22:04 INFO >> All links are good for ENS2  ...
2026-07-22 14:22:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:22:04 INFO Hourly dataset computed and listing created
2026-07-22 14:22:05 INFO Hourly dataset computed
2026-07-22 14:22:05 INFO Asked to restart from control run ...
2026-07-22 14:22:05 INFO Copying EMIS ...
2026-07-22 14:22:05 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens03.nc
2026-07-22 14:22:05 INFO Linking first END ...
2026-07-22 14:22:05 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020500_24_ENS3.nc
2026-07-22 14:22:05 INFO >> Checking links...
2026-07-22 14:22:06 INFO >> All links are good for ENS3  ...
2026-07-22 14:22:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:22:07 INFO Hourly dataset computed and listing created
2026-07-22 14:22:07 INFO Hourly dataset computed
2026-07-22 14:22:07 INFO Asked to restart from control run ...
2026-07-22 14:22:07 INFO Copying EMIS ...
2026-07-22 14:22:08 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens04.nc
2026-07-22 14:22:08 INFO Linking first END ...
2026-07-22 14:22:08 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020500_24_ENS4.nc
2026-07-22 14:22:08 INFO >> Checking links...
2026-07-22 14:22:08 INFO >> All links are good for ENS4  ...
2026-07-22 14:22:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:22:09 INFO Hourly dataset computed and listing created
2026-07-22 14:22:09 INFO Hourly dataset computed
2026-07-22 14:22:09 INFO Asked to restart from control run ...
2026-07-22 14:22:09 INFO Copying EMIS ...
2026-07-22 14:22:10 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens05.nc
2026-07-22 14:22:10 INFO Linking first END ...
2026-07-22 14:22:10 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020500_24_ENS5.nc
2026-07-22 14:22:10 INFO >> Checking links...
2026-07-22 14:22:11 INFO >> All links are good for ENS5  ...
2026-07-22 14:22:11 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:22:11 INFO Hourly dataset computed and listing created
2026-07-22 14:22:12 INFO Hourly dataset computed
2026-07-22 14:22:12 INFO Asked to restart from control run ...
2026-07-22 14:22:12 INFO Copying EMIS ...
2026-07-22 14:22:13 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens06.nc
2026-07-22 14:22:13 INFO Linking first END ...
2026-07-22 14:22:13 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020500_24_ENS6.nc
2026-07-22 14:22:13 INFO >> Checking links...
2026-07-22 14:22:13 INFO >> All links are good for ENS6  ...
2026-07-22 14:22:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:22:14 INFO Hourly dataset computed and listing created
2026-07-22 14:22:14 INFO Hourly dataset computed
2026-07-22 14:22:14 INFO Asked to restart from control run ...
2026-07-22 14:22:14 INFO Copying EMIS ...
2026-07-22 14:22:15 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens07.nc
2026-07-22 14:22:15 INFO Linking first END ...
2026-07-22 14:22:15 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020500_24_ENS7.nc
2026-07-22 14:22:15 INFO >> Checking links...
2026-07-22 14:22:15 INFO >> All links are good for ENS7  ...
2026-07-22 14:22:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:22:16 INFO Hourly dataset computed and listing created
2026-07-22 14:22:17 INFO Hourly dataset computed
2026-07-22 14:22:17 INFO Asked to restart from control run ...
2026-07-22 14:22:17 INFO Copying EMIS ...
2026-07-22 14:22:17 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens08.nc
2026-07-22 14:22:17 INFO Linking first END ...
2026-07-22 14:22:17 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020500_24_ENS8.nc
2026-07-22 14:22:17 INFO >> Checking links...
2026-07-22 14:22:17 INFO >> All links are good for ENS8  ...
2026-07-22 14:22:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:22:19 INFO Hourly dataset computed and listing created
2026-07-22 14:22:24 INFO Hourly dataset computed
2026-07-22 14:22:24 INFO Asked to restart from control run ...
2026-07-22 14:22:24 INFO Copying EMIS ...
2026-07-22 14:22:24 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens09.nc
2026-07-22 14:22:24 INFO Linking first END ...
2026-07-22 14:22:24 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020500_24_ENS9.nc
2026-07-22 14:22:24 INFO >> Checking links...
2026-07-22 14:22:24 INFO >> All links are good for ENS9  ...
2026-07-22 14:22:24 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:22:26 INFO Hourly dataset computed and listing created
2026-07-22 14:22:26 INFO Hourly dataset computed
2026-07-22 14:22:27 INFO Asked to restart from control run ...
2026-07-22 14:22:27 INFO Copying EMIS ...
2026-07-22 14:22:27 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens10.nc
2026-07-22 14:22:27 INFO Linking first END ...
2026-07-22 14:22:27 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020500_24_ENS10.nc
2026-07-22 14:22:27 INFO >> Checking links...
2026-07-22 14:22:27 INFO >> All links are good for ENS10  ...
2026-07-22 14:22:27 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:22:28 INFO Hourly dataset computed and listing created
2026-07-22 14:22:29 INFO Hourly dataset computed
2026-07-22 14:22:29 INFO Asked to restart from control run ...
2026-07-22 14:22:29 INFO Copying EMIS ...
2026-07-22 14:22:29 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens11.nc
2026-07-22 14:22:29 INFO Linking first END ...
2026-07-22 14:22:29 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020500_24_ENS11.nc
2026-07-22 14:22:29 INFO >> Checking links...
2026-07-22 14:22:30 INFO >> All links are good for ENS11  ...
2026-07-22 14:22:30 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:22:30 INFO Hourly dataset computed and listing created
2026-07-22 14:22:31 INFO Hourly dataset computed
2026-07-22 14:22:31 INFO Asked to restart from control run ...
2026-07-22 14:22:31 INFO Copying EMIS ...
2026-07-22 14:22:31 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens12.nc
2026-07-22 14:22:31 INFO Linking first END ...
2026-07-22 14:22:31 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020500_24_ENS12.nc
2026-07-22 14:22:31 INFO >> Checking links...
2026-07-22 14:22:32 INFO >> All links are good for ENS12  ...
2026-07-22 14:22:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:22:33 INFO Hourly dataset computed and listing created
2026-07-22 14:22:33 INFO Hourly dataset computed
2026-07-22 14:22:33 INFO Asked to restart from control run ...
2026-07-22 14:22:33 INFO Copying EMIS ...
2026-07-22 14:22:34 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens13.nc
2026-07-22 14:22:34 INFO Linking first END ...
2026-07-22 14:22:34 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020500_24_ENS13.nc
2026-07-22 14:22:34 INFO >> Checking links...
2026-07-22 14:22:34 INFO >> All links are good for ENS13  ...
2026-07-22 14:22:34 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:22:35 INFO Hourly dataset computed and listing created
2026-07-22 14:22:35 INFO Hourly dataset computed
2026-07-22 14:22:36 INFO Asked to restart from control run ...
2026-07-22 14:22:36 INFO Copying EMIS ...
2026-07-22 14:22:36 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens14.nc
2026-07-22 14:22:36 INFO Linking first END ...
2026-07-22 14:22:36 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020500_24_ENS14.nc
2026-07-22 14:22:36 INFO >> Checking links...
2026-07-22 14:22:36 INFO >> All links are good for ENS14  ...
2026-07-22 14:22:36 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:22:37 INFO Hourly dataset computed and listing created
2026-07-22 14:22:38 INFO Hourly dataset computed
2026-07-22 14:22:38 INFO Asked to restart from control run ...
2026-07-22 14:22:38 INFO Copying EMIS ...
2026-07-22 14:22:38 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens15.nc
2026-07-22 14:22:38 INFO Linking first END ...
2026-07-22 14:22:38 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020500_24_ENS15.nc
2026-07-22 14:22:38 INFO >> Checking links...
2026-07-22 14:22:38 INFO >> All links are good for ENS15  ...
2026-07-22 14:22:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:22:39 INFO Hourly dataset computed and listing created
2026-07-22 14:22:40 INFO Hourly dataset computed
2026-07-22 14:22:40 INFO ---------->>> Running CHIMERE model from 2020-02-06 00:00:00 to 2020-02-06 01:00:00
2026-07-22 14:22:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:22:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1
2026-07-22 14:22:40 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020500_24_ENS1.nc
2026-07-22 14:22:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-22 14:22:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:22:40 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-22 14:22:40 INFO Queuing job for member 1...
2026-07-22 14:22:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:22:40 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-22 14:22:41 INFO Found: ['5251627']
2026-07-22 14:22:46 INFO [TGCC-IRENE] Submitted job with ID:['5251627']
2026-07-22 14:22:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:22:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2
2026-07-22 14:22:46 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020500_24_ENS2.nc
2026-07-22 14:22:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-22 14:22:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:22:46 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-22 14:22:46 INFO Queuing job for member 2...
2026-07-22 14:22:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:22:46 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-22 14:22:49 INFO Found: ['5251631']
2026-07-22 14:22:54 INFO [TGCC-IRENE] Submitted job with ID:['5251631']
2026-07-22 14:22:54 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:22:54 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3
2026-07-22 14:22:54 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020500_24_ENS3.nc
2026-07-22 14:22:54 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-22 14:22:54 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:22:54 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-22 14:22:54 INFO Queuing job for member 3...
2026-07-22 14:22:54 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:22:54 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-22 14:22:56 INFO Found: ['5251632']
2026-07-22 14:23:01 INFO [TGCC-IRENE] Submitted job with ID:['5251632']
2026-07-22 14:23:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:23:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4
2026-07-22 14:23:01 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020500_24_ENS4.nc
2026-07-22 14:23:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-22 14:23:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:23:01 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-22 14:23:01 INFO Queuing job for member 4...
2026-07-22 14:23:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:23:01 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-22 14:23:04 INFO Found: ['5251634']
2026-07-22 14:23:09 INFO [TGCC-IRENE] Submitted job with ID:['5251634']
2026-07-22 14:23:09 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:23:09 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5
2026-07-22 14:23:09 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020500_24_ENS5.nc
2026-07-22 14:23:09 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-22 14:23:09 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:23:09 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-22 14:23:09 INFO Queuing job for member 5...
2026-07-22 14:23:09 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:23:09 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-22 14:23:10 INFO Found: ['5251639']
2026-07-22 14:23:15 INFO [TGCC-IRENE] Submitted job with ID:['5251639']
2026-07-22 14:23:15 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:23:15 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6
2026-07-22 14:23:15 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020500_24_ENS6.nc
2026-07-22 14:23:15 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-22 14:23:15 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:23:15 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-22 14:23:15 INFO Queuing job for member 6...
2026-07-22 14:23:15 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:23:15 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-22 14:23:15 INFO Found: ['5251641']
2026-07-22 14:23:20 INFO [TGCC-IRENE] Submitted job with ID:['5251641']
2026-07-22 14:23:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:23:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7
2026-07-22 14:23:20 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020500_24_ENS7.nc
2026-07-22 14:23:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-22 14:23:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:23:20 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-22 14:23:20 INFO Queuing job for member 7...
2026-07-22 14:23:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:23:20 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-22 14:23:21 INFO Found: ['5251644']
2026-07-22 14:23:26 INFO [TGCC-IRENE] Submitted job with ID:['5251644']
2026-07-22 14:23:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:23:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8
2026-07-22 14:23:26 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020500_24_ENS8.nc
2026-07-22 14:23:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-22 14:23:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:23:26 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-22 14:23:26 INFO Queuing job for member 8...
2026-07-22 14:23:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:23:26 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-22 14:23:27 INFO Found: ['5251651']
2026-07-22 14:23:32 INFO [TGCC-IRENE] Submitted job with ID:['5251651']
2026-07-22 14:23:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:23:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9
2026-07-22 14:23:32 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020500_24_ENS9.nc
2026-07-22 14:23:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-22 14:23:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:23:32 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-22 14:23:32 INFO Queuing job for member 9...
2026-07-22 14:23:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:23:32 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-22 14:23:33 INFO Found: ['5251656']
2026-07-22 14:23:38 INFO [TGCC-IRENE] Submitted job with ID:['5251656']
2026-07-22 14:23:38 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:23:38 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10
2026-07-22 14:23:38 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020500_24_ENS10.nc
2026-07-22 14:23:38 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-22 14:23:38 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:23:38 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-22 14:23:38 INFO Queuing job for member 10...
2026-07-22 14:23:38 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:23:38 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-22 14:23:38 INFO Found: ['5251659']
2026-07-22 14:23:43 INFO [TGCC-IRENE] Submitted job with ID:['5251659']
2026-07-22 14:23:43 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:23:43 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11
2026-07-22 14:23:43 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020500_24_ENS11.nc
2026-07-22 14:23:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-22 14:23:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:23:43 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-22 14:23:43 INFO Queuing job for member 11...
2026-07-22 14:23:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:23:43 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-22 14:23:44 INFO Found: ['5251668']
2026-07-22 14:23:49 INFO [TGCC-IRENE] Submitted job with ID:['5251668']
2026-07-22 14:23:49 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:23:49 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12
2026-07-22 14:23:49 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020500_24_ENS12.nc
2026-07-22 14:23:49 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-22 14:23:49 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:23:49 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-22 14:23:49 INFO Queuing job for member 12...
2026-07-22 14:23:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:23:49 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-22 14:23:50 INFO Found: ['5251677']
2026-07-22 14:23:55 INFO [TGCC-IRENE] Submitted job with ID:['5251677']
2026-07-22 14:23:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:23:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13
2026-07-22 14:23:55 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020500_24_ENS13.nc
2026-07-22 14:23:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-22 14:23:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:23:55 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-22 14:23:55 INFO Queuing job for member 13...
2026-07-22 14:23:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:23:55 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-22 14:23:57 INFO Found: ['5251685']
2026-07-22 14:24:02 INFO [TGCC-IRENE] Submitted job with ID:['5251685']
2026-07-22 14:24:02 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:24:02 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14
2026-07-22 14:24:02 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020500_24_ENS14.nc
2026-07-22 14:24:02 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-22 14:24:02 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:24:02 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-22 14:24:02 INFO Queuing job for member 14...
2026-07-22 14:24:02 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:24:02 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-22 14:24:04 INFO Found: ['5251691']
2026-07-22 14:24:09 INFO [TGCC-IRENE] Submitted job with ID:['5251691']
2026-07-22 14:24:09 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:24:09 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15
2026-07-22 14:24:09 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020500_24_ENS15.nc
2026-07-22 14:24:09 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-22 14:24:09 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:24:09 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-22 14:24:09 INFO Queuing job for member 15...
2026-07-22 14:24:09 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:24:09 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-22 14:24:12 INFO Found: ['5251695']
2026-07-22 14:24:17 INFO [TGCC-IRENE] Submitted job with ID:['5251695']
2026-07-22 14:24:17 INFO Checking job status ...
2026-07-22 14:24:17 INFO None 5251627: status RUNNING/PENDING
2026-07-22 14:24:17 INFO None 5251631: status RUNNING/PENDING
2026-07-22 14:24:17 INFO None 5251632: status RUNNING/PENDING
2026-07-22 14:24:17 INFO None 5251634: status RUNNING/PENDING
2026-07-22 14:24:17 INFO None 5251639: status RUNNING/PENDING
2026-07-22 14:24:17 INFO None 5251641: status RUNNING/PENDING
2026-07-22 14:24:17 INFO None 5251644: status RUNNING/PENDING
2026-07-22 14:24:17 INFO None 5251651: status RUNNING/PENDING
2026-07-22 14:24:17 INFO None 5251656: status RUNNING/PENDING
2026-07-22 14:24:17 INFO None 5251659: status RUNNING/PENDING
2026-07-22 14:24:17 INFO None 5251668: status RUNNING/PENDING
2026-07-22 14:24:17 INFO None 5251677: status RUNNING/PENDING
2026-07-22 14:24:17 INFO None 5251685: status RUNNING/PENDING
2026-07-22 14:24:17 INFO None 5251691: status RUNNING/PENDING
2026-07-22 14:24:17 INFO None 5251695: status RUNNING/PENDING
2026-07-22 14:24:17 INFO Jobs still running: ['5251627', '5251631', '5251632', '5251634', '5251639', '5251641', '5251644', '5251651', '5251656', '5251659', '5251668', '5251677', '5251685', '5251691', '5251695']. Waiting...
2026-07-22 14:24:32 INFO None 5251627: status RUNNING/PENDING
2026-07-22 14:24:32 INFO None 5251631: status RUNNING/PENDING
2026-07-22 14:24:32 INFO None 5251632: status RUNNING/PENDING
2026-07-22 14:24:32 INFO None 5251634: status RUNNING/PENDING
2026-07-22 14:24:32 INFO None 5251639: status RUNNING/PENDING
2026-07-22 14:24:32 INFO None 5251641: status RUNNING/PENDING
2026-07-22 14:24:32 INFO None 5251644: status RUNNING/PENDING
2026-07-22 14:24:32 INFO None 5251651: status RUNNING/PENDING
2026-07-22 14:24:32 INFO None 5251656: status RUNNING/PENDING
2026-07-22 14:24:32 INFO None 5251659: status RUNNING/PENDING
2026-07-22 14:24:32 INFO None 5251668: status RUNNING/PENDING
2026-07-22 14:24:32 INFO None 5251677: status RUNNING/PENDING
2026-07-22 14:24:32 INFO None 5251685: status RUNNING/PENDING
2026-07-22 14:24:32 INFO None 5251691: status RUNNING/PENDING
2026-07-22 14:24:32 INFO None 5251695: status RUNNING/PENDING
2026-07-22 14:24:32 INFO Jobs still running: ['5251627', '5251631', '5251632', '5251634', '5251639', '5251641', '5251644', '5251651', '5251656', '5251659', '5251668', '5251677', '5251685', '5251691', '5251695']. Waiting...
2026-07-22 14:24:47 INFO None 5251627: status RUNNING/PENDING
2026-07-22 14:24:47 INFO None 5251631: status RUNNING/PENDING
2026-07-22 14:24:47 INFO None 5251632: status RUNNING/PENDING
2026-07-22 14:24:47 INFO None 5251634: status RUNNING/PENDING
2026-07-22 14:24:47 INFO None 5251639: status RUNNING/PENDING
2026-07-22 14:24:47 INFO None 5251641: status RUNNING/PENDING
2026-07-22 14:24:47 INFO None 5251644: status RUNNING/PENDING
2026-07-22 14:24:47 INFO None 5251651: status RUNNING/PENDING
2026-07-22 14:24:47 INFO None 5251656: status RUNNING/PENDING
2026-07-22 14:24:48 INFO None 5251659: status RUNNING/PENDING
2026-07-22 14:24:48 INFO None 5251668: status RUNNING/PENDING
2026-07-22 14:24:48 INFO None 5251677: status RUNNING/PENDING
2026-07-22 14:24:48 INFO None 5251685: status RUNNING/PENDING
2026-07-22 14:24:48 INFO None 5251691: status RUNNING/PENDING
2026-07-22 14:24:48 INFO None 5251695: status RUNNING/PENDING
2026-07-22 14:24:48 INFO Jobs still running: ['5251627', '5251631', '5251632', '5251634', '5251639', '5251641', '5251644', '5251651', '5251656', '5251659', '5251668', '5251677', '5251685', '5251691', '5251695']. Waiting...
2026-07-22 14:25:03 INFO None 5251627: status RUNNING/PENDING
2026-07-22 14:25:03 INFO None 5251631: status RUNNING/PENDING
2026-07-22 14:25:03 INFO None 5251632: status RUNNING/PENDING
2026-07-22 14:25:03 INFO None 5251634: status RUNNING/PENDING
2026-07-22 14:25:03 INFO None 5251639: status RUNNING/PENDING
2026-07-22 14:25:03 INFO None 5251641: status RUNNING/PENDING
2026-07-22 14:25:03 INFO None 5251644: status RUNNING/PENDING
2026-07-22 14:25:03 INFO None 5251651: status RUNNING/PENDING
2026-07-22 14:25:03 INFO None 5251656: status RUNNING/PENDING
2026-07-22 14:25:03 INFO None 5251659: status RUNNING/PENDING
2026-07-22 14:25:03 INFO None 5251668: status RUNNING/PENDING
2026-07-22 14:25:03 INFO None 5251677: status RUNNING/PENDING
2026-07-22 14:25:03 INFO None 5251685: status RUNNING/PENDING
2026-07-22 14:25:03 INFO None 5251691: status RUNNING/PENDING
2026-07-22 14:25:03 INFO None 5251695: status RUNNING/PENDING
2026-07-22 14:25:03 INFO Jobs still running: ['5251627', '5251631', '5251632', '5251634', '5251639', '5251641', '5251644', '5251651', '5251656', '5251659', '5251668', '5251677', '5251685', '5251691', '5251695']. Waiting...
2026-07-22 14:25:20 INFO None 5251627: status RUNNING/PENDING
2026-07-22 14:25:20 INFO None 5251631: status RUNNING/PENDING
2026-07-22 14:25:20 INFO None 5251632: status RUNNING/PENDING
2026-07-22 14:25:20 INFO None 5251634: status RUNNING/PENDING
2026-07-22 14:25:20 INFO None 5251639: status RUNNING/PENDING
2026-07-22 14:25:20 INFO None 5251641: status RUNNING/PENDING
2026-07-22 14:25:20 INFO None 5251644: status RUNNING/PENDING
2026-07-22 14:25:20 INFO None 5251651: status RUNNING/PENDING
2026-07-22 14:25:20 INFO None 5251656: status RUNNING/PENDING
2026-07-22 14:25:20 INFO None 5251659: status RUNNING/PENDING
2026-07-22 14:25:20 INFO None 5251668: status RUNNING/PENDING
2026-07-22 14:25:20 INFO None 5251677: status RUNNING/PENDING
2026-07-22 14:25:20 INFO None 5251685: status RUNNING/PENDING
2026-07-22 14:25:20 INFO None 5251691: status RUNNING/PENDING
2026-07-22 14:25:20 INFO None 5251695: status RUNNING/PENDING
2026-07-22 14:25:20 INFO Jobs still running: ['5251627', '5251631', '5251632', '5251634', '5251639', '5251641', '5251644', '5251651', '5251656', '5251659', '5251668', '5251677', '5251685', '5251691', '5251695']. Waiting...
2026-07-22 14:25:35 INFO None 5251627: status RUNNING/PENDING
2026-07-22 14:25:35 INFO None 5251631: status RUNNING/PENDING
2026-07-22 14:25:35 INFO None 5251632: status RUNNING/PENDING
2026-07-22 14:25:35 INFO None 5251634: status RUNNING/PENDING
2026-07-22 14:25:35 INFO None 5251639: status RUNNING/PENDING
2026-07-22 14:25:35 INFO None 5251641: status RUNNING/PENDING
2026-07-22 14:25:37 INFO None 5251644: status RUNNING/PENDING
2026-07-22 14:25:37 INFO None 5251651: status RUNNING/PENDING
2026-07-22 14:25:38 INFO None 5251656: status RUNNING/PENDING
2026-07-22 14:25:38 INFO None 5251659: status RUNNING/PENDING
2026-07-22 14:25:38 INFO None 5251668: status RUNNING/PENDING
2026-07-22 14:25:38 INFO None 5251677: status RUNNING/PENDING
2026-07-22 14:25:38 INFO None 5251685: status RUNNING/PENDING
2026-07-22 14:25:38 INFO None 5251691: status RUNNING/PENDING
2026-07-22 14:25:38 INFO None 5251695: status RUNNING/PENDING
2026-07-22 14:25:38 INFO Jobs still running: ['5251627', '5251631', '5251632', '5251634', '5251639', '5251641', '5251644', '5251651', '5251656', '5251659', '5251668', '5251677', '5251685', '5251691', '5251695']. Waiting...
2026-07-22 14:25:53 INFO None 5251627: status RUNNING/PENDING
2026-07-22 14:25:53 INFO None 5251631: status RUNNING/PENDING
2026-07-22 14:25:53 INFO None 5251632: status RUNNING/PENDING
2026-07-22 14:25:53 INFO None 5251634: status RUNNING/PENDING
2026-07-22 14:25:53 INFO None 5251639: status RUNNING/PENDING
2026-07-22 14:25:53 INFO None 5251641: status RUNNING/PENDING
2026-07-22 14:25:53 INFO None 5251644: status RUNNING/PENDING
2026-07-22 14:25:53 INFO None 5251651: status RUNNING/PENDING
2026-07-22 14:25:53 INFO None 5251656: status RUNNING/PENDING
2026-07-22 14:25:53 INFO None 5251659: status RUNNING/PENDING
2026-07-22 14:25:53 INFO None 5251668: status RUNNING/PENDING
2026-07-22 14:25:53 INFO None 5251677: status RUNNING/PENDING
2026-07-22 14:25:53 INFO None 5251685: status RUNNING/PENDING
2026-07-22 14:25:53 INFO None 5251691: status RUNNING/PENDING
2026-07-22 14:25:53 INFO None 5251695: status RUNNING/PENDING
2026-07-22 14:25:53 INFO Jobs still running: ['5251627', '5251631', '5251632', '5251634', '5251639', '5251641', '5251644', '5251651', '5251656', '5251659', '5251668', '5251677', '5251685', '5251691', '5251695']. Waiting...
2026-07-22 14:26:08 INFO None 5251627: status RUNNING/PENDING
2026-07-22 14:26:08 INFO None 5251631: status RUNNING/PENDING
2026-07-22 14:26:08 INFO None 5251632: status RUNNING/PENDING
2026-07-22 14:26:08 INFO None 5251634: status RUNNING/PENDING
2026-07-22 14:26:08 INFO None 5251639: status RUNNING/PENDING
2026-07-22 14:26:08 INFO None 5251641: status RUNNING/PENDING
2026-07-22 14:26:09 INFO None 5251644: status RUNNING/PENDING
2026-07-22 14:26:09 INFO None 5251651: status RUNNING/PENDING
2026-07-22 14:26:09 INFO None 5251656: status RUNNING/PENDING
2026-07-22 14:26:09 INFO None 5251659: status RUNNING/PENDING
2026-07-22 14:26:09 INFO None 5251668: status RUNNING/PENDING
2026-07-22 14:26:09 INFO None 5251677: status RUNNING/PENDING
2026-07-22 14:26:09 INFO None 5251685: status RUNNING/PENDING
2026-07-22 14:26:09 INFO None 5251691: status RUNNING/PENDING
2026-07-22 14:26:09 INFO None 5251695: status RUNNING/PENDING
2026-07-22 14:26:09 INFO Jobs still running: ['5251627', '5251631', '5251632', '5251634', '5251639', '5251641', '5251644', '5251651', '5251656', '5251659', '5251668', '5251677', '5251685', '5251691', '5251695']. Waiting...
2026-07-22 14:26:24 INFO None 5251627: status RUNNING/PENDING
2026-07-22 14:26:24 INFO None 5251631: status RUNNING/PENDING
2026-07-22 14:26:24 INFO None 5251632: status RUNNING/PENDING
2026-07-22 14:26:24 INFO None 5251634: status RUNNING/PENDING
2026-07-22 14:26:24 INFO None 5251639: status RUNNING/PENDING
2026-07-22 14:26:24 INFO None 5251641: status RUNNING/PENDING
2026-07-22 14:26:24 INFO None 5251644: status RUNNING/PENDING
2026-07-22 14:26:24 INFO None 5251651: status RUNNING/PENDING
2026-07-22 14:26:24 INFO None 5251656: status RUNNING/PENDING
2026-07-22 14:26:24 INFO None 5251659: status RUNNING/PENDING
2026-07-22 14:26:24 INFO None 5251668: status RUNNING/PENDING
2026-07-22 14:26:24 INFO None 5251677: status RUNNING/PENDING
2026-07-22 14:26:24 INFO None 5251685: status RUNNING/PENDING
2026-07-22 14:26:24 INFO None 5251691: status RUNNING/PENDING
2026-07-22 14:26:24 INFO None 5251695: status RUNNING/PENDING
2026-07-22 14:26:24 INFO Jobs still running: ['5251627', '5251631', '5251632', '5251634', '5251639', '5251641', '5251644', '5251651', '5251656', '5251659', '5251668', '5251677', '5251685', '5251691', '5251695']. Waiting...
2026-07-22 14:26:41 INFO None 5251627: status RUNNING/PENDING
2026-07-22 14:26:41 INFO None 5251631: status RUNNING/PENDING
2026-07-22 14:26:41 INFO None 5251632: status RUNNING/PENDING
2026-07-22 14:26:41 INFO None 5251634: status RUNNING/PENDING
2026-07-22 14:26:41 INFO None 5251639: status RUNNING/PENDING
2026-07-22 14:26:41 INFO None 5251641: status FINISHED
2026-07-22 14:26:41 INFO None 5251644: status RUNNING/PENDING
2026-07-22 14:26:41 INFO None 5251651: status FINISHED
2026-07-22 14:26:41 INFO None 5251656: status RUNNING/PENDING
2026-07-22 14:26:42 INFO None 5251659: status RUNNING/PENDING
2026-07-22 14:26:42 INFO None 5251668: status RUNNING/PENDING
2026-07-22 14:26:42 INFO None 5251677: status RUNNING/PENDING
2026-07-22 14:26:42 INFO None 5251685: status RUNNING/PENDING
2026-07-22 14:26:42 INFO None 5251691: status RUNNING/PENDING
2026-07-22 14:26:42 INFO None 5251695: status RUNNING/PENDING
2026-07-22 14:26:42 INFO Jobs still running: ['5251627', '5251631', '5251632', '5251634', '5251639', '5251644', '5251656', '5251659', '5251668', '5251677', '5251685', '5251691', '5251695']. Waiting...
2026-07-22 14:26:57 INFO None 5251627: status RUNNING/PENDING
2026-07-22 14:26:57 INFO None 5251631: status RUNNING/PENDING
2026-07-22 14:26:57 INFO None 5251632: status RUNNING/PENDING
2026-07-22 14:26:57 INFO None 5251634: status RUNNING/PENDING
2026-07-22 14:26:59 INFO None 5251639: status RUNNING/PENDING
2026-07-22 14:26:59 INFO None 5251641: status FINISHED
2026-07-22 14:26:59 INFO None 5251644: status RUNNING/PENDING
2026-07-22 14:26:59 INFO None 5251651: status FINISHED
2026-07-22 14:26:59 INFO None 5251656: status RUNNING/PENDING
2026-07-22 14:26:59 INFO None 5251659: status RUNNING/PENDING
2026-07-22 14:26:59 INFO None 5251668: status RUNNING/PENDING
2026-07-22 14:26:59 INFO None 5251677: status RUNNING/PENDING
2026-07-22 14:26:59 INFO None 5251685: status RUNNING/PENDING
2026-07-22 14:26:59 INFO None 5251691: status RUNNING/PENDING
2026-07-22 14:26:59 INFO None 5251695: status RUNNING/PENDING
2026-07-22 14:26:59 INFO Jobs still running: ['5251627', '5251631', '5251632', '5251634', '5251639', '5251644', '5251656', '5251659', '5251668', '5251677', '5251685', '5251691', '5251695']. Waiting...
2026-07-22 14:27:15 INFO None 5251627: status RUNNING/PENDING
2026-07-22 14:27:15 INFO None 5251631: status RUNNING/PENDING
2026-07-22 14:27:15 INFO None 5251632: status RUNNING/PENDING
2026-07-22 14:27:15 INFO None 5251634: status RUNNING/PENDING
2026-07-22 14:27:15 INFO None 5251639: status RUNNING/PENDING
2026-07-22 14:27:15 INFO None 5251641: status FINISHED
2026-07-22 14:27:15 INFO None 5251644: status RUNNING/PENDING
2026-07-22 14:27:15 INFO None 5251651: status FINISHED
2026-07-22 14:27:15 INFO None 5251656: status RUNNING/PENDING
2026-07-22 14:27:15 INFO None 5251659: status RUNNING/PENDING
2026-07-22 14:27:15 INFO None 5251668: status RUNNING/PENDING
2026-07-22 14:27:15 INFO None 5251677: status RUNNING/PENDING
2026-07-22 14:27:15 INFO None 5251685: status RUNNING/PENDING
2026-07-22 14:27:15 INFO None 5251691: status RUNNING/PENDING
2026-07-22 14:27:15 INFO None 5251695: status RUNNING/PENDING
2026-07-22 14:27:15 INFO Jobs still running: ['5251627', '5251631', '5251632', '5251634', '5251639', '5251644', '5251656', '5251659', '5251668', '5251677', '5251685', '5251691', '5251695']. Waiting...
2026-07-22 14:27:30 INFO None 5251627: status RUNNING/PENDING
2026-07-22 14:27:30 INFO None 5251631: status RUNNING/PENDING
2026-07-22 14:27:30 INFO None 5251632: status RUNNING/PENDING
2026-07-22 14:27:30 INFO None 5251634: status RUNNING/PENDING
2026-07-22 14:27:30 INFO None 5251639: status RUNNING/PENDING
2026-07-22 14:27:30 INFO None 5251641: status FINISHED
2026-07-22 14:27:31 INFO None 5251644: status RUNNING/PENDING
2026-07-22 14:27:31 INFO None 5251651: status FINISHED
2026-07-22 14:27:31 INFO None 5251656: status RUNNING/PENDING
2026-07-22 14:27:31 INFO None 5251659: status RUNNING/PENDING
2026-07-22 14:27:31 INFO None 5251668: status RUNNING/PENDING
2026-07-22 14:27:31 INFO None 5251677: status FINISHED
2026-07-22 14:27:31 INFO None 5251685: status RUNNING/PENDING
2026-07-22 14:27:31 INFO None 5251691: status RUNNING/PENDING
2026-07-22 14:27:31 INFO None 5251695: status RUNNING/PENDING
2026-07-22 14:27:31 INFO Jobs still running: ['5251627', '5251631', '5251632', '5251634', '5251639', '5251644', '5251656', '5251659', '5251668', '5251685', '5251691', '5251695']. Waiting...
2026-07-22 14:27:46 INFO None 5251627: status RUNNING/PENDING
2026-07-22 14:27:46 INFO None 5251631: status RUNNING/PENDING
2026-07-22 14:27:46 INFO None 5251632: status RUNNING/PENDING
2026-07-22 14:27:46 INFO None 5251634: status RUNNING/PENDING
2026-07-22 14:27:46 INFO None 5251639: status RUNNING/PENDING
2026-07-22 14:27:46 INFO None 5251641: status FINISHED
2026-07-22 14:27:46 INFO None 5251644: status FINISHED
2026-07-22 14:27:46 INFO None 5251651: status FINISHED
2026-07-22 14:27:46 INFO None 5251656: status FINISHED
2026-07-22 14:27:46 INFO None 5251659: status FINISHED
2026-07-22 14:27:46 INFO None 5251668: status FINISHED
2026-07-22 14:27:46 INFO None 5251677: status FINISHED
2026-07-22 14:27:46 INFO None 5251685: status FINISHED
2026-07-22 14:27:46 INFO None 5251691: status FINISHED
2026-07-22 14:27:46 INFO None 5251695: status RUNNING/PENDING
2026-07-22 14:27:46 INFO Jobs still running: ['5251627', '5251631', '5251632', '5251634', '5251639', '5251695']. Waiting...
2026-07-22 14:28:01 INFO None 5251627: status FINISHED
2026-07-22 14:28:01 INFO None 5251631: status FINISHED
2026-07-22 14:28:01 INFO None 5251632: status FINISHED
2026-07-22 14:28:02 INFO None 5251634: status FINISHED
2026-07-22 14:28:02 INFO None 5251639: status FINISHED
2026-07-22 14:28:02 INFO None 5251641: status FINISHED
2026-07-22 14:28:02 INFO None 5251644: status FINISHED
2026-07-22 14:28:02 INFO None 5251651: status FINISHED
2026-07-22 14:28:02 INFO None 5251656: status FINISHED
2026-07-22 14:28:02 INFO None 5251659: status FINISHED
2026-07-22 14:28:02 INFO None 5251668: status FINISHED
2026-07-22 14:28:02 INFO None 5251677: status FINISHED
2026-07-22 14:28:02 INFO None 5251685: status FINISHED
2026-07-22 14:28:02 INFO None 5251691: status FINISHED
2026-07-22 14:28:02 INFO None 5251695: status FINISHED
2026-07-22 14:28:02 INFO Jobs ['5251627', '5251631', '5251632', '5251634', '5251639', '5251641', '5251644', '5251651', '5251656', '5251659', '5251668', '5251677', '5251685', '5251691', '5251695'] have finished
2026-07-22 14:28:02 INFO Checking restart files were created ...
2026-07-22 14:28:02 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020600_1_ENS1.nc(668832435 bytes)
2026-07-22 14:28:02 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020600_1_ENS2.nc(668832435 bytes)
2026-07-22 14:28:02 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020600_1_ENS3.nc(668832435 bytes)
2026-07-22 14:28:02 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020600_1_ENS4.nc(668832435 bytes)
2026-07-22 14:28:02 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020600_1_ENS5.nc(668832435 bytes)
2026-07-22 14:28:02 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020600_1_ENS6.nc(668832435 bytes)
2026-07-22 14:28:02 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020600_1_ENS7.nc(668832435 bytes)
2026-07-22 14:28:02 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020600_1_ENS8.nc(668832435 bytes)
2026-07-22 14:28:02 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020600_1_ENS9.nc(668832435 bytes)
2026-07-22 14:28:02 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020600_1_ENS10.nc(668832435 bytes)
2026-07-22 14:28:02 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020600_1_ENS11.nc(668832435 bytes)
2026-07-22 14:28:02 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020600_1_ENS12.nc(668832435 bytes)
2026-07-22 14:28:02 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020600_1_ENS13.nc(668832435 bytes)
2026-07-22 14:28:02 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020600_1_ENS14.nc(668832435 bytes)
2026-07-22 14:28:02 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020600_1_ENS15.nc(668832435 bytes)
2026-07-22 14:28:02 INFO  Run_model() completed successfully.
2026-07-22 14:28:02 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 00:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 14:28:02 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 01:00:00 days=153072 seconds=3600
2026-07-22 14:28:02 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-22 14:28:02 INFO [TIME] increment current_time 2020-02-06 00:00:00 -> 2020-02-06 01:00:00
2026-07-22 14:28:02 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 14:28:02 INFO ---------->>> Running process_satellite_data()
2026-07-22 14:28:02 INFO [DART] No satellite data found, skipping assimilation
2026-07-22 14:28:02 INFO after_assimilation() skipped
2026-07-22 14:28:02 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-22 14:28:02 INFO [TIME] step_end current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 14:28:02 INFO [TIME] step_start current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 14:28:02 INFO [TIME] window start=2020-02-06 01:00:00 end=2020-02-06 09:00:00 run_hours=8 has_assimilation=True
2026-07-22 14:28:02 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:28:04 INFO Hourly dataset computed and listing created
2026-07-22 14:28:20 INFO Hourly dataset computed
2026-07-22 14:28:20 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:28:22 INFO Hourly dataset computed and listing created
2026-07-22 14:28:31 INFO Hourly dataset computed
2026-07-22 14:28:31 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:28:33 INFO Hourly dataset computed and listing created
2026-07-22 14:28:47 INFO Hourly dataset computed
2026-07-22 14:28:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:28:49 INFO Hourly dataset computed and listing created
2026-07-22 14:28:55 INFO Hourly dataset computed
2026-07-22 14:28:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:28:56 INFO Hourly dataset computed and listing created
2026-07-22 14:28:58 INFO Hourly dataset computed
2026-07-22 14:28:58 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:28:59 INFO Hourly dataset computed and listing created
2026-07-22 14:29:01 INFO Hourly dataset computed
2026-07-22 14:29:02 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:29:03 INFO Hourly dataset computed and listing created
2026-07-22 14:29:05 INFO Hourly dataset computed
2026-07-22 14:29:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:29:06 INFO Hourly dataset computed and listing created
2026-07-22 14:29:08 INFO Hourly dataset computed
2026-07-22 14:29:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:29:10 INFO Hourly dataset computed and listing created
2026-07-22 14:29:12 INFO Hourly dataset computed
2026-07-22 14:29:12 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:29:13 INFO Hourly dataset computed and listing created
2026-07-22 14:29:15 INFO Hourly dataset computed
2026-07-22 14:29:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:29:16 INFO Hourly dataset computed and listing created
2026-07-22 14:29:18 INFO Hourly dataset computed
2026-07-22 14:29:18 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:29:19 INFO Hourly dataset computed and listing created
2026-07-22 14:29:21 INFO Hourly dataset computed
2026-07-22 14:29:21 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:29:22 INFO Hourly dataset computed and listing created
2026-07-22 14:29:24 INFO Hourly dataset computed
2026-07-22 14:29:24 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:29:26 INFO Hourly dataset computed and listing created
2026-07-22 14:29:28 INFO Hourly dataset computed
2026-07-22 14:29:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:29:29 INFO Hourly dataset computed and listing created
2026-07-22 14:29:31 INFO Hourly dataset computed
2026-07-22 14:29:31 INFO ---------->>> Running CHIMERE model from 2020-02-06 01:00:00 to 2020-02-06 09:00:00
2026-07-22 14:29:31 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:29:31 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1
2026-07-22 14:29:31 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020600_1_ENS1.nc
2026-07-22 14:29:31 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-22 14:29:31 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:29:31 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-22 14:29:31 INFO Queuing job for member 1...
2026-07-22 14:29:31 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:29:31 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-22 14:29:32 INFO Found: ['5251892']
2026-07-22 14:29:37 INFO [TGCC-IRENE] Submitted job with ID:['5251892']
2026-07-22 14:29:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:29:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2
2026-07-22 14:29:37 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020600_1_ENS2.nc
2026-07-22 14:29:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-22 14:29:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:29:37 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-22 14:29:37 INFO Queuing job for member 2...
2026-07-22 14:29:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:29:37 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-22 14:29:40 INFO Found: ['5251895']
2026-07-22 14:29:45 INFO [TGCC-IRENE] Submitted job with ID:['5251895']
2026-07-22 14:29:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:29:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3
2026-07-22 14:29:45 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020600_1_ENS3.nc
2026-07-22 14:29:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-22 14:29:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:29:45 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-22 14:29:45 INFO Queuing job for member 3...
2026-07-22 14:29:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:29:45 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-22 14:29:47 INFO Found: ['5251896']
2026-07-22 14:29:52 INFO [TGCC-IRENE] Submitted job with ID:['5251896']
2026-07-22 14:29:52 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:29:52 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4
2026-07-22 14:29:52 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020600_1_ENS4.nc
2026-07-22 14:29:52 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-22 14:29:52 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:29:52 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-22 14:29:52 INFO Queuing job for member 4...
2026-07-22 14:29:52 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:29:52 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-22 14:29:55 INFO Found: ['5251900']
2026-07-22 14:30:00 INFO [TGCC-IRENE] Submitted job with ID:['5251900']
2026-07-22 14:30:00 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:30:00 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5
2026-07-22 14:30:00 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020600_1_ENS5.nc
2026-07-22 14:30:00 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-22 14:30:00 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:30:00 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-22 14:30:00 INFO Queuing job for member 5...
2026-07-22 14:30:00 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:30:00 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-22 14:30:01 INFO Found: ['5251905']
2026-07-22 14:30:06 INFO [TGCC-IRENE] Submitted job with ID:['5251905']
2026-07-22 14:30:06 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:30:06 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6
2026-07-22 14:30:06 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020600_1_ENS6.nc
2026-07-22 14:30:06 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-22 14:30:06 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:30:06 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-22 14:30:06 INFO Queuing job for member 6...
2026-07-22 14:30:06 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:30:06 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-22 14:30:06 INFO Found: ['5251908']
2026-07-22 14:30:11 INFO [TGCC-IRENE] Submitted job with ID:['5251908']
2026-07-22 14:30:11 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:30:11 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7
2026-07-22 14:30:11 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020600_1_ENS7.nc
2026-07-22 14:30:11 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-22 14:30:11 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:30:11 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-22 14:30:11 INFO Queuing job for member 7...
2026-07-22 14:30:11 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:30:11 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-22 14:30:12 INFO Found: ['5251914']
2026-07-22 14:30:17 INFO [TGCC-IRENE] Submitted job with ID:['5251914']
2026-07-22 14:30:17 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:30:17 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8
2026-07-22 14:30:17 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020600_1_ENS8.nc
2026-07-22 14:30:17 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-22 14:30:17 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:30:17 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-22 14:30:17 INFO Queuing job for member 8...
2026-07-22 14:30:17 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:30:17 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-22 14:30:18 INFO Found: ['5251917']
2026-07-22 14:30:23 INFO [TGCC-IRENE] Submitted job with ID:['5251917']
2026-07-22 14:30:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:30:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9
2026-07-22 14:30:23 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020600_1_ENS9.nc
2026-07-22 14:30:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-22 14:30:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:30:23 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-22 14:30:23 INFO Queuing job for member 9...
2026-07-22 14:30:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:30:23 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-22 14:30:24 INFO Found: ['5251918']
2026-07-22 14:30:29 INFO [TGCC-IRENE] Submitted job with ID:['5251918']
2026-07-22 14:30:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:30:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10
2026-07-22 14:30:29 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020600_1_ENS10.nc
2026-07-22 14:30:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-22 14:30:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:30:29 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-22 14:30:29 INFO Queuing job for member 10...
2026-07-22 14:30:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:30:29 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-22 14:30:29 INFO Found: ['5251924']
2026-07-22 14:30:34 INFO [TGCC-IRENE] Submitted job with ID:['5251924']
2026-07-22 14:30:34 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:30:34 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11
2026-07-22 14:30:34 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020600_1_ENS11.nc
2026-07-22 14:30:34 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-22 14:30:34 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:30:34 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-22 14:30:35 INFO Queuing job for member 11...
2026-07-22 14:30:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:30:35 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-22 14:30:35 INFO Found: ['5251936']
2026-07-22 14:30:40 INFO [TGCC-IRENE] Submitted job with ID:['5251936']
2026-07-22 14:30:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:30:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12
2026-07-22 14:30:40 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020600_1_ENS12.nc
2026-07-22 14:30:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-22 14:30:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:30:40 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-22 14:30:40 INFO Queuing job for member 12...
2026-07-22 14:30:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:30:40 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-22 14:30:41 INFO Found: ['5251938']
2026-07-22 14:30:46 INFO [TGCC-IRENE] Submitted job with ID:['5251938']
2026-07-22 14:30:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:30:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13
2026-07-22 14:30:46 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020600_1_ENS13.nc
2026-07-22 14:30:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-22 14:30:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:30:46 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-22 14:30:46 INFO Queuing job for member 13...
2026-07-22 14:30:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:30:46 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-22 14:30:49 INFO Found: ['5251939']
2026-07-22 14:30:54 INFO [TGCC-IRENE] Submitted job with ID:['5251939']
2026-07-22 14:30:54 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:30:54 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14
2026-07-22 14:30:54 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020600_1_ENS14.nc
2026-07-22 14:30:54 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-22 14:30:54 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:30:54 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-22 14:30:54 INFO Queuing job for member 14...
2026-07-22 14:30:54 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:30:54 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-22 14:30:56 INFO Found: ['5251941']
2026-07-22 14:31:01 INFO [TGCC-IRENE] Submitted job with ID:['5251941']
2026-07-22 14:31:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:31:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15
2026-07-22 14:31:01 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020600_1_ENS15.nc
2026-07-22 14:31:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-22 14:31:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:31:01 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-22 14:31:01 INFO Queuing job for member 15...
2026-07-22 14:31:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:31:01 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-22 14:31:04 INFO Found: ['5251943']
2026-07-22 14:31:09 INFO [TGCC-IRENE] Submitted job with ID:['5251943']
2026-07-22 14:31:09 INFO Checking job status ...
2026-07-22 14:31:09 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:31:09 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:31:09 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:31:09 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:31:09 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:31:09 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:31:09 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:31:09 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:31:09 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:31:09 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:31:09 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:31:09 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:31:09 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:31:09 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:31:11 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:31:11 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:31:26 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:31:26 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:31:26 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:31:26 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:31:26 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:31:26 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:31:26 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:31:26 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:31:26 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:31:26 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:31:26 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:31:26 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:31:26 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:31:26 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:31:26 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:31:26 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:31:42 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:31:42 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:31:42 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:31:42 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:31:42 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:31:42 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:31:42 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:31:42 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:31:42 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:31:42 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:31:42 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:31:42 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:31:42 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:31:42 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:31:42 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:31:42 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:31:57 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:31:57 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:31:57 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:31:57 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:31:57 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:31:57 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:31:57 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:31:57 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:31:57 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:31:57 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:31:57 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:31:57 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:31:57 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:31:57 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:31:57 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:31:57 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:32:12 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:32:12 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:32:12 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:32:12 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:32:12 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:32:12 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:32:12 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:32:12 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:32:12 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:32:12 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:32:15 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:32:15 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:32:15 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:32:15 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:32:15 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:32:15 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:32:30 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:32:30 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:32:30 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:32:30 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:32:30 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:32:30 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:32:30 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:32:30 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:32:30 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:32:30 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:32:30 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:32:30 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:32:30 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:32:30 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:32:30 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:32:30 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:32:45 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:32:45 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:32:45 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:32:45 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:32:47 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:32:47 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:32:47 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:32:47 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:32:47 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:32:47 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:32:47 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:32:47 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:32:47 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:32:47 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:32:47 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:32:47 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:33:02 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:33:02 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:33:02 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:33:02 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:33:03 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:33:03 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:33:03 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:33:03 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:33:03 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:33:03 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:33:03 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:33:03 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:33:03 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:33:03 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:33:03 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:33:03 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:33:18 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:33:18 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:33:18 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:33:18 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:33:18 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:33:18 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:33:18 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:33:18 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:33:18 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:33:18 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:33:18 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:33:18 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:33:18 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:33:18 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:33:18 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:33:18 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:33:35 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:33:35 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:33:35 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:33:35 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:33:35 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:33:35 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:33:35 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:33:35 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:33:35 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:33:35 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:33:35 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:33:35 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:33:35 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:33:35 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:33:35 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:33:35 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:33:50 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:33:50 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:33:50 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:33:50 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:33:50 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:33:50 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:33:52 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:33:52 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:33:52 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:33:52 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:33:52 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:33:52 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:33:52 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:33:53 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:33:53 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:33:53 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:34:08 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:34:08 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:34:08 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:34:08 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:34:08 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:34:08 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:34:08 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:34:08 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:34:08 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:34:08 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:34:08 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:34:08 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:34:08 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:34:08 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:34:08 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:34:08 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:34:23 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:34:23 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:34:23 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:34:23 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:34:23 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:34:23 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:34:23 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:34:23 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:34:23 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:34:23 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:34:23 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:34:23 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:34:23 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:34:23 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:34:23 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:34:23 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:34:38 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:34:38 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:34:38 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:34:38 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:34:38 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:34:38 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:34:38 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:34:38 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:34:38 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:34:38 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:34:38 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:34:39 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:34:39 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:34:39 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:34:39 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:34:39 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:34:54 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:34:54 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:34:54 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:34:54 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:34:54 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:34:54 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:34:54 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:34:54 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:34:54 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:34:54 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:34:54 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:34:54 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:34:54 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:34:54 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:34:54 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:34:54 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:35:09 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:35:09 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:35:09 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:35:09 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:35:09 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:35:09 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:35:09 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:35:11 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:35:11 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:35:11 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:35:11 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:35:11 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:35:11 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:35:11 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:35:11 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:35:11 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:35:26 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:35:26 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:35:26 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:35:26 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:35:26 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:35:26 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:35:26 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:35:26 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:35:26 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:35:27 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:35:27 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:35:27 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:35:27 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:35:27 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:35:27 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:35:27 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:35:42 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:35:42 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:35:42 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:35:42 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:35:42 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:35:42 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:35:42 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:35:42 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:35:42 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:35:42 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:35:42 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:35:42 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:35:42 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:35:42 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:35:42 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:35:42 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:35:57 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:35:57 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:35:57 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:35:57 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:35:57 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:35:57 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:35:57 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:35:57 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:35:57 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:35:57 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:35:57 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:35:57 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:35:57 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:35:57 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:35:57 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:35:57 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:36:12 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:36:12 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:36:12 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:36:12 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:36:12 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:36:12 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:36:12 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:36:13 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:36:13 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:36:13 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:36:13 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:36:13 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:36:13 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:36:13 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:36:13 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:36:13 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:36:30 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:36:30 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:36:30 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:36:30 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:36:30 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:36:30 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:36:30 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:36:30 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:36:30 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:36:30 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:36:30 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:36:30 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:36:30 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:36:30 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:36:30 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:36:30 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:36:45 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:36:45 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:36:45 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:36:45 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:36:45 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:36:45 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:36:45 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:36:45 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:36:47 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:36:47 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:36:47 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:36:47 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:36:47 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:36:47 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:36:47 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:36:47 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:37:02 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:37:02 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:37:02 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:37:02 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:37:02 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:37:02 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:37:02 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:37:02 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:37:02 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:37:03 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:37:03 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:37:03 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:37:03 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:37:03 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:37:03 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:37:03 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:37:18 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:37:18 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:37:18 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:37:18 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:37:18 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:37:18 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:37:18 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:37:18 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:37:18 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:37:18 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:37:18 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:37:18 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:37:18 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:37:18 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:37:18 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:37:18 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:37:33 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:37:33 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:37:33 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:37:33 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:37:33 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:37:33 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:37:33 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:37:33 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:37:33 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:37:33 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:37:33 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:37:33 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:37:33 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:37:33 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:37:33 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:37:33 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:37:48 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:37:48 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:37:48 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:37:48 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:37:48 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:37:49 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:37:49 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:37:49 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:37:49 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:37:49 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:37:49 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:37:49 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:37:49 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:37:49 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:37:51 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:37:51 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:38:06 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:38:06 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:38:06 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:38:06 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:38:06 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:38:06 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:38:06 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:38:06 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:38:06 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:38:06 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:38:06 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:38:06 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:38:06 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:38:06 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:38:06 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:38:06 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:38:21 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:38:21 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:38:21 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:38:21 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:38:21 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:38:21 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:38:21 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:38:21 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:38:21 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:38:21 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:38:21 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:38:21 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:38:21 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:38:21 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:38:21 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:38:21 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:38:36 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:38:37 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:38:37 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:38:37 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:38:37 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:38:37 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:38:37 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:38:37 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:38:37 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:38:37 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:38:37 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:38:37 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:38:37 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:38:37 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:38:37 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:38:37 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:38:52 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:38:52 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:38:52 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:38:52 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:38:52 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:38:52 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:38:52 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:38:52 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:38:52 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:38:52 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:38:52 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:38:52 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:38:52 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:38:52 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:38:52 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:38:52 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:39:07 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:39:07 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:39:07 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:39:07 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:39:07 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:39:07 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:39:07 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:39:07 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:39:07 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:39:07 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:39:09 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:39:09 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:39:09 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:39:09 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:39:09 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:39:09 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:39:25 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:39:25 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:39:25 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:39:25 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:39:25 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:39:25 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:39:25 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:39:25 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:39:25 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:39:25 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:39:25 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:39:25 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:39:25 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:39:25 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:39:25 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:39:25 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:39:40 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:39:40 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:39:40 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:39:40 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:39:40 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:39:40 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:39:40 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:39:40 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:39:42 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:39:42 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:39:42 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:39:42 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:39:42 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:39:42 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:39:42 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:39:42 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:39:57 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:39:57 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:39:57 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:39:57 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:39:57 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:39:57 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:39:57 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:39:57 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:39:57 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:39:57 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:39:57 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:39:57 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:39:57 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:39:57 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:39:57 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:39:57 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:40:12 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:40:12 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:40:12 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:40:12 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:40:12 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:40:12 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:40:12 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:40:12 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:40:12 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:40:12 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:40:12 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:40:12 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:40:13 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:40:13 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:40:13 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:40:13 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:40:29 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:40:29 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:40:29 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:40:29 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:40:29 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:40:29 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:40:29 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:40:29 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:40:29 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:40:29 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:40:29 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:40:29 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:40:31 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:40:31 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:40:31 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:40:31 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:40:46 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:40:46 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:40:46 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:40:46 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:40:46 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:40:46 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:40:46 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:40:47 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:40:47 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:40:47 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:40:47 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:40:47 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:40:47 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:40:47 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:40:47 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:40:47 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:41:02 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:41:02 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:41:02 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:41:02 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:41:02 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:41:02 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:41:04 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:41:04 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:41:04 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:41:04 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:41:04 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:41:04 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:41:04 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:41:04 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:41:04 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:41:04 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:41:19 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:41:19 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:41:19 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:41:19 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:41:19 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:41:19 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:41:19 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:41:19 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:41:19 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:41:19 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:41:19 INFO None 5251936: status RUNNING/PENDING
2026-07-22 14:41:19 INFO None 5251938: status RUNNING/PENDING
2026-07-22 14:41:19 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:41:19 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:41:19 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:41:19 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:41:34 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:41:34 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:41:34 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:41:34 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:41:35 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:41:35 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:41:35 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:41:35 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:41:35 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:41:35 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:41:35 INFO None 5251936: status FINISHED
2026-07-22 14:41:35 INFO None 5251938: status FINISHED
2026-07-22 14:41:35 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:41:35 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:41:35 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:41:35 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:41:50 INFO None 5251892: status RUNNING/PENDING
2026-07-22 14:41:50 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:41:50 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:41:50 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:41:50 INFO None 5251905: status RUNNING/PENDING
2026-07-22 14:41:50 INFO None 5251908: status RUNNING/PENDING
2026-07-22 14:41:50 INFO None 5251914: status RUNNING/PENDING
2026-07-22 14:41:50 INFO None 5251917: status RUNNING/PENDING
2026-07-22 14:41:50 INFO None 5251918: status RUNNING/PENDING
2026-07-22 14:41:50 INFO None 5251924: status RUNNING/PENDING
2026-07-22 14:41:50 INFO None 5251936: status FINISHED
2026-07-22 14:41:50 INFO None 5251938: status FINISHED
2026-07-22 14:41:50 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:41:50 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:41:50 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:41:50 INFO Jobs still running: ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:42:06 INFO None 5251892: status FINISHED
2026-07-22 14:42:06 INFO None 5251895: status RUNNING/PENDING
2026-07-22 14:42:06 INFO None 5251896: status RUNNING/PENDING
2026-07-22 14:42:06 INFO None 5251900: status RUNNING/PENDING
2026-07-22 14:42:06 INFO None 5251905: status FINISHED
2026-07-22 14:42:06 INFO None 5251908: status FINISHED
2026-07-22 14:42:06 INFO None 5251914: status FINISHED
2026-07-22 14:42:06 INFO None 5251917: status FINISHED
2026-07-22 14:42:06 INFO None 5251918: status FINISHED
2026-07-22 14:42:06 INFO None 5251924: status FINISHED
2026-07-22 14:42:06 INFO None 5251936: status FINISHED
2026-07-22 14:42:06 INFO None 5251938: status FINISHED
2026-07-22 14:42:06 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:42:06 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:42:06 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:42:06 INFO Jobs still running: ['5251895', '5251896', '5251900', '5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:42:21 INFO None 5251892: status FINISHED
2026-07-22 14:42:21 INFO None 5251895: status FINISHED
2026-07-22 14:42:21 INFO None 5251896: status FINISHED
2026-07-22 14:42:21 INFO None 5251900: status FINISHED
2026-07-22 14:42:21 INFO None 5251905: status FINISHED
2026-07-22 14:42:21 INFO None 5251908: status FINISHED
2026-07-22 14:42:21 INFO None 5251914: status FINISHED
2026-07-22 14:42:21 INFO None 5251917: status FINISHED
2026-07-22 14:42:21 INFO None 5251918: status FINISHED
2026-07-22 14:42:21 INFO None 5251924: status FINISHED
2026-07-22 14:42:21 INFO None 5251936: status FINISHED
2026-07-22 14:42:21 INFO None 5251938: status FINISHED
2026-07-22 14:42:21 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:42:21 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:42:23 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:42:23 INFO Jobs still running: ['5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:42:38 INFO None 5251892: status FINISHED
2026-07-22 14:42:38 INFO None 5251895: status FINISHED
2026-07-22 14:42:38 INFO None 5251896: status FINISHED
2026-07-22 14:42:38 INFO None 5251900: status FINISHED
2026-07-22 14:42:39 INFO None 5251905: status FINISHED
2026-07-22 14:42:39 INFO None 5251908: status FINISHED
2026-07-22 14:42:39 INFO None 5251914: status FINISHED
2026-07-22 14:42:39 INFO None 5251917: status FINISHED
2026-07-22 14:42:39 INFO None 5251918: status FINISHED
2026-07-22 14:42:39 INFO None 5251924: status FINISHED
2026-07-22 14:42:39 INFO None 5251936: status FINISHED
2026-07-22 14:42:39 INFO None 5251938: status FINISHED
2026-07-22 14:42:39 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:42:39 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:42:39 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:42:39 INFO Jobs still running: ['5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:42:54 INFO None 5251892: status FINISHED
2026-07-22 14:42:54 INFO None 5251895: status FINISHED
2026-07-22 14:42:54 INFO None 5251896: status FINISHED
2026-07-22 14:42:54 INFO None 5251900: status FINISHED
2026-07-22 14:42:54 INFO None 5251905: status FINISHED
2026-07-22 14:42:54 INFO None 5251908: status FINISHED
2026-07-22 14:42:54 INFO None 5251914: status FINISHED
2026-07-22 14:42:54 INFO None 5251917: status FINISHED
2026-07-22 14:42:54 INFO None 5251918: status FINISHED
2026-07-22 14:42:54 INFO None 5251924: status FINISHED
2026-07-22 14:42:54 INFO None 5251936: status FINISHED
2026-07-22 14:42:54 INFO None 5251938: status FINISHED
2026-07-22 14:42:54 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:42:54 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:42:54 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:42:54 INFO Jobs still running: ['5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:43:09 INFO None 5251892: status FINISHED
2026-07-22 14:43:09 INFO None 5251895: status FINISHED
2026-07-22 14:43:09 INFO None 5251896: status FINISHED
2026-07-22 14:43:09 INFO None 5251900: status FINISHED
2026-07-22 14:43:09 INFO None 5251905: status FINISHED
2026-07-22 14:43:09 INFO None 5251908: status FINISHED
2026-07-22 14:43:09 INFO None 5251914: status FINISHED
2026-07-22 14:43:09 INFO None 5251917: status FINISHED
2026-07-22 14:43:09 INFO None 5251918: status FINISHED
2026-07-22 14:43:09 INFO None 5251924: status FINISHED
2026-07-22 14:43:09 INFO None 5251936: status FINISHED
2026-07-22 14:43:09 INFO None 5251938: status FINISHED
2026-07-22 14:43:09 INFO None 5251939: status RUNNING/PENDING
2026-07-22 14:43:09 INFO None 5251941: status RUNNING/PENDING
2026-07-22 14:43:09 INFO None 5251943: status RUNNING/PENDING
2026-07-22 14:43:09 INFO Jobs still running: ['5251939', '5251941', '5251943']. Waiting...
2026-07-22 14:43:24 INFO None 5251892: status FINISHED
2026-07-22 14:43:24 INFO None 5251895: status FINISHED
2026-07-22 14:43:24 INFO None 5251896: status FINISHED
2026-07-22 14:43:24 INFO None 5251900: status FINISHED
2026-07-22 14:43:24 INFO None 5251905: status FINISHED
2026-07-22 14:43:24 INFO None 5251908: status FINISHED
2026-07-22 14:43:24 INFO None 5251914: status FINISHED
2026-07-22 14:43:25 INFO None 5251917: status FINISHED
2026-07-22 14:43:25 INFO None 5251918: status FINISHED
2026-07-22 14:43:25 INFO None 5251924: status FINISHED
2026-07-22 14:43:25 INFO None 5251936: status FINISHED
2026-07-22 14:43:25 INFO None 5251938: status FINISHED
2026-07-22 14:43:25 INFO None 5251939: status FINISHED
2026-07-22 14:43:25 INFO None 5251941: status FINISHED
2026-07-22 14:43:25 INFO None 5251943: status FINISHED
2026-07-22 14:43:25 INFO Jobs ['5251892', '5251895', '5251896', '5251900', '5251905', '5251908', '5251914', '5251917', '5251918', '5251924', '5251936', '5251938', '5251939', '5251941', '5251943'] have finished
2026-07-22 14:43:25 INFO Checking restart files were created ...
2026-07-22 14:43:25 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020601_8_ENS1.nc(3005806795 bytes)
2026-07-22 14:43:25 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020601_8_ENS2.nc(3005806795 bytes)
2026-07-22 14:43:25 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020601_8_ENS3.nc(3005806795 bytes)
2026-07-22 14:43:25 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020601_8_ENS4.nc(3005806795 bytes)
2026-07-22 14:43:25 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020601_8_ENS5.nc(3005806795 bytes)
2026-07-22 14:43:25 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020601_8_ENS6.nc(3005806795 bytes)
2026-07-22 14:43:25 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020601_8_ENS7.nc(3005806795 bytes)
2026-07-22 14:43:25 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020601_8_ENS8.nc(3005806795 bytes)
2026-07-22 14:43:25 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020601_8_ENS9.nc(3005806795 bytes)
2026-07-22 14:43:25 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020601_8_ENS10.nc(3005806795 bytes)
2026-07-22 14:43:25 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020601_8_ENS11.nc(3005806795 bytes)
2026-07-22 14:43:25 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020601_8_ENS12.nc(3005806795 bytes)
2026-07-22 14:43:25 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020601_8_ENS13.nc(3005806795 bytes)
2026-07-22 14:43:25 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020601_8_ENS14.nc(3005806795 bytes)
2026-07-22 14:43:25 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020601_8_ENS15.nc(3005806795 bytes)
2026-07-22 14:43:25 INFO  Run_model() completed successfully.
2026-07-22 14:43:25 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 14:43:25 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 09:00:00 days=153072 seconds=32400
2026-07-22 14:43:25 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-22 14:43:25 INFO [TIME] increment current_time 2020-02-06 01:00:00 -> 2020-02-06 09:00:00
2026-07-22 14:43:25 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 14:43:25 INFO ---------->>> Running process_satellite_data()
2026-07-22 14:43:25 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12002.nc
2026-07-22 14:43:25 INFO ---------->>> Running run_obs_converter()
2026-07-22 14:43:25 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_33067_153072.out
2026-07-22 14:43:25 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_33067_153072.out
2026-07-22 14:43:25 INFO ---------->>> Running DART
2026-07-22 14:43:25 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-22 14:43:25 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/chim_ENS1_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/chim_ENS1_2020020609_1_out_toDART.nc
2026-07-22 14:43:25 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/chim_ENS2_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/chim_ENS2_2020020609_1_out_toDART.nc
2026-07-22 14:43:26 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/chim_ENS3_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/chim_ENS3_2020020609_1_out_toDART.nc
2026-07-22 14:43:26 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/chim_ENS4_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/chim_ENS4_2020020609_1_out_toDART.nc
2026-07-22 14:43:26 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/chim_ENS5_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/chim_ENS5_2020020609_1_out_toDART.nc
2026-07-22 14:43:27 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/chim_ENS6_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/chim_ENS6_2020020609_1_out_toDART.nc
2026-07-22 14:43:27 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/chim_ENS7_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/chim_ENS7_2020020609_1_out_toDART.nc
2026-07-22 14:43:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/chim_ENS8_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/chim_ENS8_2020020609_1_out_toDART.nc
2026-07-22 14:43:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/chim_ENS9_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/chim_ENS9_2020020609_1_out_toDART.nc
2026-07-22 14:43:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/chim_ENS10_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/chim_ENS10_2020020609_1_out_toDART.nc
2026-07-22 14:43:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/chim_ENS11_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/chim_ENS11_2020020609_1_out_toDART.nc
2026-07-22 14:43:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/chim_ENS12_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/chim_ENS12_2020020609_1_out_toDART.nc
2026-07-22 14:43:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/chim_ENS13_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/chim_ENS13_2020020609_1_out_toDART.nc
2026-07-22 14:43:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/chim_ENS14_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/chim_ENS14_2020020609_1_out_toDART.nc
2026-07-22 14:43:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/chim_ENS15_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/chim_ENS15_2020020609_1_out_toDART.nc
2026-07-22 14:43:31 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-22 14:43:31 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-22 14:43:31 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-22 14:43:31 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-22 14:43:31 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-22 14:43:31 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-22 14:43:38 INFO Found: []
2026-07-22 14:43:38 INFO No job id returned by command ./run_filter.bsh
2026-07-22 14:43:38 INFO No monitoring will be performed
2026-07-22 14:43:38 INFO Moving DART output files to analysis and preassim directories for date 2020020609 if present ...
2026-07-22 14:43:38 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 14:43:38 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 14:43:38 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 14:43:38 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 14:43:38 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 14:43:38 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 14:43:38 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 14:43:38 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 14:43:38 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 14:43:38 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 14:43:38 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 14:43:38 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 14:43:38 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 14:43:38 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 14:43:38 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 14:43:38 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 14:43:38 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 14:43:38 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 14:43:38 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 14:43:38 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 14:43:38 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 14:43:38 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 14:43:38 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 14:43:38 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 14:43:38 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 14:43:38 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 14:43:38 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 14:43:38 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 14:43:38 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 14:43:38 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 14:43:38 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 14:43:38 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 14:43:38 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 14:43:38 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 14:43:38 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-22 14:43:41 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-22 14:43:41 INFO run_dart() is DONE.
2026-07-22 14:43:41 INFO ---------->>> Running update_pollutant_in_end()
2026-07-22 14:43:42 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens01.nc
2026-07-22 14:43:42 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:43:42 INFO No previous orbit memory found.
2026-07-22 14:43:43 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:43:43 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS1.nc
2026-07-22 14:43:43 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:43:44 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:43:44 INFO No previous orbit memory found.
2026-07-22 14:43:44 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:43:45 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS1.nc
2026-07-22 14:43:45 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:43:45 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 14:43:45 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens02.nc
2026-07-22 14:43:46 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:43:46 INFO No previous orbit memory found.
2026-07-22 14:43:46 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:43:47 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS2.nc
2026-07-22 14:43:47 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:43:47 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:43:47 INFO No previous orbit memory found.
2026-07-22 14:43:48 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:43:48 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS2.nc
2026-07-22 14:43:48 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:43:48 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 14:43:49 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens03.nc
2026-07-22 14:43:49 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:43:49 INFO No previous orbit memory found.
2026-07-22 14:43:50 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:43:50 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS3.nc
2026-07-22 14:43:50 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:43:51 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:43:51 INFO No previous orbit memory found.
2026-07-22 14:43:52 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:43:52 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS3.nc
2026-07-22 14:43:52 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:43:52 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 14:43:53 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens04.nc
2026-07-22 14:43:53 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:43:53 INFO No previous orbit memory found.
2026-07-22 14:43:54 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:43:54 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS4.nc
2026-07-22 14:43:54 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:43:55 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:43:55 INFO No previous orbit memory found.
2026-07-22 14:43:55 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:43:56 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS4.nc
2026-07-22 14:43:56 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:43:56 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 14:43:56 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens05.nc
2026-07-22 14:43:57 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:43:57 INFO No previous orbit memory found.
2026-07-22 14:43:57 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:43:58 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS5.nc
2026-07-22 14:43:58 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:43:58 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:43:58 INFO No previous orbit memory found.
2026-07-22 14:43:59 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:43:59 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS5.nc
2026-07-22 14:43:59 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:43:59 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 14:44:00 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens06.nc
2026-07-22 14:44:00 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:44:00 INFO No previous orbit memory found.
2026-07-22 14:44:01 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:44:01 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS6.nc
2026-07-22 14:44:01 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:44:01 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:44:01 INFO No previous orbit memory found.
2026-07-22 14:44:02 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:44:02 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS6.nc
2026-07-22 14:44:03 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:44:03 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 14:44:03 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens07.nc
2026-07-22 14:44:03 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:44:03 INFO No previous orbit memory found.
2026-07-22 14:44:04 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:44:04 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS7.nc
2026-07-22 14:44:04 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:44:05 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:44:05 INFO No previous orbit memory found.
2026-07-22 14:44:05 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:44:06 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS7.nc
2026-07-22 14:44:06 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:44:06 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 14:44:06 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens08.nc
2026-07-22 14:44:07 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:44:07 INFO No previous orbit memory found.
2026-07-22 14:44:07 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:44:08 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS8.nc
2026-07-22 14:44:08 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:44:08 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:44:08 INFO No previous orbit memory found.
2026-07-22 14:44:09 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:44:09 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS8.nc
2026-07-22 14:44:09 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:44:09 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 14:44:10 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens09.nc
2026-07-22 14:44:10 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:44:10 INFO No previous orbit memory found.
2026-07-22 14:44:11 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:44:11 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS9.nc
2026-07-22 14:44:11 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:44:11 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:44:11 INFO No previous orbit memory found.
2026-07-22 14:44:12 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:44:13 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS9.nc
2026-07-22 14:44:13 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:44:13 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 14:44:13 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens10.nc
2026-07-22 14:44:13 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:44:13 INFO No previous orbit memory found.
2026-07-22 14:44:14 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:44:15 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS10.nc
2026-07-22 14:44:15 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:44:15 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:44:15 INFO No previous orbit memory found.
2026-07-22 14:44:15 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:44:16 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS10.nc
2026-07-22 14:44:16 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:44:16 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 14:44:17 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens11.nc
2026-07-22 14:44:17 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:44:17 INFO No previous orbit memory found.
2026-07-22 14:44:18 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:44:18 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS11.nc
2026-07-22 14:44:18 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:44:19 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:44:19 INFO No previous orbit memory found.
2026-07-22 14:44:19 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:44:20 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS11.nc
2026-07-22 14:44:20 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:44:20 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 14:44:20 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens12.nc
2026-07-22 14:44:20 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:44:20 INFO No previous orbit memory found.
2026-07-22 14:44:21 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:44:22 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS12.nc
2026-07-22 14:44:22 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:44:22 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:44:22 INFO No previous orbit memory found.
2026-07-22 14:44:23 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:44:23 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS12.nc
2026-07-22 14:44:23 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:44:23 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 14:44:24 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens13.nc
2026-07-22 14:44:24 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:44:24 INFO No previous orbit memory found.
2026-07-22 14:44:24 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:44:25 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS13.nc
2026-07-22 14:44:25 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:44:25 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:44:25 INFO No previous orbit memory found.
2026-07-22 14:44:26 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:44:26 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS13.nc
2026-07-22 14:44:26 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:44:27 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 14:44:27 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens14.nc
2026-07-22 14:44:27 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:44:27 INFO No previous orbit memory found.
2026-07-22 14:44:28 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:44:28 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS14.nc
2026-07-22 14:44:28 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:44:29 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:44:29 INFO No previous orbit memory found.
2026-07-22 14:44:29 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:44:30 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS14.nc
2026-07-22 14:44:30 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:44:30 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 14:44:30 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens15.nc
2026-07-22 14:44:30 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:44:30 INFO No previous orbit memory found.
2026-07-22 14:44:31 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:44:32 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS15.nc
2026-07-22 14:44:32 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:44:32 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:44:32 INFO No previous orbit memory found.
2026-07-22 14:44:33 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:44:33 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS15.nc
2026-07-22 14:44:33 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:44:33 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 14:44:33 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-22 14:44:33 INFO [TIME] step_end current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 14:44:33 INFO [TIME] step_start current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 14:44:33 INFO [TIME] window start=2020-02-06 09:00:00 end=2020-02-06 11:00:00 run_hours=2 has_assimilation=True
2026-07-22 14:44:33 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:44:35 INFO Hourly dataset computed and listing created
2026-07-22 14:44:42 INFO Hourly dataset computed
2026-07-22 14:44:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:44:43 INFO Hourly dataset computed and listing created
2026-07-22 14:44:43 INFO Hourly dataset computed
2026-07-22 14:44:43 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:44:44 INFO Hourly dataset computed and listing created
2026-07-22 14:44:45 INFO Hourly dataset computed
2026-07-22 14:44:45 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:44:46 INFO Hourly dataset computed and listing created
2026-07-22 14:44:47 INFO Hourly dataset computed
2026-07-22 14:44:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:44:48 INFO Hourly dataset computed and listing created
2026-07-22 14:44:49 INFO Hourly dataset computed
2026-07-22 14:44:49 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:44:50 INFO Hourly dataset computed and listing created
2026-07-22 14:44:51 INFO Hourly dataset computed
2026-07-22 14:44:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:44:51 INFO Hourly dataset computed and listing created
2026-07-22 14:44:52 INFO Hourly dataset computed
2026-07-22 14:44:52 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:44:53 INFO Hourly dataset computed and listing created
2026-07-22 14:44:54 INFO Hourly dataset computed
2026-07-22 14:44:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:44:55 INFO Hourly dataset computed and listing created
2026-07-22 14:44:56 INFO Hourly dataset computed
2026-07-22 14:44:56 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:44:57 INFO Hourly dataset computed and listing created
2026-07-22 14:44:58 INFO Hourly dataset computed
2026-07-22 14:44:58 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:44:59 INFO Hourly dataset computed and listing created
2026-07-22 14:45:00 INFO Hourly dataset computed
2026-07-22 14:45:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:45:01 INFO Hourly dataset computed and listing created
2026-07-22 14:45:02 INFO Hourly dataset computed
2026-07-22 14:45:02 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:45:03 INFO Hourly dataset computed and listing created
2026-07-22 14:45:03 INFO Hourly dataset computed
2026-07-22 14:45:03 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:45:05 INFO Hourly dataset computed and listing created
2026-07-22 14:45:05 INFO Hourly dataset computed
2026-07-22 14:45:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:45:06 INFO Hourly dataset computed and listing created
2026-07-22 14:45:07 INFO Hourly dataset computed
2026-07-22 14:45:07 INFO ---------->>> Running CHIMERE model from 2020-02-06 09:00:00 to 2020-02-06 11:00:00
2026-07-22 14:45:07 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:45:07 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1
2026-07-22 14:45:07 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020601_8_ENS1.nc
2026-07-22 14:45:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-22 14:45:07 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:45:07 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-22 14:45:07 INFO Queuing job for member 1...
2026-07-22 14:45:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:45:07 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-22 14:45:08 INFO Found: ['5252329']
2026-07-22 14:45:13 INFO [TGCC-IRENE] Submitted job with ID:['5252329']
2026-07-22 14:45:13 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:45:13 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2
2026-07-22 14:45:13 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020601_8_ENS2.nc
2026-07-22 14:45:13 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-22 14:45:13 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:45:13 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-22 14:45:13 INFO Queuing job for member 2...
2026-07-22 14:45:13 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:45:13 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-22 14:45:15 INFO Found: ['5252332']
2026-07-22 14:45:20 INFO [TGCC-IRENE] Submitted job with ID:['5252332']
2026-07-22 14:45:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:45:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3
2026-07-22 14:45:20 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020601_8_ENS3.nc
2026-07-22 14:45:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-22 14:45:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:45:20 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-22 14:45:20 INFO Queuing job for member 3...
2026-07-22 14:45:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:45:20 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-22 14:45:21 INFO Found: ['5252333']
2026-07-22 14:45:26 INFO [TGCC-IRENE] Submitted job with ID:['5252333']
2026-07-22 14:45:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:45:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4
2026-07-22 14:45:26 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020601_8_ENS4.nc
2026-07-22 14:45:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-22 14:45:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:45:26 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-22 14:45:26 INFO Queuing job for member 4...
2026-07-22 14:45:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:45:26 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-22 14:45:27 INFO Found: ['5252335']
2026-07-22 14:45:32 INFO [TGCC-IRENE] Submitted job with ID:['5252335']
2026-07-22 14:45:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:45:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5
2026-07-22 14:45:32 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020601_8_ENS5.nc
2026-07-22 14:45:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-22 14:45:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:45:32 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-22 14:45:32 INFO Queuing job for member 5...
2026-07-22 14:45:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:45:32 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-22 14:45:32 INFO Found: ['5252340']
2026-07-22 14:45:37 INFO [TGCC-IRENE] Submitted job with ID:['5252340']
2026-07-22 14:45:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:45:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6
2026-07-22 14:45:37 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020601_8_ENS6.nc
2026-07-22 14:45:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-22 14:45:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:45:37 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-22 14:45:37 INFO Queuing job for member 6...
2026-07-22 14:45:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:45:37 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-22 14:45:38 INFO Found: ['5252343']
2026-07-22 14:45:43 INFO [TGCC-IRENE] Submitted job with ID:['5252343']
2026-07-22 14:45:43 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:45:43 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7
2026-07-22 14:45:43 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020601_8_ENS7.nc
2026-07-22 14:45:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-22 14:45:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:45:43 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-22 14:45:43 INFO Queuing job for member 7...
2026-07-22 14:45:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:45:43 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-22 14:45:44 INFO Found: ['5252345']
2026-07-22 14:45:49 INFO [TGCC-IRENE] Submitted job with ID:['5252345']
2026-07-22 14:45:49 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:45:49 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8
2026-07-22 14:45:49 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020601_8_ENS8.nc
2026-07-22 14:45:49 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-22 14:45:49 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:45:49 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-22 14:45:49 INFO Queuing job for member 8...
2026-07-22 14:45:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:45:49 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-22 14:45:50 INFO Found: ['5252348']
2026-07-22 14:45:55 INFO [TGCC-IRENE] Submitted job with ID:['5252348']
2026-07-22 14:45:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:45:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9
2026-07-22 14:45:55 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020601_8_ENS9.nc
2026-07-22 14:45:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-22 14:45:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:45:55 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-22 14:45:55 INFO Queuing job for member 9...
2026-07-22 14:45:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:45:55 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-22 14:45:55 INFO Found: ['5252350']
2026-07-22 14:46:00 INFO [TGCC-IRENE] Submitted job with ID:['5252350']
2026-07-22 14:46:00 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:46:00 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10
2026-07-22 14:46:00 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020601_8_ENS10.nc
2026-07-22 14:46:00 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-22 14:46:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:46:01 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-22 14:46:01 INFO Queuing job for member 10...
2026-07-22 14:46:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:46:01 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-22 14:46:02 INFO Found: ['5252352']
2026-07-22 14:46:07 INFO [TGCC-IRENE] Submitted job with ID:['5252352']
2026-07-22 14:46:07 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:46:07 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11
2026-07-22 14:46:07 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020601_8_ENS11.nc
2026-07-22 14:46:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-22 14:46:07 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:46:07 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-22 14:46:07 INFO Queuing job for member 11...
2026-07-22 14:46:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:46:07 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-22 14:46:09 INFO Found: ['5252354']
2026-07-22 14:46:14 INFO [TGCC-IRENE] Submitted job with ID:['5252354']
2026-07-22 14:46:14 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:46:14 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12
2026-07-22 14:46:14 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020601_8_ENS12.nc
2026-07-22 14:46:14 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-22 14:46:14 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:46:14 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-22 14:46:14 INFO Queuing job for member 12...
2026-07-22 14:46:14 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:46:14 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-22 14:46:17 INFO Found: ['5252356']
2026-07-22 14:46:22 INFO [TGCC-IRENE] Submitted job with ID:['5252356']
2026-07-22 14:46:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:46:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13
2026-07-22 14:46:22 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020601_8_ENS13.nc
2026-07-22 14:46:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-22 14:46:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:46:22 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-22 14:46:22 INFO Queuing job for member 13...
2026-07-22 14:46:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:46:22 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-22 14:46:24 INFO Found: ['5252359']
2026-07-22 14:46:29 INFO [TGCC-IRENE] Submitted job with ID:['5252359']
2026-07-22 14:46:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:46:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14
2026-07-22 14:46:29 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020601_8_ENS14.nc
2026-07-22 14:46:30 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-22 14:46:30 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:46:30 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-22 14:46:30 INFO Queuing job for member 14...
2026-07-22 14:46:30 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:46:30 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-22 14:46:32 INFO Found: ['5252362']
2026-07-22 14:46:37 INFO [TGCC-IRENE] Submitted job with ID:['5252362']
2026-07-22 14:46:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:46:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15
2026-07-22 14:46:37 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020601_8_ENS15.nc
2026-07-22 14:46:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-22 14:46:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:46:37 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-22 14:46:37 INFO Queuing job for member 15...
2026-07-22 14:46:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:46:37 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-22 14:46:39 INFO Found: ['5252363']
2026-07-22 14:46:44 INFO [TGCC-IRENE] Submitted job with ID:['5252363']
2026-07-22 14:46:44 INFO Checking job status ...
2026-07-22 14:46:44 INFO None 5252329: status RUNNING/PENDING
2026-07-22 14:46:44 INFO None 5252332: status RUNNING/PENDING
2026-07-22 14:46:44 INFO None 5252333: status RUNNING/PENDING
2026-07-22 14:46:44 INFO None 5252335: status RUNNING/PENDING
2026-07-22 14:46:44 INFO None 5252340: status RUNNING/PENDING
2026-07-22 14:46:45 INFO None 5252343: status RUNNING/PENDING
2026-07-22 14:46:45 INFO None 5252345: status RUNNING/PENDING
2026-07-22 14:46:45 INFO None 5252348: status RUNNING/PENDING
2026-07-22 14:46:45 INFO None 5252350: status RUNNING/PENDING
2026-07-22 14:46:45 INFO None 5252352: status RUNNING/PENDING
2026-07-22 14:46:45 INFO None 5252354: status RUNNING/PENDING
2026-07-22 14:46:45 INFO None 5252356: status RUNNING/PENDING
2026-07-22 14:46:45 INFO None 5252359: status RUNNING/PENDING
2026-07-22 14:46:45 INFO None 5252362: status RUNNING/PENDING
2026-07-22 14:46:45 INFO None 5252363: status RUNNING/PENDING
2026-07-22 14:46:45 INFO Jobs still running: ['5252329', '5252332', '5252333', '5252335', '5252340', '5252343', '5252345', '5252348', '5252350', '5252352', '5252354', '5252356', '5252359', '5252362', '5252363']. Waiting...
2026-07-22 14:47:00 INFO None 5252329: status RUNNING/PENDING
2026-07-22 14:47:00 INFO None 5252332: status RUNNING/PENDING
2026-07-22 14:47:00 INFO None 5252333: status RUNNING/PENDING
2026-07-22 14:47:00 INFO None 5252335: status RUNNING/PENDING
2026-07-22 14:47:00 INFO None 5252340: status RUNNING/PENDING
2026-07-22 14:47:00 INFO None 5252343: status RUNNING/PENDING
2026-07-22 14:47:00 INFO None 5252345: status RUNNING/PENDING
2026-07-22 14:47:00 INFO None 5252348: status RUNNING/PENDING
2026-07-22 14:47:00 INFO None 5252350: status RUNNING/PENDING
2026-07-22 14:47:00 INFO None 5252352: status RUNNING/PENDING
2026-07-22 14:47:00 INFO None 5252354: status RUNNING/PENDING
2026-07-22 14:47:00 INFO None 5252356: status RUNNING/PENDING
2026-07-22 14:47:00 INFO None 5252359: status RUNNING/PENDING
2026-07-22 14:47:00 INFO None 5252362: status RUNNING/PENDING
2026-07-22 14:47:00 INFO None 5252363: status RUNNING/PENDING
2026-07-22 14:47:00 INFO Jobs still running: ['5252329', '5252332', '5252333', '5252335', '5252340', '5252343', '5252345', '5252348', '5252350', '5252352', '5252354', '5252356', '5252359', '5252362', '5252363']. Waiting...
2026-07-22 14:47:15 INFO None 5252329: status RUNNING/PENDING
2026-07-22 14:47:15 INFO None 5252332: status RUNNING/PENDING
2026-07-22 14:47:15 INFO None 5252333: status RUNNING/PENDING
2026-07-22 14:47:15 INFO None 5252335: status RUNNING/PENDING
2026-07-22 14:47:15 INFO None 5252340: status RUNNING/PENDING
2026-07-22 14:47:15 INFO None 5252343: status RUNNING/PENDING
2026-07-22 14:47:15 INFO None 5252345: status RUNNING/PENDING
2026-07-22 14:47:15 INFO None 5252348: status RUNNING/PENDING
2026-07-22 14:47:15 INFO None 5252350: status RUNNING/PENDING
2026-07-22 14:47:15 INFO None 5252352: status RUNNING/PENDING
2026-07-22 14:47:15 INFO None 5252354: status RUNNING/PENDING
2026-07-22 14:47:15 INFO None 5252356: status RUNNING/PENDING
2026-07-22 14:47:15 INFO None 5252359: status RUNNING/PENDING
2026-07-22 14:47:15 INFO None 5252362: status RUNNING/PENDING
2026-07-22 14:47:15 INFO None 5252363: status RUNNING/PENDING
2026-07-22 14:47:15 INFO Jobs still running: ['5252329', '5252332', '5252333', '5252335', '5252340', '5252343', '5252345', '5252348', '5252350', '5252352', '5252354', '5252356', '5252359', '5252362', '5252363']. Waiting...
2026-07-22 14:47:31 INFO None 5252329: status RUNNING/PENDING
2026-07-22 14:47:32 INFO None 5252332: status RUNNING/PENDING
2026-07-22 14:47:32 INFO None 5252333: status RUNNING/PENDING
2026-07-22 14:47:32 INFO None 5252335: status RUNNING/PENDING
2026-07-22 14:47:32 INFO None 5252340: status RUNNING/PENDING
2026-07-22 14:47:32 INFO None 5252343: status RUNNING/PENDING
2026-07-22 14:47:34 INFO None 5252345: status RUNNING/PENDING
2026-07-22 14:47:34 INFO None 5252348: status RUNNING/PENDING
2026-07-22 14:47:34 INFO None 5252350: status RUNNING/PENDING
2026-07-22 14:47:34 INFO None 5252352: status RUNNING/PENDING
2026-07-22 14:47:34 INFO None 5252354: status RUNNING/PENDING
2026-07-22 14:47:34 INFO None 5252356: status RUNNING/PENDING
2026-07-22 14:47:34 INFO None 5252359: status RUNNING/PENDING
2026-07-22 14:47:34 INFO None 5252362: status RUNNING/PENDING
2026-07-22 14:47:34 INFO None 5252363: status RUNNING/PENDING
2026-07-22 14:47:34 INFO Jobs still running: ['5252329', '5252332', '5252333', '5252335', '5252340', '5252343', '5252345', '5252348', '5252350', '5252352', '5252354', '5252356', '5252359', '5252362', '5252363']. Waiting...
2026-07-22 14:47:49 INFO None 5252329: status RUNNING/PENDING
2026-07-22 14:47:49 INFO None 5252332: status RUNNING/PENDING
2026-07-22 14:47:49 INFO None 5252333: status RUNNING/PENDING
2026-07-22 14:47:49 INFO None 5252335: status RUNNING/PENDING
2026-07-22 14:47:49 INFO None 5252340: status RUNNING/PENDING
2026-07-22 14:47:49 INFO None 5252343: status RUNNING/PENDING
2026-07-22 14:47:49 INFO None 5252345: status RUNNING/PENDING
2026-07-22 14:47:49 INFO None 5252348: status RUNNING/PENDING
2026-07-22 14:47:49 INFO None 5252350: status RUNNING/PENDING
2026-07-22 14:47:49 INFO None 5252352: status RUNNING/PENDING
2026-07-22 14:47:49 INFO None 5252354: status RUNNING/PENDING
2026-07-22 14:47:49 INFO None 5252356: status RUNNING/PENDING
2026-07-22 14:47:49 INFO None 5252359: status RUNNING/PENDING
2026-07-22 14:47:49 INFO None 5252362: status RUNNING/PENDING
2026-07-22 14:47:51 INFO None 5252363: status RUNNING/PENDING
2026-07-22 14:47:51 INFO Jobs still running: ['5252329', '5252332', '5252333', '5252335', '5252340', '5252343', '5252345', '5252348', '5252350', '5252352', '5252354', '5252356', '5252359', '5252362', '5252363']. Waiting...
2026-07-22 14:48:06 INFO None 5252329: status RUNNING/PENDING
2026-07-22 14:48:06 INFO None 5252332: status RUNNING/PENDING
2026-07-22 14:48:06 INFO None 5252333: status RUNNING/PENDING
2026-07-22 14:48:07 INFO None 5252335: status RUNNING/PENDING
2026-07-22 14:48:07 INFO None 5252340: status RUNNING/PENDING
2026-07-22 14:48:07 INFO None 5252343: status RUNNING/PENDING
2026-07-22 14:48:07 INFO None 5252345: status RUNNING/PENDING
2026-07-22 14:48:07 INFO None 5252348: status RUNNING/PENDING
2026-07-22 14:48:07 INFO None 5252350: status RUNNING/PENDING
2026-07-22 14:48:07 INFO None 5252352: status RUNNING/PENDING
2026-07-22 14:48:07 INFO None 5252354: status RUNNING/PENDING
2026-07-22 14:48:07 INFO None 5252356: status RUNNING/PENDING
2026-07-22 14:48:07 INFO None 5252359: status RUNNING/PENDING
2026-07-22 14:48:07 INFO None 5252362: status RUNNING/PENDING
2026-07-22 14:48:07 INFO None 5252363: status RUNNING/PENDING
2026-07-22 14:48:07 INFO Jobs still running: ['5252329', '5252332', '5252333', '5252335', '5252340', '5252343', '5252345', '5252348', '5252350', '5252352', '5252354', '5252356', '5252359', '5252362', '5252363']. Waiting...
2026-07-22 14:48:22 INFO None 5252329: status RUNNING/PENDING
2026-07-22 14:48:22 INFO None 5252332: status RUNNING/PENDING
2026-07-22 14:48:22 INFO None 5252333: status RUNNING/PENDING
2026-07-22 14:48:22 INFO None 5252335: status RUNNING/PENDING
2026-07-22 14:48:22 INFO None 5252340: status RUNNING/PENDING
2026-07-22 14:48:22 INFO None 5252343: status RUNNING/PENDING
2026-07-22 14:48:22 INFO None 5252345: status RUNNING/PENDING
2026-07-22 14:48:22 INFO None 5252348: status RUNNING/PENDING
2026-07-22 14:48:22 INFO None 5252350: status RUNNING/PENDING
2026-07-22 14:48:22 INFO None 5252352: status RUNNING/PENDING
2026-07-22 14:48:22 INFO None 5252354: status RUNNING/PENDING
2026-07-22 14:48:22 INFO None 5252356: status RUNNING/PENDING
2026-07-22 14:48:22 INFO None 5252359: status RUNNING/PENDING
2026-07-22 14:48:22 INFO None 5252362: status RUNNING/PENDING
2026-07-22 14:48:22 INFO None 5252363: status RUNNING/PENDING
2026-07-22 14:48:22 INFO Jobs still running: ['5252329', '5252332', '5252333', '5252335', '5252340', '5252343', '5252345', '5252348', '5252350', '5252352', '5252354', '5252356', '5252359', '5252362', '5252363']. Waiting...
2026-07-22 14:48:37 INFO None 5252329: status RUNNING/PENDING
2026-07-22 14:48:37 INFO None 5252332: status RUNNING/PENDING
2026-07-22 14:48:37 INFO None 5252333: status RUNNING/PENDING
2026-07-22 14:48:37 INFO None 5252335: status RUNNING/PENDING
2026-07-22 14:48:37 INFO None 5252340: status RUNNING/PENDING
2026-07-22 14:48:37 INFO None 5252343: status RUNNING/PENDING
2026-07-22 14:48:37 INFO None 5252345: status RUNNING/PENDING
2026-07-22 14:48:37 INFO None 5252348: status RUNNING/PENDING
2026-07-22 14:48:37 INFO None 5252350: status RUNNING/PENDING
2026-07-22 14:48:37 INFO None 5252352: status RUNNING/PENDING
2026-07-22 14:48:37 INFO None 5252354: status RUNNING/PENDING
2026-07-22 14:48:37 INFO None 5252356: status RUNNING/PENDING
2026-07-22 14:48:37 INFO None 5252359: status RUNNING/PENDING
2026-07-22 14:48:38 INFO None 5252362: status RUNNING/PENDING
2026-07-22 14:48:38 INFO None 5252363: status RUNNING/PENDING
2026-07-22 14:48:38 INFO Jobs still running: ['5252329', '5252332', '5252333', '5252335', '5252340', '5252343', '5252345', '5252348', '5252350', '5252352', '5252354', '5252356', '5252359', '5252362', '5252363']. Waiting...
2026-07-22 14:48:54 INFO None 5252329: status RUNNING/PENDING
2026-07-22 14:48:54 INFO None 5252332: status RUNNING/PENDING
2026-07-22 14:48:54 INFO None 5252333: status RUNNING/PENDING
2026-07-22 14:48:54 INFO None 5252335: status RUNNING/PENDING
2026-07-22 14:48:54 INFO None 5252340: status RUNNING/PENDING
2026-07-22 14:48:54 INFO None 5252343: status RUNNING/PENDING
2026-07-22 14:48:54 INFO None 5252345: status RUNNING/PENDING
2026-07-22 14:48:54 INFO None 5252348: status RUNNING/PENDING
2026-07-22 14:48:54 INFO None 5252350: status RUNNING/PENDING
2026-07-22 14:48:54 INFO None 5252352: status RUNNING/PENDING
2026-07-22 14:48:54 INFO None 5252354: status RUNNING/PENDING
2026-07-22 14:48:54 INFO None 5252356: status RUNNING/PENDING
2026-07-22 14:48:54 INFO None 5252359: status RUNNING/PENDING
2026-07-22 14:48:54 INFO None 5252362: status RUNNING/PENDING
2026-07-22 14:48:54 INFO None 5252363: status RUNNING/PENDING
2026-07-22 14:48:54 INFO Jobs still running: ['5252329', '5252332', '5252333', '5252335', '5252340', '5252343', '5252345', '5252348', '5252350', '5252352', '5252354', '5252356', '5252359', '5252362', '5252363']. Waiting...
2026-07-22 14:49:09 INFO None 5252329: status RUNNING/PENDING
2026-07-22 14:49:09 INFO None 5252332: status RUNNING/PENDING
2026-07-22 14:49:09 INFO None 5252333: status RUNNING/PENDING
2026-07-22 14:49:09 INFO None 5252335: status RUNNING/PENDING
2026-07-22 14:49:09 INFO None 5252340: status RUNNING/PENDING
2026-07-22 14:49:09 INFO None 5252343: status RUNNING/PENDING
2026-07-22 14:49:09 INFO None 5252345: status RUNNING/PENDING
2026-07-22 14:49:09 INFO None 5252348: status RUNNING/PENDING
2026-07-22 14:49:10 INFO None 5252350: status FINISHED
2026-07-22 14:49:12 INFO None 5252352: status RUNNING/PENDING
2026-07-22 14:49:12 INFO None 5252354: status RUNNING/PENDING
2026-07-22 14:49:12 INFO None 5252356: status RUNNING/PENDING
2026-07-22 14:49:12 INFO None 5252359: status RUNNING/PENDING
2026-07-22 14:49:12 INFO None 5252362: status RUNNING/PENDING
2026-07-22 14:49:12 INFO None 5252363: status RUNNING/PENDING
2026-07-22 14:49:12 INFO Jobs still running: ['5252329', '5252332', '5252333', '5252335', '5252340', '5252343', '5252345', '5252348', '5252352', '5252354', '5252356', '5252359', '5252362', '5252363']. Waiting...
2026-07-22 14:49:27 INFO None 5252329: status RUNNING/PENDING
2026-07-22 14:49:27 INFO None 5252332: status RUNNING/PENDING
2026-07-22 14:49:27 INFO None 5252333: status RUNNING/PENDING
2026-07-22 14:49:27 INFO None 5252335: status RUNNING/PENDING
2026-07-22 14:49:27 INFO None 5252340: status RUNNING/PENDING
2026-07-22 14:49:27 INFO None 5252343: status RUNNING/PENDING
2026-07-22 14:49:27 INFO None 5252345: status RUNNING/PENDING
2026-07-22 14:49:27 INFO None 5252348: status RUNNING/PENDING
2026-07-22 14:49:27 INFO None 5252350: status FINISHED
2026-07-22 14:49:27 INFO None 5252352: status RUNNING/PENDING
2026-07-22 14:49:27 INFO None 5252354: status RUNNING/PENDING
2026-07-22 14:49:27 INFO None 5252356: status RUNNING/PENDING
2026-07-22 14:49:27 INFO None 5252359: status RUNNING/PENDING
2026-07-22 14:49:27 INFO None 5252362: status RUNNING/PENDING
2026-07-22 14:49:27 INFO None 5252363: status RUNNING/PENDING
2026-07-22 14:49:27 INFO Jobs still running: ['5252329', '5252332', '5252333', '5252335', '5252340', '5252343', '5252345', '5252348', '5252352', '5252354', '5252356', '5252359', '5252362', '5252363']. Waiting...
2026-07-22 14:49:42 INFO None 5252329: status RUNNING/PENDING
2026-07-22 14:49:42 INFO None 5252332: status RUNNING/PENDING
2026-07-22 14:49:42 INFO None 5252333: status RUNNING/PENDING
2026-07-22 14:49:42 INFO None 5252335: status RUNNING/PENDING
2026-07-22 14:49:42 INFO None 5252340: status RUNNING/PENDING
2026-07-22 14:49:42 INFO None 5252343: status RUNNING/PENDING
2026-07-22 14:49:42 INFO None 5252345: status RUNNING/PENDING
2026-07-22 14:49:42 INFO None 5252348: status RUNNING/PENDING
2026-07-22 14:49:42 INFO None 5252350: status FINISHED
2026-07-22 14:49:42 INFO None 5252352: status RUNNING/PENDING
2026-07-22 14:49:42 INFO None 5252354: status RUNNING/PENDING
2026-07-22 14:49:42 INFO None 5252356: status RUNNING/PENDING
2026-07-22 14:49:42 INFO None 5252359: status RUNNING/PENDING
2026-07-22 14:49:42 INFO None 5252362: status RUNNING/PENDING
2026-07-22 14:49:42 INFO None 5252363: status RUNNING/PENDING
2026-07-22 14:49:42 INFO Jobs still running: ['5252329', '5252332', '5252333', '5252335', '5252340', '5252343', '5252345', '5252348', '5252352', '5252354', '5252356', '5252359', '5252362', '5252363']. Waiting...
2026-07-22 14:49:57 INFO None 5252329: status RUNNING/PENDING
2026-07-22 14:49:58 INFO None 5252332: status RUNNING/PENDING
2026-07-22 14:49:58 INFO None 5252333: status RUNNING/PENDING
2026-07-22 14:49:58 INFO None 5252335: status RUNNING/PENDING
2026-07-22 14:49:58 INFO None 5252340: status RUNNING/PENDING
2026-07-22 14:49:58 INFO None 5252343: status RUNNING/PENDING
2026-07-22 14:49:58 INFO None 5252345: status RUNNING/PENDING
2026-07-22 14:49:58 INFO None 5252348: status RUNNING/PENDING
2026-07-22 14:49:58 INFO None 5252350: status FINISHED
2026-07-22 14:49:58 INFO None 5252352: status RUNNING/PENDING
2026-07-22 14:49:58 INFO None 5252354: status RUNNING/PENDING
2026-07-22 14:49:58 INFO None 5252356: status RUNNING/PENDING
2026-07-22 14:49:58 INFO None 5252359: status RUNNING/PENDING
2026-07-22 14:49:58 INFO None 5252362: status RUNNING/PENDING
2026-07-22 14:49:58 INFO None 5252363: status RUNNING/PENDING
2026-07-22 14:49:58 INFO Jobs still running: ['5252329', '5252332', '5252333', '5252335', '5252340', '5252343', '5252345', '5252348', '5252352', '5252354', '5252356', '5252359', '5252362', '5252363']. Waiting...
2026-07-22 14:50:15 INFO None 5252329: status RUNNING/PENDING
2026-07-22 14:50:15 INFO None 5252332: status RUNNING/PENDING
2026-07-22 14:50:15 INFO None 5252333: status RUNNING/PENDING
2026-07-22 14:50:15 INFO None 5252335: status RUNNING/PENDING
2026-07-22 14:50:15 INFO None 5252340: status FINISHED
2026-07-22 14:50:15 INFO None 5252343: status FINISHED
2026-07-22 14:50:15 INFO None 5252345: status FINISHED
2026-07-22 14:50:17 INFO None 5252348: status FINISHED
2026-07-22 14:50:17 INFO None 5252350: status FINISHED
2026-07-22 14:50:17 INFO None 5252352: status RUNNING/PENDING
2026-07-22 14:50:17 INFO None 5252354: status RUNNING/PENDING
2026-07-22 14:50:17 INFO None 5252356: status RUNNING/PENDING
2026-07-22 14:50:17 INFO None 5252359: status RUNNING/PENDING
2026-07-22 14:50:18 INFO None 5252362: status RUNNING/PENDING
2026-07-22 14:50:18 INFO None 5252363: status RUNNING/PENDING
2026-07-22 14:50:18 INFO Jobs still running: ['5252329', '5252332', '5252333', '5252335', '5252352', '5252354', '5252356', '5252359', '5252362', '5252363']. Waiting...
2026-07-22 14:50:33 INFO None 5252329: status RUNNING/PENDING
2026-07-22 14:50:33 INFO None 5252332: status FINISHED
2026-07-22 14:50:33 INFO None 5252333: status FINISHED
2026-07-22 14:50:33 INFO None 5252335: status RUNNING/PENDING
2026-07-22 14:50:33 INFO None 5252340: status FINISHED
2026-07-22 14:50:33 INFO None 5252343: status FINISHED
2026-07-22 14:50:33 INFO None 5252345: status FINISHED
2026-07-22 14:50:33 INFO None 5252348: status FINISHED
2026-07-22 14:50:33 INFO None 5252350: status FINISHED
2026-07-22 14:50:33 INFO None 5252352: status RUNNING/PENDING
2026-07-22 14:50:33 INFO None 5252354: status RUNNING/PENDING
2026-07-22 14:50:33 INFO None 5252356: status RUNNING/PENDING
2026-07-22 14:50:33 INFO None 5252359: status RUNNING/PENDING
2026-07-22 14:50:33 INFO None 5252362: status RUNNING/PENDING
2026-07-22 14:50:33 INFO None 5252363: status RUNNING/PENDING
2026-07-22 14:50:33 INFO Jobs still running: ['5252329', '5252335', '5252352', '5252354', '5252356', '5252359', '5252362', '5252363']. Waiting...
2026-07-22 14:50:48 INFO None 5252329: status FINISHED
2026-07-22 14:50:48 INFO None 5252332: status FINISHED
2026-07-22 14:50:48 INFO None 5252333: status FINISHED
2026-07-22 14:50:48 INFO None 5252335: status FINISHED
2026-07-22 14:50:49 INFO None 5252340: status FINISHED
2026-07-22 14:50:49 INFO None 5252343: status FINISHED
2026-07-22 14:50:49 INFO None 5252345: status FINISHED
2026-07-22 14:50:49 INFO None 5252348: status FINISHED
2026-07-22 14:50:49 INFO None 5252350: status FINISHED
2026-07-22 14:50:49 INFO None 5252352: status FINISHED
2026-07-22 14:50:49 INFO None 5252354: status FINISHED
2026-07-22 14:50:49 INFO None 5252356: status RUNNING/PENDING
2026-07-22 14:50:49 INFO None 5252359: status RUNNING/PENDING
2026-07-22 14:50:49 INFO None 5252362: status RUNNING/PENDING
2026-07-22 14:50:49 INFO None 5252363: status RUNNING/PENDING
2026-07-22 14:50:49 INFO Jobs still running: ['5252356', '5252359', '5252362', '5252363']. Waiting...
2026-07-22 14:51:04 INFO None 5252329: status FINISHED
2026-07-22 14:51:04 INFO None 5252332: status FINISHED
2026-07-22 14:51:04 INFO None 5252333: status FINISHED
2026-07-22 14:51:04 INFO None 5252335: status FINISHED
2026-07-22 14:51:04 INFO None 5252340: status FINISHED
2026-07-22 14:51:04 INFO None 5252343: status FINISHED
2026-07-22 14:51:04 INFO None 5252345: status FINISHED
2026-07-22 14:51:04 INFO None 5252348: status FINISHED
2026-07-22 14:51:04 INFO None 5252350: status FINISHED
2026-07-22 14:51:04 INFO None 5252352: status FINISHED
2026-07-22 14:51:04 INFO None 5252354: status FINISHED
2026-07-22 14:51:04 INFO None 5252356: status RUNNING/PENDING
2026-07-22 14:51:04 INFO None 5252359: status RUNNING/PENDING
2026-07-22 14:51:05 INFO None 5252362: status RUNNING/PENDING
2026-07-22 14:51:05 INFO None 5252363: status RUNNING/PENDING
2026-07-22 14:51:05 INFO Jobs still running: ['5252356', '5252359', '5252362', '5252363']. Waiting...
2026-07-22 14:51:20 INFO None 5252329: status FINISHED
2026-07-22 14:51:20 INFO None 5252332: status FINISHED
2026-07-22 14:51:20 INFO None 5252333: status FINISHED
2026-07-22 14:51:20 INFO None 5252335: status FINISHED
2026-07-22 14:51:20 INFO None 5252340: status FINISHED
2026-07-22 14:51:20 INFO None 5252343: status FINISHED
2026-07-22 14:51:20 INFO None 5252345: status FINISHED
2026-07-22 14:51:20 INFO None 5252348: status FINISHED
2026-07-22 14:51:20 INFO None 5252350: status FINISHED
2026-07-22 14:51:20 INFO None 5252352: status FINISHED
2026-07-22 14:51:20 INFO None 5252354: status FINISHED
2026-07-22 14:51:20 INFO None 5252356: status RUNNING/PENDING
2026-07-22 14:51:20 INFO None 5252359: status RUNNING/PENDING
2026-07-22 14:51:20 INFO None 5252362: status RUNNING/PENDING
2026-07-22 14:51:20 INFO None 5252363: status RUNNING/PENDING
2026-07-22 14:51:20 INFO Jobs still running: ['5252356', '5252359', '5252362', '5252363']. Waiting...
2026-07-22 14:51:36 INFO None 5252329: status FINISHED
2026-07-22 14:51:36 INFO None 5252332: status FINISHED
2026-07-22 14:51:36 INFO None 5252333: status FINISHED
2026-07-22 14:51:36 INFO None 5252335: status FINISHED
2026-07-22 14:51:36 INFO None 5252340: status FINISHED
2026-07-22 14:51:36 INFO None 5252343: status FINISHED
2026-07-22 14:51:36 INFO None 5252345: status FINISHED
2026-07-22 14:51:36 INFO None 5252348: status FINISHED
2026-07-22 14:51:36 INFO None 5252350: status FINISHED
2026-07-22 14:51:36 INFO None 5252352: status FINISHED
2026-07-22 14:51:36 INFO None 5252354: status FINISHED
2026-07-22 14:51:36 INFO None 5252356: status RUNNING/PENDING
2026-07-22 14:51:36 INFO None 5252359: status RUNNING/PENDING
2026-07-22 14:51:36 INFO None 5252362: status RUNNING/PENDING
2026-07-22 14:51:36 INFO None 5252363: status RUNNING/PENDING
2026-07-22 14:51:36 INFO Jobs still running: ['5252356', '5252359', '5252362', '5252363']. Waiting...
2026-07-22 14:51:51 INFO None 5252329: status FINISHED
2026-07-22 14:51:51 INFO None 5252332: status FINISHED
2026-07-22 14:51:51 INFO None 5252333: status FINISHED
2026-07-22 14:51:51 INFO None 5252335: status FINISHED
2026-07-22 14:51:51 INFO None 5252340: status FINISHED
2026-07-22 14:51:51 INFO None 5252343: status FINISHED
2026-07-22 14:51:51 INFO None 5252345: status FINISHED
2026-07-22 14:51:51 INFO None 5252348: status FINISHED
2026-07-22 14:51:51 INFO None 5252350: status FINISHED
2026-07-22 14:51:51 INFO None 5252352: status FINISHED
2026-07-22 14:51:53 INFO None 5252354: status FINISHED
2026-07-22 14:51:53 INFO None 5252356: status FINISHED
2026-07-22 14:51:54 INFO None 5252359: status RUNNING/PENDING
2026-07-22 14:51:54 INFO None 5252362: status RUNNING/PENDING
2026-07-22 14:51:54 INFO None 5252363: status RUNNING/PENDING
2026-07-22 14:51:54 INFO Jobs still running: ['5252359', '5252362', '5252363']. Waiting...
2026-07-22 14:52:09 INFO None 5252329: status FINISHED
2026-07-22 14:52:09 INFO None 5252332: status FINISHED
2026-07-22 14:52:09 INFO None 5252333: status FINISHED
2026-07-22 14:52:09 INFO None 5252335: status FINISHED
2026-07-22 14:52:09 INFO None 5252340: status FINISHED
2026-07-22 14:52:09 INFO None 5252343: status FINISHED
2026-07-22 14:52:09 INFO None 5252345: status FINISHED
2026-07-22 14:52:09 INFO None 5252348: status FINISHED
2026-07-22 14:52:09 INFO None 5252350: status FINISHED
2026-07-22 14:52:09 INFO None 5252352: status FINISHED
2026-07-22 14:52:09 INFO None 5252354: status FINISHED
2026-07-22 14:52:09 INFO None 5252356: status FINISHED
2026-07-22 14:52:09 INFO None 5252359: status RUNNING/PENDING
2026-07-22 14:52:09 INFO None 5252362: status FINISHED
2026-07-22 14:52:11 INFO None 5252363: status FINISHED
2026-07-22 14:52:11 INFO Jobs still running: ['5252359']. Waiting...
2026-07-22 14:52:26 INFO None 5252329: status FINISHED
2026-07-22 14:52:26 INFO None 5252332: status FINISHED
2026-07-22 14:52:26 INFO None 5252333: status FINISHED
2026-07-22 14:52:26 INFO None 5252335: status FINISHED
2026-07-22 14:52:26 INFO None 5252340: status FINISHED
2026-07-22 14:52:26 INFO None 5252343: status FINISHED
2026-07-22 14:52:26 INFO None 5252345: status FINISHED
2026-07-22 14:52:26 INFO None 5252348: status FINISHED
2026-07-22 14:52:26 INFO None 5252350: status FINISHED
2026-07-22 14:52:26 INFO None 5252352: status FINISHED
2026-07-22 14:52:26 INFO None 5252354: status FINISHED
2026-07-22 14:52:26 INFO None 5252356: status FINISHED
2026-07-22 14:52:26 INFO None 5252359: status FINISHED
2026-07-22 14:52:26 INFO None 5252362: status FINISHED
2026-07-22 14:52:26 INFO None 5252363: status FINISHED
2026-07-22 14:52:26 INFO Jobs ['5252329', '5252332', '5252333', '5252335', '5252340', '5252343', '5252345', '5252348', '5252350', '5252352', '5252354', '5252356', '5252359', '5252362', '5252363'] have finished
2026-07-22 14:52:26 INFO Checking restart files were created ...
2026-07-22 14:52:26 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc(1002685915 bytes)
2026-07-22 14:52:26 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc(1002685915 bytes)
2026-07-22 14:52:26 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc(1002685915 bytes)
2026-07-22 14:52:26 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc(1002685915 bytes)
2026-07-22 14:52:26 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc(1002685915 bytes)
2026-07-22 14:52:26 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc(1002685915 bytes)
2026-07-22 14:52:26 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc(1002685915 bytes)
2026-07-22 14:52:26 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc(1002685915 bytes)
2026-07-22 14:52:26 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc(1002685915 bytes)
2026-07-22 14:52:26 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc(1002685915 bytes)
2026-07-22 14:52:26 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc(1002685915 bytes)
2026-07-22 14:52:26 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc(1002685915 bytes)
2026-07-22 14:52:26 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc(1002685915 bytes)
2026-07-22 14:52:26 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc(1002685915 bytes)
2026-07-22 14:52:26 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc(1002685915 bytes)
2026-07-22 14:52:26 INFO  Run_model() completed successfully.
2026-07-22 14:52:26 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 14:52:26 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 11:00:00 days=153072 seconds=39600
2026-07-22 14:52:26 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-22 14:52:26 INFO [TIME] increment current_time 2020-02-06 09:00:00 -> 2020-02-06 11:00:00
2026-07-22 14:52:26 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 14:52:26 INFO ---------->>> Running process_satellite_data()
2026-07-22 14:52:27 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12003.nc
2026-07-22 14:52:27 INFO ---------->>> Running run_obs_converter()
2026-07-22 14:52:27 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_39045_153072.out
2026-07-22 14:52:27 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_39045_153072.out
2026-07-22 14:52:27 INFO ---------->>> Running DART
2026-07-22 14:52:27 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-22 14:52:27 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/chim_ENS1_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/chim_ENS1_2020020611_1_out_toDART.nc
2026-07-22 14:52:27 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/chim_ENS2_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/chim_ENS2_2020020611_1_out_toDART.nc
2026-07-22 14:52:27 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/chim_ENS3_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/chim_ENS3_2020020611_1_out_toDART.nc
2026-07-22 14:52:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/chim_ENS4_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/chim_ENS4_2020020611_1_out_toDART.nc
2026-07-22 14:52:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/chim_ENS5_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/chim_ENS5_2020020611_1_out_toDART.nc
2026-07-22 14:52:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/chim_ENS6_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/chim_ENS6_2020020611_1_out_toDART.nc
2026-07-22 14:52:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/chim_ENS7_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/chim_ENS7_2020020611_1_out_toDART.nc
2026-07-22 14:52:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/chim_ENS8_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/chim_ENS8_2020020611_1_out_toDART.nc
2026-07-22 14:52:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/chim_ENS9_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/chim_ENS9_2020020611_1_out_toDART.nc
2026-07-22 14:52:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/chim_ENS10_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/chim_ENS10_2020020611_1_out_toDART.nc
2026-07-22 14:52:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/chim_ENS11_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/chim_ENS11_2020020611_1_out_toDART.nc
2026-07-22 14:52:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/chim_ENS12_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/chim_ENS12_2020020611_1_out_toDART.nc
2026-07-22 14:52:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/chim_ENS13_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/chim_ENS13_2020020611_1_out_toDART.nc
2026-07-22 14:52:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/chim_ENS14_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/chim_ENS14_2020020611_1_out_toDART.nc
2026-07-22 14:52:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/chim_ENS15_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/chim_ENS15_2020020611_1_out_toDART.nc
2026-07-22 14:52:32 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-22 14:52:32 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-22 14:52:32 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-22 14:52:32 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-22 14:52:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-22 14:52:32 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-22 14:52:45 INFO Found: []
2026-07-22 14:52:45 INFO No job id returned by command ./run_filter.bsh
2026-07-22 14:52:45 INFO No monitoring will be performed
2026-07-22 14:52:45 INFO Moving DART output files to analysis and preassim directories for date 2020020611 if present ...
2026-07-22 14:52:45 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 14:52:45 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 14:52:45 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 14:52:45 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 14:52:45 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 14:52:45 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 14:52:45 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 14:52:45 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 14:52:45 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 14:52:45 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 14:52:45 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 14:52:45 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 14:52:45 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 14:52:45 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 14:52:45 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 14:52:45 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 14:52:45 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 14:52:45 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 14:52:45 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 14:52:45 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 14:52:45 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 14:52:45 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 14:52:45 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 14:52:45 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 14:52:45 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 14:52:45 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 14:52:45 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 14:52:45 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 14:52:45 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 14:52:45 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 14:52:45 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 14:52:45 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 14:52:45 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 14:52:45 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 14:52:45 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-22 14:52:45 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-22 14:52:45 INFO run_dart() is DONE.
2026-07-22 14:52:45 INFO ---------->>> Running update_pollutant_in_end()
2026-07-22 14:52:46 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:52:46 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS1.nc
2026-07-22 14:52:46 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:52:47 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS1.nc
2026-07-22 14:52:47 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:52:47 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:52:47 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS1.nc
2026-07-22 14:52:48 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:52:49 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS1.nc
2026-07-22 14:52:49 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:52:49 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 14:52:49 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:52:49 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS2.nc
2026-07-22 14:52:50 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:52:51 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS2.nc
2026-07-22 14:52:51 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:52:51 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:52:51 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS2.nc
2026-07-22 14:52:52 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:52:52 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS2.nc
2026-07-22 14:52:52 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:52:52 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 14:52:53 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:52:53 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS3.nc
2026-07-22 14:52:53 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:52:54 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS3.nc
2026-07-22 14:52:54 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:52:54 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:52:54 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS3.nc
2026-07-22 14:52:55 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:52:55 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS3.nc
2026-07-22 14:52:55 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:52:55 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 14:52:56 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:52:56 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS4.nc
2026-07-22 14:52:56 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:52:57 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS4.nc
2026-07-22 14:52:57 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:52:58 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:52:58 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS4.nc
2026-07-22 14:52:58 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:52:59 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS4.nc
2026-07-22 14:52:59 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:52:59 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 14:52:59 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:52:59 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS5.nc
2026-07-22 14:53:00 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:53:00 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS5.nc
2026-07-22 14:53:00 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:53:01 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:53:01 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS5.nc
2026-07-22 14:53:01 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:53:02 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS5.nc
2026-07-22 14:53:02 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:53:02 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 14:53:02 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:53:02 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS6.nc
2026-07-22 14:53:03 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:53:03 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS6.nc
2026-07-22 14:53:03 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:53:04 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:53:04 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS6.nc
2026-07-22 14:53:04 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:53:05 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS6.nc
2026-07-22 14:53:05 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:53:05 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 14:53:05 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:53:05 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS7.nc
2026-07-22 14:53:06 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:53:07 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS7.nc
2026-07-22 14:53:07 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:53:07 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:53:07 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS7.nc
2026-07-22 14:53:08 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:53:08 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS7.nc
2026-07-22 14:53:08 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:53:08 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 14:53:08 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:53:08 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS8.nc
2026-07-22 14:53:09 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:53:10 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS8.nc
2026-07-22 14:53:10 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:53:10 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:53:10 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS8.nc
2026-07-22 14:53:11 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:53:12 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS8.nc
2026-07-22 14:53:12 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:53:12 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 14:53:12 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:53:12 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS9.nc
2026-07-22 14:53:13 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:53:13 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS9.nc
2026-07-22 14:53:13 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:53:14 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:53:14 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS9.nc
2026-07-22 14:53:14 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:53:15 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS9.nc
2026-07-22 14:53:15 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:53:15 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 14:53:15 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:53:15 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS10.nc
2026-07-22 14:53:16 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:53:17 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS10.nc
2026-07-22 14:53:17 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:53:17 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:53:17 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS10.nc
2026-07-22 14:53:17 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:53:18 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS10.nc
2026-07-22 14:53:18 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:53:18 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 14:53:18 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:53:18 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS11.nc
2026-07-22 14:53:19 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:53:20 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS11.nc
2026-07-22 14:53:20 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:53:20 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:53:20 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS11.nc
2026-07-22 14:53:21 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:53:22 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS11.nc
2026-07-22 14:53:22 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:53:22 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 14:53:22 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:53:22 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS12.nc
2026-07-22 14:53:23 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:53:23 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS12.nc
2026-07-22 14:53:23 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:53:24 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:53:24 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS12.nc
2026-07-22 14:53:24 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:53:25 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS12.nc
2026-07-22 14:53:25 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:53:25 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 14:53:25 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:53:25 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS13.nc
2026-07-22 14:53:26 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:53:27 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS13.nc
2026-07-22 14:53:27 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:53:27 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:53:27 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS13.nc
2026-07-22 14:53:28 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:53:28 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS13.nc
2026-07-22 14:53:28 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:53:28 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 14:53:29 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:53:29 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS14.nc
2026-07-22 14:53:29 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:53:30 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS14.nc
2026-07-22 14:53:30 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:53:30 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:53:30 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS14.nc
2026-07-22 14:53:31 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:53:32 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS14.nc
2026-07-22 14:53:32 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:53:32 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 14:53:32 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:53:32 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS15.nc
2026-07-22 14:53:33 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:53:34 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS15.nc
2026-07-22 14:53:34 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:53:34 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 14:53:34 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS15.nc
2026-07-22 14:53:34 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 14:53:35 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS15.nc
2026-07-22 14:53:35 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 14:53:35 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 14:53:35 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-22 14:53:35 INFO [TIME] step_end current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 14:53:35 INFO [TIME] step_start current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 14:53:35 INFO [TIME] window start=2020-02-06 11:00:00 end=2020-02-06 13:00:00 run_hours=2 has_assimilation=True
2026-07-22 14:53:35 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:53:36 INFO Hourly dataset computed and listing created
2026-07-22 14:53:43 INFO Hourly dataset computed
2026-07-22 14:53:43 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:53:44 INFO Hourly dataset computed and listing created
2026-07-22 14:53:45 INFO Hourly dataset computed
2026-07-22 14:53:45 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:53:46 INFO Hourly dataset computed and listing created
2026-07-22 14:53:47 INFO Hourly dataset computed
2026-07-22 14:53:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:53:48 INFO Hourly dataset computed and listing created
2026-07-22 14:53:49 INFO Hourly dataset computed
2026-07-22 14:53:49 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:53:50 INFO Hourly dataset computed and listing created
2026-07-22 14:53:50 INFO Hourly dataset computed
2026-07-22 14:53:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:53:51 INFO Hourly dataset computed and listing created
2026-07-22 14:53:52 INFO Hourly dataset computed
2026-07-22 14:53:52 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:53:53 INFO Hourly dataset computed and listing created
2026-07-22 14:53:54 INFO Hourly dataset computed
2026-07-22 14:53:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:53:55 INFO Hourly dataset computed and listing created
2026-07-22 14:53:56 INFO Hourly dataset computed
2026-07-22 14:53:56 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:53:57 INFO Hourly dataset computed and listing created
2026-07-22 14:53:57 INFO Hourly dataset computed
2026-07-22 14:53:58 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:53:58 INFO Hourly dataset computed and listing created
2026-07-22 14:53:59 INFO Hourly dataset computed
2026-07-22 14:53:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:54:00 INFO Hourly dataset computed and listing created
2026-07-22 14:54:01 INFO Hourly dataset computed
2026-07-22 14:54:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:54:02 INFO Hourly dataset computed and listing created
2026-07-22 14:54:03 INFO Hourly dataset computed
2026-07-22 14:54:03 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:54:04 INFO Hourly dataset computed and listing created
2026-07-22 14:54:04 INFO Hourly dataset computed
2026-07-22 14:54:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:54:05 INFO Hourly dataset computed and listing created
2026-07-22 14:54:06 INFO Hourly dataset computed
2026-07-22 14:54:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 14:54:07 INFO Hourly dataset computed and listing created
2026-07-22 14:54:08 INFO Hourly dataset computed
2026-07-22 14:54:08 INFO ---------->>> Running CHIMERE model from 2020-02-06 11:00:00 to 2020-02-06 13:00:00
2026-07-22 14:54:08 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:54:08 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1
2026-07-22 14:54:08 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc
2026-07-22 14:54:08 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-22 14:54:08 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:54:08 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-22 14:54:08 INFO Queuing job for member 1...
2026-07-22 14:54:08 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:54:08 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-22 14:54:09 INFO Found: ['5252897']
2026-07-22 14:54:14 INFO [TGCC-IRENE] Submitted job with ID:['5252897']
2026-07-22 14:54:14 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:54:14 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2
2026-07-22 14:54:14 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc
2026-07-22 14:54:14 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-22 14:54:14 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:54:14 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-22 14:54:14 INFO Queuing job for member 2...
2026-07-22 14:54:14 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:54:14 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-22 14:54:15 INFO Found: ['5252908']
2026-07-22 14:54:20 INFO [TGCC-IRENE] Submitted job with ID:['5252908']
2026-07-22 14:54:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:54:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3
2026-07-22 14:54:20 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc
2026-07-22 14:54:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-22 14:54:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:54:20 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-22 14:54:20 INFO Queuing job for member 3...
2026-07-22 14:54:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:54:20 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-22 14:54:21 INFO Found: ['5252921']
2026-07-22 14:54:26 INFO [TGCC-IRENE] Submitted job with ID:['5252921']
2026-07-22 14:54:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:54:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4
2026-07-22 14:54:26 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc
2026-07-22 14:54:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-22 14:54:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:54:26 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-22 14:54:26 INFO Queuing job for member 4...
2026-07-22 14:54:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:54:26 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-22 14:54:27 INFO Found: ['5252933']
2026-07-22 14:54:32 INFO [TGCC-IRENE] Submitted job with ID:['5252933']
2026-07-22 14:54:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:54:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5
2026-07-22 14:54:32 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc
2026-07-22 14:54:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-22 14:54:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:54:32 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-22 14:54:32 INFO Queuing job for member 5...
2026-07-22 14:54:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:54:32 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-22 14:54:35 INFO Found: ['5252938']
2026-07-22 14:54:40 INFO [TGCC-IRENE] Submitted job with ID:['5252938']
2026-07-22 14:54:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:54:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6
2026-07-22 14:54:40 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc
2026-07-22 14:54:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-22 14:54:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:54:40 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-22 14:54:40 INFO Queuing job for member 6...
2026-07-22 14:54:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:54:40 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-22 14:54:42 INFO Found: ['5252946']
2026-07-22 14:54:47 INFO [TGCC-IRENE] Submitted job with ID:['5252946']
2026-07-22 14:54:47 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:54:47 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7
2026-07-22 14:54:47 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc
2026-07-22 14:54:47 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-22 14:54:47 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:54:47 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-22 14:54:47 INFO Queuing job for member 7...
2026-07-22 14:54:47 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:54:47 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-22 14:54:50 INFO Found: ['5252950']
2026-07-22 14:54:55 INFO [TGCC-IRENE] Submitted job with ID:['5252950']
2026-07-22 14:54:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:54:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8
2026-07-22 14:54:55 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc
2026-07-22 14:54:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-22 14:54:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:54:55 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-22 14:54:55 INFO Queuing job for member 8...
2026-07-22 14:54:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:54:55 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-22 14:54:57 INFO Found: ['5252956']
2026-07-22 14:55:02 INFO [TGCC-IRENE] Submitted job with ID:['5252956']
2026-07-22 14:55:02 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:55:02 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9
2026-07-22 14:55:02 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc
2026-07-22 14:55:02 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-22 14:55:02 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:55:02 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-22 14:55:02 INFO Queuing job for member 9...
2026-07-22 14:55:02 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:55:02 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-22 14:55:04 INFO Found: ['5252959']
2026-07-22 14:55:09 INFO [TGCC-IRENE] Submitted job with ID:['5252959']
2026-07-22 14:55:09 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:55:09 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10
2026-07-22 14:55:09 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc
2026-07-22 14:55:09 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-22 14:55:09 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:55:09 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-22 14:55:09 INFO Queuing job for member 10...
2026-07-22 14:55:09 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:55:09 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-22 14:55:10 INFO Found: ['5252963']
2026-07-22 14:55:15 INFO [TGCC-IRENE] Submitted job with ID:['5252963']
2026-07-22 14:55:15 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:55:15 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11
2026-07-22 14:55:15 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc
2026-07-22 14:55:15 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-22 14:55:15 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:55:15 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-22 14:55:15 INFO Queuing job for member 11...
2026-07-22 14:55:15 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:55:15 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-22 14:55:16 INFO Found: ['5252973']
2026-07-22 14:55:21 INFO [TGCC-IRENE] Submitted job with ID:['5252973']
2026-07-22 14:55:21 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:55:21 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12
2026-07-22 14:55:21 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc
2026-07-22 14:55:21 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-22 14:55:21 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:55:21 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-22 14:55:21 INFO Queuing job for member 12...
2026-07-22 14:55:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:55:21 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-22 14:55:22 INFO Found: ['5252977']
2026-07-22 14:55:27 INFO [TGCC-IRENE] Submitted job with ID:['5252977']
2026-07-22 14:55:27 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:55:27 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13
2026-07-22 14:55:27 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc
2026-07-22 14:55:27 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-22 14:55:27 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:55:27 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-22 14:55:27 INFO Queuing job for member 13...
2026-07-22 14:55:27 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:55:27 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-22 14:55:27 INFO Found: ['5252980']
2026-07-22 14:55:32 INFO [TGCC-IRENE] Submitted job with ID:['5252980']
2026-07-22 14:55:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:55:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14
2026-07-22 14:55:32 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc
2026-07-22 14:55:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-22 14:55:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:55:32 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-22 14:55:33 INFO Queuing job for member 14...
2026-07-22 14:55:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:55:33 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-22 14:55:33 INFO Found: ['5252988']
2026-07-22 14:55:38 INFO [TGCC-IRENE] Submitted job with ID:['5252988']
2026-07-22 14:55:38 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 14:55:38 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15
2026-07-22 14:55:38 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc
2026-07-22 14:55:38 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-22 14:55:38 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 14:55:38 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-22 14:55:38 INFO Queuing job for member 15...
2026-07-22 14:55:38 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 14:55:38 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-22 14:55:39 INFO Found: ['5252995']
2026-07-22 14:55:44 INFO [TGCC-IRENE] Submitted job with ID:['5252995']
2026-07-22 14:55:44 INFO Checking job status ...
2026-07-22 14:55:44 INFO None 5252897: status RUNNING/PENDING
2026-07-22 14:55:44 INFO None 5252908: status RUNNING/PENDING
2026-07-22 14:55:44 INFO None 5252921: status RUNNING/PENDING
2026-07-22 14:55:44 INFO None 5252933: status RUNNING/PENDING
2026-07-22 14:55:44 INFO None 5252938: status RUNNING/PENDING
2026-07-22 14:55:44 INFO None 5252946: status RUNNING/PENDING
2026-07-22 14:55:44 INFO None 5252950: status RUNNING/PENDING
2026-07-22 14:55:44 INFO None 5252956: status RUNNING/PENDING
2026-07-22 14:55:44 INFO None 5252959: status RUNNING/PENDING
2026-07-22 14:55:44 INFO None 5252963: status RUNNING/PENDING
2026-07-22 14:55:44 INFO None 5252973: status RUNNING/PENDING
2026-07-22 14:55:44 INFO None 5252977: status RUNNING/PENDING
2026-07-22 14:55:44 INFO None 5252980: status RUNNING/PENDING
2026-07-22 14:55:44 INFO None 5252988: status RUNNING/PENDING
2026-07-22 14:55:44 INFO None 5252995: status RUNNING/PENDING
2026-07-22 14:55:44 INFO Jobs still running: ['5252897', '5252908', '5252921', '5252933', '5252938', '5252946', '5252950', '5252956', '5252959', '5252963', '5252973', '5252977', '5252980', '5252988', '5252995']. Waiting...
2026-07-22 14:56:01 INFO None 5252897: status RUNNING/PENDING
2026-07-22 14:56:01 INFO None 5252908: status RUNNING/PENDING
2026-07-22 14:56:01 INFO None 5252921: status RUNNING/PENDING
2026-07-22 14:56:01 INFO None 5252933: status RUNNING/PENDING
2026-07-22 14:56:01 INFO None 5252938: status RUNNING/PENDING
2026-07-22 14:56:01 INFO None 5252946: status RUNNING/PENDING
2026-07-22 14:56:01 INFO None 5252950: status RUNNING/PENDING
2026-07-22 14:56:01 INFO None 5252956: status RUNNING/PENDING
2026-07-22 14:56:01 INFO None 5252959: status RUNNING/PENDING
2026-07-22 14:56:01 INFO None 5252963: status RUNNING/PENDING
2026-07-22 14:56:01 INFO None 5252973: status RUNNING/PENDING
2026-07-22 14:56:01 INFO None 5252977: status RUNNING/PENDING
2026-07-22 14:56:02 INFO None 5252980: status RUNNING/PENDING
2026-07-22 14:56:02 INFO None 5252988: status RUNNING/PENDING
2026-07-22 14:56:02 INFO None 5252995: status RUNNING/PENDING
2026-07-22 14:56:02 INFO Jobs still running: ['5252897', '5252908', '5252921', '5252933', '5252938', '5252946', '5252950', '5252956', '5252959', '5252963', '5252973', '5252977', '5252980', '5252988', '5252995']. Waiting...
2026-07-22 14:56:17 INFO None 5252897: status RUNNING/PENDING
2026-07-22 14:56:17 INFO None 5252908: status RUNNING/PENDING
2026-07-22 14:56:17 INFO None 5252921: status RUNNING/PENDING
2026-07-22 14:56:17 INFO None 5252933: status RUNNING/PENDING
2026-07-22 14:56:17 INFO None 5252938: status RUNNING/PENDING
2026-07-22 14:56:17 INFO None 5252946: status RUNNING/PENDING
2026-07-22 14:56:17 INFO None 5252950: status RUNNING/PENDING
2026-07-22 14:56:17 INFO None 5252956: status RUNNING/PENDING
2026-07-22 14:56:17 INFO None 5252959: status RUNNING/PENDING
2026-07-22 14:56:19 INFO None 5252963: status RUNNING/PENDING
2026-07-22 14:56:19 INFO None 5252973: status RUNNING/PENDING
2026-07-22 14:56:19 INFO None 5252977: status RUNNING/PENDING
2026-07-22 14:56:21 INFO None 5252980: status RUNNING/PENDING
2026-07-22 14:56:22 INFO None 5252988: status RUNNING/PENDING
2026-07-22 14:56:22 INFO None 5252995: status RUNNING/PENDING
2026-07-22 14:56:22 INFO Jobs still running: ['5252897', '5252908', '5252921', '5252933', '5252938', '5252946', '5252950', '5252956', '5252959', '5252963', '5252973', '5252977', '5252980', '5252988', '5252995']. Waiting...
2026-07-22 14:56:37 INFO None 5252897: status RUNNING/PENDING
2026-07-22 14:56:37 INFO None 5252908: status RUNNING/PENDING
2026-07-22 14:56:37 INFO None 5252921: status RUNNING/PENDING
2026-07-22 14:56:37 INFO None 5252933: status RUNNING/PENDING
2026-07-22 14:56:37 INFO None 5252938: status RUNNING/PENDING
2026-07-22 14:56:37 INFO None 5252946: status RUNNING/PENDING
2026-07-22 14:56:37 INFO None 5252950: status RUNNING/PENDING
2026-07-22 14:56:37 INFO None 5252956: status RUNNING/PENDING
2026-07-22 14:56:37 INFO None 5252959: status RUNNING/PENDING
2026-07-22 14:56:37 INFO None 5252963: status RUNNING/PENDING
2026-07-22 14:56:37 INFO None 5252973: status RUNNING/PENDING
2026-07-22 14:56:37 INFO None 5252977: status RUNNING/PENDING
2026-07-22 14:56:37 INFO None 5252980: status RUNNING/PENDING
2026-07-22 14:56:37 INFO None 5252988: status RUNNING/PENDING
2026-07-22 14:56:37 INFO None 5252995: status RUNNING/PENDING
2026-07-22 14:56:37 INFO Jobs still running: ['5252897', '5252908', '5252921', '5252933', '5252938', '5252946', '5252950', '5252956', '5252959', '5252963', '5252973', '5252977', '5252980', '5252988', '5252995']. Waiting...
2026-07-22 14:56:52 INFO None 5252897: status RUNNING/PENDING
2026-07-22 14:56:52 INFO None 5252908: status RUNNING/PENDING
2026-07-22 14:56:52 INFO None 5252921: status RUNNING/PENDING
2026-07-22 14:56:52 INFO None 5252933: status RUNNING/PENDING
2026-07-22 14:56:52 INFO None 5252938: status RUNNING/PENDING
2026-07-22 14:56:52 INFO None 5252946: status RUNNING/PENDING
2026-07-22 14:56:52 INFO None 5252950: status RUNNING/PENDING
2026-07-22 14:56:52 INFO None 5252956: status RUNNING/PENDING
2026-07-22 14:56:52 INFO None 5252959: status RUNNING/PENDING
2026-07-22 14:56:52 INFO None 5252963: status RUNNING/PENDING
2026-07-22 14:56:52 INFO None 5252973: status RUNNING/PENDING
2026-07-22 14:56:52 INFO None 5252977: status RUNNING/PENDING
2026-07-22 14:56:52 INFO None 5252980: status RUNNING/PENDING
2026-07-22 14:56:52 INFO None 5252988: status RUNNING/PENDING
2026-07-22 14:56:52 INFO None 5252995: status RUNNING/PENDING
2026-07-22 14:56:52 INFO Jobs still running: ['5252897', '5252908', '5252921', '5252933', '5252938', '5252946', '5252950', '5252956', '5252959', '5252963', '5252973', '5252977', '5252980', '5252988', '5252995']. Waiting...
2026-07-22 14:57:07 INFO None 5252897: status RUNNING/PENDING
2026-07-22 14:57:07 INFO None 5252908: status RUNNING/PENDING
2026-07-22 14:57:07 INFO None 5252921: status RUNNING/PENDING
2026-07-22 14:57:07 INFO None 5252933: status RUNNING/PENDING
2026-07-22 14:57:07 INFO None 5252938: status RUNNING/PENDING
2026-07-22 14:57:07 INFO None 5252946: status RUNNING/PENDING
2026-07-22 14:57:07 INFO None 5252950: status RUNNING/PENDING
2026-07-22 14:57:07 INFO None 5252956: status RUNNING/PENDING
2026-07-22 14:57:07 INFO None 5252959: status RUNNING/PENDING
2026-07-22 14:57:07 INFO None 5252963: status RUNNING/PENDING
2026-07-22 14:57:07 INFO None 5252973: status RUNNING/PENDING
2026-07-22 14:57:07 INFO None 5252977: status RUNNING/PENDING
2026-07-22 14:57:07 INFO None 5252980: status RUNNING/PENDING
2026-07-22 14:57:07 INFO None 5252988: status RUNNING/PENDING
2026-07-22 14:57:07 INFO None 5252995: status RUNNING/PENDING
2026-07-22 14:57:07 INFO Jobs still running: ['5252897', '5252908', '5252921', '5252933', '5252938', '5252946', '5252950', '5252956', '5252959', '5252963', '5252973', '5252977', '5252980', '5252988', '5252995']. Waiting...
2026-07-22 14:57:24 INFO None 5252897: status RUNNING/PENDING
2026-07-22 14:57:24 INFO None 5252908: status RUNNING/PENDING
2026-07-22 14:57:24 INFO None 5252921: status RUNNING/PENDING
2026-07-22 14:57:24 INFO None 5252933: status RUNNING/PENDING
2026-07-22 14:57:24 INFO None 5252938: status RUNNING/PENDING
2026-07-22 14:57:24 INFO None 5252946: status RUNNING/PENDING
2026-07-22 14:57:24 INFO None 5252950: status RUNNING/PENDING
2026-07-22 14:57:24 INFO None 5252956: status RUNNING/PENDING
2026-07-22 14:57:24 INFO None 5252959: status RUNNING/PENDING
2026-07-22 14:57:24 INFO None 5252963: status RUNNING/PENDING
2026-07-22 14:57:24 INFO None 5252973: status RUNNING/PENDING
2026-07-22 14:57:24 INFO None 5252977: status RUNNING/PENDING
2026-07-22 14:57:24 INFO None 5252980: status RUNNING/PENDING
2026-07-22 14:57:24 INFO None 5252988: status RUNNING/PENDING
2026-07-22 14:57:24 INFO None 5252995: status RUNNING/PENDING
2026-07-22 14:57:24 INFO Jobs still running: ['5252897', '5252908', '5252921', '5252933', '5252938', '5252946', '5252950', '5252956', '5252959', '5252963', '5252973', '5252977', '5252980', '5252988', '5252995']. Waiting...
2026-07-22 14:57:39 INFO None 5252897: status FINISHED
2026-07-22 14:57:39 INFO None 5252908: status RUNNING/PENDING
2026-07-22 14:57:42 INFO None 5252921: status RUNNING/PENDING
2026-07-22 14:57:42 INFO None 5252933: status RUNNING/PENDING
2026-07-22 14:57:42 INFO None 5252938: status RUNNING/PENDING
2026-07-22 14:57:42 INFO None 5252946: status RUNNING/PENDING
2026-07-22 14:57:42 INFO None 5252950: status RUNNING/PENDING
2026-07-22 14:57:42 INFO None 5252956: status RUNNING/PENDING
2026-07-22 14:57:42 INFO None 5252959: status RUNNING/PENDING
2026-07-22 14:57:42 INFO None 5252963: status RUNNING/PENDING
2026-07-22 14:57:42 INFO None 5252973: status RUNNING/PENDING
2026-07-22 14:57:42 INFO None 5252977: status RUNNING/PENDING
2026-07-22 14:57:42 INFO None 5252980: status RUNNING/PENDING
2026-07-22 14:57:42 INFO None 5252988: status RUNNING/PENDING
2026-07-22 14:57:42 INFO None 5252995: status RUNNING/PENDING
2026-07-22 14:57:42 INFO Jobs still running: ['5252908', '5252921', '5252933', '5252938', '5252946', '5252950', '5252956', '5252959', '5252963', '5252973', '5252977', '5252980', '5252988', '5252995']. Waiting...
2026-07-22 14:57:57 INFO None 5252897: status FINISHED
2026-07-22 14:57:57 INFO None 5252908: status RUNNING/PENDING
2026-07-22 14:57:57 INFO None 5252921: status RUNNING/PENDING
2026-07-22 14:57:57 INFO None 5252933: status RUNNING/PENDING
2026-07-22 14:57:57 INFO None 5252938: status RUNNING/PENDING
2026-07-22 14:57:57 INFO None 5252946: status RUNNING/PENDING
2026-07-22 14:57:57 INFO None 5252950: status RUNNING/PENDING
2026-07-22 14:57:57 INFO None 5252956: status RUNNING/PENDING
2026-07-22 14:57:58 INFO None 5252959: status RUNNING/PENDING
2026-07-22 14:57:58 INFO None 5252963: status RUNNING/PENDING
2026-07-22 14:57:58 INFO None 5252973: status RUNNING/PENDING
2026-07-22 14:57:58 INFO None 5252977: status RUNNING/PENDING
2026-07-22 14:57:58 INFO None 5252980: status RUNNING/PENDING
2026-07-22 14:57:58 INFO None 5252988: status RUNNING/PENDING
2026-07-22 14:57:58 INFO None 5252995: status RUNNING/PENDING
2026-07-22 14:57:58 INFO Jobs still running: ['5252908', '5252921', '5252933', '5252938', '5252946', '5252950', '5252956', '5252959', '5252963', '5252973', '5252977', '5252980', '5252988', '5252995']. Waiting...
2026-07-22 14:58:13 INFO None 5252897: status FINISHED
2026-07-22 14:58:13 INFO None 5252908: status RUNNING/PENDING
2026-07-22 14:58:13 INFO None 5252921: status FINISHED
2026-07-22 14:58:13 INFO None 5252933: status RUNNING/PENDING
2026-07-22 14:58:13 INFO None 5252938: status RUNNING/PENDING
2026-07-22 14:58:13 INFO None 5252946: status RUNNING/PENDING
2026-07-22 14:58:13 INFO None 5252950: status RUNNING/PENDING
2026-07-22 14:58:13 INFO None 5252956: status RUNNING/PENDING
2026-07-22 14:58:13 INFO None 5252959: status RUNNING/PENDING
2026-07-22 14:58:13 INFO None 5252963: status RUNNING/PENDING
2026-07-22 14:58:13 INFO None 5252973: status RUNNING/PENDING
2026-07-22 14:58:13 INFO None 5252977: status RUNNING/PENDING
2026-07-22 14:58:13 INFO None 5252980: status RUNNING/PENDING
2026-07-22 14:58:14 INFO None 5252988: status RUNNING/PENDING
2026-07-22 14:58:14 INFO None 5252995: status RUNNING/PENDING
2026-07-22 14:58:14 INFO Jobs still running: ['5252908', '5252933', '5252938', '5252946', '5252950', '5252956', '5252959', '5252963', '5252973', '5252977', '5252980', '5252988', '5252995']. Waiting...
2026-07-22 14:58:29 INFO None 5252897: status FINISHED
2026-07-22 14:58:29 INFO None 5252908: status RUNNING/PENDING
2026-07-22 14:58:29 INFO None 5252921: status FINISHED
2026-07-22 14:58:29 INFO None 5252933: status RUNNING/PENDING
2026-07-22 14:58:29 INFO None 5252938: status RUNNING/PENDING
2026-07-22 14:58:29 INFO None 5252946: status RUNNING/PENDING
2026-07-22 14:58:29 INFO None 5252950: status RUNNING/PENDING
2026-07-22 14:58:29 INFO None 5252956: status RUNNING/PENDING
2026-07-22 14:58:29 INFO None 5252959: status RUNNING/PENDING
2026-07-22 14:58:29 INFO None 5252963: status RUNNING/PENDING
2026-07-22 14:58:29 INFO None 5252973: status RUNNING/PENDING
2026-07-22 14:58:29 INFO None 5252977: status RUNNING/PENDING
2026-07-22 14:58:29 INFO None 5252980: status RUNNING/PENDING
2026-07-22 14:58:29 INFO None 5252988: status RUNNING/PENDING
2026-07-22 14:58:29 INFO None 5252995: status RUNNING/PENDING
2026-07-22 14:58:29 INFO Jobs still running: ['5252908', '5252933', '5252938', '5252946', '5252950', '5252956', '5252959', '5252963', '5252973', '5252977', '5252980', '5252988', '5252995']. Waiting...
2026-07-22 14:58:45 INFO None 5252897: status FINISHED
2026-07-22 14:58:45 INFO None 5252908: status RUNNING/PENDING
2026-07-22 14:58:45 INFO None 5252921: status FINISHED
2026-07-22 14:58:45 INFO None 5252933: status RUNNING/PENDING
2026-07-22 14:58:45 INFO None 5252938: status RUNNING/PENDING
2026-07-22 14:58:45 INFO None 5252946: status RUNNING/PENDING
2026-07-22 14:58:45 INFO None 5252950: status RUNNING/PENDING
2026-07-22 14:58:45 INFO None 5252956: status RUNNING/PENDING
2026-07-22 14:58:45 INFO None 5252959: status RUNNING/PENDING
2026-07-22 14:58:45 INFO None 5252963: status RUNNING/PENDING
2026-07-22 14:58:45 INFO None 5252973: status RUNNING/PENDING
2026-07-22 14:58:46 INFO None 5252977: status RUNNING/PENDING
2026-07-22 14:58:48 INFO None 5252980: status RUNNING/PENDING
2026-07-22 14:58:48 INFO None 5252988: status RUNNING/PENDING
2026-07-22 14:58:48 INFO None 5252995: status RUNNING/PENDING
2026-07-22 14:58:48 INFO Jobs still running: ['5252908', '5252933', '5252938', '5252946', '5252950', '5252956', '5252959', '5252963', '5252973', '5252977', '5252980', '5252988', '5252995']. Waiting...
2026-07-22 14:59:03 INFO None 5252897: status FINISHED
2026-07-22 14:59:03 INFO None 5252908: status FINISHED
2026-07-22 14:59:03 INFO None 5252921: status FINISHED
2026-07-22 14:59:03 INFO None 5252933: status RUNNING/PENDING
2026-07-22 14:59:03 INFO None 5252938: status RUNNING/PENDING
2026-07-22 14:59:03 INFO None 5252946: status FINISHED
2026-07-22 14:59:03 INFO None 5252950: status RUNNING/PENDING
2026-07-22 14:59:03 INFO None 5252956: status RUNNING/PENDING
2026-07-22 14:59:03 INFO None 5252959: status RUNNING/PENDING
2026-07-22 14:59:03 INFO None 5252963: status RUNNING/PENDING
2026-07-22 14:59:03 INFO None 5252973: status RUNNING/PENDING
2026-07-22 14:59:03 INFO None 5252977: status RUNNING/PENDING
2026-07-22 14:59:03 INFO None 5252980: status RUNNING/PENDING
2026-07-22 14:59:05 INFO None 5252988: status RUNNING/PENDING
2026-07-22 14:59:05 INFO None 5252995: status RUNNING/PENDING
2026-07-22 14:59:05 INFO Jobs still running: ['5252933', '5252938', '5252950', '5252956', '5252959', '5252963', '5252973', '5252977', '5252980', '5252988', '5252995']. Waiting...
2026-07-22 14:59:20 INFO None 5252897: status FINISHED
2026-07-22 14:59:20 INFO None 5252908: status FINISHED
2026-07-22 14:59:20 INFO None 5252921: status FINISHED
2026-07-22 14:59:20 INFO None 5252933: status RUNNING/PENDING
2026-07-22 14:59:20 INFO None 5252938: status RUNNING/PENDING
2026-07-22 14:59:21 INFO None 5252946: status FINISHED
2026-07-22 14:59:21 INFO None 5252950: status RUNNING/PENDING
2026-07-22 14:59:21 INFO None 5252956: status RUNNING/PENDING
2026-07-22 14:59:21 INFO None 5252959: status RUNNING/PENDING
2026-07-22 14:59:21 INFO None 5252963: status RUNNING/PENDING
2026-07-22 14:59:21 INFO None 5252973: status RUNNING/PENDING
2026-07-22 14:59:21 INFO None 5252977: status RUNNING/PENDING
2026-07-22 14:59:21 INFO None 5252980: status RUNNING/PENDING
2026-07-22 14:59:21 INFO None 5252988: status RUNNING/PENDING
2026-07-22 14:59:21 INFO None 5252995: status RUNNING/PENDING
2026-07-22 14:59:21 INFO Jobs still running: ['5252933', '5252938', '5252950', '5252956', '5252959', '5252963', '5252973', '5252977', '5252980', '5252988', '5252995']. Waiting...
2026-07-22 14:59:36 INFO None 5252897: status FINISHED
2026-07-22 14:59:36 INFO None 5252908: status FINISHED
2026-07-22 14:59:36 INFO None 5252921: status FINISHED
2026-07-22 14:59:36 INFO None 5252933: status RUNNING/PENDING
2026-07-22 14:59:36 INFO None 5252938: status RUNNING/PENDING
2026-07-22 14:59:36 INFO None 5252946: status FINISHED
2026-07-22 14:59:36 INFO None 5252950: status RUNNING/PENDING
2026-07-22 14:59:36 INFO None 5252956: status RUNNING/PENDING
2026-07-22 14:59:36 INFO None 5252959: status RUNNING/PENDING
2026-07-22 14:59:36 INFO None 5252963: status RUNNING/PENDING
2026-07-22 14:59:36 INFO None 5252973: status RUNNING/PENDING
2026-07-22 14:59:36 INFO None 5252977: status RUNNING/PENDING
2026-07-22 14:59:36 INFO None 5252980: status RUNNING/PENDING
2026-07-22 14:59:36 INFO None 5252988: status RUNNING/PENDING
2026-07-22 14:59:36 INFO None 5252995: status RUNNING/PENDING
2026-07-22 14:59:36 INFO Jobs still running: ['5252933', '5252938', '5252950', '5252956', '5252959', '5252963', '5252973', '5252977', '5252980', '5252988', '5252995']. Waiting...
2026-07-22 14:59:51 INFO None 5252897: status FINISHED
2026-07-22 14:59:51 INFO None 5252908: status FINISHED
2026-07-22 14:59:51 INFO None 5252921: status FINISHED
2026-07-22 14:59:51 INFO None 5252933: status RUNNING/PENDING
2026-07-22 14:59:51 INFO None 5252938: status FINISHED
2026-07-22 14:59:51 INFO None 5252946: status FINISHED
2026-07-22 14:59:51 INFO None 5252950: status RUNNING/PENDING
2026-07-22 14:59:51 INFO None 5252956: status RUNNING/PENDING
2026-07-22 14:59:51 INFO None 5252959: status RUNNING/PENDING
2026-07-22 14:59:51 INFO None 5252963: status RUNNING/PENDING
2026-07-22 14:59:51 INFO None 5252973: status RUNNING/PENDING
2026-07-22 14:59:51 INFO None 5252977: status RUNNING/PENDING
2026-07-22 14:59:51 INFO None 5252980: status RUNNING/PENDING
2026-07-22 14:59:51 INFO None 5252988: status RUNNING/PENDING
2026-07-22 14:59:51 INFO None 5252995: status RUNNING/PENDING
2026-07-22 14:59:51 INFO Jobs still running: ['5252933', '5252950', '5252956', '5252959', '5252963', '5252973', '5252977', '5252980', '5252988', '5252995']. Waiting...
2026-07-22 15:00:07 INFO None 5252897: status FINISHED
2026-07-22 15:00:07 INFO None 5252908: status FINISHED
2026-07-22 15:00:07 INFO None 5252921: status FINISHED
2026-07-22 15:00:07 INFO None 5252933: status RUNNING/PENDING
2026-07-22 15:00:07 INFO None 5252938: status FINISHED
2026-07-22 15:00:07 INFO None 5252946: status FINISHED
2026-07-22 15:00:07 INFO None 5252950: status FINISHED
2026-07-22 15:00:07 INFO None 5252956: status FINISHED
2026-07-22 15:00:07 INFO None 5252959: status FINISHED
2026-07-22 15:00:07 INFO None 5252963: status FINISHED
2026-07-22 15:00:07 INFO None 5252973: status FINISHED
2026-07-22 15:00:07 INFO None 5252977: status FINISHED
2026-07-22 15:00:07 INFO None 5252980: status RUNNING/PENDING
2026-07-22 15:00:08 INFO None 5252988: status FINISHED
2026-07-22 15:00:08 INFO None 5252995: status RUNNING/PENDING
2026-07-22 15:00:08 INFO Jobs still running: ['5252933', '5252980', '5252995']. Waiting...
2026-07-22 15:00:23 INFO None 5252897: status FINISHED
2026-07-22 15:00:23 INFO None 5252908: status FINISHED
2026-07-22 15:00:23 INFO None 5252921: status FINISHED
2026-07-22 15:00:23 INFO None 5252933: status RUNNING/PENDING
2026-07-22 15:00:23 INFO None 5252938: status FINISHED
2026-07-22 15:00:23 INFO None 5252946: status FINISHED
2026-07-22 15:00:23 INFO None 5252950: status FINISHED
2026-07-22 15:00:23 INFO None 5252956: status FINISHED
2026-07-22 15:00:23 INFO None 5252959: status FINISHED
2026-07-22 15:00:23 INFO None 5252963: status FINISHED
2026-07-22 15:00:23 INFO None 5252973: status FINISHED
2026-07-22 15:00:23 INFO None 5252977: status FINISHED
2026-07-22 15:00:23 INFO None 5252980: status RUNNING/PENDING
2026-07-22 15:00:23 INFO None 5252988: status FINISHED
2026-07-22 15:00:23 INFO None 5252995: status FINISHED
2026-07-22 15:00:23 INFO Jobs still running: ['5252933', '5252980']. Waiting...
2026-07-22 15:00:38 INFO None 5252897: status FINISHED
2026-07-22 15:00:38 INFO None 5252908: status FINISHED
2026-07-22 15:00:40 INFO None 5252921: status FINISHED
2026-07-22 15:00:40 INFO None 5252933: status RUNNING/PENDING
2026-07-22 15:00:40 INFO None 5252938: status FINISHED
2026-07-22 15:00:40 INFO None 5252946: status FINISHED
2026-07-22 15:00:40 INFO None 5252950: status FINISHED
2026-07-22 15:00:40 INFO None 5252956: status FINISHED
2026-07-22 15:00:40 INFO None 5252959: status FINISHED
2026-07-22 15:00:40 INFO None 5252963: status FINISHED
2026-07-22 15:00:40 INFO None 5252973: status FINISHED
2026-07-22 15:00:40 INFO None 5252977: status FINISHED
2026-07-22 15:00:40 INFO None 5252980: status RUNNING/PENDING
2026-07-22 15:00:40 INFO None 5252988: status FINISHED
2026-07-22 15:00:40 INFO None 5252995: status FINISHED
2026-07-22 15:00:40 INFO Jobs still running: ['5252933', '5252980']. Waiting...
2026-07-22 15:00:55 INFO None 5252897: status FINISHED
2026-07-22 15:00:55 INFO None 5252908: status FINISHED
2026-07-22 15:00:55 INFO None 5252921: status FINISHED
2026-07-22 15:00:55 INFO None 5252933: status RUNNING/PENDING
2026-07-22 15:00:55 INFO None 5252938: status FINISHED
2026-07-22 15:00:55 INFO None 5252946: status FINISHED
2026-07-22 15:00:55 INFO None 5252950: status FINISHED
2026-07-22 15:00:55 INFO None 5252956: status FINISHED
2026-07-22 15:00:55 INFO None 5252959: status FINISHED
2026-07-22 15:00:55 INFO None 5252963: status FINISHED
2026-07-22 15:00:55 INFO None 5252973: status FINISHED
2026-07-22 15:00:55 INFO None 5252977: status FINISHED
2026-07-22 15:00:56 INFO None 5252980: status RUNNING/PENDING
2026-07-22 15:00:56 INFO None 5252988: status FINISHED
2026-07-22 15:00:56 INFO None 5252995: status FINISHED
2026-07-22 15:00:56 INFO Jobs still running: ['5252933', '5252980']. Waiting...
2026-07-22 15:01:11 INFO None 5252897: status FINISHED
2026-07-22 15:01:11 INFO None 5252908: status FINISHED
2026-07-22 15:01:11 INFO None 5252921: status FINISHED
2026-07-22 15:01:11 INFO None 5252933: status FINISHED
2026-07-22 15:01:11 INFO None 5252938: status FINISHED
2026-07-22 15:01:11 INFO None 5252946: status FINISHED
2026-07-22 15:01:11 INFO None 5252950: status FINISHED
2026-07-22 15:01:11 INFO None 5252956: status FINISHED
2026-07-22 15:01:11 INFO None 5252959: status FINISHED
2026-07-22 15:01:11 INFO None 5252963: status FINISHED
2026-07-22 15:01:11 INFO None 5252973: status FINISHED
2026-07-22 15:01:11 INFO None 5252977: status FINISHED
2026-07-22 15:01:11 INFO None 5252980: status RUNNING/PENDING
2026-07-22 15:01:11 INFO None 5252988: status FINISHED
2026-07-22 15:01:11 INFO None 5252995: status FINISHED
2026-07-22 15:01:11 INFO Jobs still running: ['5252980']. Waiting...
2026-07-22 15:01:28 INFO None 5252897: status FINISHED
2026-07-22 15:01:28 INFO None 5252908: status FINISHED
2026-07-22 15:01:28 INFO None 5252921: status FINISHED
2026-07-22 15:01:28 INFO None 5252933: status FINISHED
2026-07-22 15:01:28 INFO None 5252938: status FINISHED
2026-07-22 15:01:28 INFO None 5252946: status FINISHED
2026-07-22 15:01:28 INFO None 5252950: status FINISHED
2026-07-22 15:01:28 INFO None 5252956: status FINISHED
2026-07-22 15:01:28 INFO None 5252959: status FINISHED
2026-07-22 15:01:28 INFO None 5252963: status FINISHED
2026-07-22 15:01:28 INFO None 5252973: status FINISHED
2026-07-22 15:01:28 INFO None 5252977: status FINISHED
2026-07-22 15:01:28 INFO None 5252980: status FINISHED
2026-07-22 15:01:28 INFO None 5252988: status FINISHED
2026-07-22 15:01:28 INFO None 5252995: status FINISHED
2026-07-22 15:01:28 INFO Jobs ['5252897', '5252908', '5252921', '5252933', '5252938', '5252946', '5252950', '5252956', '5252959', '5252963', '5252973', '5252977', '5252980', '5252988', '5252995'] have finished
2026-07-22 15:01:28 INFO Checking restart files were created ...
2026-07-22 15:01:28 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020611_2_ENS1.nc(1002685915 bytes)
2026-07-22 15:01:28 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020611_2_ENS2.nc(1002685915 bytes)
2026-07-22 15:01:28 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020611_2_ENS3.nc(1002685915 bytes)
2026-07-22 15:01:28 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020611_2_ENS4.nc(1002685915 bytes)
2026-07-22 15:01:28 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020611_2_ENS5.nc(1002685915 bytes)
2026-07-22 15:01:28 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020611_2_ENS6.nc(1002685915 bytes)
2026-07-22 15:01:28 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020611_2_ENS7.nc(1002685915 bytes)
2026-07-22 15:01:28 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020611_2_ENS8.nc(1002685915 bytes)
2026-07-22 15:01:28 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020611_2_ENS9.nc(1002685915 bytes)
2026-07-22 15:01:28 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020611_2_ENS10.nc(1002685915 bytes)
2026-07-22 15:01:28 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020611_2_ENS11.nc(1002685915 bytes)
2026-07-22 15:01:28 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020611_2_ENS12.nc(1002685915 bytes)
2026-07-22 15:01:28 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020611_2_ENS13.nc(1002685915 bytes)
2026-07-22 15:01:28 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020611_2_ENS14.nc(1002685915 bytes)
2026-07-22 15:01:28 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020611_2_ENS15.nc(1002685915 bytes)
2026-07-22 15:01:28 INFO  Run_model() completed successfully.
2026-07-22 15:01:28 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 15:01:28 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 13:00:00 days=153072 seconds=46800
2026-07-22 15:01:28 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-22 15:01:28 INFO [TIME] increment current_time 2020-02-06 11:00:00 -> 2020-02-06 13:00:00
2026-07-22 15:01:28 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 15:01:28 INFO ---------->>> Running process_satellite_data()
2026-07-22 15:01:28 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12004.nc
2026-07-22 15:01:28 INFO ---------->>> Running run_obs_converter()
2026-07-22 15:01:28 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_45135_153072.out
2026-07-22 15:01:28 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_45135_153072.out
2026-07-22 15:01:28 INFO ---------->>> Running DART
2026-07-22 15:01:28 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-22 15:01:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/chim_ENS1_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/chim_ENS1_2020020613_1_out_toDART.nc
2026-07-22 15:01:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/chim_ENS2_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/chim_ENS2_2020020613_1_out_toDART.nc
2026-07-22 15:01:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/chim_ENS3_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/chim_ENS3_2020020613_1_out_toDART.nc
2026-07-22 15:01:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/chim_ENS4_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/chim_ENS4_2020020613_1_out_toDART.nc
2026-07-22 15:01:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/chim_ENS5_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/chim_ENS5_2020020613_1_out_toDART.nc
2026-07-22 15:01:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/chim_ENS6_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/chim_ENS6_2020020613_1_out_toDART.nc
2026-07-22 15:01:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/chim_ENS7_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/chim_ENS7_2020020613_1_out_toDART.nc
2026-07-22 15:01:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/chim_ENS8_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/chim_ENS8_2020020613_1_out_toDART.nc
2026-07-22 15:01:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/chim_ENS9_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/chim_ENS9_2020020613_1_out_toDART.nc
2026-07-22 15:01:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/chim_ENS10_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/chim_ENS10_2020020613_1_out_toDART.nc
2026-07-22 15:01:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/chim_ENS11_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/chim_ENS11_2020020613_1_out_toDART.nc
2026-07-22 15:01:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/chim_ENS12_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/chim_ENS12_2020020613_1_out_toDART.nc
2026-07-22 15:01:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/chim_ENS13_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/chim_ENS13_2020020613_1_out_toDART.nc
2026-07-22 15:01:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/chim_ENS14_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/chim_ENS14_2020020613_1_out_toDART.nc
2026-07-22 15:01:33 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/chim_ENS15_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/chim_ENS15_2020020613_1_out_toDART.nc
2026-07-22 15:01:33 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-22 15:01:33 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-22 15:01:33 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-22 15:01:33 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-22 15:01:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-22 15:01:33 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-22 15:01:47 INFO Found: []
2026-07-22 15:01:47 INFO No job id returned by command ./run_filter.bsh
2026-07-22 15:01:47 INFO No monitoring will be performed
2026-07-22 15:01:47 INFO Moving DART output files to analysis and preassim directories for date 2020020613 if present ...
2026-07-22 15:01:47 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 15:01:47 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 15:01:47 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 15:01:47 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 15:01:47 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 15:01:47 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 15:01:47 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 15:01:47 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 15:01:47 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 15:01:47 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 15:01:48 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 15:01:48 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 15:01:48 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 15:01:48 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 15:01:48 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 15:01:48 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 15:01:48 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 15:01:48 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 15:01:48 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 15:01:48 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 15:01:48 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 15:01:48 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 15:01:48 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 15:01:48 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 15:01:48 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 15:01:48 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 15:01:48 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 15:01:48 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 15:01:48 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 15:01:48 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 15:01:48 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 15:01:48 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 15:01:48 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 15:01:48 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 15:01:48 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-22 15:01:48 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-22 15:01:48 INFO run_dart() is DONE.
2026-07-22 15:01:48 INFO ---------->>> Running update_pollutant_in_end()
2026-07-22 15:01:48 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:01:48 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS1.nc
2026-07-22 15:01:49 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:01:50 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS1.nc
2026-07-22 15:01:50 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:01:50 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:01:50 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS1.nc
2026-07-22 15:01:50 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:01:51 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS1.nc
2026-07-22 15:01:51 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:01:51 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 15:01:51 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:01:51 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS2.nc
2026-07-22 15:01:52 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:01:53 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS2.nc
2026-07-22 15:01:53 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:01:53 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:01:53 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS2.nc
2026-07-22 15:01:54 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:01:54 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS2.nc
2026-07-22 15:01:54 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:01:54 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 15:01:55 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:01:55 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS3.nc
2026-07-22 15:01:55 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:01:56 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS3.nc
2026-07-22 15:01:56 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:01:56 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:01:56 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS3.nc
2026-07-22 15:01:57 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:01:57 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS3.nc
2026-07-22 15:01:57 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:01:58 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 15:01:58 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:01:58 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS4.nc
2026-07-22 15:01:58 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:02:00 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS4.nc
2026-07-22 15:02:00 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:02:00 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:02:00 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS4.nc
2026-07-22 15:02:01 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:02:01 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS4.nc
2026-07-22 15:02:01 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:02:01 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 15:02:02 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:02:02 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS5.nc
2026-07-22 15:02:02 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:02:03 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS5.nc
2026-07-22 15:02:03 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:02:03 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:02:03 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS5.nc
2026-07-22 15:02:04 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:02:04 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS5.nc
2026-07-22 15:02:04 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:02:04 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 15:02:05 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:02:05 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS6.nc
2026-07-22 15:02:06 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:02:06 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS6.nc
2026-07-22 15:02:06 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:02:07 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:02:07 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS6.nc
2026-07-22 15:02:07 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:02:08 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS6.nc
2026-07-22 15:02:08 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:02:08 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 15:02:08 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:02:08 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS7.nc
2026-07-22 15:02:09 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:02:09 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS7.nc
2026-07-22 15:02:09 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:02:10 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:02:10 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS7.nc
2026-07-22 15:02:10 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:02:11 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS7.nc
2026-07-22 15:02:11 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:02:11 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 15:02:11 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:02:11 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS8.nc
2026-07-22 15:02:12 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:02:13 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS8.nc
2026-07-22 15:02:13 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:02:13 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:02:13 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS8.nc
2026-07-22 15:02:14 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:02:14 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS8.nc
2026-07-22 15:02:14 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:02:14 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 15:02:15 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:02:15 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS9.nc
2026-07-22 15:02:15 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:02:16 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS9.nc
2026-07-22 15:02:16 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:02:17 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:02:17 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS9.nc
2026-07-22 15:02:17 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:02:18 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS9.nc
2026-07-22 15:02:18 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:02:18 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 15:02:18 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:02:18 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS10.nc
2026-07-22 15:02:19 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:02:19 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS10.nc
2026-07-22 15:02:19 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:02:20 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:02:20 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS10.nc
2026-07-22 15:02:20 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:02:21 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS10.nc
2026-07-22 15:02:21 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:02:21 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 15:02:21 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:02:21 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS11.nc
2026-07-22 15:02:22 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:02:23 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS11.nc
2026-07-22 15:02:23 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:02:23 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:02:23 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS11.nc
2026-07-22 15:02:24 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:02:24 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS11.nc
2026-07-22 15:02:24 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:02:24 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 15:02:25 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:02:25 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS12.nc
2026-07-22 15:02:25 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:02:26 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS12.nc
2026-07-22 15:02:26 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:02:26 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:02:26 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS12.nc
2026-07-22 15:02:27 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:02:28 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS12.nc
2026-07-22 15:02:28 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:02:28 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 15:02:28 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:02:28 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS13.nc
2026-07-22 15:02:29 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:02:29 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS13.nc
2026-07-22 15:02:29 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:02:30 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:02:30 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS13.nc
2026-07-22 15:02:30 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:02:31 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS13.nc
2026-07-22 15:02:31 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:02:31 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 15:02:31 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:02:31 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS14.nc
2026-07-22 15:02:32 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:02:33 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS14.nc
2026-07-22 15:02:33 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:02:33 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:02:33 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS14.nc
2026-07-22 15:02:33 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:02:34 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS14.nc
2026-07-22 15:02:34 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:02:34 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 15:02:34 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:02:34 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS15.nc
2026-07-22 15:02:35 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:02:36 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS15.nc
2026-07-22 15:02:36 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:02:36 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:02:36 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS15.nc
2026-07-22 15:02:37 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:02:37 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS15.nc
2026-07-22 15:02:37 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:02:37 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 15:02:37 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-22 15:02:37 INFO [TIME] step_end current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 15:02:37 INFO [TIME] step_start current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 15:02:37 INFO [TIME] window start=2020-02-06 13:00:00 end=2020-02-06 14:00:00 run_hours=1 has_assimilation=True
2026-07-22 15:02:37 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 15:02:39 INFO Hourly dataset computed and listing created
2026-07-22 15:02:41 INFO Hourly dataset computed
2026-07-22 15:02:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 15:02:42 INFO Hourly dataset computed and listing created
2026-07-22 15:02:42 INFO Hourly dataset computed
2026-07-22 15:02:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 15:02:43 INFO Hourly dataset computed and listing created
2026-07-22 15:02:44 INFO Hourly dataset computed
2026-07-22 15:02:44 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 15:02:45 INFO Hourly dataset computed and listing created
2026-07-22 15:02:45 INFO Hourly dataset computed
2026-07-22 15:02:45 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 15:02:46 INFO Hourly dataset computed and listing created
2026-07-22 15:02:47 INFO Hourly dataset computed
2026-07-22 15:02:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 15:02:48 INFO Hourly dataset computed and listing created
2026-07-22 15:02:48 INFO Hourly dataset computed
2026-07-22 15:02:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 15:02:49 INFO Hourly dataset computed and listing created
2026-07-22 15:02:50 INFO Hourly dataset computed
2026-07-22 15:02:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 15:02:51 INFO Hourly dataset computed and listing created
2026-07-22 15:02:51 INFO Hourly dataset computed
2026-07-22 15:02:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 15:02:52 INFO Hourly dataset computed and listing created
2026-07-22 15:02:53 INFO Hourly dataset computed
2026-07-22 15:02:53 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 15:02:54 INFO Hourly dataset computed and listing created
2026-07-22 15:02:54 INFO Hourly dataset computed
2026-07-22 15:02:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 15:02:55 INFO Hourly dataset computed and listing created
2026-07-22 15:02:56 INFO Hourly dataset computed
2026-07-22 15:02:56 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 15:02:57 INFO Hourly dataset computed and listing created
2026-07-22 15:02:57 INFO Hourly dataset computed
2026-07-22 15:02:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 15:02:58 INFO Hourly dataset computed and listing created
2026-07-22 15:02:59 INFO Hourly dataset computed
2026-07-22 15:02:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 15:03:00 INFO Hourly dataset computed and listing created
2026-07-22 15:03:00 INFO Hourly dataset computed
2026-07-22 15:03:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 15:03:01 INFO Hourly dataset computed and listing created
2026-07-22 15:03:02 INFO Hourly dataset computed
2026-07-22 15:03:02 INFO ---------->>> Running CHIMERE model from 2020-02-06 13:00:00 to 2020-02-06 14:00:00
2026-07-22 15:03:02 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 15:03:02 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1
2026-07-22 15:03:02 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020611_2_ENS1.nc
2026-07-22 15:03:02 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-22 15:03:02 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 15:03:02 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-22 15:03:02 INFO Queuing job for member 1...
2026-07-22 15:03:02 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 15:03:02 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-22 15:03:04 INFO Found: ['5253212']
2026-07-22 15:03:09 INFO [TGCC-IRENE] Submitted job with ID:['5253212']
2026-07-22 15:03:09 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 15:03:09 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2
2026-07-22 15:03:09 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020611_2_ENS2.nc
2026-07-22 15:03:09 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-22 15:03:09 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 15:03:09 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-22 15:03:09 INFO Queuing job for member 2...
2026-07-22 15:03:09 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 15:03:09 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-22 15:03:11 INFO Found: ['5253216']
2026-07-22 15:03:16 INFO [TGCC-IRENE] Submitted job with ID:['5253216']
2026-07-22 15:03:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 15:03:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3
2026-07-22 15:03:16 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020611_2_ENS3.nc
2026-07-22 15:03:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-22 15:03:17 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 15:03:17 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-22 15:03:17 INFO Queuing job for member 3...
2026-07-22 15:03:17 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 15:03:17 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-22 15:03:18 INFO Found: ['5253220']
2026-07-22 15:03:23 INFO [TGCC-IRENE] Submitted job with ID:['5253220']
2026-07-22 15:03:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 15:03:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4
2026-07-22 15:03:23 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020611_2_ENS4.nc
2026-07-22 15:03:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-22 15:03:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 15:03:23 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-22 15:03:23 INFO Queuing job for member 4...
2026-07-22 15:03:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 15:03:23 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-22 15:03:24 INFO Found: ['5253231']
2026-07-22 15:03:29 INFO [TGCC-IRENE] Submitted job with ID:['5253231']
2026-07-22 15:03:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 15:03:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5
2026-07-22 15:03:29 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020611_2_ENS5.nc
2026-07-22 15:03:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-22 15:03:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 15:03:29 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-22 15:03:29 INFO Queuing job for member 5...
2026-07-22 15:03:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 15:03:29 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-22 15:03:30 INFO Found: ['5253235']
2026-07-22 15:03:35 INFO [TGCC-IRENE] Submitted job with ID:['5253235']
2026-07-22 15:03:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 15:03:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6
2026-07-22 15:03:35 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020611_2_ENS6.nc
2026-07-22 15:03:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-22 15:03:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 15:03:35 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-22 15:03:35 INFO Queuing job for member 6...
2026-07-22 15:03:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 15:03:35 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-22 15:03:36 INFO Found: ['5253237']
2026-07-22 15:03:41 INFO [TGCC-IRENE] Submitted job with ID:['5253237']
2026-07-22 15:03:41 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 15:03:41 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7
2026-07-22 15:03:41 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020611_2_ENS7.nc
2026-07-22 15:03:41 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-22 15:03:41 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 15:03:41 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-22 15:03:41 INFO Queuing job for member 7...
2026-07-22 15:03:41 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 15:03:41 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-22 15:03:41 INFO Found: ['5253241']
2026-07-22 15:03:46 INFO [TGCC-IRENE] Submitted job with ID:['5253241']
2026-07-22 15:03:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 15:03:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8
2026-07-22 15:03:46 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020611_2_ENS8.nc
2026-07-22 15:03:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-22 15:03:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 15:03:46 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-22 15:03:46 INFO Queuing job for member 8...
2026-07-22 15:03:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 15:03:46 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-22 15:03:47 INFO Found: ['5253246']
2026-07-22 15:03:52 INFO [TGCC-IRENE] Submitted job with ID:['5253246']
2026-07-22 15:03:52 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 15:03:52 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9
2026-07-22 15:03:52 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020611_2_ENS9.nc
2026-07-22 15:03:52 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-22 15:03:52 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 15:03:52 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-22 15:03:52 INFO Queuing job for member 9...
2026-07-22 15:03:52 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 15:03:52 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-22 15:03:53 INFO Found: ['5253254']
2026-07-22 15:03:58 INFO [TGCC-IRENE] Submitted job with ID:['5253254']
2026-07-22 15:03:58 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 15:03:58 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10
2026-07-22 15:03:58 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020611_2_ENS10.nc
2026-07-22 15:03:58 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-22 15:03:58 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 15:03:58 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-22 15:03:58 INFO Queuing job for member 10...
2026-07-22 15:03:58 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 15:03:58 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-22 15:03:59 INFO Found: ['5253264']
2026-07-22 15:04:04 INFO [TGCC-IRENE] Submitted job with ID:['5253264']
2026-07-22 15:04:04 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 15:04:04 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11
2026-07-22 15:04:04 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020611_2_ENS11.nc
2026-07-22 15:04:04 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-22 15:04:04 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 15:04:04 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-22 15:04:04 INFO Queuing job for member 11...
2026-07-22 15:04:04 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 15:04:04 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-22 15:04:05 INFO Found: ['5253272']
2026-07-22 15:04:10 INFO [TGCC-IRENE] Submitted job with ID:['5253272']
2026-07-22 15:04:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 15:04:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12
2026-07-22 15:04:10 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020611_2_ENS12.nc
2026-07-22 15:04:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-22 15:04:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 15:04:10 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-22 15:04:10 INFO Queuing job for member 12...
2026-07-22 15:04:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 15:04:10 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-22 15:04:13 INFO Found: ['5253276']
2026-07-22 15:04:18 INFO [TGCC-IRENE] Submitted job with ID:['5253276']
2026-07-22 15:04:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 15:04:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13
2026-07-22 15:04:18 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020611_2_ENS13.nc
2026-07-22 15:04:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-22 15:04:18 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 15:04:18 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-22 15:04:18 INFO Queuing job for member 13...
2026-07-22 15:04:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 15:04:18 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-22 15:04:20 INFO Found: ['5253280']
2026-07-22 15:04:25 INFO [TGCC-IRENE] Submitted job with ID:['5253280']
2026-07-22 15:04:25 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 15:04:25 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14
2026-07-22 15:04:25 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020611_2_ENS14.nc
2026-07-22 15:04:25 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-22 15:04:25 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 15:04:25 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-22 15:04:25 INFO Queuing job for member 14...
2026-07-22 15:04:25 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 15:04:25 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-22 15:04:28 INFO Found: ['5253282']
2026-07-22 15:04:33 INFO [TGCC-IRENE] Submitted job with ID:['5253282']
2026-07-22 15:04:33 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 15:04:33 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15
2026-07-22 15:04:33 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020611_2_ENS15.nc
2026-07-22 15:04:33 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-22 15:04:33 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 15:04:33 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-22 15:04:33 INFO Queuing job for member 15...
2026-07-22 15:04:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 15:04:33 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-22 15:04:35 INFO Found: ['5253284']
2026-07-22 15:04:40 INFO [TGCC-IRENE] Submitted job with ID:['5253284']
2026-07-22 15:04:40 INFO Checking job status ...
2026-07-22 15:04:40 INFO None 5253212: status RUNNING/PENDING
2026-07-22 15:04:40 INFO None 5253216: status RUNNING/PENDING
2026-07-22 15:04:40 INFO None 5253220: status RUNNING/PENDING
2026-07-22 15:04:40 INFO None 5253231: status RUNNING/PENDING
2026-07-22 15:04:41 INFO None 5253235: status RUNNING/PENDING
2026-07-22 15:04:41 INFO None 5253237: status RUNNING/PENDING
2026-07-22 15:04:41 INFO None 5253241: status RUNNING/PENDING
2026-07-22 15:04:41 INFO None 5253246: status RUNNING/PENDING
2026-07-22 15:04:41 INFO None 5253254: status RUNNING/PENDING
2026-07-22 15:04:41 INFO None 5253264: status RUNNING/PENDING
2026-07-22 15:04:41 INFO None 5253272: status RUNNING/PENDING
2026-07-22 15:04:41 INFO None 5253276: status RUNNING/PENDING
2026-07-22 15:04:41 INFO None 5253280: status RUNNING/PENDING
2026-07-22 15:04:41 INFO None 5253282: status RUNNING/PENDING
2026-07-22 15:04:41 INFO None 5253284: status RUNNING/PENDING
2026-07-22 15:04:41 INFO Jobs still running: ['5253212', '5253216', '5253220', '5253231', '5253235', '5253237', '5253241', '5253246', '5253254', '5253264', '5253272', '5253276', '5253280', '5253282', '5253284']. Waiting...
2026-07-22 15:04:56 INFO None 5253212: status RUNNING/PENDING
2026-07-22 15:04:56 INFO None 5253216: status RUNNING/PENDING
2026-07-22 15:04:56 INFO None 5253220: status RUNNING/PENDING
2026-07-22 15:04:56 INFO None 5253231: status RUNNING/PENDING
2026-07-22 15:04:56 INFO None 5253235: status RUNNING/PENDING
2026-07-22 15:04:56 INFO None 5253237: status RUNNING/PENDING
2026-07-22 15:04:56 INFO None 5253241: status RUNNING/PENDING
2026-07-22 15:04:56 INFO None 5253246: status RUNNING/PENDING
2026-07-22 15:04:56 INFO None 5253254: status RUNNING/PENDING
2026-07-22 15:04:56 INFO None 5253264: status RUNNING/PENDING
2026-07-22 15:04:56 INFO None 5253272: status RUNNING/PENDING
2026-07-22 15:04:56 INFO None 5253276: status RUNNING/PENDING
2026-07-22 15:04:56 INFO None 5253280: status RUNNING/PENDING
2026-07-22 15:04:56 INFO None 5253282: status RUNNING/PENDING
2026-07-22 15:04:56 INFO None 5253284: status RUNNING/PENDING
2026-07-22 15:04:56 INFO Jobs still running: ['5253212', '5253216', '5253220', '5253231', '5253235', '5253237', '5253241', '5253246', '5253254', '5253264', '5253272', '5253276', '5253280', '5253282', '5253284']. Waiting...
2026-07-22 15:05:11 INFO None 5253212: status RUNNING/PENDING
2026-07-22 15:05:11 INFO None 5253216: status RUNNING/PENDING
2026-07-22 15:05:11 INFO None 5253220: status RUNNING/PENDING
2026-07-22 15:05:11 INFO None 5253231: status RUNNING/PENDING
2026-07-22 15:05:11 INFO None 5253235: status RUNNING/PENDING
2026-07-22 15:05:11 INFO None 5253237: status RUNNING/PENDING
2026-07-22 15:05:11 INFO None 5253241: status RUNNING/PENDING
2026-07-22 15:05:11 INFO None 5253246: status RUNNING/PENDING
2026-07-22 15:05:11 INFO None 5253254: status RUNNING/PENDING
2026-07-22 15:05:11 INFO None 5253264: status RUNNING/PENDING
2026-07-22 15:05:11 INFO None 5253272: status RUNNING/PENDING
2026-07-22 15:05:11 INFO None 5253276: status RUNNING/PENDING
2026-07-22 15:05:11 INFO None 5253280: status RUNNING/PENDING
2026-07-22 15:05:11 INFO None 5253282: status RUNNING/PENDING
2026-07-22 15:05:11 INFO None 5253284: status RUNNING/PENDING
2026-07-22 15:05:11 INFO Jobs still running: ['5253212', '5253216', '5253220', '5253231', '5253235', '5253237', '5253241', '5253246', '5253254', '5253264', '5253272', '5253276', '5253280', '5253282', '5253284']. Waiting...
2026-07-22 15:05:28 INFO None 5253212: status RUNNING/PENDING
2026-07-22 15:05:28 INFO None 5253216: status RUNNING/PENDING
2026-07-22 15:05:28 INFO None 5253220: status RUNNING/PENDING
2026-07-22 15:05:28 INFO None 5253231: status RUNNING/PENDING
2026-07-22 15:05:28 INFO None 5253235: status RUNNING/PENDING
2026-07-22 15:05:28 INFO None 5253237: status RUNNING/PENDING
2026-07-22 15:05:29 INFO None 5253241: status RUNNING/PENDING
2026-07-22 15:05:29 INFO None 5253246: status RUNNING/PENDING
2026-07-22 15:05:29 INFO None 5253254: status RUNNING/PENDING
2026-07-22 15:05:29 INFO None 5253264: status RUNNING/PENDING
2026-07-22 15:05:29 INFO None 5253272: status RUNNING/PENDING
2026-07-22 15:05:29 INFO None 5253276: status RUNNING/PENDING
2026-07-22 15:05:29 INFO None 5253280: status RUNNING/PENDING
2026-07-22 15:05:29 INFO None 5253282: status RUNNING/PENDING
2026-07-22 15:05:29 INFO None 5253284: status RUNNING/PENDING
2026-07-22 15:05:29 INFO Jobs still running: ['5253212', '5253216', '5253220', '5253231', '5253235', '5253237', '5253241', '5253246', '5253254', '5253264', '5253272', '5253276', '5253280', '5253282', '5253284']. Waiting...
2026-07-22 15:05:44 INFO None 5253212: status RUNNING/PENDING
2026-07-22 15:05:44 INFO None 5253216: status RUNNING/PENDING
2026-07-22 15:05:44 INFO None 5253220: status RUNNING/PENDING
2026-07-22 15:05:44 INFO None 5253231: status RUNNING/PENDING
2026-07-22 15:05:44 INFO None 5253235: status RUNNING/PENDING
2026-07-22 15:05:44 INFO None 5253237: status RUNNING/PENDING
2026-07-22 15:05:44 INFO None 5253241: status RUNNING/PENDING
2026-07-22 15:05:44 INFO None 5253246: status RUNNING/PENDING
2026-07-22 15:05:44 INFO None 5253254: status RUNNING/PENDING
2026-07-22 15:05:46 INFO None 5253264: status RUNNING/PENDING
2026-07-22 15:05:46 INFO None 5253272: status RUNNING/PENDING
2026-07-22 15:05:46 INFO None 5253276: status RUNNING/PENDING
2026-07-22 15:05:46 INFO None 5253280: status RUNNING/PENDING
2026-07-22 15:05:46 INFO None 5253282: status RUNNING/PENDING
2026-07-22 15:05:46 INFO None 5253284: status RUNNING/PENDING
2026-07-22 15:05:46 INFO Jobs still running: ['5253212', '5253216', '5253220', '5253231', '5253235', '5253237', '5253241', '5253246', '5253254', '5253264', '5253272', '5253276', '5253280', '5253282', '5253284']. Waiting...
2026-07-22 15:06:01 INFO None 5253212: status RUNNING/PENDING
2026-07-22 15:06:01 INFO None 5253216: status RUNNING/PENDING
2026-07-22 15:06:01 INFO None 5253220: status RUNNING/PENDING
2026-07-22 15:06:01 INFO None 5253231: status RUNNING/PENDING
2026-07-22 15:06:01 INFO None 5253235: status RUNNING/PENDING
2026-07-22 15:06:01 INFO None 5253237: status RUNNING/PENDING
2026-07-22 15:06:01 INFO None 5253241: status RUNNING/PENDING
2026-07-22 15:06:01 INFO None 5253246: status RUNNING/PENDING
2026-07-22 15:06:01 INFO None 5253254: status RUNNING/PENDING
2026-07-22 15:06:01 INFO None 5253264: status RUNNING/PENDING
2026-07-22 15:06:02 INFO None 5253272: status RUNNING/PENDING
2026-07-22 15:06:02 INFO None 5253276: status RUNNING/PENDING
2026-07-22 15:06:02 INFO None 5253280: status RUNNING/PENDING
2026-07-22 15:06:02 INFO None 5253282: status RUNNING/PENDING
2026-07-22 15:06:02 INFO None 5253284: status RUNNING/PENDING
2026-07-22 15:06:02 INFO Jobs still running: ['5253212', '5253216', '5253220', '5253231', '5253235', '5253237', '5253241', '5253246', '5253254', '5253264', '5253272', '5253276', '5253280', '5253282', '5253284']. Waiting...
2026-07-22 15:06:17 INFO None 5253212: status FINISHED
2026-07-22 15:06:17 INFO None 5253216: status RUNNING/PENDING
2026-07-22 15:06:17 INFO None 5253220: status RUNNING/PENDING
2026-07-22 15:06:17 INFO None 5253231: status FINISHED
2026-07-22 15:06:17 INFO None 5253235: status RUNNING/PENDING
2026-07-22 15:06:17 INFO None 5253237: status RUNNING/PENDING
2026-07-22 15:06:17 INFO None 5253241: status RUNNING/PENDING
2026-07-22 15:06:17 INFO None 5253246: status RUNNING/PENDING
2026-07-22 15:06:17 INFO None 5253254: status FINISHED
2026-07-22 15:06:17 INFO None 5253264: status FINISHED
2026-07-22 15:06:17 INFO None 5253272: status RUNNING/PENDING
2026-07-22 15:06:17 INFO None 5253276: status RUNNING/PENDING
2026-07-22 15:06:17 INFO None 5253280: status RUNNING/PENDING
2026-07-22 15:06:17 INFO None 5253282: status RUNNING/PENDING
2026-07-22 15:06:17 INFO None 5253284: status RUNNING/PENDING
2026-07-22 15:06:17 INFO Jobs still running: ['5253216', '5253220', '5253235', '5253237', '5253241', '5253246', '5253272', '5253276', '5253280', '5253282', '5253284']. Waiting...
2026-07-22 15:06:32 INFO None 5253212: status FINISHED
2026-07-22 15:06:32 INFO None 5253216: status RUNNING/PENDING
2026-07-22 15:06:32 INFO None 5253220: status RUNNING/PENDING
2026-07-22 15:06:32 INFO None 5253231: status FINISHED
2026-07-22 15:06:32 INFO None 5253235: status RUNNING/PENDING
2026-07-22 15:06:32 INFO None 5253237: status RUNNING/PENDING
2026-07-22 15:06:32 INFO None 5253241: status RUNNING/PENDING
2026-07-22 15:06:32 INFO None 5253246: status FINISHED
2026-07-22 15:06:32 INFO None 5253254: status FINISHED
2026-07-22 15:06:32 INFO None 5253264: status FINISHED
2026-07-22 15:06:32 INFO None 5253272: status RUNNING/PENDING
2026-07-22 15:06:32 INFO None 5253276: status RUNNING/PENDING
2026-07-22 15:06:32 INFO None 5253280: status RUNNING/PENDING
2026-07-22 15:06:32 INFO None 5253282: status RUNNING/PENDING
2026-07-22 15:06:32 INFO None 5253284: status RUNNING/PENDING
2026-07-22 15:06:32 INFO Jobs still running: ['5253216', '5253220', '5253235', '5253237', '5253241', '5253272', '5253276', '5253280', '5253282', '5253284']. Waiting...
2026-07-22 15:06:49 INFO None 5253212: status FINISHED
2026-07-22 15:06:49 INFO None 5253216: status RUNNING/PENDING
2026-07-22 15:06:49 INFO None 5253220: status RUNNING/PENDING
2026-07-22 15:06:49 INFO None 5253231: status FINISHED
2026-07-22 15:06:49 INFO None 5253235: status RUNNING/PENDING
2026-07-22 15:06:49 INFO None 5253237: status RUNNING/PENDING
2026-07-22 15:06:49 INFO None 5253241: status RUNNING/PENDING
2026-07-22 15:06:49 INFO None 5253246: status FINISHED
2026-07-22 15:06:49 INFO None 5253254: status FINISHED
2026-07-22 15:06:49 INFO None 5253264: status FINISHED
2026-07-22 15:06:49 INFO None 5253272: status RUNNING/PENDING
2026-07-22 15:06:49 INFO None 5253276: status RUNNING/PENDING
2026-07-22 15:06:51 INFO None 5253280: status RUNNING/PENDING
2026-07-22 15:06:51 INFO None 5253282: status RUNNING/PENDING
2026-07-22 15:06:51 INFO None 5253284: status RUNNING/PENDING
2026-07-22 15:06:51 INFO Jobs still running: ['5253216', '5253220', '5253235', '5253237', '5253241', '5253272', '5253276', '5253280', '5253282', '5253284']. Waiting...
2026-07-22 15:07:07 INFO None 5253212: status FINISHED
2026-07-22 15:07:07 INFO None 5253216: status RUNNING/PENDING
2026-07-22 15:07:07 INFO None 5253220: status FINISHED
2026-07-22 15:07:07 INFO None 5253231: status FINISHED
2026-07-22 15:07:07 INFO None 5253235: status FINISHED
2026-07-22 15:07:07 INFO None 5253237: status FINISHED
2026-07-22 15:07:07 INFO None 5253241: status FINISHED
2026-07-22 15:07:07 INFO None 5253246: status FINISHED
2026-07-22 15:07:07 INFO None 5253254: status FINISHED
2026-07-22 15:07:07 INFO None 5253264: status FINISHED
2026-07-22 15:07:07 INFO None 5253272: status RUNNING/PENDING
2026-07-22 15:07:07 INFO None 5253276: status RUNNING/PENDING
2026-07-22 15:07:07 INFO None 5253280: status RUNNING/PENDING
2026-07-22 15:07:07 INFO None 5253282: status RUNNING/PENDING
2026-07-22 15:07:07 INFO None 5253284: status RUNNING/PENDING
2026-07-22 15:07:07 INFO Jobs still running: ['5253216', '5253272', '5253276', '5253280', '5253282', '5253284']. Waiting...
2026-07-22 15:07:22 INFO None 5253212: status FINISHED
2026-07-22 15:07:22 INFO None 5253216: status RUNNING/PENDING
2026-07-22 15:07:22 INFO None 5253220: status FINISHED
2026-07-22 15:07:22 INFO None 5253231: status FINISHED
2026-07-22 15:07:22 INFO None 5253235: status FINISHED
2026-07-22 15:07:22 INFO None 5253237: status FINISHED
2026-07-22 15:07:22 INFO None 5253241: status FINISHED
2026-07-22 15:07:22 INFO None 5253246: status FINISHED
2026-07-22 15:07:22 INFO None 5253254: status FINISHED
2026-07-22 15:07:23 INFO None 5253264: status FINISHED
2026-07-22 15:07:23 INFO None 5253272: status RUNNING/PENDING
2026-07-22 15:07:23 INFO None 5253276: status RUNNING/PENDING
2026-07-22 15:07:23 INFO None 5253280: status RUNNING/PENDING
2026-07-22 15:07:23 INFO None 5253282: status RUNNING/PENDING
2026-07-22 15:07:23 INFO None 5253284: status RUNNING/PENDING
2026-07-22 15:07:23 INFO Jobs still running: ['5253216', '5253272', '5253276', '5253280', '5253282', '5253284']. Waiting...
2026-07-22 15:07:38 INFO None 5253212: status FINISHED
2026-07-22 15:07:38 INFO None 5253216: status FINISHED
2026-07-22 15:07:38 INFO None 5253220: status FINISHED
2026-07-22 15:07:38 INFO None 5253231: status FINISHED
2026-07-22 15:07:38 INFO None 5253235: status FINISHED
2026-07-22 15:07:38 INFO None 5253237: status FINISHED
2026-07-22 15:07:38 INFO None 5253241: status FINISHED
2026-07-22 15:07:38 INFO None 5253246: status FINISHED
2026-07-22 15:07:38 INFO None 5253254: status FINISHED
2026-07-22 15:07:38 INFO None 5253264: status FINISHED
2026-07-22 15:07:38 INFO None 5253272: status RUNNING/PENDING
2026-07-22 15:07:38 INFO None 5253276: status RUNNING/PENDING
2026-07-22 15:07:38 INFO None 5253280: status RUNNING/PENDING
2026-07-22 15:07:39 INFO None 5253282: status RUNNING/PENDING
2026-07-22 15:07:39 INFO None 5253284: status RUNNING/PENDING
2026-07-22 15:07:39 INFO Jobs still running: ['5253272', '5253276', '5253280', '5253282', '5253284']. Waiting...
2026-07-22 15:07:54 INFO None 5253212: status FINISHED
2026-07-22 15:07:54 INFO None 5253216: status FINISHED
2026-07-22 15:07:54 INFO None 5253220: status FINISHED
2026-07-22 15:07:54 INFO None 5253231: status FINISHED
2026-07-22 15:07:54 INFO None 5253235: status FINISHED
2026-07-22 15:07:54 INFO None 5253237: status FINISHED
2026-07-22 15:07:54 INFO None 5253241: status FINISHED
2026-07-22 15:07:54 INFO None 5253246: status FINISHED
2026-07-22 15:07:54 INFO None 5253254: status FINISHED
2026-07-22 15:07:54 INFO None 5253264: status FINISHED
2026-07-22 15:07:54 INFO None 5253272: status RUNNING/PENDING
2026-07-22 15:07:54 INFO None 5253276: status RUNNING/PENDING
2026-07-22 15:07:54 INFO None 5253280: status RUNNING/PENDING
2026-07-22 15:07:54 INFO None 5253282: status RUNNING/PENDING
2026-07-22 15:07:54 INFO None 5253284: status FINISHED
2026-07-22 15:07:54 INFO Jobs still running: ['5253272', '5253276', '5253280', '5253282']. Waiting...
2026-07-22 15:08:10 INFO None 5253212: status FINISHED
2026-07-22 15:08:10 INFO None 5253216: status FINISHED
2026-07-22 15:08:10 INFO None 5253220: status FINISHED
2026-07-22 15:08:10 INFO None 5253231: status FINISHED
2026-07-22 15:08:10 INFO None 5253235: status FINISHED
2026-07-22 15:08:10 INFO None 5253237: status FINISHED
2026-07-22 15:08:10 INFO None 5253241: status FINISHED
2026-07-22 15:08:10 INFO None 5253246: status FINISHED
2026-07-22 15:08:10 INFO None 5253254: status FINISHED
2026-07-22 15:08:10 INFO None 5253264: status FINISHED
2026-07-22 15:08:10 INFO None 5253272: status FINISHED
2026-07-22 15:08:10 INFO None 5253276: status RUNNING/PENDING
2026-07-22 15:08:10 INFO None 5253280: status RUNNING/PENDING
2026-07-22 15:08:10 INFO None 5253282: status RUNNING/PENDING
2026-07-22 15:08:10 INFO None 5253284: status FINISHED
2026-07-22 15:08:10 INFO Jobs still running: ['5253276', '5253280', '5253282']. Waiting...
2026-07-22 15:08:25 INFO None 5253212: status FINISHED
2026-07-22 15:08:25 INFO None 5253216: status FINISHED
2026-07-22 15:08:25 INFO None 5253220: status FINISHED
2026-07-22 15:08:26 INFO None 5253231: status FINISHED
2026-07-22 15:08:26 INFO None 5253235: status FINISHED
2026-07-22 15:08:26 INFO None 5253237: status FINISHED
2026-07-22 15:08:26 INFO None 5253241: status FINISHED
2026-07-22 15:08:26 INFO None 5253246: status FINISHED
2026-07-22 15:08:26 INFO None 5253254: status FINISHED
2026-07-22 15:08:28 INFO None 5253264: status FINISHED
2026-07-22 15:08:28 INFO None 5253272: status FINISHED
2026-07-22 15:08:28 INFO None 5253276: status RUNNING/PENDING
2026-07-22 15:08:28 INFO None 5253280: status RUNNING/PENDING
2026-07-22 15:08:28 INFO None 5253282: status RUNNING/PENDING
2026-07-22 15:08:28 INFO None 5253284: status FINISHED
2026-07-22 15:08:28 INFO Jobs still running: ['5253276', '5253280', '5253282']. Waiting...
2026-07-22 15:08:43 INFO None 5253212: status FINISHED
2026-07-22 15:08:43 INFO None 5253216: status FINISHED
2026-07-22 15:08:43 INFO None 5253220: status FINISHED
2026-07-22 15:08:43 INFO None 5253231: status FINISHED
2026-07-22 15:08:43 INFO None 5253235: status FINISHED
2026-07-22 15:08:43 INFO None 5253237: status FINISHED
2026-07-22 15:08:43 INFO None 5253241: status FINISHED
2026-07-22 15:08:43 INFO None 5253246: status FINISHED
2026-07-22 15:08:43 INFO None 5253254: status FINISHED
2026-07-22 15:08:43 INFO None 5253264: status FINISHED
2026-07-22 15:08:43 INFO None 5253272: status FINISHED
2026-07-22 15:08:43 INFO None 5253276: status RUNNING/PENDING
2026-07-22 15:08:43 INFO None 5253280: status RUNNING/PENDING
2026-07-22 15:08:43 INFO None 5253282: status RUNNING/PENDING
2026-07-22 15:08:45 INFO None 5253284: status FINISHED
2026-07-22 15:08:45 INFO Jobs still running: ['5253276', '5253280', '5253282']. Waiting...
2026-07-22 15:09:00 INFO None 5253212: status FINISHED
2026-07-22 15:09:00 INFO None 5253216: status FINISHED
2026-07-22 15:09:00 INFO None 5253220: status FINISHED
2026-07-22 15:09:01 INFO None 5253231: status FINISHED
2026-07-22 15:09:01 INFO None 5253235: status FINISHED
2026-07-22 15:09:01 INFO None 5253237: status FINISHED
2026-07-22 15:09:01 INFO None 5253241: status FINISHED
2026-07-22 15:09:01 INFO None 5253246: status FINISHED
2026-07-22 15:09:01 INFO None 5253254: status FINISHED
2026-07-22 15:09:01 INFO None 5253264: status FINISHED
2026-07-22 15:09:01 INFO None 5253272: status FINISHED
2026-07-22 15:09:01 INFO None 5253276: status RUNNING/PENDING
2026-07-22 15:09:01 INFO None 5253280: status FINISHED
2026-07-22 15:09:01 INFO None 5253282: status FINISHED
2026-07-22 15:09:01 INFO None 5253284: status FINISHED
2026-07-22 15:09:01 INFO Jobs still running: ['5253276']. Waiting...
2026-07-22 15:09:16 INFO None 5253212: status FINISHED
2026-07-22 15:09:16 INFO None 5253216: status FINISHED
2026-07-22 15:09:16 INFO None 5253220: status FINISHED
2026-07-22 15:09:16 INFO None 5253231: status FINISHED
2026-07-22 15:09:16 INFO None 5253235: status FINISHED
2026-07-22 15:09:16 INFO None 5253237: status FINISHED
2026-07-22 15:09:16 INFO None 5253241: status FINISHED
2026-07-22 15:09:16 INFO None 5253246: status FINISHED
2026-07-22 15:09:16 INFO None 5253254: status FINISHED
2026-07-22 15:09:16 INFO None 5253264: status FINISHED
2026-07-22 15:09:16 INFO None 5253272: status FINISHED
2026-07-22 15:09:16 INFO None 5253276: status RUNNING/PENDING
2026-07-22 15:09:16 INFO None 5253280: status FINISHED
2026-07-22 15:09:16 INFO None 5253282: status FINISHED
2026-07-22 15:09:16 INFO None 5253284: status FINISHED
2026-07-22 15:09:16 INFO Jobs still running: ['5253276']. Waiting...
2026-07-22 15:09:31 INFO None 5253212: status FINISHED
2026-07-22 15:09:31 INFO None 5253216: status FINISHED
2026-07-22 15:09:31 INFO None 5253220: status FINISHED
2026-07-22 15:09:31 INFO None 5253231: status FINISHED
2026-07-22 15:09:31 INFO None 5253235: status FINISHED
2026-07-22 15:09:31 INFO None 5253237: status FINISHED
2026-07-22 15:09:31 INFO None 5253241: status FINISHED
2026-07-22 15:09:31 INFO None 5253246: status FINISHED
2026-07-22 15:09:31 INFO None 5253254: status FINISHED
2026-07-22 15:09:31 INFO None 5253264: status FINISHED
2026-07-22 15:09:31 INFO None 5253272: status FINISHED
2026-07-22 15:09:31 INFO None 5253276: status RUNNING/PENDING
2026-07-22 15:09:33 INFO None 5253280: status FINISHED
2026-07-22 15:09:33 INFO None 5253282: status FINISHED
2026-07-22 15:09:33 INFO None 5253284: status FINISHED
2026-07-22 15:09:33 INFO Jobs still running: ['5253276']. Waiting...
2026-07-22 15:09:48 INFO None 5253212: status FINISHED
2026-07-22 15:09:48 INFO None 5253216: status FINISHED
2026-07-22 15:09:49 INFO None 5253220: status FINISHED
2026-07-22 15:09:49 INFO None 5253231: status FINISHED
2026-07-22 15:09:49 INFO None 5253235: status FINISHED
2026-07-22 15:09:49 INFO None 5253237: status FINISHED
2026-07-22 15:09:49 INFO None 5253241: status FINISHED
2026-07-22 15:09:49 INFO None 5253246: status FINISHED
2026-07-22 15:09:49 INFO None 5253254: status FINISHED
2026-07-22 15:09:49 INFO None 5253264: status FINISHED
2026-07-22 15:09:49 INFO None 5253272: status FINISHED
2026-07-22 15:09:49 INFO None 5253276: status RUNNING/PENDING
2026-07-22 15:09:49 INFO None 5253280: status FINISHED
2026-07-22 15:09:49 INFO None 5253282: status FINISHED
2026-07-22 15:09:49 INFO None 5253284: status FINISHED
2026-07-22 15:09:49 INFO Jobs still running: ['5253276']. Waiting...
2026-07-22 15:10:04 INFO None 5253212: status FINISHED
2026-07-22 15:10:04 INFO None 5253216: status FINISHED
2026-07-22 15:10:04 INFO None 5253220: status FINISHED
2026-07-22 15:10:04 INFO None 5253231: status FINISHED
2026-07-22 15:10:04 INFO None 5253235: status FINISHED
2026-07-22 15:10:04 INFO None 5253237: status FINISHED
2026-07-22 15:10:06 INFO None 5253241: status FINISHED
2026-07-22 15:10:06 INFO None 5253246: status FINISHED
2026-07-22 15:10:06 INFO None 5253254: status FINISHED
2026-07-22 15:10:06 INFO None 5253264: status FINISHED
2026-07-22 15:10:06 INFO None 5253272: status FINISHED
2026-07-22 15:10:06 INFO None 5253276: status FINISHED
2026-07-22 15:10:06 INFO None 5253280: status FINISHED
2026-07-22 15:10:06 INFO None 5253282: status FINISHED
2026-07-22 15:10:06 INFO None 5253284: status FINISHED
2026-07-22 15:10:06 INFO Jobs ['5253212', '5253216', '5253220', '5253231', '5253235', '5253237', '5253241', '5253246', '5253254', '5253264', '5253272', '5253276', '5253280', '5253282', '5253284'] have finished
2026-07-22 15:10:06 INFO Checking restart files were created ...
2026-07-22 15:10:06 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020613_1_ENS1.nc(668832435 bytes)
2026-07-22 15:10:06 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020613_1_ENS2.nc(668832435 bytes)
2026-07-22 15:10:06 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020613_1_ENS3.nc(668832435 bytes)
2026-07-22 15:10:06 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020613_1_ENS4.nc(668832435 bytes)
2026-07-22 15:10:06 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020613_1_ENS5.nc(668832435 bytes)
2026-07-22 15:10:06 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020613_1_ENS6.nc(668832435 bytes)
2026-07-22 15:10:06 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020613_1_ENS7.nc(668832435 bytes)
2026-07-22 15:10:06 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020613_1_ENS8.nc(668832435 bytes)
2026-07-22 15:10:06 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020613_1_ENS9.nc(668832435 bytes)
2026-07-22 15:10:06 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020613_1_ENS10.nc(668832435 bytes)
2026-07-22 15:10:06 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020613_1_ENS11.nc(668832435 bytes)
2026-07-22 15:10:06 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020613_1_ENS12.nc(668832435 bytes)
2026-07-22 15:10:06 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020613_1_ENS13.nc(668832435 bytes)
2026-07-22 15:10:06 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020613_1_ENS14.nc(668832435 bytes)
2026-07-22 15:10:06 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020613_1_ENS15.nc(668832435 bytes)
2026-07-22 15:10:06 INFO  Run_model() completed successfully.
2026-07-22 15:10:06 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 15:10:06 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 14:00:00 days=153072 seconds=50400
2026-07-22 15:10:06 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-22 15:10:06 INFO [TIME] increment current_time 2020-02-06 13:00:00 -> 2020-02-06 14:00:00
2026-07-22 15:10:06 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 15:10:06 INFO ---------->>> Running process_satellite_data()
2026-07-22 15:10:06 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12005.nc
2026-07-22 15:10:06 INFO ---------->>> Running run_obs_converter()
2026-07-22 15:10:06 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_51225_153072.out
2026-07-22 15:10:06 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_51225_153072.out
2026-07-22 15:10:06 INFO ---------->>> Running DART
2026-07-22 15:10:06 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-22 15:10:06 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/chim_ENS1_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/chim_ENS1_2020020614_1_out_toDART.nc
2026-07-22 15:10:07 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/chim_ENS2_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/chim_ENS2_2020020614_1_out_toDART.nc
2026-07-22 15:10:07 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/chim_ENS3_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/chim_ENS3_2020020614_1_out_toDART.nc
2026-07-22 15:10:07 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/chim_ENS4_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/chim_ENS4_2020020614_1_out_toDART.nc
2026-07-22 15:10:07 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/chim_ENS5_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/chim_ENS5_2020020614_1_out_toDART.nc
2026-07-22 15:10:08 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/chim_ENS6_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/chim_ENS6_2020020614_1_out_toDART.nc
2026-07-22 15:10:08 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/chim_ENS7_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/chim_ENS7_2020020614_1_out_toDART.nc
2026-07-22 15:10:08 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/chim_ENS8_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/chim_ENS8_2020020614_1_out_toDART.nc
2026-07-22 15:10:09 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/chim_ENS9_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/chim_ENS9_2020020614_1_out_toDART.nc
2026-07-22 15:10:09 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/chim_ENS10_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/chim_ENS10_2020020614_1_out_toDART.nc
2026-07-22 15:10:09 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/chim_ENS11_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/chim_ENS11_2020020614_1_out_toDART.nc
2026-07-22 15:10:10 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/chim_ENS12_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/chim_ENS12_2020020614_1_out_toDART.nc
2026-07-22 15:10:10 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/chim_ENS13_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/chim_ENS13_2020020614_1_out_toDART.nc
2026-07-22 15:10:10 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/chim_ENS14_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/chim_ENS14_2020020614_1_out_toDART.nc
2026-07-22 15:10:11 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/chim_ENS15_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/chim_ENS15_2020020614_1_out_toDART.nc
2026-07-22 15:10:11 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-22 15:10:11 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-22 15:10:11 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-22 15:10:11 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-22 15:10:11 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-22 15:10:11 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-22 15:10:19 INFO Found: []
2026-07-22 15:10:19 INFO No job id returned by command ./run_filter.bsh
2026-07-22 15:10:19 INFO No monitoring will be performed
2026-07-22 15:10:19 INFO Moving DART output files to analysis and preassim directories for date 2020020614 if present ...
2026-07-22 15:10:19 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 15:10:19 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 15:10:19 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 15:10:19 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 15:10:19 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 15:10:19 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 15:10:19 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 15:10:19 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 15:10:19 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 15:10:19 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 15:10:19 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 15:10:19 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 15:10:19 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 15:10:19 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 15:10:19 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 15:10:19 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 15:10:19 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 15:10:19 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 15:10:19 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 15:10:19 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 15:10:19 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 15:10:19 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 15:10:19 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 15:10:19 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 15:10:20 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 15:10:20 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 15:10:20 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 15:10:20 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 15:10:20 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 15:10:20 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 15:10:20 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 15:10:20 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 15:10:20 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 15:10:20 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 15:10:20 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-22 15:10:20 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-22 15:10:20 INFO run_dart() is DONE.
2026-07-22 15:10:20 INFO ---------->>> Running update_pollutant_in_end()
2026-07-22 15:10:20 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:10:20 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS1.nc
2026-07-22 15:10:21 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:10:21 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS1.nc
2026-07-22 15:10:21 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:10:22 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:10:22 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS1.nc
2026-07-22 15:10:22 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:10:23 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS1.nc
2026-07-22 15:10:23 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:10:23 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 15:10:23 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:10:23 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS2.nc
2026-07-22 15:10:24 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:10:25 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS2.nc
2026-07-22 15:10:25 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:10:25 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:10:25 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS2.nc
2026-07-22 15:10:26 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:10:26 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS2.nc
2026-07-22 15:10:26 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:10:26 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 15:10:26 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:10:26 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS3.nc
2026-07-22 15:10:27 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:10:28 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS3.nc
2026-07-22 15:10:28 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:10:28 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:10:28 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS3.nc
2026-07-22 15:10:29 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:10:29 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS3.nc
2026-07-22 15:10:29 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:10:29 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 15:10:30 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:10:30 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS4.nc
2026-07-22 15:10:30 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:10:31 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS4.nc
2026-07-22 15:10:31 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:10:31 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:10:31 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS4.nc
2026-07-22 15:10:32 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:10:32 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS4.nc
2026-07-22 15:10:32 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:10:33 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 15:10:33 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:10:33 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS5.nc
2026-07-22 15:10:33 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:10:35 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS5.nc
2026-07-22 15:10:35 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:10:35 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:10:35 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS5.nc
2026-07-22 15:10:35 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:10:36 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS5.nc
2026-07-22 15:10:36 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:10:36 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 15:10:36 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:10:36 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS6.nc
2026-07-22 15:10:37 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:10:38 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS6.nc
2026-07-22 15:10:38 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:10:38 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:10:38 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS6.nc
2026-07-22 15:10:39 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:10:39 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS6.nc
2026-07-22 15:10:39 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:10:39 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 15:10:39 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:10:39 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS7.nc
2026-07-22 15:10:40 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:10:41 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS7.nc
2026-07-22 15:10:41 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:10:41 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:10:41 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS7.nc
2026-07-22 15:10:42 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:10:42 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS7.nc
2026-07-22 15:10:42 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:10:42 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 15:10:43 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:10:43 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS8.nc
2026-07-22 15:10:43 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:10:44 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS8.nc
2026-07-22 15:10:44 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:10:44 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:10:44 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS8.nc
2026-07-22 15:10:45 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:10:46 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS8.nc
2026-07-22 15:10:46 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:10:46 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 15:10:46 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:10:46 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS9.nc
2026-07-22 15:10:47 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:10:47 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS9.nc
2026-07-22 15:10:47 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:10:48 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:10:48 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS9.nc
2026-07-22 15:10:48 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:10:49 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS9.nc
2026-07-22 15:10:49 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:10:49 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 15:10:49 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:10:49 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS10.nc
2026-07-22 15:10:50 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:10:50 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS10.nc
2026-07-22 15:10:50 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:10:51 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:10:51 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS10.nc
2026-07-22 15:10:51 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:10:52 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS10.nc
2026-07-22 15:10:52 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:10:52 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 15:10:52 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:10:52 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS11.nc
2026-07-22 15:10:53 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:10:54 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS11.nc
2026-07-22 15:10:54 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:10:54 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:10:54 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS11.nc
2026-07-22 15:10:54 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:10:55 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS11.nc
2026-07-22 15:10:55 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:10:55 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 15:10:56 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:10:56 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS12.nc
2026-07-22 15:10:56 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:10:57 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS12.nc
2026-07-22 15:10:57 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:10:57 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:10:57 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS12.nc
2026-07-22 15:10:58 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:10:58 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS12.nc
2026-07-22 15:10:58 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:10:58 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 15:10:59 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:10:59 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS13.nc
2026-07-22 15:10:59 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:11:00 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS13.nc
2026-07-22 15:11:00 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:11:00 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:11:00 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS13.nc
2026-07-22 15:11:01 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:11:02 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS13.nc
2026-07-22 15:11:02 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:11:02 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 15:11:02 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:11:02 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS14.nc
2026-07-22 15:11:03 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:11:03 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS14.nc
2026-07-22 15:11:03 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:11:04 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:11:04 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS14.nc
2026-07-22 15:11:04 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:11:05 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS14.nc
2026-07-22 15:11:05 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:11:05 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 15:11:05 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:11:05 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS15.nc
2026-07-22 15:11:06 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:11:07 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS15.nc
2026-07-22 15:11:07 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:11:07 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 15:11:07 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS15.nc
2026-07-22 15:11:07 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 15:11:08 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS15.nc
2026-07-22 15:11:08 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 15:11:08 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 15:11:08 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-22 15:11:08 INFO [TIME] step_end current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 15:11:08 INFO [TIME] step_start current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 15:11:08 INFO [TIME] window start=2020-02-06 14:00:00 end=2020-02-07 00:00:00 run_hours=10 has_assimilation=False
Traceback (most recent call last):
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/main.py", line 103, in <module>
    pipeline.run_pipeline()
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/pipelines/base_pipeline.py", line 135, in run_pipeline
    self.before_step()
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/pipelines/chimere2023/pipeline.py", line 220, in before_step
    if self.time_manager.slot_time.strftime("%H") == '00' and not persistence_until_next_day:
                                                                  ^^^^^^^^^^^^^^^^^^^^^^^^^^
NameError: name 'persistence_until_next_day' is not defined. Did you mean: 'self.persistence_until_next_day'?
+ exit 0
