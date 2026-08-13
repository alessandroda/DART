+ SCRIPT_PID=3798525
+ /bin/bash -x /tmp/tmp.tTGMBxEK46
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
2026-07-22 16:28:56 INFO 
███    ███ ██ ███    ███ ███████ ███████ ██
████  ████ ██ ████  ████ ██      ██      ██
██ ████ ██ ██ ██ ████ ██ █████   ███████ ██
██  ██  ██ ██ ██  ██  ██ ██           ██ ██
██      ██ ██ ██      ██ ███████ ███████ ██



2026-07-22 16:28:56 INFO [PIPELINE] =======================================
2026-07-22 16:28:56 INFO [PIPELINE] Starting chimere–DART orchestrator
2026-07-22 16:28:56 INFO [PIPELINE] Config file: config/config_irene_IM.yaml
2026-07-22 16:28:56 INFO [PIPELINE] Run dir: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart
2026-07-22 16:28:56 INFO [PIPELINE] Log file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/mimesi_orchestrator_logs/chimere_DART_20260722_162856.log
2026-07-22 16:28:56 INFO [PIPELINE] =======================================
2026-07-22 16:28:56 INFO Running assimilation with model_type=ModelType.CHIMERE
2026-07-22 16:28:56 INFO Using scheduler=Scheduler.SLURM, queue=rome
2026-07-22 16:28:56 INFO [STEP] ---- TIME LOOP START ----
2026-07-22 16:28:56 INFO [TIME] step_start current_time=2020-02-06 00:00:00 simulated_time=None dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 16:28:56 INFO [TIME] window start=2020-02-06 00:00:00 end=2020-02-06 01:00:00 run_hours=1 has_assimilation=False
2026-07-22 16:28:56 INFO Asked to restart from control run ...
2026-07-22 16:28:56 INFO Copying EMIS ...
2026-07-22 16:28:57 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens01.nc
2026-07-22 16:28:57 INFO Linking first END ...
2026-07-22 16:28:57 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020500_24_ENS1.nc
2026-07-22 16:28:57 INFO >> Checking links...
2026-07-22 16:28:58 INFO >> All links are good for ENS1  ...
2026-07-22 16:28:58 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:29:06 INFO Hourly dataset computed and listing created
2026-07-22 16:29:11 INFO Hourly dataset computed
2026-07-22 16:29:11 INFO Asked to restart from control run ...
2026-07-22 16:29:11 INFO Copying EMIS ...
2026-07-22 16:29:11 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens02.nc
2026-07-22 16:29:11 INFO Linking first END ...
2026-07-22 16:29:11 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020500_24_ENS2.nc
2026-07-22 16:29:11 INFO >> Checking links...
2026-07-22 16:29:12 INFO >> All links are good for ENS2  ...
2026-07-22 16:29:12 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:29:13 INFO Hourly dataset computed and listing created
2026-07-22 16:29:14 INFO Hourly dataset computed
2026-07-22 16:29:14 INFO Asked to restart from control run ...
2026-07-22 16:29:14 INFO Copying EMIS ...
2026-07-22 16:29:14 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens03.nc
2026-07-22 16:29:14 INFO Linking first END ...
2026-07-22 16:29:14 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020500_24_ENS3.nc
2026-07-22 16:29:14 INFO >> Checking links...
2026-07-22 16:29:15 INFO >> All links are good for ENS3  ...
2026-07-22 16:29:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:29:16 INFO Hourly dataset computed and listing created
2026-07-22 16:29:17 INFO Hourly dataset computed
2026-07-22 16:29:17 INFO Asked to restart from control run ...
2026-07-22 16:29:17 INFO Copying EMIS ...
2026-07-22 16:29:18 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens04.nc
2026-07-22 16:29:18 INFO Linking first END ...
2026-07-22 16:29:18 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020500_24_ENS4.nc
2026-07-22 16:29:18 INFO >> Checking links...
2026-07-22 16:29:19 INFO >> All links are good for ENS4  ...
2026-07-22 16:29:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:29:19 INFO Hourly dataset computed and listing created
2026-07-22 16:29:20 INFO Hourly dataset computed
2026-07-22 16:29:20 INFO Asked to restart from control run ...
2026-07-22 16:29:20 INFO Copying EMIS ...
2026-07-22 16:29:21 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens05.nc
2026-07-22 16:29:21 INFO Linking first END ...
2026-07-22 16:29:21 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020500_24_ENS5.nc
2026-07-22 16:29:21 INFO >> Checking links...
2026-07-22 16:29:22 INFO >> All links are good for ENS5  ...
2026-07-22 16:29:22 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:29:22 INFO Hourly dataset computed and listing created
2026-07-22 16:29:23 INFO Hourly dataset computed
2026-07-22 16:29:23 INFO Asked to restart from control run ...
2026-07-22 16:29:23 INFO Copying EMIS ...
2026-07-22 16:29:24 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens06.nc
2026-07-22 16:29:24 INFO Linking first END ...
2026-07-22 16:29:24 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020500_24_ENS6.nc
2026-07-22 16:29:24 INFO >> Checking links...
2026-07-22 16:29:25 INFO >> All links are good for ENS6  ...
2026-07-22 16:29:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:29:26 INFO Hourly dataset computed and listing created
2026-07-22 16:29:26 INFO Hourly dataset computed
2026-07-22 16:29:26 INFO Asked to restart from control run ...
2026-07-22 16:29:26 INFO Copying EMIS ...
2026-07-22 16:29:27 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens07.nc
2026-07-22 16:29:27 INFO Linking first END ...
2026-07-22 16:29:27 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020500_24_ENS7.nc
2026-07-22 16:29:27 INFO >> Checking links...
2026-07-22 16:29:28 INFO >> All links are good for ENS7  ...
2026-07-22 16:29:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:29:29 INFO Hourly dataset computed and listing created
2026-07-22 16:29:29 INFO Hourly dataset computed
2026-07-22 16:29:29 INFO Asked to restart from control run ...
2026-07-22 16:29:29 INFO Copying EMIS ...
2026-07-22 16:29:30 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens08.nc
2026-07-22 16:29:30 INFO Linking first END ...
2026-07-22 16:29:30 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020500_24_ENS8.nc
2026-07-22 16:29:30 INFO >> Checking links...
2026-07-22 16:29:31 INFO >> All links are good for ENS8  ...
2026-07-22 16:29:31 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:29:32 INFO Hourly dataset computed and listing created
2026-07-22 16:29:32 INFO Hourly dataset computed
2026-07-22 16:29:32 INFO Asked to restart from control run ...
2026-07-22 16:29:32 INFO Copying EMIS ...
2026-07-22 16:29:33 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens09.nc
2026-07-22 16:29:33 INFO Linking first END ...
2026-07-22 16:29:33 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020500_24_ENS9.nc
2026-07-22 16:29:33 INFO >> Checking links...
2026-07-22 16:29:34 INFO >> All links are good for ENS9  ...
2026-07-22 16:29:34 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:29:35 INFO Hourly dataset computed and listing created
2026-07-22 16:29:35 INFO Hourly dataset computed
2026-07-22 16:29:35 INFO Asked to restart from control run ...
2026-07-22 16:29:35 INFO Copying EMIS ...
2026-07-22 16:29:36 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens10.nc
2026-07-22 16:29:36 INFO Linking first END ...
2026-07-22 16:29:36 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020500_24_ENS10.nc
2026-07-22 16:29:36 INFO >> Checking links...
2026-07-22 16:29:37 INFO >> All links are good for ENS10  ...
2026-07-22 16:29:37 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:29:38 INFO Hourly dataset computed and listing created
2026-07-22 16:29:38 INFO Hourly dataset computed
2026-07-22 16:29:38 INFO Asked to restart from control run ...
2026-07-22 16:29:38 INFO Copying EMIS ...
2026-07-22 16:29:39 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens11.nc
2026-07-22 16:29:39 INFO Linking first END ...
2026-07-22 16:29:39 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020500_24_ENS11.nc
2026-07-22 16:29:39 INFO >> Checking links...
2026-07-22 16:29:40 INFO >> All links are good for ENS11  ...
2026-07-22 16:29:40 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:29:40 INFO Hourly dataset computed and listing created
2026-07-22 16:29:41 INFO Hourly dataset computed
2026-07-22 16:29:41 INFO Asked to restart from control run ...
2026-07-22 16:29:41 INFO Copying EMIS ...
2026-07-22 16:29:42 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens12.nc
2026-07-22 16:29:42 INFO Linking first END ...
2026-07-22 16:29:42 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020500_24_ENS12.nc
2026-07-22 16:29:42 INFO >> Checking links...
2026-07-22 16:29:43 INFO >> All links are good for ENS12  ...
2026-07-22 16:29:43 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:29:43 INFO Hourly dataset computed and listing created
2026-07-22 16:29:44 INFO Hourly dataset computed
2026-07-22 16:29:44 INFO Asked to restart from control run ...
2026-07-22 16:29:44 INFO Copying EMIS ...
2026-07-22 16:29:44 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens13.nc
2026-07-22 16:29:45 INFO Linking first END ...
2026-07-22 16:29:45 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020500_24_ENS13.nc
2026-07-22 16:29:45 INFO >> Checking links...
2026-07-22 16:29:46 INFO >> All links are good for ENS13  ...
2026-07-22 16:29:46 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:29:46 INFO Hourly dataset computed and listing created
2026-07-22 16:29:47 INFO Hourly dataset computed
2026-07-22 16:29:47 INFO Asked to restart from control run ...
2026-07-22 16:29:47 INFO Copying EMIS ...
2026-07-22 16:29:47 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens14.nc
2026-07-22 16:29:47 INFO Linking first END ...
2026-07-22 16:29:47 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020500_24_ENS14.nc
2026-07-22 16:29:47 INFO >> Checking links...
2026-07-22 16:29:48 INFO >> All links are good for ENS14  ...
2026-07-22 16:29:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:29:49 INFO Hourly dataset computed and listing created
2026-07-22 16:29:50 INFO Hourly dataset computed
2026-07-22 16:29:50 INFO Asked to restart from control run ...
2026-07-22 16:29:50 INFO Copying EMIS ...
2026-07-22 16:29:50 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens15.nc
2026-07-22 16:29:50 INFO Linking first END ...
2026-07-22 16:29:50 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020500_24_ENS15.nc
2026-07-22 16:29:50 INFO >> Checking links...
2026-07-22 16:29:51 INFO >> All links are good for ENS15  ...
2026-07-22 16:29:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:29:52 INFO Hourly dataset computed and listing created
2026-07-22 16:29:53 INFO Hourly dataset computed
2026-07-22 16:29:53 INFO ---------->>> Running CHIMERE model from 2020-02-06 00:00:00 to 2020-02-06 01:00:00
2026-07-22 16:29:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:29:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1
2026-07-22 16:29:53 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020500_24_ENS1.nc
2026-07-22 16:29:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-22 16:29:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:29:53 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-22 16:29:53 INFO Queuing job for member 1...
2026-07-22 16:29:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:29:53 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-22 16:29:56 INFO Found: ['5255871']
2026-07-22 16:30:01 INFO [TGCC-IRENE] Submitted job with ID:['5255871']
2026-07-22 16:30:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:30:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2
2026-07-22 16:30:01 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020500_24_ENS2.nc
2026-07-22 16:30:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-22 16:30:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:30:01 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-22 16:30:01 INFO Queuing job for member 2...
2026-07-22 16:30:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:30:01 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-22 16:30:04 INFO Found: ['5255874']
2026-07-22 16:30:09 INFO [TGCC-IRENE] Submitted job with ID:['5255874']
2026-07-22 16:30:09 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:30:09 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3
2026-07-22 16:30:09 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020500_24_ENS3.nc
2026-07-22 16:30:09 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-22 16:30:09 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:30:09 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-22 16:30:09 INFO Queuing job for member 3...
2026-07-22 16:30:09 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:30:09 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-22 16:30:11 INFO Found: ['5255880']
2026-07-22 16:30:16 INFO [TGCC-IRENE] Submitted job with ID:['5255880']
2026-07-22 16:30:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:30:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4
2026-07-22 16:30:16 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020500_24_ENS4.nc
2026-07-22 16:30:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-22 16:30:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:30:16 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-22 16:30:16 INFO Queuing job for member 4...
2026-07-22 16:30:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:30:16 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-22 16:30:19 INFO Found: ['5255885']
2026-07-22 16:30:24 INFO [TGCC-IRENE] Submitted job with ID:['5255885']
2026-07-22 16:30:24 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:30:24 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5
2026-07-22 16:30:24 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020500_24_ENS5.nc
2026-07-22 16:30:24 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-22 16:30:24 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:30:24 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-22 16:30:24 INFO Queuing job for member 5...
2026-07-22 16:30:24 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:30:24 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-22 16:30:26 INFO Found: ['5255888']
2026-07-22 16:30:31 INFO [TGCC-IRENE] Submitted job with ID:['5255888']
2026-07-22 16:30:31 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:30:31 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6
2026-07-22 16:30:31 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020500_24_ENS6.nc
2026-07-22 16:30:31 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-22 16:30:31 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:30:31 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-22 16:30:31 INFO Queuing job for member 6...
2026-07-22 16:30:31 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:30:31 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-22 16:30:34 INFO Found: ['5255893']
2026-07-22 16:30:39 INFO [TGCC-IRENE] Submitted job with ID:['5255893']
2026-07-22 16:30:39 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:30:39 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7
2026-07-22 16:30:39 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020500_24_ENS7.nc
2026-07-22 16:30:39 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-22 16:30:39 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:30:39 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-22 16:30:39 INFO Queuing job for member 7...
2026-07-22 16:30:39 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:30:39 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-22 16:30:39 INFO Found: ['5255899']
2026-07-22 16:30:44 INFO [TGCC-IRENE] Submitted job with ID:['5255899']
2026-07-22 16:30:44 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:30:44 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8
2026-07-22 16:30:44 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020500_24_ENS8.nc
2026-07-22 16:30:44 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-22 16:30:44 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:30:44 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-22 16:30:44 INFO Queuing job for member 8...
2026-07-22 16:30:44 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:30:44 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-22 16:30:45 INFO Found: ['5255917']
2026-07-22 16:30:50 INFO [TGCC-IRENE] Submitted job with ID:['5255917']
2026-07-22 16:30:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:30:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9
2026-07-22 16:30:50 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020500_24_ENS9.nc
2026-07-22 16:30:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-22 16:30:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:30:50 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-22 16:30:50 INFO Queuing job for member 9...
2026-07-22 16:30:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:30:50 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-22 16:30:51 INFO Found: ['5255927']
2026-07-22 16:30:56 INFO [TGCC-IRENE] Submitted job with ID:['5255927']
2026-07-22 16:30:56 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:30:56 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10
2026-07-22 16:30:56 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020500_24_ENS10.nc
2026-07-22 16:30:56 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-22 16:30:56 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:30:56 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-22 16:30:56 INFO Queuing job for member 10...
2026-07-22 16:30:56 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:30:56 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-22 16:30:57 INFO Found: ['5255930']
2026-07-22 16:31:02 INFO [TGCC-IRENE] Submitted job with ID:['5255930']
2026-07-22 16:31:02 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:31:02 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11
2026-07-22 16:31:02 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020500_24_ENS11.nc
2026-07-22 16:31:02 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-22 16:31:02 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:31:02 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-22 16:31:02 INFO Queuing job for member 11...
2026-07-22 16:31:02 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:31:02 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-22 16:31:02 INFO Found: ['5255934']
2026-07-22 16:31:07 INFO [TGCC-IRENE] Submitted job with ID:['5255934']
2026-07-22 16:31:07 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:31:07 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12
2026-07-22 16:31:07 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020500_24_ENS12.nc
2026-07-22 16:31:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-22 16:31:08 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:31:08 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-22 16:31:08 INFO Queuing job for member 12...
2026-07-22 16:31:08 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:31:08 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-22 16:31:08 INFO Found: ['5255936']
2026-07-22 16:31:13 INFO [TGCC-IRENE] Submitted job with ID:['5255936']
2026-07-22 16:31:13 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:31:13 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13
2026-07-22 16:31:13 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020500_24_ENS13.nc
2026-07-22 16:31:13 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-22 16:31:13 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:31:13 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-22 16:31:13 INFO Queuing job for member 13...
2026-07-22 16:31:13 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:31:13 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-22 16:31:14 INFO Found: ['5255941']
2026-07-22 16:31:19 INFO [TGCC-IRENE] Submitted job with ID:['5255941']
2026-07-22 16:31:19 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:31:19 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14
2026-07-22 16:31:19 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020500_24_ENS14.nc
2026-07-22 16:31:19 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-22 16:31:19 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:31:19 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-22 16:31:19 INFO Queuing job for member 14...
2026-07-22 16:31:19 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:31:19 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-22 16:31:20 INFO Found: ['5255944']
2026-07-22 16:31:25 INFO [TGCC-IRENE] Submitted job with ID:['5255944']
2026-07-22 16:31:25 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:31:25 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15
2026-07-22 16:31:25 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020500_24_ENS15.nc
2026-07-22 16:31:25 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-22 16:31:25 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:31:25 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-22 16:31:25 INFO Queuing job for member 15...
2026-07-22 16:31:25 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:31:25 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-22 16:31:26 INFO Found: ['5255946']
2026-07-22 16:31:31 INFO [TGCC-IRENE] Submitted job with ID:['5255946']
2026-07-22 16:31:31 INFO Checking job status ...
2026-07-22 16:31:31 INFO None 5255871: status RUNNING/PENDING
2026-07-22 16:31:31 INFO None 5255874: status RUNNING/PENDING
2026-07-22 16:31:31 INFO None 5255880: status RUNNING/PENDING
2026-07-22 16:31:31 INFO None 5255885: status RUNNING/PENDING
2026-07-22 16:31:31 INFO None 5255888: status RUNNING/PENDING
2026-07-22 16:31:31 INFO None 5255893: status RUNNING/PENDING
2026-07-22 16:31:31 INFO None 5255899: status RUNNING/PENDING
2026-07-22 16:31:31 INFO None 5255917: status RUNNING/PENDING
2026-07-22 16:31:31 INFO None 5255927: status RUNNING/PENDING
2026-07-22 16:31:31 INFO None 5255930: status RUNNING/PENDING
2026-07-22 16:31:33 INFO None 5255934: status RUNNING/PENDING
2026-07-22 16:31:33 INFO None 5255936: status RUNNING/PENDING
2026-07-22 16:31:33 INFO None 5255941: status RUNNING/PENDING
2026-07-22 16:31:33 INFO None 5255944: status RUNNING/PENDING
2026-07-22 16:31:33 INFO None 5255946: status RUNNING/PENDING
2026-07-22 16:31:33 INFO Jobs still running: ['5255871', '5255874', '5255880', '5255885', '5255888', '5255893', '5255899', '5255917', '5255927', '5255930', '5255934', '5255936', '5255941', '5255944', '5255946']. Waiting...
2026-07-22 16:31:48 INFO None 5255871: status RUNNING/PENDING
2026-07-22 16:31:48 INFO None 5255874: status RUNNING/PENDING
2026-07-22 16:31:48 INFO None 5255880: status RUNNING/PENDING
2026-07-22 16:31:48 INFO None 5255885: status RUNNING/PENDING
2026-07-22 16:31:48 INFO None 5255888: status RUNNING/PENDING
2026-07-22 16:31:48 INFO None 5255893: status RUNNING/PENDING
2026-07-22 16:31:48 INFO None 5255899: status RUNNING/PENDING
2026-07-22 16:31:48 INFO None 5255917: status RUNNING/PENDING
2026-07-22 16:31:48 INFO None 5255927: status RUNNING/PENDING
2026-07-22 16:31:48 INFO None 5255930: status RUNNING/PENDING
2026-07-22 16:31:48 INFO None 5255934: status RUNNING/PENDING
2026-07-22 16:31:51 INFO None 5255936: status RUNNING/PENDING
2026-07-22 16:31:51 INFO None 5255941: status RUNNING/PENDING
2026-07-22 16:31:51 INFO None 5255944: status RUNNING/PENDING
2026-07-22 16:31:51 INFO None 5255946: status RUNNING/PENDING
2026-07-22 16:31:51 INFO Jobs still running: ['5255871', '5255874', '5255880', '5255885', '5255888', '5255893', '5255899', '5255917', '5255927', '5255930', '5255934', '5255936', '5255941', '5255944', '5255946']. Waiting...
2026-07-22 16:32:06 INFO None 5255871: status RUNNING/PENDING
2026-07-22 16:32:06 INFO None 5255874: status RUNNING/PENDING
2026-07-22 16:32:06 INFO None 5255880: status RUNNING/PENDING
2026-07-22 16:32:06 INFO None 5255885: status RUNNING/PENDING
2026-07-22 16:32:06 INFO None 5255888: status RUNNING/PENDING
2026-07-22 16:32:06 INFO None 5255893: status RUNNING/PENDING
2026-07-22 16:32:06 INFO None 5255899: status RUNNING/PENDING
2026-07-22 16:32:06 INFO None 5255917: status RUNNING/PENDING
2026-07-22 16:32:06 INFO None 5255927: status RUNNING/PENDING
2026-07-22 16:32:06 INFO None 5255930: status RUNNING/PENDING
2026-07-22 16:32:06 INFO None 5255934: status RUNNING/PENDING
2026-07-22 16:32:06 INFO None 5255936: status RUNNING/PENDING
2026-07-22 16:32:06 INFO None 5255941: status RUNNING/PENDING
2026-07-22 16:32:08 INFO None 5255944: status RUNNING/PENDING
2026-07-22 16:32:08 INFO None 5255946: status RUNNING/PENDING
2026-07-22 16:32:08 INFO Jobs still running: ['5255871', '5255874', '5255880', '5255885', '5255888', '5255893', '5255899', '5255917', '5255927', '5255930', '5255934', '5255936', '5255941', '5255944', '5255946']. Waiting...
2026-07-22 16:32:23 INFO None 5255871: status RUNNING/PENDING
2026-07-22 16:32:23 INFO None 5255874: status RUNNING/PENDING
2026-07-22 16:32:23 INFO None 5255880: status RUNNING/PENDING
2026-07-22 16:32:23 INFO None 5255885: status RUNNING/PENDING
2026-07-22 16:32:23 INFO None 5255888: status RUNNING/PENDING
2026-07-22 16:32:23 INFO None 5255893: status RUNNING/PENDING
2026-07-22 16:32:23 INFO None 5255899: status RUNNING/PENDING
2026-07-22 16:32:23 INFO None 5255917: status RUNNING/PENDING
2026-07-22 16:32:23 INFO None 5255927: status RUNNING/PENDING
2026-07-22 16:32:23 INFO None 5255930: status RUNNING/PENDING
2026-07-22 16:32:23 INFO None 5255934: status RUNNING/PENDING
2026-07-22 16:32:23 INFO None 5255936: status RUNNING/PENDING
2026-07-22 16:32:23 INFO None 5255941: status RUNNING/PENDING
2026-07-22 16:32:24 INFO None 5255944: status RUNNING/PENDING
2026-07-22 16:32:24 INFO None 5255946: status RUNNING/PENDING
2026-07-22 16:32:24 INFO Jobs still running: ['5255871', '5255874', '5255880', '5255885', '5255888', '5255893', '5255899', '5255917', '5255927', '5255930', '5255934', '5255936', '5255941', '5255944', '5255946']. Waiting...
2026-07-22 16:32:39 INFO None 5255871: status RUNNING/PENDING
2026-07-22 16:32:39 INFO None 5255874: status RUNNING/PENDING
2026-07-22 16:32:39 INFO None 5255880: status RUNNING/PENDING
2026-07-22 16:32:39 INFO None 5255885: status RUNNING/PENDING
2026-07-22 16:32:39 INFO None 5255888: status RUNNING/PENDING
2026-07-22 16:32:39 INFO None 5255893: status RUNNING/PENDING
2026-07-22 16:32:39 INFO None 5255899: status RUNNING/PENDING
2026-07-22 16:32:39 INFO None 5255917: status RUNNING/PENDING
2026-07-22 16:32:39 INFO None 5255927: status RUNNING/PENDING
2026-07-22 16:32:39 INFO None 5255930: status RUNNING/PENDING
2026-07-22 16:32:39 INFO None 5255934: status RUNNING/PENDING
2026-07-22 16:32:39 INFO None 5255936: status RUNNING/PENDING
2026-07-22 16:32:39 INFO None 5255941: status RUNNING/PENDING
2026-07-22 16:32:39 INFO None 5255944: status RUNNING/PENDING
2026-07-22 16:32:39 INFO None 5255946: status RUNNING/PENDING
2026-07-22 16:32:39 INFO Jobs still running: ['5255871', '5255874', '5255880', '5255885', '5255888', '5255893', '5255899', '5255917', '5255927', '5255930', '5255934', '5255936', '5255941', '5255944', '5255946']. Waiting...
2026-07-22 16:32:54 INFO None 5255871: status RUNNING/PENDING
2026-07-22 16:32:54 INFO None 5255874: status RUNNING/PENDING
2026-07-22 16:32:54 INFO None 5255880: status RUNNING/PENDING
2026-07-22 16:32:54 INFO None 5255885: status RUNNING/PENDING
2026-07-22 16:32:54 INFO None 5255888: status RUNNING/PENDING
2026-07-22 16:32:54 INFO None 5255893: status RUNNING/PENDING
2026-07-22 16:32:54 INFO None 5255899: status RUNNING/PENDING
2026-07-22 16:32:54 INFO None 5255917: status RUNNING/PENDING
2026-07-22 16:32:54 INFO None 5255927: status RUNNING/PENDING
2026-07-22 16:32:54 INFO None 5255930: status RUNNING/PENDING
2026-07-22 16:32:54 INFO None 5255934: status RUNNING/PENDING
2026-07-22 16:32:54 INFO None 5255936: status RUNNING/PENDING
2026-07-22 16:32:54 INFO None 5255941: status RUNNING/PENDING
2026-07-22 16:32:54 INFO None 5255944: status RUNNING/PENDING
2026-07-22 16:32:54 INFO None 5255946: status RUNNING/PENDING
2026-07-22 16:32:54 INFO Jobs still running: ['5255871', '5255874', '5255880', '5255885', '5255888', '5255893', '5255899', '5255917', '5255927', '5255930', '5255934', '5255936', '5255941', '5255944', '5255946']. Waiting...
2026-07-22 16:33:12 INFO None 5255871: status RUNNING/PENDING
2026-07-22 16:33:12 INFO None 5255874: status RUNNING/PENDING
2026-07-22 16:33:12 INFO None 5255880: status RUNNING/PENDING
2026-07-22 16:33:12 INFO None 5255885: status RUNNING/PENDING
2026-07-22 16:33:12 INFO None 5255888: status RUNNING/PENDING
2026-07-22 16:33:12 INFO None 5255893: status RUNNING/PENDING
2026-07-22 16:33:12 INFO None 5255899: status RUNNING/PENDING
2026-07-22 16:33:12 INFO None 5255917: status RUNNING/PENDING
2026-07-22 16:33:12 INFO None 5255927: status RUNNING/PENDING
2026-07-22 16:33:12 INFO None 5255930: status RUNNING/PENDING
2026-07-22 16:33:12 INFO None 5255934: status RUNNING/PENDING
2026-07-22 16:33:12 INFO None 5255936: status RUNNING/PENDING
2026-07-22 16:33:12 INFO None 5255941: status RUNNING/PENDING
2026-07-22 16:33:12 INFO None 5255944: status RUNNING/PENDING
2026-07-22 16:33:12 INFO None 5255946: status RUNNING/PENDING
2026-07-22 16:33:12 INFO Jobs still running: ['5255871', '5255874', '5255880', '5255885', '5255888', '5255893', '5255899', '5255917', '5255927', '5255930', '5255934', '5255936', '5255941', '5255944', '5255946']. Waiting...
2026-07-22 16:33:27 INFO None 5255871: status RUNNING/PENDING
2026-07-22 16:33:27 INFO None 5255874: status RUNNING/PENDING
2026-07-22 16:33:27 INFO None 5255880: status RUNNING/PENDING
2026-07-22 16:33:27 INFO None 5255885: status RUNNING/PENDING
2026-07-22 16:33:27 INFO None 5255888: status RUNNING/PENDING
2026-07-22 16:33:27 INFO None 5255893: status RUNNING/PENDING
2026-07-22 16:33:27 INFO None 5255899: status RUNNING/PENDING
2026-07-22 16:33:27 INFO None 5255917: status RUNNING/PENDING
2026-07-22 16:33:27 INFO None 5255927: status RUNNING/PENDING
2026-07-22 16:33:29 INFO None 5255930: status RUNNING/PENDING
2026-07-22 16:33:29 INFO None 5255934: status RUNNING/PENDING
2026-07-22 16:33:29 INFO None 5255936: status RUNNING/PENDING
2026-07-22 16:33:29 INFO None 5255941: status RUNNING/PENDING
2026-07-22 16:33:29 INFO None 5255944: status RUNNING/PENDING
2026-07-22 16:33:29 INFO None 5255946: status RUNNING/PENDING
2026-07-22 16:33:29 INFO Jobs still running: ['5255871', '5255874', '5255880', '5255885', '5255888', '5255893', '5255899', '5255917', '5255927', '5255930', '5255934', '5255936', '5255941', '5255944', '5255946']. Waiting...
2026-07-22 16:33:44 INFO None 5255871: status RUNNING/PENDING
2026-07-22 16:33:44 INFO None 5255874: status RUNNING/PENDING
2026-07-22 16:33:44 INFO None 5255880: status RUNNING/PENDING
2026-07-22 16:33:44 INFO None 5255885: status RUNNING/PENDING
2026-07-22 16:33:44 INFO None 5255888: status RUNNING/PENDING
2026-07-22 16:33:44 INFO None 5255893: status RUNNING/PENDING
2026-07-22 16:33:45 INFO None 5255899: status RUNNING/PENDING
2026-07-22 16:33:45 INFO None 5255917: status RUNNING/PENDING
2026-07-22 16:33:45 INFO None 5255927: status RUNNING/PENDING
2026-07-22 16:33:45 INFO None 5255930: status RUNNING/PENDING
2026-07-22 16:33:45 INFO None 5255934: status RUNNING/PENDING
2026-07-22 16:33:45 INFO None 5255936: status RUNNING/PENDING
2026-07-22 16:33:45 INFO None 5255941: status RUNNING/PENDING
2026-07-22 16:33:45 INFO None 5255944: status RUNNING/PENDING
2026-07-22 16:33:45 INFO None 5255946: status RUNNING/PENDING
2026-07-22 16:33:45 INFO Jobs still running: ['5255871', '5255874', '5255880', '5255885', '5255888', '5255893', '5255899', '5255917', '5255927', '5255930', '5255934', '5255936', '5255941', '5255944', '5255946']. Waiting...
2026-07-22 16:34:00 INFO None 5255871: status RUNNING/PENDING
2026-07-22 16:34:00 INFO None 5255874: status RUNNING/PENDING
2026-07-22 16:34:00 INFO None 5255880: status RUNNING/PENDING
2026-07-22 16:34:00 INFO None 5255885: status RUNNING/PENDING
2026-07-22 16:34:00 INFO None 5255888: status RUNNING/PENDING
2026-07-22 16:34:00 INFO None 5255893: status RUNNING/PENDING
2026-07-22 16:34:00 INFO None 5255899: status RUNNING/PENDING
2026-07-22 16:34:00 INFO None 5255917: status RUNNING/PENDING
2026-07-22 16:34:00 INFO None 5255927: status RUNNING/PENDING
2026-07-22 16:34:00 INFO None 5255930: status RUNNING/PENDING
2026-07-22 16:34:00 INFO None 5255934: status RUNNING/PENDING
2026-07-22 16:34:00 INFO None 5255936: status RUNNING/PENDING
2026-07-22 16:34:00 INFO None 5255941: status RUNNING/PENDING
2026-07-22 16:34:00 INFO None 5255944: status RUNNING/PENDING
2026-07-22 16:34:00 INFO None 5255946: status RUNNING/PENDING
2026-07-22 16:34:00 INFO Jobs still running: ['5255871', '5255874', '5255880', '5255885', '5255888', '5255893', '5255899', '5255917', '5255927', '5255930', '5255934', '5255936', '5255941', '5255944', '5255946']. Waiting...
2026-07-22 16:34:15 INFO None 5255871: status RUNNING/PENDING
2026-07-22 16:34:15 INFO None 5255874: status RUNNING/PENDING
2026-07-22 16:34:15 INFO None 5255880: status RUNNING/PENDING
2026-07-22 16:34:15 INFO None 5255885: status RUNNING/PENDING
2026-07-22 16:34:15 INFO None 5255888: status RUNNING/PENDING
2026-07-22 16:34:15 INFO None 5255893: status RUNNING/PENDING
2026-07-22 16:34:15 INFO None 5255899: status FINISHED
2026-07-22 16:34:15 INFO None 5255917: status RUNNING/PENDING
2026-07-22 16:34:15 INFO None 5255927: status RUNNING/PENDING
2026-07-22 16:34:15 INFO None 5255930: status RUNNING/PENDING
2026-07-22 16:34:15 INFO None 5255934: status RUNNING/PENDING
2026-07-22 16:34:15 INFO None 5255936: status RUNNING/PENDING
2026-07-22 16:34:15 INFO None 5255941: status RUNNING/PENDING
2026-07-22 16:34:15 INFO None 5255944: status RUNNING/PENDING
2026-07-22 16:34:15 INFO None 5255946: status RUNNING/PENDING
2026-07-22 16:34:15 INFO Jobs still running: ['5255871', '5255874', '5255880', '5255885', '5255888', '5255893', '5255917', '5255927', '5255930', '5255934', '5255936', '5255941', '5255944', '5255946']. Waiting...
2026-07-22 16:34:30 INFO None 5255871: status RUNNING/PENDING
2026-07-22 16:34:30 INFO None 5255874: status RUNNING/PENDING
2026-07-22 16:34:30 INFO None 5255880: status RUNNING/PENDING
2026-07-22 16:34:30 INFO None 5255885: status RUNNING/PENDING
2026-07-22 16:34:30 INFO None 5255888: status FINISHED
2026-07-22 16:34:30 INFO None 5255893: status FINISHED
2026-07-22 16:34:31 INFO None 5255899: status FINISHED
2026-07-22 16:34:31 INFO None 5255917: status FINISHED
2026-07-22 16:34:31 INFO None 5255927: status FINISHED
2026-07-22 16:34:31 INFO None 5255930: status RUNNING/PENDING
2026-07-22 16:34:31 INFO None 5255934: status RUNNING/PENDING
2026-07-22 16:34:31 INFO None 5255936: status RUNNING/PENDING
2026-07-22 16:34:31 INFO None 5255941: status RUNNING/PENDING
2026-07-22 16:34:31 INFO None 5255944: status RUNNING/PENDING
2026-07-22 16:34:31 INFO None 5255946: status RUNNING/PENDING
2026-07-22 16:34:31 INFO Jobs still running: ['5255871', '5255874', '5255880', '5255885', '5255930', '5255934', '5255936', '5255941', '5255944', '5255946']. Waiting...
2026-07-22 16:34:46 INFO None 5255871: status RUNNING/PENDING
2026-07-22 16:34:46 INFO None 5255874: status RUNNING/PENDING
2026-07-22 16:34:46 INFO None 5255880: status RUNNING/PENDING
2026-07-22 16:34:46 INFO None 5255885: status RUNNING/PENDING
2026-07-22 16:34:46 INFO None 5255888: status FINISHED
2026-07-22 16:34:46 INFO None 5255893: status FINISHED
2026-07-22 16:34:46 INFO None 5255899: status FINISHED
2026-07-22 16:34:46 INFO None 5255917: status FINISHED
2026-07-22 16:34:46 INFO None 5255927: status FINISHED
2026-07-22 16:34:46 INFO None 5255930: status RUNNING/PENDING
2026-07-22 16:34:46 INFO None 5255934: status RUNNING/PENDING
2026-07-22 16:34:46 INFO None 5255936: status RUNNING/PENDING
2026-07-22 16:34:46 INFO None 5255941: status RUNNING/PENDING
2026-07-22 16:34:46 INFO None 5255944: status RUNNING/PENDING
2026-07-22 16:34:46 INFO None 5255946: status RUNNING/PENDING
2026-07-22 16:34:46 INFO Jobs still running: ['5255871', '5255874', '5255880', '5255885', '5255930', '5255934', '5255936', '5255941', '5255944', '5255946']. Waiting...
2026-07-22 16:35:01 INFO None 5255871: status FINISHED
2026-07-22 16:35:01 INFO None 5255874: status RUNNING/PENDING
2026-07-22 16:35:01 INFO None 5255880: status RUNNING/PENDING
2026-07-22 16:35:02 INFO None 5255885: status RUNNING/PENDING
2026-07-22 16:35:02 INFO None 5255888: status FINISHED
2026-07-22 16:35:02 INFO None 5255893: status FINISHED
2026-07-22 16:35:02 INFO None 5255899: status FINISHED
2026-07-22 16:35:02 INFO None 5255917: status FINISHED
2026-07-22 16:35:02 INFO None 5255927: status FINISHED
2026-07-22 16:35:02 INFO None 5255930: status FINISHED
2026-07-22 16:35:02 INFO None 5255934: status FINISHED
2026-07-22 16:35:02 INFO None 5255936: status FINISHED
2026-07-22 16:35:04 INFO None 5255941: status RUNNING/PENDING
2026-07-22 16:35:04 INFO None 5255944: status RUNNING/PENDING
2026-07-22 16:35:04 INFO None 5255946: status RUNNING/PENDING
2026-07-22 16:35:04 INFO Jobs still running: ['5255874', '5255880', '5255885', '5255941', '5255944', '5255946']. Waiting...
2026-07-22 16:35:19 INFO None 5255871: status FINISHED
2026-07-22 16:35:19 INFO None 5255874: status FINISHED
2026-07-22 16:35:19 INFO None 5255880: status RUNNING/PENDING
2026-07-22 16:35:19 INFO None 5255885: status RUNNING/PENDING
2026-07-22 16:35:19 INFO None 5255888: status FINISHED
2026-07-22 16:35:19 INFO None 5255893: status FINISHED
2026-07-22 16:35:19 INFO None 5255899: status FINISHED
2026-07-22 16:35:19 INFO None 5255917: status FINISHED
2026-07-22 16:35:19 INFO None 5255927: status FINISHED
2026-07-22 16:35:19 INFO None 5255930: status FINISHED
2026-07-22 16:35:19 INFO None 5255934: status FINISHED
2026-07-22 16:35:19 INFO None 5255936: status FINISHED
2026-07-22 16:35:19 INFO None 5255941: status RUNNING/PENDING
2026-07-22 16:35:19 INFO None 5255944: status RUNNING/PENDING
2026-07-22 16:35:19 INFO None 5255946: status RUNNING/PENDING
2026-07-22 16:35:19 INFO Jobs still running: ['5255880', '5255885', '5255941', '5255944', '5255946']. Waiting...
2026-07-22 16:35:34 INFO None 5255871: status FINISHED
2026-07-22 16:35:34 INFO None 5255874: status FINISHED
2026-07-22 16:35:34 INFO None 5255880: status FINISHED
2026-07-22 16:35:34 INFO None 5255885: status FINISHED
2026-07-22 16:35:34 INFO None 5255888: status FINISHED
2026-07-22 16:35:34 INFO None 5255893: status FINISHED
2026-07-22 16:35:34 INFO None 5255899: status FINISHED
2026-07-22 16:35:34 INFO None 5255917: status FINISHED
2026-07-22 16:35:34 INFO None 5255927: status FINISHED
2026-07-22 16:35:34 INFO None 5255930: status FINISHED
2026-07-22 16:35:34 INFO None 5255934: status FINISHED
2026-07-22 16:35:34 INFO None 5255936: status FINISHED
2026-07-22 16:35:34 INFO None 5255941: status FINISHED
2026-07-22 16:35:34 INFO None 5255944: status FINISHED
2026-07-22 16:35:34 INFO None 5255946: status RUNNING/PENDING
2026-07-22 16:35:34 INFO Jobs still running: ['5255946']. Waiting...
2026-07-22 16:35:49 INFO None 5255871: status FINISHED
2026-07-22 16:35:49 INFO None 5255874: status FINISHED
2026-07-22 16:35:50 INFO None 5255880: status FINISHED
2026-07-22 16:35:50 INFO None 5255885: status FINISHED
2026-07-22 16:35:50 INFO None 5255888: status FINISHED
2026-07-22 16:35:50 INFO None 5255893: status FINISHED
2026-07-22 16:35:50 INFO None 5255899: status FINISHED
2026-07-22 16:35:50 INFO None 5255917: status FINISHED
2026-07-22 16:35:50 INFO None 5255927: status FINISHED
2026-07-22 16:35:50 INFO None 5255930: status FINISHED
2026-07-22 16:35:50 INFO None 5255934: status FINISHED
2026-07-22 16:35:50 INFO None 5255936: status FINISHED
2026-07-22 16:35:50 INFO None 5255941: status FINISHED
2026-07-22 16:35:50 INFO None 5255944: status FINISHED
2026-07-22 16:35:50 INFO None 5255946: status RUNNING/PENDING
2026-07-22 16:35:50 INFO Jobs still running: ['5255946']. Waiting...
2026-07-22 16:36:05 INFO None 5255871: status FINISHED
2026-07-22 16:36:05 INFO None 5255874: status FINISHED
2026-07-22 16:36:05 INFO None 5255880: status FINISHED
2026-07-22 16:36:05 INFO None 5255885: status FINISHED
2026-07-22 16:36:05 INFO None 5255888: status FINISHED
2026-07-22 16:36:05 INFO None 5255893: status FINISHED
2026-07-22 16:36:05 INFO None 5255899: status FINISHED
2026-07-22 16:36:05 INFO None 5255917: status FINISHED
2026-07-22 16:36:05 INFO None 5255927: status FINISHED
2026-07-22 16:36:05 INFO None 5255930: status FINISHED
2026-07-22 16:36:05 INFO None 5255934: status FINISHED
2026-07-22 16:36:05 INFO None 5255936: status FINISHED
2026-07-22 16:36:05 INFO None 5255941: status FINISHED
2026-07-22 16:36:05 INFO None 5255944: status FINISHED
2026-07-22 16:36:05 INFO None 5255946: status RUNNING/PENDING
2026-07-22 16:36:05 INFO Jobs still running: ['5255946']. Waiting...
2026-07-22 16:36:20 INFO None 5255871: status FINISHED
2026-07-22 16:36:20 INFO None 5255874: status FINISHED
2026-07-22 16:36:20 INFO None 5255880: status FINISHED
2026-07-22 16:36:20 INFO None 5255885: status FINISHED
2026-07-22 16:36:22 INFO None 5255888: status FINISHED
2026-07-22 16:36:22 INFO None 5255893: status FINISHED
2026-07-22 16:36:22 INFO None 5255899: status FINISHED
2026-07-22 16:36:22 INFO None 5255917: status FINISHED
2026-07-22 16:36:22 INFO None 5255927: status FINISHED
2026-07-22 16:36:22 INFO None 5255930: status FINISHED
2026-07-22 16:36:23 INFO None 5255934: status FINISHED
2026-07-22 16:36:23 INFO None 5255936: status FINISHED
2026-07-22 16:36:23 INFO None 5255941: status FINISHED
2026-07-22 16:36:23 INFO None 5255944: status FINISHED
2026-07-22 16:36:23 INFO None 5255946: status RUNNING/PENDING
2026-07-22 16:36:23 INFO Jobs still running: ['5255946']. Waiting...
2026-07-22 16:36:38 INFO None 5255871: status FINISHED
2026-07-22 16:36:38 INFO None 5255874: status FINISHED
2026-07-22 16:36:38 INFO None 5255880: status FINISHED
2026-07-22 16:36:38 INFO None 5255885: status FINISHED
2026-07-22 16:36:38 INFO None 5255888: status FINISHED
2026-07-22 16:36:38 INFO None 5255893: status FINISHED
2026-07-22 16:36:38 INFO None 5255899: status FINISHED
2026-07-22 16:36:38 INFO None 5255917: status FINISHED
2026-07-22 16:36:38 INFO None 5255927: status FINISHED
2026-07-22 16:36:38 INFO None 5255930: status FINISHED
2026-07-22 16:36:38 INFO None 5255934: status FINISHED
2026-07-22 16:36:38 INFO None 5255936: status FINISHED
2026-07-22 16:36:38 INFO None 5255941: status FINISHED
2026-07-22 16:36:38 INFO None 5255944: status FINISHED
2026-07-22 16:36:40 INFO None 5255946: status RUNNING/PENDING
2026-07-22 16:36:40 INFO Jobs still running: ['5255946']. Waiting...
2026-07-22 16:36:55 INFO None 5255871: status FINISHED
2026-07-22 16:36:55 INFO None 5255874: status FINISHED
2026-07-22 16:36:55 INFO None 5255880: status FINISHED
2026-07-22 16:36:55 INFO None 5255885: status FINISHED
2026-07-22 16:36:55 INFO None 5255888: status FINISHED
2026-07-22 16:36:55 INFO None 5255893: status FINISHED
2026-07-22 16:36:55 INFO None 5255899: status FINISHED
2026-07-22 16:36:55 INFO None 5255917: status FINISHED
2026-07-22 16:36:55 INFO None 5255927: status FINISHED
2026-07-22 16:36:55 INFO None 5255930: status FINISHED
2026-07-22 16:36:55 INFO None 5255934: status FINISHED
2026-07-22 16:36:55 INFO None 5255936: status FINISHED
2026-07-22 16:36:55 INFO None 5255941: status FINISHED
2026-07-22 16:36:55 INFO None 5255944: status FINISHED
2026-07-22 16:36:55 INFO None 5255946: status RUNNING/PENDING
2026-07-22 16:36:55 INFO Jobs still running: ['5255946']. Waiting...
2026-07-22 16:37:10 INFO None 5255871: status FINISHED
2026-07-22 16:37:10 INFO None 5255874: status FINISHED
2026-07-22 16:37:10 INFO None 5255880: status FINISHED
2026-07-22 16:37:10 INFO None 5255885: status FINISHED
2026-07-22 16:37:10 INFO None 5255888: status FINISHED
2026-07-22 16:37:10 INFO None 5255893: status FINISHED
2026-07-22 16:37:10 INFO None 5255899: status FINISHED
2026-07-22 16:37:10 INFO None 5255917: status FINISHED
2026-07-22 16:37:10 INFO None 5255927: status FINISHED
2026-07-22 16:37:10 INFO None 5255930: status FINISHED
2026-07-22 16:37:10 INFO None 5255934: status FINISHED
2026-07-22 16:37:10 INFO None 5255936: status FINISHED
2026-07-22 16:37:10 INFO None 5255941: status FINISHED
2026-07-22 16:37:11 INFO None 5255944: status FINISHED
2026-07-22 16:37:11 INFO None 5255946: status FINISHED
2026-07-22 16:37:11 INFO Jobs ['5255871', '5255874', '5255880', '5255885', '5255888', '5255893', '5255899', '5255917', '5255927', '5255930', '5255934', '5255936', '5255941', '5255944', '5255946'] have finished
2026-07-22 16:37:11 INFO Checking restart files were created ...
2026-07-22 16:37:11 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020600_1_ENS1.nc(668832435 bytes)
2026-07-22 16:37:11 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020600_1_ENS2.nc(668832435 bytes)
2026-07-22 16:37:11 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020600_1_ENS3.nc(668832435 bytes)
2026-07-22 16:37:11 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020600_1_ENS4.nc(668832435 bytes)
2026-07-22 16:37:11 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020600_1_ENS5.nc(668832435 bytes)
2026-07-22 16:37:11 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020600_1_ENS6.nc(668832435 bytes)
2026-07-22 16:37:11 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020600_1_ENS7.nc(668832435 bytes)
2026-07-22 16:37:11 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020600_1_ENS8.nc(668832435 bytes)
2026-07-22 16:37:11 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020600_1_ENS9.nc(668832435 bytes)
2026-07-22 16:37:11 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020600_1_ENS10.nc(668832435 bytes)
2026-07-22 16:37:11 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020600_1_ENS11.nc(668832435 bytes)
2026-07-22 16:37:11 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020600_1_ENS12.nc(668832435 bytes)
2026-07-22 16:37:11 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020600_1_ENS13.nc(668832435 bytes)
2026-07-22 16:37:11 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020600_1_ENS14.nc(668832435 bytes)
2026-07-22 16:37:11 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020600_1_ENS15.nc(668832435 bytes)
2026-07-22 16:37:11 INFO  Run_model() completed successfully.
2026-07-22 16:37:11 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 00:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 16:37:11 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 01:00:00 days=153072 seconds=3600
2026-07-22 16:37:11 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-22 16:37:11 INFO [TIME] increment current_time 2020-02-06 00:00:00 -> 2020-02-06 01:00:00
2026-07-22 16:37:11 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 16:37:11 INFO ---------->>> Running process_satellite_data()
2026-07-22 16:37:11 INFO [DART] No satellite data found, skipping assimilation
2026-07-22 16:37:11 INFO after_assimilation() skipped
2026-07-22 16:37:11 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-22 16:37:11 INFO [TIME] step_end current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 16:37:11 INFO [TIME] step_start current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 16:37:11 INFO [TIME] window start=2020-02-06 01:00:00 end=2020-02-06 09:00:00 run_hours=8 has_assimilation=True
2026-07-22 16:37:11 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:37:12 INFO Hourly dataset computed and listing created
2026-07-22 16:37:25 INFO Hourly dataset computed
2026-07-22 16:37:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:37:26 INFO Hourly dataset computed and listing created
2026-07-22 16:37:38 INFO Hourly dataset computed
2026-07-22 16:37:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:37:40 INFO Hourly dataset computed and listing created
2026-07-22 16:37:50 INFO Hourly dataset computed
2026-07-22 16:37:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:37:52 INFO Hourly dataset computed and listing created
2026-07-22 16:38:04 INFO Hourly dataset computed
2026-07-22 16:38:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:38:06 INFO Hourly dataset computed and listing created
2026-07-22 16:38:17 INFO Hourly dataset computed
2026-07-22 16:38:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:38:19 INFO Hourly dataset computed and listing created
2026-07-22 16:38:31 INFO Hourly dataset computed
2026-07-22 16:38:31 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:38:33 INFO Hourly dataset computed and listing created
2026-07-22 16:38:44 INFO Hourly dataset computed
2026-07-22 16:38:44 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:38:46 INFO Hourly dataset computed and listing created
2026-07-22 16:38:58 INFO Hourly dataset computed
2026-07-22 16:38:58 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:39:01 INFO Hourly dataset computed and listing created
2026-07-22 16:39:11 INFO Hourly dataset computed
2026-07-22 16:39:11 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:39:13 INFO Hourly dataset computed and listing created
2026-07-22 16:39:27 INFO Hourly dataset computed
2026-07-22 16:39:27 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:39:28 INFO Hourly dataset computed and listing created
2026-07-22 16:39:41 INFO Hourly dataset computed
2026-07-22 16:39:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:39:42 INFO Hourly dataset computed and listing created
2026-07-22 16:39:55 INFO Hourly dataset computed
2026-07-22 16:39:56 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:39:57 INFO Hourly dataset computed and listing created
2026-07-22 16:40:08 INFO Hourly dataset computed
2026-07-22 16:40:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:40:10 INFO Hourly dataset computed and listing created
2026-07-22 16:40:36 INFO Hourly dataset computed
2026-07-22 16:40:36 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:40:38 INFO Hourly dataset computed and listing created
2026-07-22 16:40:56 INFO Hourly dataset computed
2026-07-22 16:40:56 INFO ---------->>> Running CHIMERE model from 2020-02-06 01:00:00 to 2020-02-06 09:00:00
2026-07-22 16:40:56 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:40:56 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1
2026-07-22 16:40:56 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020600_1_ENS1.nc
2026-07-22 16:40:56 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-22 16:40:56 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:40:56 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-22 16:40:56 INFO Queuing job for member 1...
2026-07-22 16:40:56 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:40:56 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-22 16:40:57 INFO Found: ['5256224']
2026-07-22 16:41:02 INFO [TGCC-IRENE] Submitted job with ID:['5256224']
2026-07-22 16:41:02 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:41:02 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2
2026-07-22 16:41:02 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020600_1_ENS2.nc
2026-07-22 16:41:02 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-22 16:41:02 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:41:02 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-22 16:41:02 INFO Queuing job for member 2...
2026-07-22 16:41:02 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:41:02 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-22 16:41:02 INFO Found: ['5256236']
2026-07-22 16:41:07 INFO [TGCC-IRENE] Submitted job with ID:['5256236']
2026-07-22 16:41:07 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:41:07 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3
2026-07-22 16:41:07 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020600_1_ENS3.nc
2026-07-22 16:41:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-22 16:41:07 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:41:07 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-22 16:41:07 INFO Queuing job for member 3...
2026-07-22 16:41:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:41:07 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-22 16:41:08 INFO Found: ['5256246']
2026-07-22 16:41:13 INFO [TGCC-IRENE] Submitted job with ID:['5256246']
2026-07-22 16:41:13 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:41:13 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4
2026-07-22 16:41:13 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020600_1_ENS4.nc
2026-07-22 16:41:13 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-22 16:41:13 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:41:13 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-22 16:41:13 INFO Queuing job for member 4...
2026-07-22 16:41:13 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:41:13 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-22 16:41:14 INFO Found: ['5256255']
2026-07-22 16:41:19 INFO [TGCC-IRENE] Submitted job with ID:['5256255']
2026-07-22 16:41:19 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:41:19 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5
2026-07-22 16:41:19 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020600_1_ENS5.nc
2026-07-22 16:41:19 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-22 16:41:19 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:41:19 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-22 16:41:19 INFO Queuing job for member 5...
2026-07-22 16:41:19 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:41:19 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-22 16:41:20 INFO Found: ['5256259']
2026-07-22 16:41:25 INFO [TGCC-IRENE] Submitted job with ID:['5256259']
2026-07-22 16:41:25 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:41:25 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6
2026-07-22 16:41:25 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020600_1_ENS6.nc
2026-07-22 16:41:25 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-22 16:41:25 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:41:25 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-22 16:41:25 INFO Queuing job for member 6...
2026-07-22 16:41:25 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:41:25 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-22 16:41:27 INFO Found: ['5256262']
2026-07-22 16:41:32 INFO [TGCC-IRENE] Submitted job with ID:['5256262']
2026-07-22 16:41:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:41:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7
2026-07-22 16:41:32 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020600_1_ENS7.nc
2026-07-22 16:41:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-22 16:41:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:41:32 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-22 16:41:32 INFO Queuing job for member 7...
2026-07-22 16:41:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:41:32 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-22 16:41:34 INFO Found: ['5256266']
2026-07-22 16:41:39 INFO [TGCC-IRENE] Submitted job with ID:['5256266']
2026-07-22 16:41:39 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:41:39 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8
2026-07-22 16:41:39 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020600_1_ENS8.nc
2026-07-22 16:41:39 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-22 16:41:39 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:41:39 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-22 16:41:39 INFO Queuing job for member 8...
2026-07-22 16:41:39 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:41:39 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-22 16:41:42 INFO Found: ['5256270']
2026-07-22 16:41:47 INFO [TGCC-IRENE] Submitted job with ID:['5256270']
2026-07-22 16:41:47 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:41:47 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9
2026-07-22 16:41:47 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020600_1_ENS9.nc
2026-07-22 16:41:47 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-22 16:41:47 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:41:47 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-22 16:41:47 INFO Queuing job for member 9...
2026-07-22 16:41:47 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:41:47 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-22 16:41:50 INFO Found: ['5256275']
2026-07-22 16:41:55 INFO [TGCC-IRENE] Submitted job with ID:['5256275']
2026-07-22 16:41:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:41:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10
2026-07-22 16:41:55 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020600_1_ENS10.nc
2026-07-22 16:41:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-22 16:41:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:41:55 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-22 16:41:55 INFO Queuing job for member 10...
2026-07-22 16:41:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:41:55 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-22 16:41:57 INFO Found: ['5256276']
2026-07-22 16:42:02 INFO [TGCC-IRENE] Submitted job with ID:['5256276']
2026-07-22 16:42:02 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:42:02 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11
2026-07-22 16:42:02 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020600_1_ENS11.nc
2026-07-22 16:42:02 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-22 16:42:02 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:42:02 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-22 16:42:02 INFO Queuing job for member 11...
2026-07-22 16:42:02 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:42:02 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-22 16:42:05 INFO Found: ['5256280']
2026-07-22 16:42:10 INFO [TGCC-IRENE] Submitted job with ID:['5256280']
2026-07-22 16:42:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:42:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12
2026-07-22 16:42:10 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020600_1_ENS12.nc
2026-07-22 16:42:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-22 16:42:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:42:10 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-22 16:42:10 INFO Queuing job for member 12...
2026-07-22 16:42:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:42:10 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-22 16:42:10 INFO Found: ['5256281']
2026-07-22 16:42:15 INFO [TGCC-IRENE] Submitted job with ID:['5256281']
2026-07-22 16:42:15 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:42:15 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13
2026-07-22 16:42:15 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020600_1_ENS13.nc
2026-07-22 16:42:15 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-22 16:42:15 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:42:15 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-22 16:42:15 INFO Queuing job for member 13...
2026-07-22 16:42:15 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:42:15 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-22 16:42:16 INFO Found: ['5256284']
2026-07-22 16:42:21 INFO [TGCC-IRENE] Submitted job with ID:['5256284']
2026-07-22 16:42:21 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:42:21 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14
2026-07-22 16:42:21 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020600_1_ENS14.nc
2026-07-22 16:42:21 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-22 16:42:21 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:42:21 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-22 16:42:21 INFO Queuing job for member 14...
2026-07-22 16:42:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:42:21 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-22 16:42:22 INFO Found: ['5256294']
2026-07-22 16:42:27 INFO [TGCC-IRENE] Submitted job with ID:['5256294']
2026-07-22 16:42:27 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:42:27 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15
2026-07-22 16:42:27 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020600_1_ENS15.nc
2026-07-22 16:42:27 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-22 16:42:27 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:42:27 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-22 16:42:27 INFO Queuing job for member 15...
2026-07-22 16:42:27 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:42:27 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-22 16:42:28 INFO Found: ['5256304']
2026-07-22 16:42:33 INFO [TGCC-IRENE] Submitted job with ID:['5256304']
2026-07-22 16:42:33 INFO Checking job status ...
2026-07-22 16:42:33 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:42:33 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:42:33 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:42:33 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:42:33 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:42:33 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:42:33 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:42:33 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:42:33 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:42:33 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:42:33 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:42:33 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:42:33 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:42:33 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:42:33 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:42:33 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:42:48 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:42:48 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:42:48 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:42:48 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:42:48 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:42:48 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:42:48 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:42:48 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:42:48 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:42:48 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:42:48 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:42:48 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:42:48 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:42:48 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:42:48 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:42:48 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:43:04 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:43:04 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:43:04 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:43:04 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:43:04 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:43:04 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:43:04 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:43:04 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:43:04 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:43:04 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:43:04 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:43:04 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:43:04 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:43:04 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:43:04 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:43:04 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:43:19 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:43:19 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:43:19 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:43:19 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:43:19 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:43:19 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:43:19 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:43:21 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:43:21 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:43:22 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:43:22 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:43:22 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:43:22 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:43:22 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:43:22 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:43:22 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:43:37 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:43:37 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:43:37 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:43:37 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:43:37 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:43:37 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:43:37 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:43:37 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:43:37 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:43:37 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:43:37 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:43:37 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:43:37 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:43:37 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:43:37 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:43:37 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:43:52 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:43:52 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:43:52 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:43:52 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:43:52 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:43:52 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:43:52 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:43:52 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:43:52 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:43:52 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:43:52 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:43:52 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:43:52 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:43:52 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:43:52 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:43:52 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:44:07 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:44:07 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:44:07 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:44:07 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:44:07 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:44:07 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:44:07 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:44:07 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:44:07 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:44:07 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:44:07 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:44:08 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:44:08 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:44:08 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:44:08 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:44:08 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:44:23 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:44:23 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:44:23 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:44:23 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:44:23 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:44:23 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:44:23 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:44:23 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:44:23 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:44:23 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:44:23 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:44:23 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:44:23 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:44:23 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:44:23 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:44:23 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:44:38 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:44:38 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:44:38 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:44:38 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:44:38 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:44:38 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:44:38 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:44:38 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:44:38 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:44:38 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:44:38 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:44:38 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:44:38 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:44:38 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:44:38 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:44:38 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:44:53 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:44:53 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:44:53 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:44:53 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:44:53 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:44:55 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:44:55 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:44:55 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:44:55 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:44:56 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:44:56 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:44:56 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:44:56 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:44:56 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:44:56 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:44:56 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:45:11 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:45:11 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:45:11 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:45:11 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:45:11 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:45:11 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:45:11 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:45:11 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:45:11 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:45:11 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:45:11 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:45:11 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:45:11 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:45:11 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:45:11 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:45:11 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:45:26 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:45:26 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:45:26 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:45:26 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:45:26 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:45:26 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:45:26 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:45:26 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:45:26 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:45:26 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:45:26 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:45:26 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:45:26 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:45:26 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:45:26 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:45:26 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:45:41 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:45:41 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:45:41 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:45:41 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:45:41 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:45:41 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:45:41 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:45:41 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:45:41 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:45:41 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:45:41 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:45:41 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:45:41 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:45:41 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:45:42 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:45:42 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:45:58 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:45:58 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:45:58 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:45:58 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:45:58 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:45:58 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:45:58 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:45:58 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:45:58 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:45:58 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:45:58 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:45:58 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:45:58 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:45:58 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:45:58 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:45:58 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:46:13 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:46:13 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:46:13 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:46:13 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:46:13 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:46:14 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:46:14 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:46:14 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:46:14 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:46:14 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:46:14 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:46:16 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:46:16 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:46:16 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:46:16 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:46:16 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:46:31 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:46:31 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:46:31 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:46:31 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:46:31 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:46:31 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:46:31 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:46:31 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:46:31 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:46:31 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:46:31 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:46:31 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:46:31 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:46:31 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:46:33 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:46:33 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:46:48 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:46:48 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:46:49 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:46:49 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:46:49 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:46:49 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:46:49 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:46:49 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:46:49 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:46:49 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:46:49 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:46:49 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:46:49 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:46:49 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:46:49 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:46:49 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:47:04 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:47:04 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:47:04 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:47:04 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:47:04 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:47:04 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:47:04 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:47:04 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:47:04 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:47:04 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:47:04 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:47:04 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:47:04 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:47:04 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:47:04 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:47:04 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:47:19 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:47:19 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:47:19 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:47:19 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:47:19 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:47:19 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:47:19 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:47:19 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:47:19 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:47:19 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:47:19 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:47:19 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:47:19 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:47:19 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:47:20 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:47:20 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:47:35 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:47:35 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:47:35 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:47:35 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:47:38 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:47:38 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:47:38 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:47:38 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:47:38 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:47:38 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:47:38 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:47:38 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:47:38 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:47:38 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:47:38 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:47:38 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:47:53 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:47:53 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:47:53 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:47:53 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:47:53 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:47:53 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:47:53 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:47:53 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:47:53 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:47:53 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:47:53 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:47:53 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:47:55 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:47:55 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:47:55 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:47:55 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:48:11 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:48:11 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:48:11 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:48:11 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:48:11 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:48:11 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:48:11 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:48:11 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:48:11 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:48:11 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:48:11 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:48:11 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:48:11 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:48:11 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:48:11 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:48:11 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:48:26 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:48:26 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:48:26 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:48:26 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:48:26 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:48:26 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:48:26 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:48:26 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:48:26 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:48:26 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:48:26 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:48:26 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:48:26 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:48:26 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:48:26 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:48:26 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:48:41 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:48:41 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:48:41 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:48:41 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:48:41 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:48:41 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:48:41 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:48:41 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:48:41 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:48:41 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:48:41 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:48:41 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:48:41 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:48:41 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:48:41 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:48:41 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:48:56 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:48:57 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:48:57 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:48:57 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:48:57 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:48:57 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:48:57 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:48:57 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:48:57 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:48:57 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:48:57 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:48:57 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:48:57 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:48:57 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:48:57 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:48:57 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:49:13 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:49:13 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:49:13 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:49:13 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:49:13 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:49:13 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:49:13 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:49:13 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:49:13 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:49:13 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:49:13 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:49:13 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:49:13 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:49:13 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:49:13 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:49:13 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:49:28 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:49:28 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:49:28 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:49:28 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:49:28 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:49:28 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:49:30 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:49:30 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:49:30 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:49:30 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:49:30 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:49:30 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:49:30 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:49:30 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:49:30 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:49:30 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:49:45 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:49:45 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:49:45 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:49:45 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:49:46 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:49:46 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:49:46 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:49:46 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:49:46 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:49:46 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:49:46 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:49:46 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:49:46 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:49:46 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:49:46 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:49:46 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:50:01 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:50:01 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:50:01 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:50:02 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:50:02 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:50:02 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:50:02 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:50:02 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:50:02 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:50:02 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:50:02 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:50:02 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:50:02 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:50:02 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:50:02 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:50:02 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:50:17 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:50:17 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:50:17 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:50:17 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:50:17 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:50:17 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:50:17 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:50:17 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:50:17 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:50:17 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:50:17 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:50:17 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:50:17 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:50:17 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:50:18 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:50:18 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:50:33 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:50:33 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:50:33 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:50:33 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:50:33 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:50:33 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:50:33 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:50:33 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:50:33 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:50:33 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:50:33 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:50:33 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:50:33 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:50:33 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:50:33 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:50:33 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:50:49 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:50:49 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:50:49 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:50:49 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:50:49 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:50:49 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:50:49 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:50:49 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:50:49 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:50:49 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:50:49 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:50:49 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:50:49 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:50:49 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:50:49 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:50:49 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:51:04 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:51:04 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:51:04 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:51:04 INFO None 5256255: status RUNNING/PENDING
2026-07-22 16:51:04 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:51:04 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:51:04 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:51:04 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:51:04 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:51:04 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:51:04 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:51:04 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:51:04 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:51:04 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:51:04 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:51:04 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:51:19 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:51:19 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:51:19 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:51:21 INFO None 5256255: status FINISHED
2026-07-22 16:51:22 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:51:22 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:51:22 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:51:22 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:51:22 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:51:22 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:51:22 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:51:22 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:51:22 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:51:22 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:51:22 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:51:22 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:51:37 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:51:37 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:51:37 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:51:37 INFO None 5256255: status FINISHED
2026-07-22 16:51:37 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:51:37 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:51:37 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:51:37 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:51:37 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:51:37 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:51:37 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:51:37 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:51:39 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:51:39 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:51:39 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:51:39 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:51:54 INFO None 5256224: status RUNNING/PENDING
2026-07-22 16:51:54 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:51:54 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:51:54 INFO None 5256255: status FINISHED
2026-07-22 16:51:54 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:51:54 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:51:54 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:51:54 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:51:54 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:51:54 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:51:54 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:51:54 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:51:54 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:51:54 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:51:54 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:51:54 INFO Jobs still running: ['5256224', '5256236', '5256246', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:52:09 INFO None 5256224: status FINISHED
2026-07-22 16:52:09 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:52:09 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:52:09 INFO None 5256255: status FINISHED
2026-07-22 16:52:10 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:52:10 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:52:10 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:52:10 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:52:10 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:52:10 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:52:10 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:52:10 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:52:10 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:52:10 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:52:10 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:52:10 INFO Jobs still running: ['5256236', '5256246', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:52:25 INFO None 5256224: status FINISHED
2026-07-22 16:52:25 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:52:25 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:52:25 INFO None 5256255: status FINISHED
2026-07-22 16:52:25 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:52:25 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:52:25 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:52:25 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:52:25 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:52:25 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:52:25 INFO None 5256280: status RUNNING/PENDING
2026-07-22 16:52:25 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:52:25 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:52:25 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:52:25 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:52:25 INFO Jobs still running: ['5256236', '5256246', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:52:41 INFO None 5256224: status FINISHED
2026-07-22 16:52:41 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:52:41 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:52:41 INFO None 5256255: status FINISHED
2026-07-22 16:52:41 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:52:41 INFO None 5256262: status RUNNING/PENDING
2026-07-22 16:52:41 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:52:41 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:52:41 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:52:41 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:52:41 INFO None 5256280: status FINISHED
2026-07-22 16:52:41 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:52:41 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:52:41 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:52:41 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:52:41 INFO Jobs still running: ['5256236', '5256246', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:52:56 INFO None 5256224: status FINISHED
2026-07-22 16:52:56 INFO None 5256236: status RUNNING/PENDING
2026-07-22 16:52:56 INFO None 5256246: status RUNNING/PENDING
2026-07-22 16:52:56 INFO None 5256255: status FINISHED
2026-07-22 16:52:56 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:52:56 INFO None 5256262: status FINISHED
2026-07-22 16:52:56 INFO None 5256266: status RUNNING/PENDING
2026-07-22 16:52:56 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:52:56 INFO None 5256275: status RUNNING/PENDING
2026-07-22 16:52:58 INFO None 5256276: status RUNNING/PENDING
2026-07-22 16:52:58 INFO None 5256280: status FINISHED
2026-07-22 16:52:58 INFO None 5256281: status RUNNING/PENDING
2026-07-22 16:52:59 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:52:59 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:52:59 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:52:59 INFO Jobs still running: ['5256236', '5256246', '5256259', '5256266', '5256270', '5256275', '5256276', '5256281', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:53:14 INFO None 5256224: status FINISHED
2026-07-22 16:53:14 INFO None 5256236: status FINISHED
2026-07-22 16:53:14 INFO None 5256246: status FINISHED
2026-07-22 16:53:14 INFO None 5256255: status FINISHED
2026-07-22 16:53:14 INFO None 5256259: status RUNNING/PENDING
2026-07-22 16:53:14 INFO None 5256262: status FINISHED
2026-07-22 16:53:14 INFO None 5256266: status FINISHED
2026-07-22 16:53:14 INFO None 5256270: status RUNNING/PENDING
2026-07-22 16:53:14 INFO None 5256275: status FINISHED
2026-07-22 16:53:14 INFO None 5256276: status FINISHED
2026-07-22 16:53:14 INFO None 5256280: status FINISHED
2026-07-22 16:53:14 INFO None 5256281: status FINISHED
2026-07-22 16:53:14 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:53:14 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:53:14 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:53:14 INFO Jobs still running: ['5256259', '5256270', '5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:53:29 INFO None 5256224: status FINISHED
2026-07-22 16:53:29 INFO None 5256236: status FINISHED
2026-07-22 16:53:29 INFO None 5256246: status FINISHED
2026-07-22 16:53:29 INFO None 5256255: status FINISHED
2026-07-22 16:53:29 INFO None 5256259: status FINISHED
2026-07-22 16:53:29 INFO None 5256262: status FINISHED
2026-07-22 16:53:29 INFO None 5256266: status FINISHED
2026-07-22 16:53:29 INFO None 5256270: status FINISHED
2026-07-22 16:53:29 INFO None 5256275: status FINISHED
2026-07-22 16:53:29 INFO None 5256276: status FINISHED
2026-07-22 16:53:29 INFO None 5256280: status FINISHED
2026-07-22 16:53:29 INFO None 5256281: status FINISHED
2026-07-22 16:53:29 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:53:29 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:53:29 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:53:29 INFO Jobs still running: ['5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:53:44 INFO None 5256224: status FINISHED
2026-07-22 16:53:44 INFO None 5256236: status FINISHED
2026-07-22 16:53:44 INFO None 5256246: status FINISHED
2026-07-22 16:53:44 INFO None 5256255: status FINISHED
2026-07-22 16:53:44 INFO None 5256259: status FINISHED
2026-07-22 16:53:44 INFO None 5256262: status FINISHED
2026-07-22 16:53:44 INFO None 5256266: status FINISHED
2026-07-22 16:53:44 INFO None 5256270: status FINISHED
2026-07-22 16:53:44 INFO None 5256275: status FINISHED
2026-07-22 16:53:44 INFO None 5256276: status FINISHED
2026-07-22 16:53:44 INFO None 5256280: status FINISHED
2026-07-22 16:53:44 INFO None 5256281: status FINISHED
2026-07-22 16:53:44 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:53:44 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:53:44 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:53:44 INFO Jobs still running: ['5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:53:59 INFO None 5256224: status FINISHED
2026-07-22 16:53:59 INFO None 5256236: status FINISHED
2026-07-22 16:53:59 INFO None 5256246: status FINISHED
2026-07-22 16:53:59 INFO None 5256255: status FINISHED
2026-07-22 16:53:59 INFO None 5256259: status FINISHED
2026-07-22 16:54:00 INFO None 5256262: status FINISHED
2026-07-22 16:54:00 INFO None 5256266: status FINISHED
2026-07-22 16:54:00 INFO None 5256270: status FINISHED
2026-07-22 16:54:00 INFO None 5256275: status FINISHED
2026-07-22 16:54:00 INFO None 5256276: status FINISHED
2026-07-22 16:54:00 INFO None 5256280: status FINISHED
2026-07-22 16:54:00 INFO None 5256281: status FINISHED
2026-07-22 16:54:00 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:54:00 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:54:00 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:54:00 INFO Jobs still running: ['5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:54:15 INFO None 5256224: status FINISHED
2026-07-22 16:54:15 INFO None 5256236: status FINISHED
2026-07-22 16:54:15 INFO None 5256246: status FINISHED
2026-07-22 16:54:15 INFO None 5256255: status FINISHED
2026-07-22 16:54:15 INFO None 5256259: status FINISHED
2026-07-22 16:54:15 INFO None 5256262: status FINISHED
2026-07-22 16:54:17 INFO None 5256266: status FINISHED
2026-07-22 16:54:17 INFO None 5256270: status FINISHED
2026-07-22 16:54:17 INFO None 5256275: status FINISHED
2026-07-22 16:54:17 INFO None 5256276: status FINISHED
2026-07-22 16:54:17 INFO None 5256280: status FINISHED
2026-07-22 16:54:17 INFO None 5256281: status FINISHED
2026-07-22 16:54:17 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:54:17 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:54:17 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:54:17 INFO Jobs still running: ['5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:54:32 INFO None 5256224: status FINISHED
2026-07-22 16:54:32 INFO None 5256236: status FINISHED
2026-07-22 16:54:32 INFO None 5256246: status FINISHED
2026-07-22 16:54:32 INFO None 5256255: status FINISHED
2026-07-22 16:54:32 INFO None 5256259: status FINISHED
2026-07-22 16:54:32 INFO None 5256262: status FINISHED
2026-07-22 16:54:32 INFO None 5256266: status FINISHED
2026-07-22 16:54:32 INFO None 5256270: status FINISHED
2026-07-22 16:54:32 INFO None 5256275: status FINISHED
2026-07-22 16:54:32 INFO None 5256276: status FINISHED
2026-07-22 16:54:32 INFO None 5256280: status FINISHED
2026-07-22 16:54:32 INFO None 5256281: status FINISHED
2026-07-22 16:54:32 INFO None 5256284: status RUNNING/PENDING
2026-07-22 16:54:32 INFO None 5256294: status RUNNING/PENDING
2026-07-22 16:54:32 INFO None 5256304: status RUNNING/PENDING
2026-07-22 16:54:32 INFO Jobs still running: ['5256284', '5256294', '5256304']. Waiting...
2026-07-22 16:54:47 INFO None 5256224: status FINISHED
2026-07-22 16:54:47 INFO None 5256236: status FINISHED
2026-07-22 16:54:50 INFO None 5256246: status FINISHED
2026-07-22 16:54:50 INFO None 5256255: status FINISHED
2026-07-22 16:54:50 INFO None 5256259: status FINISHED
2026-07-22 16:54:50 INFO None 5256262: status FINISHED
2026-07-22 16:54:50 INFO None 5256266: status FINISHED
2026-07-22 16:54:50 INFO None 5256270: status FINISHED
2026-07-22 16:54:50 INFO None 5256275: status FINISHED
2026-07-22 16:54:50 INFO None 5256276: status FINISHED
2026-07-22 16:54:50 INFO None 5256280: status FINISHED
2026-07-22 16:54:50 INFO None 5256281: status FINISHED
2026-07-22 16:54:50 INFO None 5256284: status FINISHED
2026-07-22 16:54:50 INFO None 5256294: status FINISHED
2026-07-22 16:54:50 INFO None 5256304: status FINISHED
2026-07-22 16:54:50 INFO Jobs ['5256224', '5256236', '5256246', '5256255', '5256259', '5256262', '5256266', '5256270', '5256275', '5256276', '5256280', '5256281', '5256284', '5256294', '5256304'] have finished
2026-07-22 16:54:50 INFO Checking restart files were created ...
2026-07-22 16:54:50 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020601_8_ENS1.nc(3005806795 bytes)
2026-07-22 16:54:50 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020601_8_ENS2.nc(3005806795 bytes)
2026-07-22 16:54:50 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020601_8_ENS3.nc(3005806795 bytes)
2026-07-22 16:54:50 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020601_8_ENS4.nc(3005806795 bytes)
2026-07-22 16:54:50 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020601_8_ENS5.nc(3005806795 bytes)
2026-07-22 16:54:50 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020601_8_ENS6.nc(3005806795 bytes)
2026-07-22 16:54:50 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020601_8_ENS7.nc(3005806795 bytes)
2026-07-22 16:54:50 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020601_8_ENS8.nc(3005806795 bytes)
2026-07-22 16:54:50 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020601_8_ENS9.nc(3005806795 bytes)
2026-07-22 16:54:50 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020601_8_ENS10.nc(3005806795 bytes)
2026-07-22 16:54:50 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020601_8_ENS11.nc(3005806795 bytes)
2026-07-22 16:54:50 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020601_8_ENS12.nc(3005806795 bytes)
2026-07-22 16:54:50 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020601_8_ENS13.nc(3005806795 bytes)
2026-07-22 16:54:50 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020601_8_ENS14.nc(3005806795 bytes)
2026-07-22 16:54:50 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020601_8_ENS15.nc(3005806795 bytes)
2026-07-22 16:54:50 INFO  Run_model() completed successfully.
2026-07-22 16:54:50 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 16:54:50 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 09:00:00 days=153072 seconds=32400
2026-07-22 16:54:50 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-22 16:54:50 INFO [TIME] increment current_time 2020-02-06 01:00:00 -> 2020-02-06 09:00:00
2026-07-22 16:54:50 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 16:54:50 INFO ---------->>> Running process_satellite_data()
2026-07-22 16:54:50 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12002.nc
2026-07-22 16:54:50 INFO ---------->>> Running run_obs_converter()
2026-07-22 16:54:50 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_33067_153072.out
2026-07-22 16:54:50 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_33067_153072.out
2026-07-22 16:54:50 INFO ---------->>> Running DART
2026-07-22 16:54:50 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-22 16:54:50 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/chim_ENS1_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/chim_ENS1_2020020609_1_out_toDART.nc
2026-07-22 16:54:51 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/chim_ENS2_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/chim_ENS2_2020020609_1_out_toDART.nc
2026-07-22 16:54:51 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/chim_ENS3_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/chim_ENS3_2020020609_1_out_toDART.nc
2026-07-22 16:54:51 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/chim_ENS4_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/chim_ENS4_2020020609_1_out_toDART.nc
2026-07-22 16:54:52 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/chim_ENS5_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/chim_ENS5_2020020609_1_out_toDART.nc
2026-07-22 16:54:52 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/chim_ENS6_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/chim_ENS6_2020020609_1_out_toDART.nc
2026-07-22 16:54:53 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/chim_ENS7_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/chim_ENS7_2020020609_1_out_toDART.nc
2026-07-22 16:54:53 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/chim_ENS8_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/chim_ENS8_2020020609_1_out_toDART.nc
2026-07-22 16:54:53 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/chim_ENS9_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/chim_ENS9_2020020609_1_out_toDART.nc
2026-07-22 16:54:54 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/chim_ENS10_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/chim_ENS10_2020020609_1_out_toDART.nc
2026-07-22 16:54:54 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/chim_ENS11_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/chim_ENS11_2020020609_1_out_toDART.nc
2026-07-22 16:54:54 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/chim_ENS12_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/chim_ENS12_2020020609_1_out_toDART.nc
2026-07-22 16:54:55 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/chim_ENS13_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/chim_ENS13_2020020609_1_out_toDART.nc
2026-07-22 16:54:55 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/chim_ENS14_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/chim_ENS14_2020020609_1_out_toDART.nc
2026-07-22 16:54:55 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/chim_ENS15_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/chim_ENS15_2020020609_1_out_toDART.nc
2026-07-22 16:54:56 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-22 16:54:56 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-22 16:54:56 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-22 16:54:56 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-22 16:54:56 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-22 16:54:56 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-22 16:55:03 INFO Found: []
2026-07-22 16:55:03 INFO No job id returned by command ./run_filter.bsh
2026-07-22 16:55:03 INFO No monitoring will be performed
2026-07-22 16:55:03 INFO Moving DART output files to analysis and preassim directories for date 2020020609 if present ...
2026-07-22 16:55:03 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 16:55:03 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 16:55:03 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 16:55:03 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 16:55:03 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 16:55:03 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 16:55:03 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 16:55:03 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 16:55:03 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 16:55:03 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 16:55:03 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 16:55:03 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 16:55:03 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 16:55:03 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 16:55:03 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 16:55:03 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 16:55:03 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 16:55:03 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 16:55:03 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 16:55:03 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 16:55:03 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 16:55:03 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 16:55:03 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 16:55:03 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 16:55:03 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 16:55:03 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 16:55:03 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 16:55:03 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 16:55:03 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 16:55:03 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 16:55:03 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 16:55:03 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020609'
2026-07-22 16:55:03 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 16:55:03 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020609'
2026-07-22 16:55:03 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-22 16:55:06 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-22 16:55:06 INFO run_dart() is DONE.
2026-07-22 16:55:06 INFO ---------->>> Running update_pollutant_in_end()
2026-07-22 16:55:06 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 16:55:06 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS1.nc
2026-07-22 16:55:07 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 16:55:08 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS1.nc
2026-07-22 16:55:08 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 16:55:08 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 16:55:08 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS1.nc
2026-07-22 16:55:08 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 16:55:09 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS1.nc
2026-07-22 16:55:09 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 16:55:09 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 16:55:09 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 16:55:09 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS2.nc
2026-07-22 16:55:10 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 16:55:10 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS2.nc
2026-07-22 16:55:10 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 16:55:11 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 16:55:11 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS2.nc
2026-07-22 16:55:11 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 16:55:12 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS2.nc
2026-07-22 16:55:12 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 16:55:12 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 16:55:12 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 16:55:12 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS3.nc
2026-07-22 16:55:13 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 16:55:13 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS3.nc
2026-07-22 16:55:13 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 16:55:14 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 16:55:14 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS3.nc
2026-07-22 16:55:14 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 16:55:15 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS3.nc
2026-07-22 16:55:15 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 16:55:15 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 16:55:15 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 16:55:15 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS4.nc
2026-07-22 16:55:16 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 16:55:16 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS4.nc
2026-07-22 16:55:16 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 16:55:17 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 16:55:17 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS4.nc
2026-07-22 16:55:17 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 16:55:18 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS4.nc
2026-07-22 16:55:18 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 16:55:18 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 16:55:18 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 16:55:18 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS5.nc
2026-07-22 16:55:19 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 16:55:19 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS5.nc
2026-07-22 16:55:19 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 16:55:20 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 16:55:20 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS5.nc
2026-07-22 16:55:20 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 16:55:21 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS5.nc
2026-07-22 16:55:21 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 16:55:21 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 16:55:21 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 16:55:21 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS6.nc
2026-07-22 16:55:22 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 16:55:22 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS6.nc
2026-07-22 16:55:22 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 16:55:23 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 16:55:23 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS6.nc
2026-07-22 16:55:23 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 16:55:24 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS6.nc
2026-07-22 16:55:24 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 16:55:24 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 16:55:24 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 16:55:24 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS7.nc
2026-07-22 16:55:25 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 16:55:25 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS7.nc
2026-07-22 16:55:25 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 16:55:25 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 16:55:25 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS7.nc
2026-07-22 16:55:26 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 16:55:26 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS7.nc
2026-07-22 16:55:26 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 16:55:26 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 16:55:27 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 16:55:27 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS8.nc
2026-07-22 16:55:27 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 16:55:28 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS8.nc
2026-07-22 16:55:28 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 16:55:28 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 16:55:29 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS8.nc
2026-07-22 16:55:29 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 16:55:30 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS8.nc
2026-07-22 16:55:30 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 16:55:30 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 16:55:30 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 16:55:30 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS9.nc
2026-07-22 16:55:30 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 16:55:31 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS9.nc
2026-07-22 16:55:31 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 16:55:31 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 16:55:31 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS9.nc
2026-07-22 16:55:32 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 16:55:32 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS9.nc
2026-07-22 16:55:32 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 16:55:33 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 16:55:33 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 16:55:33 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS10.nc
2026-07-22 16:55:33 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 16:55:34 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS10.nc
2026-07-22 16:55:34 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 16:55:35 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 16:55:35 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS10.nc
2026-07-22 16:55:35 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 16:55:36 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS10.nc
2026-07-22 16:55:36 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 16:55:36 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 16:55:36 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 16:55:36 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS11.nc
2026-07-22 16:55:37 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 16:55:37 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS11.nc
2026-07-22 16:55:37 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 16:55:38 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 16:55:38 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS11.nc
2026-07-22 16:55:38 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 16:55:39 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS11.nc
2026-07-22 16:55:39 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 16:55:39 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 16:55:39 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 16:55:39 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS12.nc
2026-07-22 16:55:40 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 16:55:41 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS12.nc
2026-07-22 16:55:41 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 16:55:41 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 16:55:41 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS12.nc
2026-07-22 16:55:41 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 16:55:42 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS12.nc
2026-07-22 16:55:42 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 16:55:42 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 16:55:42 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 16:55:42 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS13.nc
2026-07-22 16:55:43 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 16:55:44 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS13.nc
2026-07-22 16:55:44 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 16:55:44 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 16:55:44 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS13.nc
2026-07-22 16:55:45 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 16:55:45 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS13.nc
2026-07-22 16:55:45 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 16:55:45 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 16:55:46 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 16:55:46 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS14.nc
2026-07-22 16:55:46 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 16:55:47 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS14.nc
2026-07-22 16:55:47 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 16:55:47 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 16:55:47 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS14.nc
2026-07-22 16:55:48 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 16:55:48 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS14.nc
2026-07-22 16:55:48 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 16:55:48 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 16:55:49 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 16:55:49 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS15.nc
2026-07-22 16:55:49 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 16:55:50 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS15.nc
2026-07-22 16:55:50 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 16:55:50 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 16:55:50 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS15.nc
2026-07-22 16:55:51 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 16:55:51 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS15.nc
2026-07-22 16:55:51 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 16:55:51 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 16:55:51 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-22 16:55:51 INFO [TIME] step_end current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 16:55:51 INFO [TIME] step_start current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 16:55:51 INFO [TIME] window start=2020-02-06 09:00:00 end=2020-02-06 11:00:00 run_hours=2 has_assimilation=True
2026-07-22 16:55:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:55:53 INFO Hourly dataset computed and listing created
2026-07-22 16:55:58 INFO Hourly dataset computed
2026-07-22 16:55:58 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:55:59 INFO Hourly dataset computed and listing created
2026-07-22 16:56:00 INFO Hourly dataset computed
2026-07-22 16:56:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:56:01 INFO Hourly dataset computed and listing created
2026-07-22 16:56:01 INFO Hourly dataset computed
2026-07-22 16:56:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:56:02 INFO Hourly dataset computed and listing created
2026-07-22 16:56:03 INFO Hourly dataset computed
2026-07-22 16:56:03 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:56:04 INFO Hourly dataset computed and listing created
2026-07-22 16:56:05 INFO Hourly dataset computed
2026-07-22 16:56:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:56:06 INFO Hourly dataset computed and listing created
2026-07-22 16:56:06 INFO Hourly dataset computed
2026-07-22 16:56:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:56:07 INFO Hourly dataset computed and listing created
2026-07-22 16:56:08 INFO Hourly dataset computed
2026-07-22 16:56:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:56:09 INFO Hourly dataset computed and listing created
2026-07-22 16:56:10 INFO Hourly dataset computed
2026-07-22 16:56:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:56:11 INFO Hourly dataset computed and listing created
2026-07-22 16:56:11 INFO Hourly dataset computed
2026-07-22 16:56:11 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:56:12 INFO Hourly dataset computed and listing created
2026-07-22 16:56:13 INFO Hourly dataset computed
2026-07-22 16:56:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:56:14 INFO Hourly dataset computed and listing created
2026-07-22 16:56:15 INFO Hourly dataset computed
2026-07-22 16:56:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:56:16 INFO Hourly dataset computed and listing created
2026-07-22 16:56:17 INFO Hourly dataset computed
2026-07-22 16:56:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:56:18 INFO Hourly dataset computed and listing created
2026-07-22 16:56:18 INFO Hourly dataset computed
2026-07-22 16:56:18 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:56:19 INFO Hourly dataset computed and listing created
2026-07-22 16:56:20 INFO Hourly dataset computed
2026-07-22 16:56:20 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 16:56:21 INFO Hourly dataset computed and listing created
2026-07-22 16:56:22 INFO Hourly dataset computed
2026-07-22 16:56:22 INFO ---------->>> Running CHIMERE model from 2020-02-06 09:00:00 to 2020-02-06 11:00:00
2026-07-22 16:56:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:56:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1
2026-07-22 16:56:22 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020601_8_ENS1.nc
2026-07-22 16:56:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-22 16:56:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:56:22 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-22 16:56:22 INFO Queuing job for member 1...
2026-07-22 16:56:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:56:22 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-22 16:56:24 INFO Found: ['5256732']
2026-07-22 16:56:29 INFO [TGCC-IRENE] Submitted job with ID:['5256732']
2026-07-22 16:56:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:56:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2
2026-07-22 16:56:29 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020601_8_ENS2.nc
2026-07-22 16:56:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-22 16:56:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:56:29 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-22 16:56:29 INFO Queuing job for member 2...
2026-07-22 16:56:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:56:29 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-22 16:56:32 INFO Found: ['5256735']
2026-07-22 16:56:37 INFO [TGCC-IRENE] Submitted job with ID:['5256735']
2026-07-22 16:56:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:56:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3
2026-07-22 16:56:37 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020601_8_ENS3.nc
2026-07-22 16:56:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-22 16:56:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:56:37 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-22 16:56:37 INFO Queuing job for member 3...
2026-07-22 16:56:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:56:37 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-22 16:56:39 INFO Found: ['5256740']
2026-07-22 16:56:44 INFO [TGCC-IRENE] Submitted job with ID:['5256740']
2026-07-22 16:56:44 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:56:44 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4
2026-07-22 16:56:44 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020601_8_ENS4.nc
2026-07-22 16:56:44 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-22 16:56:44 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:56:44 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-22 16:56:44 INFO Queuing job for member 4...
2026-07-22 16:56:44 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:56:44 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-22 16:56:45 INFO Found: ['5256742']
2026-07-22 16:56:50 INFO [TGCC-IRENE] Submitted job with ID:['5256742']
2026-07-22 16:56:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:56:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5
2026-07-22 16:56:50 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020601_8_ENS5.nc
2026-07-22 16:56:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-22 16:56:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:56:50 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-22 16:56:50 INFO Queuing job for member 5...
2026-07-22 16:56:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:56:50 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-22 16:56:51 INFO Found: ['5256747']
2026-07-22 16:56:56 INFO [TGCC-IRENE] Submitted job with ID:['5256747']
2026-07-22 16:56:56 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:56:56 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6
2026-07-22 16:56:56 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020601_8_ENS6.nc
2026-07-22 16:56:56 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-22 16:56:56 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:56:56 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-22 16:56:56 INFO Queuing job for member 6...
2026-07-22 16:56:56 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:56:56 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-22 16:56:56 INFO Found: ['5256756']
2026-07-22 16:57:01 INFO [TGCC-IRENE] Submitted job with ID:['5256756']
2026-07-22 16:57:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:57:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7
2026-07-22 16:57:01 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020601_8_ENS7.nc
2026-07-22 16:57:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-22 16:57:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:57:02 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-22 16:57:02 INFO Queuing job for member 7...
2026-07-22 16:57:02 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:57:02 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-22 16:57:02 INFO Found: ['5256761']
2026-07-22 16:57:07 INFO [TGCC-IRENE] Submitted job with ID:['5256761']
2026-07-22 16:57:07 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:57:07 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8
2026-07-22 16:57:07 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020601_8_ENS8.nc
2026-07-22 16:57:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-22 16:57:07 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:57:07 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-22 16:57:07 INFO Queuing job for member 8...
2026-07-22 16:57:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:57:07 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-22 16:57:08 INFO Found: ['5256765']
2026-07-22 16:57:13 INFO [TGCC-IRENE] Submitted job with ID:['5256765']
2026-07-22 16:57:13 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:57:13 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9
2026-07-22 16:57:13 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020601_8_ENS9.nc
2026-07-22 16:57:13 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-22 16:57:13 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:57:13 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-22 16:57:13 INFO Queuing job for member 9...
2026-07-22 16:57:13 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:57:13 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-22 16:57:14 INFO Found: ['5256770']
2026-07-22 16:57:19 INFO [TGCC-IRENE] Submitted job with ID:['5256770']
2026-07-22 16:57:19 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:57:19 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10
2026-07-22 16:57:19 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020601_8_ENS10.nc
2026-07-22 16:57:19 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-22 16:57:19 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:57:19 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-22 16:57:19 INFO Queuing job for member 10...
2026-07-22 16:57:19 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:57:19 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-22 16:57:20 INFO Found: ['5256776']
2026-07-22 16:57:25 INFO [TGCC-IRENE] Submitted job with ID:['5256776']
2026-07-22 16:57:25 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:57:25 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11
2026-07-22 16:57:25 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020601_8_ENS11.nc
2026-07-22 16:57:25 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-22 16:57:25 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:57:25 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-22 16:57:25 INFO Queuing job for member 11...
2026-07-22 16:57:25 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:57:25 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-22 16:57:25 INFO Found: ['5256782']
2026-07-22 16:57:30 INFO [TGCC-IRENE] Submitted job with ID:['5256782']
2026-07-22 16:57:30 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:57:30 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12
2026-07-22 16:57:30 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020601_8_ENS12.nc
2026-07-22 16:57:30 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-22 16:57:30 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:57:30 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-22 16:57:30 INFO Queuing job for member 12...
2026-07-22 16:57:30 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:57:30 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-22 16:57:32 INFO Found: ['5256788']
2026-07-22 16:57:37 INFO [TGCC-IRENE] Submitted job with ID:['5256788']
2026-07-22 16:57:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:57:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13
2026-07-22 16:57:37 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020601_8_ENS13.nc
2026-07-22 16:57:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-22 16:57:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:57:37 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-22 16:57:37 INFO Queuing job for member 13...
2026-07-22 16:57:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:57:37 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-22 16:57:40 INFO Found: ['5256793']
2026-07-22 16:57:45 INFO [TGCC-IRENE] Submitted job with ID:['5256793']
2026-07-22 16:57:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:57:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14
2026-07-22 16:57:45 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020601_8_ENS14.nc
2026-07-22 16:57:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-22 16:57:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:57:45 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-22 16:57:45 INFO Queuing job for member 14...
2026-07-22 16:57:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:57:45 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-22 16:57:47 INFO Found: ['5256798']
2026-07-22 16:57:52 INFO [TGCC-IRENE] Submitted job with ID:['5256798']
2026-07-22 16:57:52 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 16:57:52 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15
2026-07-22 16:57:52 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020601_8_ENS15.nc
2026-07-22 16:57:52 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-22 16:57:52 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 16:57:52 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-22 16:57:52 INFO Queuing job for member 15...
2026-07-22 16:57:52 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 16:57:52 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-22 16:57:55 INFO Found: ['5256803']
2026-07-22 16:58:00 INFO [TGCC-IRENE] Submitted job with ID:['5256803']
2026-07-22 16:58:00 INFO Checking job status ...
2026-07-22 16:58:00 INFO None 5256732: status RUNNING/PENDING
2026-07-22 16:58:00 INFO None 5256735: status RUNNING/PENDING
2026-07-22 16:58:00 INFO None 5256740: status RUNNING/PENDING
2026-07-22 16:58:00 INFO None 5256742: status RUNNING/PENDING
2026-07-22 16:58:00 INFO None 5256747: status RUNNING/PENDING
2026-07-22 16:58:00 INFO None 5256756: status RUNNING/PENDING
2026-07-22 16:58:00 INFO None 5256761: status RUNNING/PENDING
2026-07-22 16:58:02 INFO None 5256765: status RUNNING/PENDING
2026-07-22 16:58:02 INFO None 5256770: status RUNNING/PENDING
2026-07-22 16:58:02 INFO None 5256776: status RUNNING/PENDING
2026-07-22 16:58:02 INFO None 5256782: status RUNNING/PENDING
2026-07-22 16:58:02 INFO None 5256788: status RUNNING/PENDING
2026-07-22 16:58:02 INFO None 5256793: status RUNNING/PENDING
2026-07-22 16:58:02 INFO None 5256798: status RUNNING/PENDING
2026-07-22 16:58:02 INFO None 5256803: status RUNNING/PENDING
2026-07-22 16:58:02 INFO Jobs still running: ['5256732', '5256735', '5256740', '5256742', '5256747', '5256756', '5256761', '5256765', '5256770', '5256776', '5256782', '5256788', '5256793', '5256798', '5256803']. Waiting...
2026-07-22 16:58:17 INFO None 5256732: status RUNNING/PENDING
2026-07-22 16:58:17 INFO None 5256735: status RUNNING/PENDING
2026-07-22 16:58:17 INFO None 5256740: status RUNNING/PENDING
2026-07-22 16:58:17 INFO None 5256742: status RUNNING/PENDING
2026-07-22 16:58:17 INFO None 5256747: status RUNNING/PENDING
2026-07-22 16:58:17 INFO None 5256756: status RUNNING/PENDING
2026-07-22 16:58:17 INFO None 5256761: status RUNNING/PENDING
2026-07-22 16:58:17 INFO None 5256765: status RUNNING/PENDING
2026-07-22 16:58:19 INFO None 5256770: status RUNNING/PENDING
2026-07-22 16:58:20 INFO None 5256776: status RUNNING/PENDING
2026-07-22 16:58:20 INFO None 5256782: status RUNNING/PENDING
2026-07-22 16:58:20 INFO None 5256788: status RUNNING/PENDING
2026-07-22 16:58:20 INFO None 5256793: status RUNNING/PENDING
2026-07-22 16:58:20 INFO None 5256798: status RUNNING/PENDING
2026-07-22 16:58:20 INFO None 5256803: status RUNNING/PENDING
2026-07-22 16:58:20 INFO Jobs still running: ['5256732', '5256735', '5256740', '5256742', '5256747', '5256756', '5256761', '5256765', '5256770', '5256776', '5256782', '5256788', '5256793', '5256798', '5256803']. Waiting...
2026-07-22 16:58:35 INFO None 5256732: status RUNNING/PENDING
2026-07-22 16:58:35 INFO None 5256735: status RUNNING/PENDING
2026-07-22 16:58:35 INFO None 5256740: status RUNNING/PENDING
2026-07-22 16:58:35 INFO None 5256742: status RUNNING/PENDING
2026-07-22 16:58:35 INFO None 5256747: status RUNNING/PENDING
2026-07-22 16:58:35 INFO None 5256756: status RUNNING/PENDING
2026-07-22 16:58:35 INFO None 5256761: status RUNNING/PENDING
2026-07-22 16:58:35 INFO None 5256765: status RUNNING/PENDING
2026-07-22 16:58:35 INFO None 5256770: status RUNNING/PENDING
2026-07-22 16:58:35 INFO None 5256776: status RUNNING/PENDING
2026-07-22 16:58:35 INFO None 5256782: status RUNNING/PENDING
2026-07-22 16:58:35 INFO None 5256788: status RUNNING/PENDING
2026-07-22 16:58:35 INFO None 5256793: status RUNNING/PENDING
2026-07-22 16:58:35 INFO None 5256798: status RUNNING/PENDING
2026-07-22 16:58:35 INFO None 5256803: status RUNNING/PENDING
2026-07-22 16:58:35 INFO Jobs still running: ['5256732', '5256735', '5256740', '5256742', '5256747', '5256756', '5256761', '5256765', '5256770', '5256776', '5256782', '5256788', '5256793', '5256798', '5256803']. Waiting...
2026-07-22 16:58:50 INFO None 5256732: status RUNNING/PENDING
2026-07-22 16:58:50 INFO None 5256735: status RUNNING/PENDING
2026-07-22 16:58:50 INFO None 5256740: status RUNNING/PENDING
2026-07-22 16:58:50 INFO None 5256742: status RUNNING/PENDING
2026-07-22 16:58:50 INFO None 5256747: status RUNNING/PENDING
2026-07-22 16:58:50 INFO None 5256756: status RUNNING/PENDING
2026-07-22 16:58:50 INFO None 5256761: status RUNNING/PENDING
2026-07-22 16:58:50 INFO None 5256765: status RUNNING/PENDING
2026-07-22 16:58:50 INFO None 5256770: status RUNNING/PENDING
2026-07-22 16:58:50 INFO None 5256776: status RUNNING/PENDING
2026-07-22 16:58:50 INFO None 5256782: status RUNNING/PENDING
2026-07-22 16:58:50 INFO None 5256788: status RUNNING/PENDING
2026-07-22 16:58:50 INFO None 5256793: status RUNNING/PENDING
2026-07-22 16:58:50 INFO None 5256798: status RUNNING/PENDING
2026-07-22 16:58:50 INFO None 5256803: status RUNNING/PENDING
2026-07-22 16:58:50 INFO Jobs still running: ['5256732', '5256735', '5256740', '5256742', '5256747', '5256756', '5256761', '5256765', '5256770', '5256776', '5256782', '5256788', '5256793', '5256798', '5256803']. Waiting...
2026-07-22 16:59:05 INFO None 5256732: status RUNNING/PENDING
2026-07-22 16:59:05 INFO None 5256735: status RUNNING/PENDING
2026-07-22 16:59:05 INFO None 5256740: status RUNNING/PENDING
2026-07-22 16:59:05 INFO None 5256742: status RUNNING/PENDING
2026-07-22 16:59:05 INFO None 5256747: status RUNNING/PENDING
2026-07-22 16:59:05 INFO None 5256756: status RUNNING/PENDING
2026-07-22 16:59:05 INFO None 5256761: status RUNNING/PENDING
2026-07-22 16:59:05 INFO None 5256765: status RUNNING/PENDING
2026-07-22 16:59:06 INFO None 5256770: status RUNNING/PENDING
2026-07-22 16:59:06 INFO None 5256776: status RUNNING/PENDING
2026-07-22 16:59:06 INFO None 5256782: status RUNNING/PENDING
2026-07-22 16:59:06 INFO None 5256788: status RUNNING/PENDING
2026-07-22 16:59:06 INFO None 5256793: status RUNNING/PENDING
2026-07-22 16:59:06 INFO None 5256798: status RUNNING/PENDING
2026-07-22 16:59:06 INFO None 5256803: status RUNNING/PENDING
2026-07-22 16:59:06 INFO Jobs still running: ['5256732', '5256735', '5256740', '5256742', '5256747', '5256756', '5256761', '5256765', '5256770', '5256776', '5256782', '5256788', '5256793', '5256798', '5256803']. Waiting...
2026-07-22 16:59:22 INFO None 5256732: status RUNNING/PENDING
2026-07-22 16:59:22 INFO None 5256735: status RUNNING/PENDING
2026-07-22 16:59:22 INFO None 5256740: status RUNNING/PENDING
2026-07-22 16:59:22 INFO None 5256742: status RUNNING/PENDING
2026-07-22 16:59:22 INFO None 5256747: status RUNNING/PENDING
2026-07-22 16:59:22 INFO None 5256756: status RUNNING/PENDING
2026-07-22 16:59:22 INFO None 5256761: status RUNNING/PENDING
2026-07-22 16:59:22 INFO None 5256765: status RUNNING/PENDING
2026-07-22 16:59:22 INFO None 5256770: status RUNNING/PENDING
2026-07-22 16:59:22 INFO None 5256776: status RUNNING/PENDING
2026-07-22 16:59:22 INFO None 5256782: status RUNNING/PENDING
2026-07-22 16:59:22 INFO None 5256788: status RUNNING/PENDING
2026-07-22 16:59:22 INFO None 5256793: status RUNNING/PENDING
2026-07-22 16:59:22 INFO None 5256798: status RUNNING/PENDING
2026-07-22 16:59:22 INFO None 5256803: status RUNNING/PENDING
2026-07-22 16:59:22 INFO Jobs still running: ['5256732', '5256735', '5256740', '5256742', '5256747', '5256756', '5256761', '5256765', '5256770', '5256776', '5256782', '5256788', '5256793', '5256798', '5256803']. Waiting...
2026-07-22 16:59:37 INFO None 5256732: status RUNNING/PENDING
2026-07-22 16:59:37 INFO None 5256735: status RUNNING/PENDING
2026-07-22 16:59:37 INFO None 5256740: status RUNNING/PENDING
2026-07-22 16:59:37 INFO None 5256742: status RUNNING/PENDING
2026-07-22 16:59:37 INFO None 5256747: status RUNNING/PENDING
2026-07-22 16:59:37 INFO None 5256756: status RUNNING/PENDING
2026-07-22 16:59:37 INFO None 5256761: status RUNNING/PENDING
2026-07-22 16:59:38 INFO None 5256765: status RUNNING/PENDING
2026-07-22 16:59:38 INFO None 5256770: status RUNNING/PENDING
2026-07-22 16:59:38 INFO None 5256776: status RUNNING/PENDING
2026-07-22 16:59:38 INFO None 5256782: status RUNNING/PENDING
2026-07-22 16:59:40 INFO None 5256788: status RUNNING/PENDING
2026-07-22 16:59:40 INFO None 5256793: status RUNNING/PENDING
2026-07-22 16:59:40 INFO None 5256798: status RUNNING/PENDING
2026-07-22 16:59:40 INFO None 5256803: status RUNNING/PENDING
2026-07-22 16:59:40 INFO Jobs still running: ['5256732', '5256735', '5256740', '5256742', '5256747', '5256756', '5256761', '5256765', '5256770', '5256776', '5256782', '5256788', '5256793', '5256798', '5256803']. Waiting...
2026-07-22 16:59:55 INFO None 5256732: status RUNNING/PENDING
2026-07-22 16:59:55 INFO None 5256735: status RUNNING/PENDING
2026-07-22 16:59:55 INFO None 5256740: status RUNNING/PENDING
2026-07-22 16:59:55 INFO None 5256742: status RUNNING/PENDING
2026-07-22 16:59:55 INFO None 5256747: status RUNNING/PENDING
2026-07-22 16:59:55 INFO None 5256756: status RUNNING/PENDING
2026-07-22 16:59:55 INFO None 5256761: status RUNNING/PENDING
2026-07-22 16:59:55 INFO None 5256765: status RUNNING/PENDING
2026-07-22 16:59:55 INFO None 5256770: status RUNNING/PENDING
2026-07-22 16:59:55 INFO None 5256776: status RUNNING/PENDING
2026-07-22 16:59:55 INFO None 5256782: status RUNNING/PENDING
2026-07-22 16:59:55 INFO None 5256788: status RUNNING/PENDING
2026-07-22 16:59:55 INFO None 5256793: status RUNNING/PENDING
2026-07-22 16:59:55 INFO None 5256798: status RUNNING/PENDING
2026-07-22 16:59:55 INFO None 5256803: status RUNNING/PENDING
2026-07-22 16:59:55 INFO Jobs still running: ['5256732', '5256735', '5256740', '5256742', '5256747', '5256756', '5256761', '5256765', '5256770', '5256776', '5256782', '5256788', '5256793', '5256798', '5256803']. Waiting...
2026-07-22 17:00:10 INFO None 5256732: status RUNNING/PENDING
2026-07-22 17:00:10 INFO None 5256735: status RUNNING/PENDING
2026-07-22 17:00:12 INFO None 5256740: status RUNNING/PENDING
2026-07-22 17:00:12 INFO None 5256742: status RUNNING/PENDING
2026-07-22 17:00:12 INFO None 5256747: status RUNNING/PENDING
2026-07-22 17:00:12 INFO None 5256756: status RUNNING/PENDING
2026-07-22 17:00:12 INFO None 5256761: status RUNNING/PENDING
2026-07-22 17:00:12 INFO None 5256765: status RUNNING/PENDING
2026-07-22 17:00:12 INFO None 5256770: status RUNNING/PENDING
2026-07-22 17:00:12 INFO None 5256776: status RUNNING/PENDING
2026-07-22 17:00:12 INFO None 5256782: status RUNNING/PENDING
2026-07-22 17:00:12 INFO None 5256788: status RUNNING/PENDING
2026-07-22 17:00:12 INFO None 5256793: status RUNNING/PENDING
2026-07-22 17:00:12 INFO None 5256798: status RUNNING/PENDING
2026-07-22 17:00:12 INFO None 5256803: status RUNNING/PENDING
2026-07-22 17:00:12 INFO Jobs still running: ['5256732', '5256735', '5256740', '5256742', '5256747', '5256756', '5256761', '5256765', '5256770', '5256776', '5256782', '5256788', '5256793', '5256798', '5256803']. Waiting...
2026-07-22 17:00:28 INFO None 5256732: status RUNNING/PENDING
2026-07-22 17:00:28 INFO None 5256735: status RUNNING/PENDING
2026-07-22 17:00:28 INFO None 5256740: status RUNNING/PENDING
2026-07-22 17:00:28 INFO None 5256742: status RUNNING/PENDING
2026-07-22 17:00:28 INFO None 5256747: status RUNNING/PENDING
2026-07-22 17:00:28 INFO None 5256756: status RUNNING/PENDING
2026-07-22 17:00:28 INFO None 5256761: status RUNNING/PENDING
2026-07-22 17:00:28 INFO None 5256765: status RUNNING/PENDING
2026-07-22 17:00:28 INFO None 5256770: status RUNNING/PENDING
2026-07-22 17:00:28 INFO None 5256776: status RUNNING/PENDING
2026-07-22 17:00:28 INFO None 5256782: status RUNNING/PENDING
2026-07-22 17:00:28 INFO None 5256788: status RUNNING/PENDING
2026-07-22 17:00:28 INFO None 5256793: status RUNNING/PENDING
2026-07-22 17:00:28 INFO None 5256798: status RUNNING/PENDING
2026-07-22 17:00:28 INFO None 5256803: status RUNNING/PENDING
2026-07-22 17:00:28 INFO Jobs still running: ['5256732', '5256735', '5256740', '5256742', '5256747', '5256756', '5256761', '5256765', '5256770', '5256776', '5256782', '5256788', '5256793', '5256798', '5256803']. Waiting...
2026-07-22 17:00:43 INFO None 5256732: status RUNNING/PENDING
2026-07-22 17:00:43 INFO None 5256735: status RUNNING/PENDING
2026-07-22 17:00:43 INFO None 5256740: status RUNNING/PENDING
2026-07-22 17:00:43 INFO None 5256742: status RUNNING/PENDING
2026-07-22 17:00:43 INFO None 5256747: status RUNNING/PENDING
2026-07-22 17:00:43 INFO None 5256756: status RUNNING/PENDING
2026-07-22 17:00:43 INFO None 5256761: status RUNNING/PENDING
2026-07-22 17:00:43 INFO None 5256765: status RUNNING/PENDING
2026-07-22 17:00:43 INFO None 5256770: status FINISHED
2026-07-22 17:00:43 INFO None 5256776: status RUNNING/PENDING
2026-07-22 17:00:43 INFO None 5256782: status FINISHED
2026-07-22 17:00:43 INFO None 5256788: status RUNNING/PENDING
2026-07-22 17:00:43 INFO None 5256793: status RUNNING/PENDING
2026-07-22 17:00:43 INFO None 5256798: status RUNNING/PENDING
2026-07-22 17:00:43 INFO None 5256803: status RUNNING/PENDING
2026-07-22 17:00:43 INFO Jobs still running: ['5256732', '5256735', '5256740', '5256742', '5256747', '5256756', '5256761', '5256765', '5256776', '5256788', '5256793', '5256798', '5256803']. Waiting...
2026-07-22 17:00:58 INFO None 5256732: status RUNNING/PENDING
2026-07-22 17:00:58 INFO None 5256735: status FINISHED
2026-07-22 17:00:58 INFO None 5256740: status RUNNING/PENDING
2026-07-22 17:00:58 INFO None 5256742: status RUNNING/PENDING
2026-07-22 17:00:58 INFO None 5256747: status RUNNING/PENDING
2026-07-22 17:00:58 INFO None 5256756: status RUNNING/PENDING
2026-07-22 17:00:58 INFO None 5256761: status RUNNING/PENDING
2026-07-22 17:00:58 INFO None 5256765: status RUNNING/PENDING
2026-07-22 17:00:58 INFO None 5256770: status FINISHED
2026-07-22 17:00:58 INFO None 5256776: status RUNNING/PENDING
2026-07-22 17:00:58 INFO None 5256782: status FINISHED
2026-07-22 17:00:58 INFO None 5256788: status RUNNING/PENDING
2026-07-22 17:00:58 INFO None 5256793: status RUNNING/PENDING
2026-07-22 17:00:58 INFO None 5256798: status RUNNING/PENDING
2026-07-22 17:00:58 INFO None 5256803: status RUNNING/PENDING
2026-07-22 17:00:58 INFO Jobs still running: ['5256732', '5256740', '5256742', '5256747', '5256756', '5256761', '5256765', '5256776', '5256788', '5256793', '5256798', '5256803']. Waiting...
2026-07-22 17:01:13 INFO None 5256732: status RUNNING/PENDING
2026-07-22 17:01:14 INFO None 5256735: status FINISHED
2026-07-22 17:01:16 INFO None 5256740: status RUNNING/PENDING
2026-07-22 17:01:16 INFO None 5256742: status RUNNING/PENDING
2026-07-22 17:01:16 INFO None 5256747: status FINISHED
2026-07-22 17:01:16 INFO None 5256756: status FINISHED
2026-07-22 17:01:16 INFO None 5256761: status RUNNING/PENDING
2026-07-22 17:01:16 INFO None 5256765: status FINISHED
2026-07-22 17:01:16 INFO None 5256770: status FINISHED
2026-07-22 17:01:16 INFO None 5256776: status RUNNING/PENDING
2026-07-22 17:01:16 INFO None 5256782: status FINISHED
2026-07-22 17:01:16 INFO None 5256788: status RUNNING/PENDING
2026-07-22 17:01:16 INFO None 5256793: status RUNNING/PENDING
2026-07-22 17:01:16 INFO None 5256798: status RUNNING/PENDING
2026-07-22 17:01:16 INFO None 5256803: status RUNNING/PENDING
2026-07-22 17:01:16 INFO Jobs still running: ['5256732', '5256740', '5256742', '5256761', '5256776', '5256788', '5256793', '5256798', '5256803']. Waiting...
2026-07-22 17:01:31 INFO None 5256732: status RUNNING/PENDING
2026-07-22 17:01:31 INFO None 5256735: status FINISHED
2026-07-22 17:01:31 INFO None 5256740: status FINISHED
2026-07-22 17:01:31 INFO None 5256742: status FINISHED
2026-07-22 17:01:31 INFO None 5256747: status FINISHED
2026-07-22 17:01:31 INFO None 5256756: status FINISHED
2026-07-22 17:01:31 INFO None 5256761: status FINISHED
2026-07-22 17:01:31 INFO None 5256765: status FINISHED
2026-07-22 17:01:31 INFO None 5256770: status FINISHED
2026-07-22 17:01:31 INFO None 5256776: status FINISHED
2026-07-22 17:01:31 INFO None 5256782: status FINISHED
2026-07-22 17:01:31 INFO None 5256788: status RUNNING/PENDING
2026-07-22 17:01:31 INFO None 5256793: status RUNNING/PENDING
2026-07-22 17:01:31 INFO None 5256798: status RUNNING/PENDING
2026-07-22 17:01:31 INFO None 5256803: status RUNNING/PENDING
2026-07-22 17:01:31 INFO Jobs still running: ['5256732', '5256788', '5256793', '5256798', '5256803']. Waiting...
2026-07-22 17:01:46 INFO None 5256732: status FINISHED
2026-07-22 17:01:46 INFO None 5256735: status FINISHED
2026-07-22 17:01:46 INFO None 5256740: status FINISHED
2026-07-22 17:01:48 INFO None 5256742: status FINISHED
2026-07-22 17:01:48 INFO None 5256747: status FINISHED
2026-07-22 17:01:48 INFO None 5256756: status FINISHED
2026-07-22 17:01:48 INFO None 5256761: status FINISHED
2026-07-22 17:01:48 INFO None 5256765: status FINISHED
2026-07-22 17:01:48 INFO None 5256770: status FINISHED
2026-07-22 17:01:48 INFO None 5256776: status FINISHED
2026-07-22 17:01:48 INFO None 5256782: status FINISHED
2026-07-22 17:01:48 INFO None 5256788: status RUNNING/PENDING
2026-07-22 17:01:48 INFO None 5256793: status RUNNING/PENDING
2026-07-22 17:01:48 INFO None 5256798: status RUNNING/PENDING
2026-07-22 17:01:49 INFO None 5256803: status RUNNING/PENDING
2026-07-22 17:01:49 INFO Jobs still running: ['5256788', '5256793', '5256798', '5256803']. Waiting...
2026-07-22 17:02:04 INFO None 5256732: status FINISHED
2026-07-22 17:02:04 INFO None 5256735: status FINISHED
2026-07-22 17:02:04 INFO None 5256740: status FINISHED
2026-07-22 17:02:04 INFO None 5256742: status FINISHED
2026-07-22 17:02:04 INFO None 5256747: status FINISHED
2026-07-22 17:02:04 INFO None 5256756: status FINISHED
2026-07-22 17:02:04 INFO None 5256761: status FINISHED
2026-07-22 17:02:04 INFO None 5256765: status FINISHED
2026-07-22 17:02:04 INFO None 5256770: status FINISHED
2026-07-22 17:02:04 INFO None 5256776: status FINISHED
2026-07-22 17:02:04 INFO None 5256782: status FINISHED
2026-07-22 17:02:04 INFO None 5256788: status RUNNING/PENDING
2026-07-22 17:02:04 INFO None 5256793: status RUNNING/PENDING
2026-07-22 17:02:04 INFO None 5256798: status RUNNING/PENDING
2026-07-22 17:02:04 INFO None 5256803: status RUNNING/PENDING
2026-07-22 17:02:04 INFO Jobs still running: ['5256788', '5256793', '5256798', '5256803']. Waiting...
2026-07-22 17:02:19 INFO None 5256732: status FINISHED
2026-07-22 17:02:19 INFO None 5256735: status FINISHED
2026-07-22 17:02:19 INFO None 5256740: status FINISHED
2026-07-22 17:02:19 INFO None 5256742: status FINISHED
2026-07-22 17:02:19 INFO None 5256747: status FINISHED
2026-07-22 17:02:19 INFO None 5256756: status FINISHED
2026-07-22 17:02:19 INFO None 5256761: status FINISHED
2026-07-22 17:02:19 INFO None 5256765: status FINISHED
2026-07-22 17:02:19 INFO None 5256770: status FINISHED
2026-07-22 17:02:19 INFO None 5256776: status FINISHED
2026-07-22 17:02:19 INFO None 5256782: status FINISHED
2026-07-22 17:02:19 INFO None 5256788: status RUNNING/PENDING
2026-07-22 17:02:19 INFO None 5256793: status RUNNING/PENDING
2026-07-22 17:02:19 INFO None 5256798: status RUNNING/PENDING
2026-07-22 17:02:19 INFO None 5256803: status RUNNING/PENDING
2026-07-22 17:02:19 INFO Jobs still running: ['5256788', '5256793', '5256798', '5256803']. Waiting...
2026-07-22 17:02:34 INFO None 5256732: status FINISHED
2026-07-22 17:02:34 INFO None 5256735: status FINISHED
2026-07-22 17:02:34 INFO None 5256740: status FINISHED
2026-07-22 17:02:34 INFO None 5256742: status FINISHED
2026-07-22 17:02:34 INFO None 5256747: status FINISHED
2026-07-22 17:02:34 INFO None 5256756: status FINISHED
2026-07-22 17:02:34 INFO None 5256761: status FINISHED
2026-07-22 17:02:34 INFO None 5256765: status FINISHED
2026-07-22 17:02:34 INFO None 5256770: status FINISHED
2026-07-22 17:02:34 INFO None 5256776: status FINISHED
2026-07-22 17:02:34 INFO None 5256782: status FINISHED
2026-07-22 17:02:34 INFO None 5256788: status RUNNING/PENDING
2026-07-22 17:02:34 INFO None 5256793: status RUNNING/PENDING
2026-07-22 17:02:34 INFO None 5256798: status RUNNING/PENDING
2026-07-22 17:02:34 INFO None 5256803: status RUNNING/PENDING
2026-07-22 17:02:34 INFO Jobs still running: ['5256788', '5256793', '5256798', '5256803']. Waiting...
2026-07-22 17:02:50 INFO None 5256732: status FINISHED
2026-07-22 17:02:50 INFO None 5256735: status FINISHED
2026-07-22 17:02:50 INFO None 5256740: status FINISHED
2026-07-22 17:02:50 INFO None 5256742: status FINISHED
2026-07-22 17:02:50 INFO None 5256747: status FINISHED
2026-07-22 17:02:50 INFO None 5256756: status FINISHED
2026-07-22 17:02:50 INFO None 5256761: status FINISHED
2026-07-22 17:02:50 INFO None 5256765: status FINISHED
2026-07-22 17:02:50 INFO None 5256770: status FINISHED
2026-07-22 17:02:50 INFO None 5256776: status FINISHED
2026-07-22 17:02:50 INFO None 5256782: status FINISHED
2026-07-22 17:02:50 INFO None 5256788: status RUNNING/PENDING
2026-07-22 17:02:50 INFO None 5256793: status RUNNING/PENDING
2026-07-22 17:02:50 INFO None 5256798: status RUNNING/PENDING
2026-07-22 17:02:50 INFO None 5256803: status FINISHED
2026-07-22 17:02:50 INFO Jobs still running: ['5256788', '5256793', '5256798']. Waiting...
2026-07-22 17:03:05 INFO None 5256732: status FINISHED
2026-07-22 17:03:05 INFO None 5256735: status FINISHED
2026-07-22 17:03:05 INFO None 5256740: status FINISHED
2026-07-22 17:03:05 INFO None 5256742: status FINISHED
2026-07-22 17:03:05 INFO None 5256747: status FINISHED
2026-07-22 17:03:05 INFO None 5256756: status FINISHED
2026-07-22 17:03:05 INFO None 5256761: status FINISHED
2026-07-22 17:03:05 INFO None 5256765: status FINISHED
2026-07-22 17:03:05 INFO None 5256770: status FINISHED
2026-07-22 17:03:05 INFO None 5256776: status FINISHED
2026-07-22 17:03:05 INFO None 5256782: status FINISHED
2026-07-22 17:03:05 INFO None 5256788: status RUNNING/PENDING
2026-07-22 17:03:08 INFO None 5256793: status RUNNING/PENDING
2026-07-22 17:03:08 INFO None 5256798: status RUNNING/PENDING
2026-07-22 17:03:08 INFO None 5256803: status FINISHED
2026-07-22 17:03:08 INFO Jobs still running: ['5256788', '5256793', '5256798']. Waiting...
2026-07-22 17:03:23 INFO None 5256732: status FINISHED
2026-07-22 17:03:23 INFO None 5256735: status FINISHED
2026-07-22 17:03:23 INFO None 5256740: status FINISHED
2026-07-22 17:03:23 INFO None 5256742: status FINISHED
2026-07-22 17:03:23 INFO None 5256747: status FINISHED
2026-07-22 17:03:23 INFO None 5256756: status FINISHED
2026-07-22 17:03:23 INFO None 5256761: status FINISHED
2026-07-22 17:03:23 INFO None 5256765: status FINISHED
2026-07-22 17:03:23 INFO None 5256770: status FINISHED
2026-07-22 17:03:23 INFO None 5256776: status FINISHED
2026-07-22 17:03:23 INFO None 5256782: status FINISHED
2026-07-22 17:03:23 INFO None 5256788: status RUNNING/PENDING
2026-07-22 17:03:23 INFO None 5256793: status RUNNING/PENDING
2026-07-22 17:03:23 INFO None 5256798: status RUNNING/PENDING
2026-07-22 17:03:23 INFO None 5256803: status FINISHED
2026-07-22 17:03:23 INFO Jobs still running: ['5256788', '5256793', '5256798']. Waiting...
2026-07-22 17:03:38 INFO None 5256732: status FINISHED
2026-07-22 17:03:38 INFO None 5256735: status FINISHED
2026-07-22 17:03:38 INFO None 5256740: status FINISHED
2026-07-22 17:03:38 INFO None 5256742: status FINISHED
2026-07-22 17:03:38 INFO None 5256747: status FINISHED
2026-07-22 17:03:38 INFO None 5256756: status FINISHED
2026-07-22 17:03:38 INFO None 5256761: status FINISHED
2026-07-22 17:03:38 INFO None 5256765: status FINISHED
2026-07-22 17:03:38 INFO None 5256770: status FINISHED
2026-07-22 17:03:40 INFO None 5256776: status FINISHED
2026-07-22 17:03:40 INFO None 5256782: status FINISHED
2026-07-22 17:03:40 INFO None 5256788: status RUNNING/PENDING
2026-07-22 17:03:40 INFO None 5256793: status RUNNING/PENDING
2026-07-22 17:03:40 INFO None 5256798: status FINISHED
2026-07-22 17:03:40 INFO None 5256803: status FINISHED
2026-07-22 17:03:40 INFO Jobs still running: ['5256788', '5256793']. Waiting...
2026-07-22 17:03:55 INFO None 5256732: status FINISHED
2026-07-22 17:03:55 INFO None 5256735: status FINISHED
2026-07-22 17:03:55 INFO None 5256740: status FINISHED
2026-07-22 17:03:55 INFO None 5256742: status FINISHED
2026-07-22 17:03:55 INFO None 5256747: status FINISHED
2026-07-22 17:03:55 INFO None 5256756: status FINISHED
2026-07-22 17:03:55 INFO None 5256761: status FINISHED
2026-07-22 17:03:55 INFO None 5256765: status FINISHED
2026-07-22 17:03:55 INFO None 5256770: status FINISHED
2026-07-22 17:03:55 INFO None 5256776: status FINISHED
2026-07-22 17:03:55 INFO None 5256782: status FINISHED
2026-07-22 17:03:55 INFO None 5256788: status FINISHED
2026-07-22 17:03:55 INFO None 5256793: status FINISHED
2026-07-22 17:03:55 INFO None 5256798: status FINISHED
2026-07-22 17:03:55 INFO None 5256803: status FINISHED
2026-07-22 17:03:55 INFO Jobs ['5256732', '5256735', '5256740', '5256742', '5256747', '5256756', '5256761', '5256765', '5256770', '5256776', '5256782', '5256788', '5256793', '5256798', '5256803'] have finished
2026-07-22 17:03:55 INFO Checking restart files were created ...
2026-07-22 17:03:55 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc(1002685915 bytes)
2026-07-22 17:03:55 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc(1002685915 bytes)
2026-07-22 17:03:55 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc(1002685915 bytes)
2026-07-22 17:03:55 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc(1002685915 bytes)
2026-07-22 17:03:55 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc(1002685915 bytes)
2026-07-22 17:03:55 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc(1002685915 bytes)
2026-07-22 17:03:55 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc(1002685915 bytes)
2026-07-22 17:03:55 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc(1002685915 bytes)
2026-07-22 17:03:55 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc(1002685915 bytes)
2026-07-22 17:03:55 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc(1002685915 bytes)
2026-07-22 17:03:55 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc(1002685915 bytes)
2026-07-22 17:03:56 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc(1002685915 bytes)
2026-07-22 17:03:56 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc(1002685915 bytes)
2026-07-22 17:03:56 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc(1002685915 bytes)
2026-07-22 17:03:56 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc(1002685915 bytes)
2026-07-22 17:03:56 INFO  Run_model() completed successfully.
2026-07-22 17:03:56 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 17:03:56 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 11:00:00 days=153072 seconds=39600
2026-07-22 17:03:56 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-22 17:03:56 INFO [TIME] increment current_time 2020-02-06 09:00:00 -> 2020-02-06 11:00:00
2026-07-22 17:03:56 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 17:03:56 INFO ---------->>> Running process_satellite_data()
2026-07-22 17:03:56 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12003.nc
2026-07-22 17:03:56 INFO ---------->>> Running run_obs_converter()
2026-07-22 17:03:56 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_39045_153072.out
2026-07-22 17:03:56 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_39045_153072.out
2026-07-22 17:03:56 INFO ---------->>> Running DART
2026-07-22 17:03:56 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-22 17:03:56 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/chim_ENS1_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/chim_ENS1_2020020611_1_out_toDART.nc
2026-07-22 17:03:56 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/chim_ENS2_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/chim_ENS2_2020020611_1_out_toDART.nc
2026-07-22 17:03:56 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/chim_ENS3_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/chim_ENS3_2020020611_1_out_toDART.nc
2026-07-22 17:03:57 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/chim_ENS4_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/chim_ENS4_2020020611_1_out_toDART.nc
2026-07-22 17:03:57 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/chim_ENS5_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/chim_ENS5_2020020611_1_out_toDART.nc
2026-07-22 17:03:57 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/chim_ENS6_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/chim_ENS6_2020020611_1_out_toDART.nc
2026-07-22 17:03:58 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/chim_ENS7_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/chim_ENS7_2020020611_1_out_toDART.nc
2026-07-22 17:03:58 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/chim_ENS8_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/chim_ENS8_2020020611_1_out_toDART.nc
2026-07-22 17:03:58 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/chim_ENS9_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/chim_ENS9_2020020611_1_out_toDART.nc
2026-07-22 17:03:59 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/chim_ENS10_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/chim_ENS10_2020020611_1_out_toDART.nc
2026-07-22 17:03:59 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/chim_ENS11_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/chim_ENS11_2020020611_1_out_toDART.nc
2026-07-22 17:03:59 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/chim_ENS12_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/chim_ENS12_2020020611_1_out_toDART.nc
2026-07-22 17:03:59 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/chim_ENS13_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/chim_ENS13_2020020611_1_out_toDART.nc
2026-07-22 17:04:00 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/chim_ENS14_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/chim_ENS14_2020020611_1_out_toDART.nc
2026-07-22 17:04:00 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/chim_ENS15_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/chim_ENS15_2020020611_1_out_toDART.nc
2026-07-22 17:04:00 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-22 17:04:00 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-22 17:04:00 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-22 17:04:00 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-22 17:04:00 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-22 17:04:00 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-22 17:04:13 INFO Found: []
2026-07-22 17:04:13 INFO No job id returned by command ./run_filter.bsh
2026-07-22 17:04:13 INFO No monitoring will be performed
2026-07-22 17:04:13 INFO Moving DART output files to analysis and preassim directories for date 2020020611 if present ...
2026-07-22 17:04:13 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 17:04:13 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 17:04:13 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 17:04:13 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 17:04:13 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 17:04:13 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 17:04:13 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 17:04:13 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 17:04:13 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 17:04:13 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 17:04:13 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 17:04:13 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 17:04:13 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 17:04:13 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 17:04:13 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 17:04:13 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 17:04:13 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 17:04:13 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 17:04:13 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 17:04:13 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 17:04:14 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 17:04:14 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 17:04:14 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 17:04:14 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 17:04:14 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 17:04:14 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 17:04:14 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 17:04:14 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 17:04:14 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 17:04:14 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 17:04:14 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 17:04:14 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020611'
2026-07-22 17:04:14 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 17:04:14 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020611'
2026-07-22 17:04:14 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-22 17:04:14 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-22 17:04:14 INFO run_dart() is DONE.
2026-07-22 17:04:14 INFO ---------->>> Running update_pollutant_in_end()
2026-07-22 17:04:14 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:04:14 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS1.nc
2026-07-22 17:04:15 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:04:15 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS1.nc
2026-07-22 17:04:15 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:04:16 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:04:16 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS1.nc
2026-07-22 17:04:16 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:04:17 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS1.nc
2026-07-22 17:04:17 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:04:17 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:04:17 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:04:17 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS2.nc
2026-07-22 17:04:18 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:04:18 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS2.nc
2026-07-22 17:04:18 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:04:19 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:04:19 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS2.nc
2026-07-22 17:04:19 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:04:20 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS2.nc
2026-07-22 17:04:20 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:04:20 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:04:20 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:04:20 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS3.nc
2026-07-22 17:04:21 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:04:21 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS3.nc
2026-07-22 17:04:21 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:04:22 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:04:22 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS3.nc
2026-07-22 17:04:22 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:04:23 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS3.nc
2026-07-22 17:04:23 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:04:23 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:04:23 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:04:23 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS4.nc
2026-07-22 17:04:24 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:04:24 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS4.nc
2026-07-22 17:04:24 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:04:25 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:04:25 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS4.nc
2026-07-22 17:04:25 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:04:26 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS4.nc
2026-07-22 17:04:26 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:04:26 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:04:26 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:04:26 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS5.nc
2026-07-22 17:04:27 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:04:28 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS5.nc
2026-07-22 17:04:28 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:04:28 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:04:28 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS5.nc
2026-07-22 17:04:28 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:04:29 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS5.nc
2026-07-22 17:04:29 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:04:29 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:04:29 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:04:29 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS6.nc
2026-07-22 17:04:30 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:04:31 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS6.nc
2026-07-22 17:04:31 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:04:31 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:04:31 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS6.nc
2026-07-22 17:04:32 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:04:32 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS6.nc
2026-07-22 17:04:32 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:04:32 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:04:32 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:04:32 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS7.nc
2026-07-22 17:04:33 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:04:34 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS7.nc
2026-07-22 17:04:34 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:04:34 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:04:34 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS7.nc
2026-07-22 17:04:35 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:04:35 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS7.nc
2026-07-22 17:04:35 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:04:35 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:04:36 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:04:36 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS8.nc
2026-07-22 17:04:36 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:04:37 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS8.nc
2026-07-22 17:04:37 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:04:37 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:04:37 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS8.nc
2026-07-22 17:04:38 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:04:38 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS8.nc
2026-07-22 17:04:38 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:04:38 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:04:39 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:04:39 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS9.nc
2026-07-22 17:04:39 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:04:40 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS9.nc
2026-07-22 17:04:40 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:04:40 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:04:40 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS9.nc
2026-07-22 17:04:41 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:04:41 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS9.nc
2026-07-22 17:04:41 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:04:41 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:04:42 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:04:42 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS10.nc
2026-07-22 17:04:42 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:04:43 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS10.nc
2026-07-22 17:04:43 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:04:43 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:04:43 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS10.nc
2026-07-22 17:04:44 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:04:44 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS10.nc
2026-07-22 17:04:44 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:04:44 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:04:45 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:04:45 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS11.nc
2026-07-22 17:04:46 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:04:46 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS11.nc
2026-07-22 17:04:46 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:04:47 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:04:47 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS11.nc
2026-07-22 17:04:47 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:04:48 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS11.nc
2026-07-22 17:04:48 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:04:48 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:04:48 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:04:48 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS12.nc
2026-07-22 17:04:49 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:04:49 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS12.nc
2026-07-22 17:04:49 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:04:50 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:04:50 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS12.nc
2026-07-22 17:04:50 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:04:51 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS12.nc
2026-07-22 17:04:51 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:04:51 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:04:51 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:04:51 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS13.nc
2026-07-22 17:04:52 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:04:52 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS13.nc
2026-07-22 17:04:52 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:04:53 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:04:53 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS13.nc
2026-07-22 17:04:53 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:04:54 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS13.nc
2026-07-22 17:04:54 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:04:54 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:04:54 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:04:54 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS14.nc
2026-07-22 17:04:55 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:04:56 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS14.nc
2026-07-22 17:04:56 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:04:56 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:04:56 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS14.nc
2026-07-22 17:04:56 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:04:57 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS14.nc
2026-07-22 17:04:57 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:04:57 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:04:57 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:04:57 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS15.nc
2026-07-22 17:04:58 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:04:59 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS15.nc
2026-07-22 17:04:59 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:04:59 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:04:59 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS15.nc
2026-07-22 17:05:00 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:05:00 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS15.nc
2026-07-22 17:05:00 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:05:00 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:05:00 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-22 17:05:00 INFO [TIME] step_end current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 17:05:00 INFO [TIME] step_start current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 17:05:00 INFO [TIME] window start=2020-02-06 11:00:00 end=2020-02-06 13:00:00 run_hours=2 has_assimilation=True
2026-07-22 17:05:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:05:02 INFO Hourly dataset computed and listing created
2026-07-22 17:05:06 INFO Hourly dataset computed
2026-07-22 17:05:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:05:07 INFO Hourly dataset computed and listing created
2026-07-22 17:05:08 INFO Hourly dataset computed
2026-07-22 17:05:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:05:08 INFO Hourly dataset computed and listing created
2026-07-22 17:05:09 INFO Hourly dataset computed
2026-07-22 17:05:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:05:10 INFO Hourly dataset computed and listing created
2026-07-22 17:05:11 INFO Hourly dataset computed
2026-07-22 17:05:11 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:05:12 INFO Hourly dataset computed and listing created
2026-07-22 17:05:13 INFO Hourly dataset computed
2026-07-22 17:05:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:05:14 INFO Hourly dataset computed and listing created
2026-07-22 17:05:14 INFO Hourly dataset computed
2026-07-22 17:05:14 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:05:15 INFO Hourly dataset computed and listing created
2026-07-22 17:05:16 INFO Hourly dataset computed
2026-07-22 17:05:16 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:05:17 INFO Hourly dataset computed and listing created
2026-07-22 17:05:18 INFO Hourly dataset computed
2026-07-22 17:05:18 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:05:19 INFO Hourly dataset computed and listing created
2026-07-22 17:05:19 INFO Hourly dataset computed
2026-07-22 17:05:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:05:20 INFO Hourly dataset computed and listing created
2026-07-22 17:05:21 INFO Hourly dataset computed
2026-07-22 17:05:21 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:05:22 INFO Hourly dataset computed and listing created
2026-07-22 17:05:23 INFO Hourly dataset computed
2026-07-22 17:05:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:05:24 INFO Hourly dataset computed and listing created
2026-07-22 17:05:24 INFO Hourly dataset computed
2026-07-22 17:05:24 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:05:25 INFO Hourly dataset computed and listing created
2026-07-22 17:05:26 INFO Hourly dataset computed
2026-07-22 17:05:26 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:05:27 INFO Hourly dataset computed and listing created
2026-07-22 17:05:28 INFO Hourly dataset computed
2026-07-22 17:05:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:05:29 INFO Hourly dataset computed and listing created
2026-07-22 17:05:30 INFO Hourly dataset computed
2026-07-22 17:05:30 INFO ---------->>> Running CHIMERE model from 2020-02-06 11:00:00 to 2020-02-06 13:00:00
2026-07-22 17:05:30 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:05:30 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1
2026-07-22 17:05:30 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc
2026-07-22 17:05:30 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-22 17:05:30 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:05:30 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-22 17:05:30 INFO Queuing job for member 1...
2026-07-22 17:05:30 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:05:30 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-22 17:05:31 INFO Found: ['5256891']
2026-07-22 17:05:36 INFO [TGCC-IRENE] Submitted job with ID:['5256891']
2026-07-22 17:05:36 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:05:36 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2
2026-07-22 17:05:36 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc
2026-07-22 17:05:36 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-22 17:05:36 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:05:36 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-22 17:05:36 INFO Queuing job for member 2...
2026-07-22 17:05:36 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:05:36 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-22 17:05:39 INFO Found: ['5256894']
2026-07-22 17:05:44 INFO [TGCC-IRENE] Submitted job with ID:['5256894']
2026-07-22 17:05:44 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:05:44 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3
2026-07-22 17:05:44 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc
2026-07-22 17:05:44 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-22 17:05:44 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:05:44 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-22 17:05:44 INFO Queuing job for member 3...
2026-07-22 17:05:44 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:05:44 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-22 17:05:45 INFO Found: ['5256895']
2026-07-22 17:05:50 INFO [TGCC-IRENE] Submitted job with ID:['5256895']
2026-07-22 17:05:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:05:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4
2026-07-22 17:05:50 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc
2026-07-22 17:05:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-22 17:05:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:05:50 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-22 17:05:50 INFO Queuing job for member 4...
2026-07-22 17:05:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:05:50 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-22 17:05:50 INFO Found: ['5256896']
2026-07-22 17:05:55 INFO [TGCC-IRENE] Submitted job with ID:['5256896']
2026-07-22 17:05:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:05:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5
2026-07-22 17:05:55 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc
2026-07-22 17:05:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-22 17:05:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:05:55 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-22 17:05:55 INFO Queuing job for member 5...
2026-07-22 17:05:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:05:55 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-22 17:05:56 INFO Found: ['5256897']
2026-07-22 17:06:01 INFO [TGCC-IRENE] Submitted job with ID:['5256897']
2026-07-22 17:06:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:06:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6
2026-07-22 17:06:01 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc
2026-07-22 17:06:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-22 17:06:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:06:01 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-22 17:06:01 INFO Queuing job for member 6...
2026-07-22 17:06:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:06:01 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-22 17:06:02 INFO Found: ['5256899']
2026-07-22 17:06:07 INFO [TGCC-IRENE] Submitted job with ID:['5256899']
2026-07-22 17:06:07 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:06:07 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7
2026-07-22 17:06:07 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc
2026-07-22 17:06:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-22 17:06:07 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:06:07 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-22 17:06:07 INFO Queuing job for member 7...
2026-07-22 17:06:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:06:07 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-22 17:06:08 INFO Found: ['5256902']
2026-07-22 17:06:13 INFO [TGCC-IRENE] Submitted job with ID:['5256902']
2026-07-22 17:06:13 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:06:13 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8
2026-07-22 17:06:13 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc
2026-07-22 17:06:13 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-22 17:06:13 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:06:13 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-22 17:06:13 INFO Queuing job for member 8...
2026-07-22 17:06:13 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:06:13 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-22 17:06:14 INFO Found: ['5256904']
2026-07-22 17:06:19 INFO [TGCC-IRENE] Submitted job with ID:['5256904']
2026-07-22 17:06:19 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:06:19 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9
2026-07-22 17:06:19 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc
2026-07-22 17:06:19 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-22 17:06:19 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:06:19 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-22 17:06:19 INFO Queuing job for member 9...
2026-07-22 17:06:19 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:06:19 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-22 17:06:19 INFO Found: ['5256905']
2026-07-22 17:06:24 INFO [TGCC-IRENE] Submitted job with ID:['5256905']
2026-07-22 17:06:24 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:06:24 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10
2026-07-22 17:06:24 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc
2026-07-22 17:06:24 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-22 17:06:24 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:06:24 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-22 17:06:24 INFO Queuing job for member 10...
2026-07-22 17:06:24 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:06:24 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-22 17:06:25 INFO Found: ['5256908']
2026-07-22 17:06:30 INFO [TGCC-IRENE] Submitted job with ID:['5256908']
2026-07-22 17:06:30 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:06:30 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11
2026-07-22 17:06:30 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc
2026-07-22 17:06:30 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-22 17:06:30 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:06:30 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-22 17:06:30 INFO Queuing job for member 11...
2026-07-22 17:06:30 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:06:30 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-22 17:06:32 INFO Found: ['5256909']
2026-07-22 17:06:37 INFO [TGCC-IRENE] Submitted job with ID:['5256909']
2026-07-22 17:06:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:06:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12
2026-07-22 17:06:37 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc
2026-07-22 17:06:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-22 17:06:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:06:37 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-22 17:06:37 INFO Queuing job for member 12...
2026-07-22 17:06:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:06:37 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-22 17:06:40 INFO Found: ['5256911']
2026-07-22 17:06:45 INFO [TGCC-IRENE] Submitted job with ID:['5256911']
2026-07-22 17:06:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:06:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13
2026-07-22 17:06:45 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc
2026-07-22 17:06:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-22 17:06:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:06:45 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-22 17:06:45 INFO Queuing job for member 13...
2026-07-22 17:06:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:06:45 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-22 17:06:47 INFO Found: ['5256913']
2026-07-22 17:06:52 INFO [TGCC-IRENE] Submitted job with ID:['5256913']
2026-07-22 17:06:52 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:06:52 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14
2026-07-22 17:06:52 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc
2026-07-22 17:06:52 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-22 17:06:52 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:06:52 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-22 17:06:52 INFO Queuing job for member 14...
2026-07-22 17:06:52 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:06:52 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-22 17:06:55 INFO Found: ['5256914']
2026-07-22 17:07:00 INFO [TGCC-IRENE] Submitted job with ID:['5256914']
2026-07-22 17:07:00 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:07:00 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15
2026-07-22 17:07:00 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc
2026-07-22 17:07:00 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-22 17:07:00 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:07:00 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-22 17:07:00 INFO Queuing job for member 15...
2026-07-22 17:07:00 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:07:00 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-22 17:07:02 INFO Found: ['5256916']
2026-07-22 17:07:07 INFO [TGCC-IRENE] Submitted job with ID:['5256916']
2026-07-22 17:07:07 INFO Checking job status ...
2026-07-22 17:07:07 INFO None 5256891: status RUNNING/PENDING
2026-07-22 17:07:07 INFO None 5256894: status RUNNING/PENDING
2026-07-22 17:07:07 INFO None 5256895: status RUNNING/PENDING
2026-07-22 17:07:07 INFO None 5256896: status RUNNING/PENDING
2026-07-22 17:07:07 INFO None 5256897: status RUNNING/PENDING
2026-07-22 17:07:07 INFO None 5256899: status RUNNING/PENDING
2026-07-22 17:07:07 INFO None 5256902: status RUNNING/PENDING
2026-07-22 17:07:07 INFO None 5256904: status RUNNING/PENDING
2026-07-22 17:07:07 INFO None 5256905: status RUNNING/PENDING
2026-07-22 17:07:07 INFO None 5256908: status RUNNING/PENDING
2026-07-22 17:07:07 INFO None 5256909: status RUNNING/PENDING
2026-07-22 17:07:07 INFO None 5256911: status RUNNING/PENDING
2026-07-22 17:07:07 INFO None 5256913: status RUNNING/PENDING
2026-07-22 17:07:08 INFO None 5256914: status RUNNING/PENDING
2026-07-22 17:07:08 INFO None 5256916: status RUNNING/PENDING
2026-07-22 17:07:08 INFO Jobs still running: ['5256891', '5256894', '5256895', '5256896', '5256897', '5256899', '5256902', '5256904', '5256905', '5256908', '5256909', '5256911', '5256913', '5256914', '5256916']. Waiting...
2026-07-22 17:07:23 INFO None 5256891: status RUNNING/PENDING
2026-07-22 17:07:23 INFO None 5256894: status RUNNING/PENDING
2026-07-22 17:07:23 INFO None 5256895: status RUNNING/PENDING
2026-07-22 17:07:23 INFO None 5256896: status RUNNING/PENDING
2026-07-22 17:07:23 INFO None 5256897: status RUNNING/PENDING
2026-07-22 17:07:23 INFO None 5256899: status RUNNING/PENDING
2026-07-22 17:07:23 INFO None 5256902: status RUNNING/PENDING
2026-07-22 17:07:23 INFO None 5256904: status RUNNING/PENDING
2026-07-22 17:07:23 INFO None 5256905: status RUNNING/PENDING
2026-07-22 17:07:23 INFO None 5256908: status RUNNING/PENDING
2026-07-22 17:07:25 INFO None 5256909: status RUNNING/PENDING
2026-07-22 17:07:25 INFO None 5256911: status RUNNING/PENDING
2026-07-22 17:07:25 INFO None 5256913: status RUNNING/PENDING
2026-07-22 17:07:25 INFO None 5256914: status RUNNING/PENDING
2026-07-22 17:07:25 INFO None 5256916: status RUNNING/PENDING
2026-07-22 17:07:25 INFO Jobs still running: ['5256891', '5256894', '5256895', '5256896', '5256897', '5256899', '5256902', '5256904', '5256905', '5256908', '5256909', '5256911', '5256913', '5256914', '5256916']. Waiting...
2026-07-22 17:07:40 INFO None 5256891: status RUNNING/PENDING
2026-07-22 17:07:40 INFO None 5256894: status RUNNING/PENDING
2026-07-22 17:07:40 INFO None 5256895: status RUNNING/PENDING
2026-07-22 17:07:40 INFO None 5256896: status RUNNING/PENDING
2026-07-22 17:07:40 INFO None 5256897: status RUNNING/PENDING
2026-07-22 17:07:40 INFO None 5256899: status RUNNING/PENDING
2026-07-22 17:07:40 INFO None 5256902: status RUNNING/PENDING
2026-07-22 17:07:40 INFO None 5256904: status RUNNING/PENDING
2026-07-22 17:07:40 INFO None 5256905: status RUNNING/PENDING
2026-07-22 17:07:40 INFO None 5256908: status RUNNING/PENDING
2026-07-22 17:07:40 INFO None 5256909: status RUNNING/PENDING
2026-07-22 17:07:40 INFO None 5256911: status RUNNING/PENDING
2026-07-22 17:07:40 INFO None 5256913: status RUNNING/PENDING
2026-07-22 17:07:40 INFO None 5256914: status RUNNING/PENDING
2026-07-22 17:07:40 INFO None 5256916: status RUNNING/PENDING
2026-07-22 17:07:40 INFO Jobs still running: ['5256891', '5256894', '5256895', '5256896', '5256897', '5256899', '5256902', '5256904', '5256905', '5256908', '5256909', '5256911', '5256913', '5256914', '5256916']. Waiting...
2026-07-22 17:07:55 INFO None 5256891: status RUNNING/PENDING
2026-07-22 17:07:55 INFO None 5256894: status RUNNING/PENDING
2026-07-22 17:07:55 INFO None 5256895: status RUNNING/PENDING
2026-07-22 17:07:55 INFO None 5256896: status RUNNING/PENDING
2026-07-22 17:07:55 INFO None 5256897: status RUNNING/PENDING
2026-07-22 17:07:55 INFO None 5256899: status RUNNING/PENDING
2026-07-22 17:07:55 INFO None 5256902: status RUNNING/PENDING
2026-07-22 17:07:55 INFO None 5256904: status RUNNING/PENDING
2026-07-22 17:07:55 INFO None 5256905: status RUNNING/PENDING
2026-07-22 17:07:55 INFO None 5256908: status RUNNING/PENDING
2026-07-22 17:07:55 INFO None 5256909: status RUNNING/PENDING
2026-07-22 17:07:55 INFO None 5256911: status RUNNING/PENDING
2026-07-22 17:07:55 INFO None 5256913: status RUNNING/PENDING
2026-07-22 17:07:55 INFO None 5256914: status RUNNING/PENDING
2026-07-22 17:07:55 INFO None 5256916: status RUNNING/PENDING
2026-07-22 17:07:55 INFO Jobs still running: ['5256891', '5256894', '5256895', '5256896', '5256897', '5256899', '5256902', '5256904', '5256905', '5256908', '5256909', '5256911', '5256913', '5256914', '5256916']. Waiting...
2026-07-22 17:08:10 INFO None 5256891: status RUNNING/PENDING
2026-07-22 17:08:11 INFO None 5256894: status RUNNING/PENDING
2026-07-22 17:08:11 INFO None 5256895: status RUNNING/PENDING
2026-07-22 17:08:11 INFO None 5256896: status RUNNING/PENDING
2026-07-22 17:08:11 INFO None 5256897: status RUNNING/PENDING
2026-07-22 17:08:11 INFO None 5256899: status RUNNING/PENDING
2026-07-22 17:08:11 INFO None 5256902: status RUNNING/PENDING
2026-07-22 17:08:11 INFO None 5256904: status RUNNING/PENDING
2026-07-22 17:08:11 INFO None 5256905: status RUNNING/PENDING
2026-07-22 17:08:11 INFO None 5256908: status RUNNING/PENDING
2026-07-22 17:08:11 INFO None 5256909: status RUNNING/PENDING
2026-07-22 17:08:11 INFO None 5256911: status RUNNING/PENDING
2026-07-22 17:08:11 INFO None 5256913: status RUNNING/PENDING
2026-07-22 17:08:11 INFO None 5256914: status RUNNING/PENDING
2026-07-22 17:08:11 INFO None 5256916: status RUNNING/PENDING
2026-07-22 17:08:11 INFO Jobs still running: ['5256891', '5256894', '5256895', '5256896', '5256897', '5256899', '5256902', '5256904', '5256905', '5256908', '5256909', '5256911', '5256913', '5256914', '5256916']. Waiting...
2026-07-22 17:08:26 INFO None 5256891: status RUNNING/PENDING
2026-07-22 17:08:26 INFO None 5256894: status RUNNING/PENDING
2026-07-22 17:08:26 INFO None 5256895: status RUNNING/PENDING
2026-07-22 17:08:26 INFO None 5256896: status RUNNING/PENDING
2026-07-22 17:08:26 INFO None 5256897: status RUNNING/PENDING
2026-07-22 17:08:26 INFO None 5256899: status RUNNING/PENDING
2026-07-22 17:08:26 INFO None 5256902: status RUNNING/PENDING
2026-07-22 17:08:26 INFO None 5256904: status RUNNING/PENDING
2026-07-22 17:08:26 INFO None 5256905: status RUNNING/PENDING
2026-07-22 17:08:26 INFO None 5256908: status RUNNING/PENDING
2026-07-22 17:08:28 INFO None 5256909: status RUNNING/PENDING
2026-07-22 17:08:28 INFO None 5256911: status RUNNING/PENDING
2026-07-22 17:08:28 INFO None 5256913: status RUNNING/PENDING
2026-07-22 17:08:28 INFO None 5256914: status RUNNING/PENDING
2026-07-22 17:08:28 INFO None 5256916: status RUNNING/PENDING
2026-07-22 17:08:28 INFO Jobs still running: ['5256891', '5256894', '5256895', '5256896', '5256897', '5256899', '5256902', '5256904', '5256905', '5256908', '5256909', '5256911', '5256913', '5256914', '5256916']. Waiting...
2026-07-22 17:08:43 INFO None 5256891: status RUNNING/PENDING
2026-07-22 17:08:43 INFO None 5256894: status RUNNING/PENDING
2026-07-22 17:08:43 INFO None 5256895: status RUNNING/PENDING
2026-07-22 17:08:43 INFO None 5256896: status RUNNING/PENDING
2026-07-22 17:08:43 INFO None 5256897: status RUNNING/PENDING
2026-07-22 17:08:43 INFO None 5256899: status RUNNING/PENDING
2026-07-22 17:08:43 INFO None 5256902: status RUNNING/PENDING
2026-07-22 17:08:43 INFO None 5256904: status RUNNING/PENDING
2026-07-22 17:08:43 INFO None 5256905: status RUNNING/PENDING
2026-07-22 17:08:43 INFO None 5256908: status RUNNING/PENDING
2026-07-22 17:08:43 INFO None 5256909: status RUNNING/PENDING
2026-07-22 17:08:43 INFO None 5256911: status RUNNING/PENDING
2026-07-22 17:08:43 INFO None 5256913: status RUNNING/PENDING
2026-07-22 17:08:43 INFO None 5256914: status RUNNING/PENDING
2026-07-22 17:08:43 INFO None 5256916: status RUNNING/PENDING
2026-07-22 17:08:43 INFO Jobs still running: ['5256891', '5256894', '5256895', '5256896', '5256897', '5256899', '5256902', '5256904', '5256905', '5256908', '5256909', '5256911', '5256913', '5256914', '5256916']. Waiting...
2026-07-22 17:08:58 INFO None 5256891: status RUNNING/PENDING
2026-07-22 17:08:58 INFO None 5256894: status RUNNING/PENDING
2026-07-22 17:08:58 INFO None 5256895: status RUNNING/PENDING
2026-07-22 17:08:58 INFO None 5256896: status RUNNING/PENDING
2026-07-22 17:08:59 INFO None 5256897: status RUNNING/PENDING
2026-07-22 17:08:59 INFO None 5256899: status RUNNING/PENDING
2026-07-22 17:08:59 INFO None 5256902: status RUNNING/PENDING
2026-07-22 17:08:59 INFO None 5256904: status RUNNING/PENDING
2026-07-22 17:08:59 INFO None 5256905: status RUNNING/PENDING
2026-07-22 17:08:59 INFO None 5256908: status RUNNING/PENDING
2026-07-22 17:08:59 INFO None 5256909: status RUNNING/PENDING
2026-07-22 17:09:01 INFO None 5256911: status RUNNING/PENDING
2026-07-22 17:09:01 INFO None 5256913: status RUNNING/PENDING
2026-07-22 17:09:01 INFO None 5256914: status RUNNING/PENDING
2026-07-22 17:09:01 INFO None 5256916: status RUNNING/PENDING
2026-07-22 17:09:01 INFO Jobs still running: ['5256891', '5256894', '5256895', '5256896', '5256897', '5256899', '5256902', '5256904', '5256905', '5256908', '5256909', '5256911', '5256913', '5256914', '5256916']. Waiting...
2026-07-22 17:09:16 INFO None 5256891: status RUNNING/PENDING
2026-07-22 17:09:16 INFO None 5256894: status RUNNING/PENDING
2026-07-22 17:09:16 INFO None 5256895: status RUNNING/PENDING
2026-07-22 17:09:16 INFO None 5256896: status RUNNING/PENDING
2026-07-22 17:09:16 INFO None 5256897: status RUNNING/PENDING
2026-07-22 17:09:16 INFO None 5256899: status RUNNING/PENDING
2026-07-22 17:09:16 INFO None 5256902: status RUNNING/PENDING
2026-07-22 17:09:16 INFO None 5256904: status RUNNING/PENDING
2026-07-22 17:09:16 INFO None 5256905: status RUNNING/PENDING
2026-07-22 17:09:16 INFO None 5256908: status RUNNING/PENDING
2026-07-22 17:09:16 INFO None 5256909: status RUNNING/PENDING
2026-07-22 17:09:16 INFO None 5256911: status RUNNING/PENDING
2026-07-22 17:09:16 INFO None 5256913: status RUNNING/PENDING
2026-07-22 17:09:16 INFO None 5256914: status RUNNING/PENDING
2026-07-22 17:09:16 INFO None 5256916: status RUNNING/PENDING
2026-07-22 17:09:16 INFO Jobs still running: ['5256891', '5256894', '5256895', '5256896', '5256897', '5256899', '5256902', '5256904', '5256905', '5256908', '5256909', '5256911', '5256913', '5256914', '5256916']. Waiting...
2026-07-22 17:09:31 INFO None 5256891: status RUNNING/PENDING
2026-07-22 17:09:31 INFO None 5256894: status RUNNING/PENDING
2026-07-22 17:09:31 INFO None 5256895: status RUNNING/PENDING
2026-07-22 17:09:31 INFO None 5256896: status RUNNING/PENDING
2026-07-22 17:09:31 INFO None 5256897: status RUNNING/PENDING
2026-07-22 17:09:31 INFO None 5256899: status RUNNING/PENDING
2026-07-22 17:09:31 INFO None 5256902: status RUNNING/PENDING
2026-07-22 17:09:31 INFO None 5256904: status RUNNING/PENDING
2026-07-22 17:09:31 INFO None 5256905: status RUNNING/PENDING
2026-07-22 17:09:31 INFO None 5256908: status RUNNING/PENDING
2026-07-22 17:09:31 INFO None 5256909: status RUNNING/PENDING
2026-07-22 17:09:31 INFO None 5256911: status RUNNING/PENDING
2026-07-22 17:09:31 INFO None 5256913: status RUNNING/PENDING
2026-07-22 17:09:31 INFO None 5256914: status RUNNING/PENDING
2026-07-22 17:09:31 INFO None 5256916: status RUNNING/PENDING
2026-07-22 17:09:31 INFO Jobs still running: ['5256891', '5256894', '5256895', '5256896', '5256897', '5256899', '5256902', '5256904', '5256905', '5256908', '5256909', '5256911', '5256913', '5256914', '5256916']. Waiting...
2026-07-22 17:09:46 INFO None 5256891: status RUNNING/PENDING
2026-07-22 17:09:46 INFO None 5256894: status RUNNING/PENDING
2026-07-22 17:09:46 INFO None 5256895: status RUNNING/PENDING
2026-07-22 17:09:46 INFO None 5256896: status RUNNING/PENDING
2026-07-22 17:09:47 INFO None 5256897: status RUNNING/PENDING
2026-07-22 17:09:47 INFO None 5256899: status RUNNING/PENDING
2026-07-22 17:09:47 INFO None 5256902: status RUNNING/PENDING
2026-07-22 17:09:47 INFO None 5256904: status FINISHED
2026-07-22 17:09:47 INFO None 5256905: status FINISHED
2026-07-22 17:09:47 INFO None 5256908: status FINISHED
2026-07-22 17:09:47 INFO None 5256909: status RUNNING/PENDING
2026-07-22 17:09:47 INFO None 5256911: status RUNNING/PENDING
2026-07-22 17:09:47 INFO None 5256913: status RUNNING/PENDING
2026-07-22 17:09:47 INFO None 5256914: status RUNNING/PENDING
2026-07-22 17:09:47 INFO None 5256916: status RUNNING/PENDING
2026-07-22 17:09:47 INFO Jobs still running: ['5256891', '5256894', '5256895', '5256896', '5256897', '5256899', '5256902', '5256909', '5256911', '5256913', '5256914', '5256916']. Waiting...
2026-07-22 17:10:02 INFO None 5256891: status RUNNING/PENDING
2026-07-22 17:10:02 INFO None 5256894: status RUNNING/PENDING
2026-07-22 17:10:02 INFO None 5256895: status RUNNING/PENDING
2026-07-22 17:10:02 INFO None 5256896: status FINISHED
2026-07-22 17:10:02 INFO None 5256897: status RUNNING/PENDING
2026-07-22 17:10:02 INFO None 5256899: status FINISHED
2026-07-22 17:10:02 INFO None 5256902: status RUNNING/PENDING
2026-07-22 17:10:02 INFO None 5256904: status FINISHED
2026-07-22 17:10:02 INFO None 5256905: status FINISHED
2026-07-22 17:10:02 INFO None 5256908: status FINISHED
2026-07-22 17:10:02 INFO None 5256909: status RUNNING/PENDING
2026-07-22 17:10:02 INFO None 5256911: status RUNNING/PENDING
2026-07-22 17:10:02 INFO None 5256913: status RUNNING/PENDING
2026-07-22 17:10:02 INFO None 5256914: status RUNNING/PENDING
2026-07-22 17:10:02 INFO None 5256916: status RUNNING/PENDING
2026-07-22 17:10:02 INFO Jobs still running: ['5256891', '5256894', '5256895', '5256897', '5256902', '5256909', '5256911', '5256913', '5256914', '5256916']. Waiting...
2026-07-22 17:10:18 INFO None 5256891: status RUNNING/PENDING
2026-07-22 17:10:18 INFO None 5256894: status RUNNING/PENDING
2026-07-22 17:10:19 INFO None 5256895: status RUNNING/PENDING
2026-07-22 17:10:19 INFO None 5256896: status FINISHED
2026-07-22 17:10:19 INFO None 5256897: status RUNNING/PENDING
2026-07-22 17:10:19 INFO None 5256899: status FINISHED
2026-07-22 17:10:19 INFO None 5256902: status FINISHED
2026-07-22 17:10:19 INFO None 5256904: status FINISHED
2026-07-22 17:10:19 INFO None 5256905: status FINISHED
2026-07-22 17:10:19 INFO None 5256908: status FINISHED
2026-07-22 17:10:19 INFO None 5256909: status RUNNING/PENDING
2026-07-22 17:10:19 INFO None 5256911: status RUNNING/PENDING
2026-07-22 17:10:19 INFO None 5256913: status RUNNING/PENDING
2026-07-22 17:10:19 INFO None 5256914: status RUNNING/PENDING
2026-07-22 17:10:19 INFO None 5256916: status RUNNING/PENDING
2026-07-22 17:10:19 INFO Jobs still running: ['5256891', '5256894', '5256895', '5256897', '5256909', '5256911', '5256913', '5256914', '5256916']. Waiting...
2026-07-22 17:10:34 INFO None 5256891: status RUNNING/PENDING
2026-07-22 17:10:34 INFO None 5256894: status RUNNING/PENDING
2026-07-22 17:10:34 INFO None 5256895: status RUNNING/PENDING
2026-07-22 17:10:34 INFO None 5256896: status FINISHED
2026-07-22 17:10:34 INFO None 5256897: status RUNNING/PENDING
2026-07-22 17:10:34 INFO None 5256899: status FINISHED
2026-07-22 17:10:34 INFO None 5256902: status FINISHED
2026-07-22 17:10:34 INFO None 5256904: status FINISHED
2026-07-22 17:10:34 INFO None 5256905: status FINISHED
2026-07-22 17:10:34 INFO None 5256908: status FINISHED
2026-07-22 17:10:34 INFO None 5256909: status RUNNING/PENDING
2026-07-22 17:10:34 INFO None 5256911: status RUNNING/PENDING
2026-07-22 17:10:34 INFO None 5256913: status RUNNING/PENDING
2026-07-22 17:10:36 INFO None 5256914: status RUNNING/PENDING
2026-07-22 17:10:36 INFO None 5256916: status RUNNING/PENDING
2026-07-22 17:10:36 INFO Jobs still running: ['5256891', '5256894', '5256895', '5256897', '5256909', '5256911', '5256913', '5256914', '5256916']. Waiting...
2026-07-22 17:10:51 INFO None 5256891: status FINISHED
2026-07-22 17:10:51 INFO None 5256894: status FINISHED
2026-07-22 17:10:51 INFO None 5256895: status FINISHED
2026-07-22 17:10:51 INFO None 5256896: status FINISHED
2026-07-22 17:10:51 INFO None 5256897: status FINISHED
2026-07-22 17:10:51 INFO None 5256899: status FINISHED
2026-07-22 17:10:51 INFO None 5256902: status FINISHED
2026-07-22 17:10:51 INFO None 5256904: status FINISHED
2026-07-22 17:10:51 INFO None 5256905: status FINISHED
2026-07-22 17:10:51 INFO None 5256908: status FINISHED
2026-07-22 17:10:51 INFO None 5256909: status RUNNING/PENDING
2026-07-22 17:10:51 INFO None 5256911: status RUNNING/PENDING
2026-07-22 17:10:51 INFO None 5256913: status RUNNING/PENDING
2026-07-22 17:10:51 INFO None 5256914: status RUNNING/PENDING
2026-07-22 17:10:51 INFO None 5256916: status RUNNING/PENDING
2026-07-22 17:10:51 INFO Jobs still running: ['5256909', '5256911', '5256913', '5256914', '5256916']. Waiting...
2026-07-22 17:11:06 INFO None 5256891: status FINISHED
2026-07-22 17:11:06 INFO None 5256894: status FINISHED
2026-07-22 17:11:06 INFO None 5256895: status FINISHED
2026-07-22 17:11:06 INFO None 5256896: status FINISHED
2026-07-22 17:11:07 INFO None 5256897: status FINISHED
2026-07-22 17:11:07 INFO None 5256899: status FINISHED
2026-07-22 17:11:07 INFO None 5256902: status FINISHED
2026-07-22 17:11:07 INFO None 5256904: status FINISHED
2026-07-22 17:11:07 INFO None 5256905: status FINISHED
2026-07-22 17:11:07 INFO None 5256908: status FINISHED
2026-07-22 17:11:09 INFO None 5256909: status RUNNING/PENDING
2026-07-22 17:11:09 INFO None 5256911: status RUNNING/PENDING
2026-07-22 17:11:09 INFO None 5256913: status RUNNING/PENDING
2026-07-22 17:11:09 INFO None 5256914: status RUNNING/PENDING
2026-07-22 17:11:09 INFO None 5256916: status RUNNING/PENDING
2026-07-22 17:11:09 INFO Jobs still running: ['5256909', '5256911', '5256913', '5256914', '5256916']. Waiting...
2026-07-22 17:11:24 INFO None 5256891: status FINISHED
2026-07-22 17:11:24 INFO None 5256894: status FINISHED
2026-07-22 17:11:24 INFO None 5256895: status FINISHED
2026-07-22 17:11:24 INFO None 5256896: status FINISHED
2026-07-22 17:11:24 INFO None 5256897: status FINISHED
2026-07-22 17:11:24 INFO None 5256899: status FINISHED
2026-07-22 17:11:24 INFO None 5256902: status FINISHED
2026-07-22 17:11:24 INFO None 5256904: status FINISHED
2026-07-22 17:11:24 INFO None 5256905: status FINISHED
2026-07-22 17:11:24 INFO None 5256908: status FINISHED
2026-07-22 17:11:24 INFO None 5256909: status RUNNING/PENDING
2026-07-22 17:11:24 INFO None 5256911: status RUNNING/PENDING
2026-07-22 17:11:24 INFO None 5256913: status RUNNING/PENDING
2026-07-22 17:11:24 INFO None 5256914: status RUNNING/PENDING
2026-07-22 17:11:24 INFO None 5256916: status RUNNING/PENDING
2026-07-22 17:11:24 INFO Jobs still running: ['5256909', '5256911', '5256913', '5256914', '5256916']. Waiting...
2026-07-22 17:11:39 INFO None 5256891: status FINISHED
2026-07-22 17:11:39 INFO None 5256894: status FINISHED
2026-07-22 17:11:39 INFO None 5256895: status FINISHED
2026-07-22 17:11:39 INFO None 5256896: status FINISHED
2026-07-22 17:11:39 INFO None 5256897: status FINISHED
2026-07-22 17:11:39 INFO None 5256899: status FINISHED
2026-07-22 17:11:39 INFO None 5256902: status FINISHED
2026-07-22 17:11:39 INFO None 5256904: status FINISHED
2026-07-22 17:11:39 INFO None 5256905: status FINISHED
2026-07-22 17:11:39 INFO None 5256908: status FINISHED
2026-07-22 17:11:39 INFO None 5256909: status RUNNING/PENDING
2026-07-22 17:11:39 INFO None 5256911: status FINISHED
2026-07-22 17:11:39 INFO None 5256913: status FINISHED
2026-07-22 17:11:39 INFO None 5256914: status RUNNING/PENDING
2026-07-22 17:11:39 INFO None 5256916: status RUNNING/PENDING
2026-07-22 17:11:39 INFO Jobs still running: ['5256909', '5256914', '5256916']. Waiting...
2026-07-22 17:11:54 INFO None 5256891: status FINISHED
2026-07-22 17:11:54 INFO None 5256894: status FINISHED
2026-07-22 17:11:54 INFO None 5256895: status FINISHED
2026-07-22 17:11:54 INFO None 5256896: status FINISHED
2026-07-22 17:11:54 INFO None 5256897: status FINISHED
2026-07-22 17:11:54 INFO None 5256899: status FINISHED
2026-07-22 17:11:54 INFO None 5256902: status FINISHED
2026-07-22 17:11:54 INFO None 5256904: status FINISHED
2026-07-22 17:11:54 INFO None 5256905: status FINISHED
2026-07-22 17:11:54 INFO None 5256908: status FINISHED
2026-07-22 17:11:55 INFO None 5256909: status RUNNING/PENDING
2026-07-22 17:11:55 INFO None 5256911: status FINISHED
2026-07-22 17:11:55 INFO None 5256913: status FINISHED
2026-07-22 17:11:55 INFO None 5256914: status RUNNING/PENDING
2026-07-22 17:11:55 INFO None 5256916: status RUNNING/PENDING
2026-07-22 17:11:55 INFO Jobs still running: ['5256909', '5256914', '5256916']. Waiting...
2026-07-22 17:12:11 INFO None 5256891: status FINISHED
2026-07-22 17:12:11 INFO None 5256894: status FINISHED
2026-07-22 17:12:11 INFO None 5256895: status FINISHED
2026-07-22 17:12:11 INFO None 5256896: status FINISHED
2026-07-22 17:12:11 INFO None 5256897: status FINISHED
2026-07-22 17:12:11 INFO None 5256899: status FINISHED
2026-07-22 17:12:11 INFO None 5256902: status FINISHED
2026-07-22 17:12:11 INFO None 5256904: status FINISHED
2026-07-22 17:12:11 INFO None 5256905: status FINISHED
2026-07-22 17:12:11 INFO None 5256908: status FINISHED
2026-07-22 17:12:11 INFO None 5256909: status FINISHED
2026-07-22 17:12:11 INFO None 5256911: status FINISHED
2026-07-22 17:12:11 INFO None 5256913: status FINISHED
2026-07-22 17:12:11 INFO None 5256914: status FINISHED
2026-07-22 17:12:11 INFO None 5256916: status RUNNING/PENDING
2026-07-22 17:12:11 INFO Jobs still running: ['5256916']. Waiting...
2026-07-22 17:12:26 INFO None 5256891: status FINISHED
2026-07-22 17:12:26 INFO None 5256894: status FINISHED
2026-07-22 17:12:26 INFO None 5256895: status FINISHED
2026-07-22 17:12:26 INFO None 5256896: status FINISHED
2026-07-22 17:12:26 INFO None 5256897: status FINISHED
2026-07-22 17:12:26 INFO None 5256899: status FINISHED
2026-07-22 17:12:26 INFO None 5256902: status FINISHED
2026-07-22 17:12:26 INFO None 5256904: status FINISHED
2026-07-22 17:12:26 INFO None 5256905: status FINISHED
2026-07-22 17:12:28 INFO None 5256908: status FINISHED
2026-07-22 17:12:28 INFO None 5256909: status FINISHED
2026-07-22 17:12:28 INFO None 5256911: status FINISHED
2026-07-22 17:12:28 INFO None 5256913: status FINISHED
2026-07-22 17:12:28 INFO None 5256914: status FINISHED
2026-07-22 17:12:28 INFO None 5256916: status FINISHED
2026-07-22 17:12:28 INFO Jobs ['5256891', '5256894', '5256895', '5256896', '5256897', '5256899', '5256902', '5256904', '5256905', '5256908', '5256909', '5256911', '5256913', '5256914', '5256916'] have finished
2026-07-22 17:12:28 INFO Checking restart files were created ...
2026-07-22 17:12:28 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020611_2_ENS1.nc(1002685915 bytes)
2026-07-22 17:12:28 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020611_2_ENS2.nc(1002685915 bytes)
2026-07-22 17:12:28 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020611_2_ENS3.nc(1002685915 bytes)
2026-07-22 17:12:28 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020611_2_ENS4.nc(1002685915 bytes)
2026-07-22 17:12:28 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020611_2_ENS5.nc(1002685915 bytes)
2026-07-22 17:12:28 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020611_2_ENS6.nc(1002685915 bytes)
2026-07-22 17:12:28 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020611_2_ENS7.nc(1002685915 bytes)
2026-07-22 17:12:28 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020611_2_ENS8.nc(1002685915 bytes)
2026-07-22 17:12:28 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020611_2_ENS9.nc(1002685915 bytes)
2026-07-22 17:12:28 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020611_2_ENS10.nc(1002685915 bytes)
2026-07-22 17:12:28 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020611_2_ENS11.nc(1002685915 bytes)
2026-07-22 17:12:28 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020611_2_ENS12.nc(1002685915 bytes)
2026-07-22 17:12:28 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020611_2_ENS13.nc(1002685915 bytes)
2026-07-22 17:12:28 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020611_2_ENS14.nc(1002685915 bytes)
2026-07-22 17:12:28 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020611_2_ENS15.nc(1002685915 bytes)
2026-07-22 17:12:28 INFO  Run_model() completed successfully.
2026-07-22 17:12:28 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 17:12:28 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 13:00:00 days=153072 seconds=46800
2026-07-22 17:12:28 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-22 17:12:28 INFO [TIME] increment current_time 2020-02-06 11:00:00 -> 2020-02-06 13:00:00
2026-07-22 17:12:28 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 17:12:28 INFO ---------->>> Running process_satellite_data()
2026-07-22 17:12:28 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12004.nc
2026-07-22 17:12:28 INFO ---------->>> Running run_obs_converter()
2026-07-22 17:12:28 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_45135_153072.out
2026-07-22 17:12:28 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_45135_153072.out
2026-07-22 17:12:28 INFO ---------->>> Running DART
2026-07-22 17:12:28 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-22 17:12:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/chim_ENS1_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/chim_ENS1_2020020613_1_out_toDART.nc
2026-07-22 17:12:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/chim_ENS2_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/chim_ENS2_2020020613_1_out_toDART.nc
2026-07-22 17:12:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/chim_ENS3_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/chim_ENS3_2020020613_1_out_toDART.nc
2026-07-22 17:12:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/chim_ENS4_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/chim_ENS4_2020020613_1_out_toDART.nc
2026-07-22 17:12:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/chim_ENS5_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/chim_ENS5_2020020613_1_out_toDART.nc
2026-07-22 17:12:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/chim_ENS6_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/chim_ENS6_2020020613_1_out_toDART.nc
2026-07-22 17:12:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/chim_ENS7_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/chim_ENS7_2020020613_1_out_toDART.nc
2026-07-22 17:12:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/chim_ENS8_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/chim_ENS8_2020020613_1_out_toDART.nc
2026-07-22 17:12:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/chim_ENS9_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/chim_ENS9_2020020613_1_out_toDART.nc
2026-07-22 17:12:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/chim_ENS10_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/chim_ENS10_2020020613_1_out_toDART.nc
2026-07-22 17:12:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/chim_ENS11_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/chim_ENS11_2020020613_1_out_toDART.nc
2026-07-22 17:12:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/chim_ENS12_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/chim_ENS12_2020020613_1_out_toDART.nc
2026-07-22 17:12:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/chim_ENS13_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/chim_ENS13_2020020613_1_out_toDART.nc
2026-07-22 17:12:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/chim_ENS14_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/chim_ENS14_2020020613_1_out_toDART.nc
2026-07-22 17:12:33 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/chim_ENS15_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/chim_ENS15_2020020613_1_out_toDART.nc
2026-07-22 17:12:33 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-22 17:12:33 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-22 17:12:33 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-22 17:12:33 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-22 17:12:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-22 17:12:33 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-22 17:12:47 INFO Found: []
2026-07-22 17:12:47 INFO No job id returned by command ./run_filter.bsh
2026-07-22 17:12:47 INFO No monitoring will be performed
2026-07-22 17:12:47 INFO Moving DART output files to analysis and preassim directories for date 2020020613 if present ...
2026-07-22 17:12:47 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 17:12:47 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 17:12:47 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 17:12:47 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 17:12:47 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 17:12:47 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 17:12:47 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 17:12:47 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 17:12:47 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 17:12:47 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 17:12:47 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 17:12:47 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 17:12:47 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 17:12:47 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 17:12:47 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 17:12:47 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 17:12:47 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 17:12:47 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 17:12:47 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 17:12:47 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 17:12:47 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 17:12:47 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 17:12:47 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 17:12:47 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 17:12:47 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 17:12:47 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 17:12:47 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 17:12:47 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 17:12:47 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 17:12:47 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 17:12:47 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 17:12:47 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020613'
2026-07-22 17:12:47 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 17:12:47 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020613'
2026-07-22 17:12:47 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-22 17:12:48 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-22 17:12:48 INFO run_dart() is DONE.
2026-07-22 17:12:48 INFO ---------->>> Running update_pollutant_in_end()
2026-07-22 17:12:48 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:12:48 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS1.nc
2026-07-22 17:12:48 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:12:49 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS1.nc
2026-07-22 17:12:49 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:12:49 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:12:49 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS1.nc
2026-07-22 17:12:50 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:12:50 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS1.nc
2026-07-22 17:12:50 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:12:50 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:12:51 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:12:51 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS2.nc
2026-07-22 17:12:51 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:12:52 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS2.nc
2026-07-22 17:12:52 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:12:52 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:12:52 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS2.nc
2026-07-22 17:12:52 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:12:53 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS2.nc
2026-07-22 17:12:53 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:12:53 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:12:53 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:12:53 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS3.nc
2026-07-22 17:12:54 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:12:55 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS3.nc
2026-07-22 17:12:55 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:12:55 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:12:55 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS3.nc
2026-07-22 17:12:56 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:12:56 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS3.nc
2026-07-22 17:12:56 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:12:56 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:12:57 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:12:57 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS4.nc
2026-07-22 17:12:57 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:12:58 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS4.nc
2026-07-22 17:12:58 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:12:58 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:12:58 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS4.nc
2026-07-22 17:12:59 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:12:59 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS4.nc
2026-07-22 17:12:59 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:12:59 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:12:59 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:12:59 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS5.nc
2026-07-22 17:13:00 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:13:01 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS5.nc
2026-07-22 17:13:01 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:13:01 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:13:01 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS5.nc
2026-07-22 17:13:02 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:13:02 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS5.nc
2026-07-22 17:13:02 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:13:02 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:13:02 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:13:02 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS6.nc
2026-07-22 17:13:03 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:13:04 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS6.nc
2026-07-22 17:13:04 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:13:04 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:13:04 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS6.nc
2026-07-22 17:13:04 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:13:05 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS6.nc
2026-07-22 17:13:05 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:13:05 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:13:05 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:13:05 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS7.nc
2026-07-22 17:13:06 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:13:07 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS7.nc
2026-07-22 17:13:07 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:13:07 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:13:07 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS7.nc
2026-07-22 17:13:07 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:13:08 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS7.nc
2026-07-22 17:13:08 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:13:08 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:13:08 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:13:08 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS8.nc
2026-07-22 17:13:09 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:13:09 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS8.nc
2026-07-22 17:13:09 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:13:10 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:13:10 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS8.nc
2026-07-22 17:13:10 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:13:11 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS8.nc
2026-07-22 17:13:11 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:13:11 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:13:11 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:13:11 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS9.nc
2026-07-22 17:13:12 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:13:12 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS9.nc
2026-07-22 17:13:12 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:13:13 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:13:13 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS9.nc
2026-07-22 17:13:13 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:13:14 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS9.nc
2026-07-22 17:13:14 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:13:14 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:13:14 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:13:14 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS10.nc
2026-07-22 17:13:14 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:13:15 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS10.nc
2026-07-22 17:13:15 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:13:16 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:13:16 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS10.nc
2026-07-22 17:13:16 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:13:17 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS10.nc
2026-07-22 17:13:17 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:13:17 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:13:17 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:13:17 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS11.nc
2026-07-22 17:13:18 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:13:18 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS11.nc
2026-07-22 17:13:18 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:13:19 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:13:19 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS11.nc
2026-07-22 17:13:19 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:13:20 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS11.nc
2026-07-22 17:13:20 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:13:20 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:13:20 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:13:20 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS12.nc
2026-07-22 17:13:21 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:13:21 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS12.nc
2026-07-22 17:13:21 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:13:22 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:13:22 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS12.nc
2026-07-22 17:13:22 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:13:23 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS12.nc
2026-07-22 17:13:23 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:13:23 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:13:23 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:13:23 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS13.nc
2026-07-22 17:13:24 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:13:24 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS13.nc
2026-07-22 17:13:24 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:13:24 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:13:24 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS13.nc
2026-07-22 17:13:25 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:13:26 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS13.nc
2026-07-22 17:13:26 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:13:26 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:13:26 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:13:26 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS14.nc
2026-07-22 17:13:27 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:13:28 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS14.nc
2026-07-22 17:13:28 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:13:28 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:13:28 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS14.nc
2026-07-22 17:13:28 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:13:29 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS14.nc
2026-07-22 17:13:29 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:13:29 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:13:29 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:13:29 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS15.nc
2026-07-22 17:13:30 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:13:31 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS15.nc
2026-07-22 17:13:31 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:13:31 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:13:31 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS15.nc
2026-07-22 17:13:32 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:13:32 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS15.nc
2026-07-22 17:13:32 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:13:32 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:13:32 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-22 17:13:32 INFO [TIME] step_end current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 17:13:32 INFO [TIME] step_start current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 17:13:32 INFO [TIME] window start=2020-02-06 13:00:00 end=2020-02-06 14:00:00 run_hours=1 has_assimilation=True
2026-07-22 17:13:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:13:34 INFO Hourly dataset computed and listing created
2026-07-22 17:13:37 INFO Hourly dataset computed
2026-07-22 17:13:37 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:13:38 INFO Hourly dataset computed and listing created
2026-07-22 17:13:38 INFO Hourly dataset computed
2026-07-22 17:13:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:13:39 INFO Hourly dataset computed and listing created
2026-07-22 17:13:39 INFO Hourly dataset computed
2026-07-22 17:13:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:13:40 INFO Hourly dataset computed and listing created
2026-07-22 17:13:41 INFO Hourly dataset computed
2026-07-22 17:13:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:13:42 INFO Hourly dataset computed and listing created
2026-07-22 17:13:42 INFO Hourly dataset computed
2026-07-22 17:13:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:13:43 INFO Hourly dataset computed and listing created
2026-07-22 17:13:44 INFO Hourly dataset computed
2026-07-22 17:13:44 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:13:45 INFO Hourly dataset computed and listing created
2026-07-22 17:13:45 INFO Hourly dataset computed
2026-07-22 17:13:45 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:13:46 INFO Hourly dataset computed and listing created
2026-07-22 17:13:46 INFO Hourly dataset computed
2026-07-22 17:13:46 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:13:47 INFO Hourly dataset computed and listing created
2026-07-22 17:13:48 INFO Hourly dataset computed
2026-07-22 17:13:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:13:49 INFO Hourly dataset computed and listing created
2026-07-22 17:13:49 INFO Hourly dataset computed
2026-07-22 17:13:49 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:13:50 INFO Hourly dataset computed and listing created
2026-07-22 17:13:51 INFO Hourly dataset computed
2026-07-22 17:13:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:13:51 INFO Hourly dataset computed and listing created
2026-07-22 17:13:52 INFO Hourly dataset computed
2026-07-22 17:13:52 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:13:53 INFO Hourly dataset computed and listing created
2026-07-22 17:13:54 INFO Hourly dataset computed
2026-07-22 17:13:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:13:54 INFO Hourly dataset computed and listing created
2026-07-22 17:13:55 INFO Hourly dataset computed
2026-07-22 17:13:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:13:56 INFO Hourly dataset computed and listing created
2026-07-22 17:13:56 INFO Hourly dataset computed
2026-07-22 17:13:56 INFO ---------->>> Running CHIMERE model from 2020-02-06 13:00:00 to 2020-02-06 14:00:00
2026-07-22 17:13:56 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:13:56 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1
2026-07-22 17:13:56 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020611_2_ENS1.nc
2026-07-22 17:13:56 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-22 17:13:56 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:13:56 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-22 17:13:56 INFO Queuing job for member 1...
2026-07-22 17:13:56 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:13:56 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-22 17:13:57 INFO Found: ['5257003']
2026-07-22 17:14:02 INFO [TGCC-IRENE] Submitted job with ID:['5257003']
2026-07-22 17:14:02 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:14:02 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2
2026-07-22 17:14:02 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020611_2_ENS2.nc
2026-07-22 17:14:02 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-22 17:14:02 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:14:02 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-22 17:14:02 INFO Queuing job for member 2...
2026-07-22 17:14:02 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:14:02 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-22 17:14:05 INFO Found: ['5257005']
2026-07-22 17:14:10 INFO [TGCC-IRENE] Submitted job with ID:['5257005']
2026-07-22 17:14:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:14:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3
2026-07-22 17:14:10 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020611_2_ENS3.nc
2026-07-22 17:14:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-22 17:14:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:14:10 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-22 17:14:10 INFO Queuing job for member 3...
2026-07-22 17:14:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:14:10 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-22 17:14:12 INFO Found: ['5257006']
2026-07-22 17:14:17 INFO [TGCC-IRENE] Submitted job with ID:['5257006']
2026-07-22 17:14:17 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:14:17 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4
2026-07-22 17:14:17 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020611_2_ENS4.nc
2026-07-22 17:14:17 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-22 17:14:17 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:14:17 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-22 17:14:17 INFO Queuing job for member 4...
2026-07-22 17:14:17 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:14:17 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-22 17:14:20 INFO Found: ['5257007']
2026-07-22 17:14:25 INFO [TGCC-IRENE] Submitted job with ID:['5257007']
2026-07-22 17:14:25 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:14:25 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5
2026-07-22 17:14:25 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020611_2_ENS5.nc
2026-07-22 17:14:25 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-22 17:14:25 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:14:25 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-22 17:14:25 INFO Queuing job for member 5...
2026-07-22 17:14:25 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:14:25 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-22 17:14:27 INFO Found: ['5257009']
2026-07-22 17:14:32 INFO [TGCC-IRENE] Submitted job with ID:['5257009']
2026-07-22 17:14:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:14:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6
2026-07-22 17:14:32 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020611_2_ENS6.nc
2026-07-22 17:14:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-22 17:14:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:14:32 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-22 17:14:32 INFO Queuing job for member 6...
2026-07-22 17:14:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:14:32 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-22 17:14:35 INFO Found: ['5257011']
2026-07-22 17:14:40 INFO [TGCC-IRENE] Submitted job with ID:['5257011']
2026-07-22 17:14:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:14:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7
2026-07-22 17:14:40 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020611_2_ENS7.nc
2026-07-22 17:14:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-22 17:14:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:14:40 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-22 17:14:40 INFO Queuing job for member 7...
2026-07-22 17:14:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:14:40 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-22 17:14:41 INFO Found: ['5257013']
2026-07-22 17:14:46 INFO [TGCC-IRENE] Submitted job with ID:['5257013']
2026-07-22 17:14:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:14:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8
2026-07-22 17:14:46 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020611_2_ENS8.nc
2026-07-22 17:14:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-22 17:14:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:14:46 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-22 17:14:46 INFO Queuing job for member 8...
2026-07-22 17:14:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:14:46 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-22 17:14:46 INFO Found: ['5257014']
2026-07-22 17:14:51 INFO [TGCC-IRENE] Submitted job with ID:['5257014']
2026-07-22 17:14:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:14:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9
2026-07-22 17:14:51 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020611_2_ENS9.nc
2026-07-22 17:14:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-22 17:14:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:14:51 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-22 17:14:51 INFO Queuing job for member 9...
2026-07-22 17:14:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:14:51 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-22 17:14:52 INFO Found: ['5257016']
2026-07-22 17:14:57 INFO [TGCC-IRENE] Submitted job with ID:['5257016']
2026-07-22 17:14:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:14:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10
2026-07-22 17:14:57 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020611_2_ENS10.nc
2026-07-22 17:14:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-22 17:14:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:14:57 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-22 17:14:57 INFO Queuing job for member 10...
2026-07-22 17:14:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:14:57 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-22 17:14:58 INFO Found: ['5257017']
2026-07-22 17:15:03 INFO [TGCC-IRENE] Submitted job with ID:['5257017']
2026-07-22 17:15:03 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:15:03 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11
2026-07-22 17:15:03 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020611_2_ENS11.nc
2026-07-22 17:15:03 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-22 17:15:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:15:03 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-22 17:15:03 INFO Queuing job for member 11...
2026-07-22 17:15:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:15:03 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-22 17:15:04 INFO Found: ['5257020']
2026-07-22 17:15:09 INFO [TGCC-IRENE] Submitted job with ID:['5257020']
2026-07-22 17:15:09 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:15:09 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12
2026-07-22 17:15:09 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020611_2_ENS12.nc
2026-07-22 17:15:09 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-22 17:15:09 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:15:09 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-22 17:15:09 INFO Queuing job for member 12...
2026-07-22 17:15:09 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:15:09 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-22 17:15:09 INFO Found: ['5257021']
2026-07-22 17:15:14 INFO [TGCC-IRENE] Submitted job with ID:['5257021']
2026-07-22 17:15:14 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:15:14 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13
2026-07-22 17:15:14 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020611_2_ENS13.nc
2026-07-22 17:15:14 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-22 17:15:14 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:15:14 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-22 17:15:14 INFO Queuing job for member 13...
2026-07-22 17:15:14 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:15:14 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-22 17:15:15 INFO Found: ['5257022']
2026-07-22 17:15:20 INFO [TGCC-IRENE] Submitted job with ID:['5257022']
2026-07-22 17:15:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:15:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14
2026-07-22 17:15:20 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020611_2_ENS14.nc
2026-07-22 17:15:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-22 17:15:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:15:20 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-22 17:15:20 INFO Queuing job for member 14...
2026-07-22 17:15:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:15:20 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-22 17:15:21 INFO Found: ['5257024']
2026-07-22 17:15:26 INFO [TGCC-IRENE] Submitted job with ID:['5257024']
2026-07-22 17:15:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:15:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15
2026-07-22 17:15:26 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020611_2_ENS15.nc
2026-07-22 17:15:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-22 17:15:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:15:26 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-22 17:15:26 INFO Queuing job for member 15...
2026-07-22 17:15:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:15:26 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-22 17:15:28 INFO Found: ['5257025']
2026-07-22 17:15:33 INFO [TGCC-IRENE] Submitted job with ID:['5257025']
2026-07-22 17:15:33 INFO Checking job status ...
2026-07-22 17:15:33 INFO None 5257003: status RUNNING/PENDING
2026-07-22 17:15:33 INFO None 5257005: status RUNNING/PENDING
2026-07-22 17:15:33 INFO None 5257006: status RUNNING/PENDING
2026-07-22 17:15:33 INFO None 5257007: status RUNNING/PENDING
2026-07-22 17:15:33 INFO None 5257009: status RUNNING/PENDING
2026-07-22 17:15:33 INFO None 5257011: status RUNNING/PENDING
2026-07-22 17:15:33 INFO None 5257013: status RUNNING/PENDING
2026-07-22 17:15:33 INFO None 5257014: status RUNNING/PENDING
2026-07-22 17:15:33 INFO None 5257016: status RUNNING/PENDING
2026-07-22 17:15:33 INFO None 5257017: status RUNNING/PENDING
2026-07-22 17:15:33 INFO None 5257020: status RUNNING/PENDING
2026-07-22 17:15:33 INFO None 5257021: status RUNNING/PENDING
2026-07-22 17:15:33 INFO None 5257022: status RUNNING/PENDING
2026-07-22 17:15:33 INFO None 5257024: status RUNNING/PENDING
2026-07-22 17:15:34 INFO None 5257025: status RUNNING/PENDING
2026-07-22 17:15:34 INFO Jobs still running: ['5257003', '5257005', '5257006', '5257007', '5257009', '5257011', '5257013', '5257014', '5257016', '5257017', '5257020', '5257021', '5257022', '5257024', '5257025']. Waiting...
2026-07-22 17:15:49 INFO None 5257003: status RUNNING/PENDING
2026-07-22 17:15:49 INFO None 5257005: status RUNNING/PENDING
2026-07-22 17:15:49 INFO None 5257006: status RUNNING/PENDING
2026-07-22 17:15:49 INFO None 5257007: status RUNNING/PENDING
2026-07-22 17:15:49 INFO None 5257009: status RUNNING/PENDING
2026-07-22 17:15:49 INFO None 5257011: status RUNNING/PENDING
2026-07-22 17:15:49 INFO None 5257013: status RUNNING/PENDING
2026-07-22 17:15:49 INFO None 5257014: status RUNNING/PENDING
2026-07-22 17:15:51 INFO None 5257016: status RUNNING/PENDING
2026-07-22 17:15:51 INFO None 5257017: status RUNNING/PENDING
2026-07-22 17:15:51 INFO None 5257020: status RUNNING/PENDING
2026-07-22 17:15:51 INFO None 5257021: status RUNNING/PENDING
2026-07-22 17:15:51 INFO None 5257022: status RUNNING/PENDING
2026-07-22 17:15:51 INFO None 5257024: status RUNNING/PENDING
2026-07-22 17:15:51 INFO None 5257025: status RUNNING/PENDING
2026-07-22 17:15:51 INFO Jobs still running: ['5257003', '5257005', '5257006', '5257007', '5257009', '5257011', '5257013', '5257014', '5257016', '5257017', '5257020', '5257021', '5257022', '5257024', '5257025']. Waiting...
2026-07-22 17:16:06 INFO None 5257003: status RUNNING/PENDING
2026-07-22 17:16:06 INFO None 5257005: status RUNNING/PENDING
2026-07-22 17:16:06 INFO None 5257006: status RUNNING/PENDING
2026-07-22 17:16:06 INFO None 5257007: status RUNNING/PENDING
2026-07-22 17:16:06 INFO None 5257009: status RUNNING/PENDING
2026-07-22 17:16:06 INFO None 5257011: status RUNNING/PENDING
2026-07-22 17:16:06 INFO None 5257013: status RUNNING/PENDING
2026-07-22 17:16:06 INFO None 5257014: status RUNNING/PENDING
2026-07-22 17:16:06 INFO None 5257016: status RUNNING/PENDING
2026-07-22 17:16:06 INFO None 5257017: status RUNNING/PENDING
2026-07-22 17:16:06 INFO None 5257020: status RUNNING/PENDING
2026-07-22 17:16:06 INFO None 5257021: status RUNNING/PENDING
2026-07-22 17:16:06 INFO None 5257022: status RUNNING/PENDING
2026-07-22 17:16:06 INFO None 5257024: status RUNNING/PENDING
2026-07-22 17:16:06 INFO None 5257025: status RUNNING/PENDING
2026-07-22 17:16:06 INFO Jobs still running: ['5257003', '5257005', '5257006', '5257007', '5257009', '5257011', '5257013', '5257014', '5257016', '5257017', '5257020', '5257021', '5257022', '5257024', '5257025']. Waiting...
2026-07-22 17:16:21 INFO None 5257003: status RUNNING/PENDING
2026-07-22 17:16:21 INFO None 5257005: status RUNNING/PENDING
2026-07-22 17:16:21 INFO None 5257006: status RUNNING/PENDING
2026-07-22 17:16:21 INFO None 5257007: status RUNNING/PENDING
2026-07-22 17:16:21 INFO None 5257009: status RUNNING/PENDING
2026-07-22 17:16:21 INFO None 5257011: status RUNNING/PENDING
2026-07-22 17:16:21 INFO None 5257013: status RUNNING/PENDING
2026-07-22 17:16:21 INFO None 5257014: status RUNNING/PENDING
2026-07-22 17:16:21 INFO None 5257016: status RUNNING/PENDING
2026-07-22 17:16:21 INFO None 5257017: status RUNNING/PENDING
2026-07-22 17:16:23 INFO None 5257020: status RUNNING/PENDING
2026-07-22 17:16:23 INFO None 5257021: status RUNNING/PENDING
2026-07-22 17:16:23 INFO None 5257022: status RUNNING/PENDING
2026-07-22 17:16:23 INFO None 5257024: status RUNNING/PENDING
2026-07-22 17:16:23 INFO None 5257025: status RUNNING/PENDING
2026-07-22 17:16:23 INFO Jobs still running: ['5257003', '5257005', '5257006', '5257007', '5257009', '5257011', '5257013', '5257014', '5257016', '5257017', '5257020', '5257021', '5257022', '5257024', '5257025']. Waiting...
2026-07-22 17:16:38 INFO None 5257003: status RUNNING/PENDING
2026-07-22 17:16:38 INFO None 5257005: status RUNNING/PENDING
2026-07-22 17:16:38 INFO None 5257006: status RUNNING/PENDING
2026-07-22 17:16:38 INFO None 5257007: status RUNNING/PENDING
2026-07-22 17:16:38 INFO None 5257009: status RUNNING/PENDING
2026-07-22 17:16:39 INFO None 5257011: status RUNNING/PENDING
2026-07-22 17:16:39 INFO None 5257013: status RUNNING/PENDING
2026-07-22 17:16:39 INFO None 5257014: status RUNNING/PENDING
2026-07-22 17:16:39 INFO None 5257016: status RUNNING/PENDING
2026-07-22 17:16:39 INFO None 5257017: status RUNNING/PENDING
2026-07-22 17:16:39 INFO None 5257020: status RUNNING/PENDING
2026-07-22 17:16:39 INFO None 5257021: status RUNNING/PENDING
2026-07-22 17:16:39 INFO None 5257022: status RUNNING/PENDING
2026-07-22 17:16:39 INFO None 5257024: status RUNNING/PENDING
2026-07-22 17:16:39 INFO None 5257025: status RUNNING/PENDING
2026-07-22 17:16:39 INFO Jobs still running: ['5257003', '5257005', '5257006', '5257007', '5257009', '5257011', '5257013', '5257014', '5257016', '5257017', '5257020', '5257021', '5257022', '5257024', '5257025']. Waiting...
2026-07-22 17:16:54 INFO None 5257003: status RUNNING/PENDING
2026-07-22 17:16:54 INFO None 5257005: status RUNNING/PENDING
2026-07-22 17:16:54 INFO None 5257006: status RUNNING/PENDING
2026-07-22 17:16:54 INFO None 5257007: status RUNNING/PENDING
2026-07-22 17:16:54 INFO None 5257009: status RUNNING/PENDING
2026-07-22 17:16:54 INFO None 5257011: status RUNNING/PENDING
2026-07-22 17:16:54 INFO None 5257013: status RUNNING/PENDING
2026-07-22 17:16:54 INFO None 5257014: status RUNNING/PENDING
2026-07-22 17:16:54 INFO None 5257016: status RUNNING/PENDING
2026-07-22 17:16:54 INFO None 5257017: status RUNNING/PENDING
2026-07-22 17:16:54 INFO None 5257020: status RUNNING/PENDING
2026-07-22 17:16:54 INFO None 5257021: status RUNNING/PENDING
2026-07-22 17:16:54 INFO None 5257022: status RUNNING/PENDING
2026-07-22 17:16:54 INFO None 5257024: status RUNNING/PENDING
2026-07-22 17:16:54 INFO None 5257025: status RUNNING/PENDING
2026-07-22 17:16:54 INFO Jobs still running: ['5257003', '5257005', '5257006', '5257007', '5257009', '5257011', '5257013', '5257014', '5257016', '5257017', '5257020', '5257021', '5257022', '5257024', '5257025']. Waiting...
2026-07-22 17:17:09 INFO None 5257003: status RUNNING/PENDING
2026-07-22 17:17:09 INFO None 5257005: status RUNNING/PENDING
2026-07-22 17:17:09 INFO None 5257006: status RUNNING/PENDING
2026-07-22 17:17:09 INFO None 5257007: status RUNNING/PENDING
2026-07-22 17:17:09 INFO None 5257009: status RUNNING/PENDING
2026-07-22 17:17:09 INFO None 5257011: status RUNNING/PENDING
2026-07-22 17:17:09 INFO None 5257013: status RUNNING/PENDING
2026-07-22 17:17:09 INFO None 5257014: status RUNNING/PENDING
2026-07-22 17:17:09 INFO None 5257016: status RUNNING/PENDING
2026-07-22 17:17:09 INFO None 5257017: status RUNNING/PENDING
2026-07-22 17:17:09 INFO None 5257020: status RUNNING/PENDING
2026-07-22 17:17:09 INFO None 5257021: status RUNNING/PENDING
2026-07-22 17:17:09 INFO None 5257022: status RUNNING/PENDING
2026-07-22 17:17:09 INFO None 5257024: status RUNNING/PENDING
2026-07-22 17:17:09 INFO None 5257025: status RUNNING/PENDING
2026-07-22 17:17:09 INFO Jobs still running: ['5257003', '5257005', '5257006', '5257007', '5257009', '5257011', '5257013', '5257014', '5257016', '5257017', '5257020', '5257021', '5257022', '5257024', '5257025']. Waiting...
2026-07-22 17:17:25 INFO None 5257003: status RUNNING/PENDING
2026-07-22 17:17:25 INFO None 5257005: status RUNNING/PENDING
2026-07-22 17:17:25 INFO None 5257006: status RUNNING/PENDING
2026-07-22 17:17:25 INFO None 5257007: status RUNNING/PENDING
2026-07-22 17:17:25 INFO None 5257009: status RUNNING/PENDING
2026-07-22 17:17:25 INFO None 5257011: status RUNNING/PENDING
2026-07-22 17:17:25 INFO None 5257013: status RUNNING/PENDING
2026-07-22 17:17:25 INFO None 5257014: status RUNNING/PENDING
2026-07-22 17:17:25 INFO None 5257016: status RUNNING/PENDING
2026-07-22 17:17:25 INFO None 5257017: status RUNNING/PENDING
2026-07-22 17:17:25 INFO None 5257020: status RUNNING/PENDING
2026-07-22 17:17:25 INFO None 5257021: status RUNNING/PENDING
2026-07-22 17:17:25 INFO None 5257022: status RUNNING/PENDING
2026-07-22 17:17:25 INFO None 5257024: status RUNNING/PENDING
2026-07-22 17:17:25 INFO None 5257025: status RUNNING/PENDING
2026-07-22 17:17:25 INFO Jobs still running: ['5257003', '5257005', '5257006', '5257007', '5257009', '5257011', '5257013', '5257014', '5257016', '5257017', '5257020', '5257021', '5257022', '5257024', '5257025']. Waiting...
2026-07-22 17:17:40 INFO None 5257003: status RUNNING/PENDING
2026-07-22 17:17:41 INFO None 5257005: status RUNNING/PENDING
2026-07-22 17:17:41 INFO None 5257006: status RUNNING/PENDING
2026-07-22 17:17:41 INFO None 5257007: status RUNNING/PENDING
2026-07-22 17:17:41 INFO None 5257009: status RUNNING/PENDING
2026-07-22 17:17:41 INFO None 5257011: status RUNNING/PENDING
2026-07-22 17:17:41 INFO None 5257013: status RUNNING/PENDING
2026-07-22 17:17:41 INFO None 5257014: status RUNNING/PENDING
2026-07-22 17:17:41 INFO None 5257016: status RUNNING/PENDING
2026-07-22 17:17:41 INFO None 5257017: status RUNNING/PENDING
2026-07-22 17:17:41 INFO None 5257020: status RUNNING/PENDING
2026-07-22 17:17:41 INFO None 5257021: status RUNNING/PENDING
2026-07-22 17:17:41 INFO None 5257022: status RUNNING/PENDING
2026-07-22 17:17:43 INFO None 5257024: status RUNNING/PENDING
2026-07-22 17:17:43 INFO None 5257025: status RUNNING/PENDING
2026-07-22 17:17:43 INFO Jobs still running: ['5257003', '5257005', '5257006', '5257007', '5257009', '5257011', '5257013', '5257014', '5257016', '5257017', '5257020', '5257021', '5257022', '5257024', '5257025']. Waiting...
2026-07-22 17:17:58 INFO None 5257003: status RUNNING/PENDING
2026-07-22 17:17:58 INFO None 5257005: status RUNNING/PENDING
2026-07-22 17:17:58 INFO None 5257006: status RUNNING/PENDING
2026-07-22 17:17:58 INFO None 5257007: status RUNNING/PENDING
2026-07-22 17:17:58 INFO None 5257009: status RUNNING/PENDING
2026-07-22 17:17:58 INFO None 5257011: status RUNNING/PENDING
2026-07-22 17:17:58 INFO None 5257013: status RUNNING/PENDING
2026-07-22 17:17:58 INFO None 5257014: status RUNNING/PENDING
2026-07-22 17:17:58 INFO None 5257016: status RUNNING/PENDING
2026-07-22 17:17:58 INFO None 5257017: status RUNNING/PENDING
2026-07-22 17:17:58 INFO None 5257020: status RUNNING/PENDING
2026-07-22 17:17:58 INFO None 5257021: status RUNNING/PENDING
2026-07-22 17:17:58 INFO None 5257022: status RUNNING/PENDING
2026-07-22 17:17:58 INFO None 5257024: status RUNNING/PENDING
2026-07-22 17:17:58 INFO None 5257025: status RUNNING/PENDING
2026-07-22 17:17:58 INFO Jobs still running: ['5257003', '5257005', '5257006', '5257007', '5257009', '5257011', '5257013', '5257014', '5257016', '5257017', '5257020', '5257021', '5257022', '5257024', '5257025']. Waiting...
2026-07-22 17:18:13 INFO None 5257003: status RUNNING/PENDING
2026-07-22 17:18:13 INFO None 5257005: status RUNNING/PENDING
2026-07-22 17:18:13 INFO None 5257006: status RUNNING/PENDING
2026-07-22 17:18:13 INFO None 5257007: status RUNNING/PENDING
2026-07-22 17:18:13 INFO None 5257009: status RUNNING/PENDING
2026-07-22 17:18:13 INFO None 5257011: status RUNNING/PENDING
2026-07-22 17:18:13 INFO None 5257013: status RUNNING/PENDING
2026-07-22 17:18:13 INFO None 5257014: status RUNNING/PENDING
2026-07-22 17:18:13 INFO None 5257016: status FINISHED
2026-07-22 17:18:13 INFO None 5257017: status RUNNING/PENDING
2026-07-22 17:18:13 INFO None 5257020: status RUNNING/PENDING
2026-07-22 17:18:13 INFO None 5257021: status RUNNING/PENDING
2026-07-22 17:18:13 INFO None 5257022: status RUNNING/PENDING
2026-07-22 17:18:15 INFO None 5257024: status FINISHED
2026-07-22 17:18:16 INFO None 5257025: status RUNNING/PENDING
2026-07-22 17:18:16 INFO Jobs still running: ['5257003', '5257005', '5257006', '5257007', '5257009', '5257011', '5257013', '5257014', '5257017', '5257020', '5257021', '5257022', '5257025']. Waiting...
2026-07-22 17:18:31 INFO None 5257003: status RUNNING/PENDING
2026-07-22 17:18:31 INFO None 5257005: status RUNNING/PENDING
2026-07-22 17:18:31 INFO None 5257006: status RUNNING/PENDING
2026-07-22 17:18:31 INFO None 5257007: status RUNNING/PENDING
2026-07-22 17:18:31 INFO None 5257009: status RUNNING/PENDING
2026-07-22 17:18:31 INFO None 5257011: status RUNNING/PENDING
2026-07-22 17:18:31 INFO None 5257013: status RUNNING/PENDING
2026-07-22 17:18:31 INFO None 5257014: status RUNNING/PENDING
2026-07-22 17:18:31 INFO None 5257016: status FINISHED
2026-07-22 17:18:31 INFO None 5257017: status RUNNING/PENDING
2026-07-22 17:18:31 INFO None 5257020: status RUNNING/PENDING
2026-07-22 17:18:31 INFO None 5257021: status FINISHED
2026-07-22 17:18:31 INFO None 5257022: status RUNNING/PENDING
2026-07-22 17:18:31 INFO None 5257024: status FINISHED
2026-07-22 17:18:31 INFO None 5257025: status RUNNING/PENDING
2026-07-22 17:18:31 INFO Jobs still running: ['5257003', '5257005', '5257006', '5257007', '5257009', '5257011', '5257013', '5257014', '5257017', '5257020', '5257022', '5257025']. Waiting...
2026-07-22 17:18:46 INFO None 5257003: status FINISHED
2026-07-22 17:18:46 INFO None 5257005: status FINISHED
2026-07-22 17:18:46 INFO None 5257006: status FINISHED
2026-07-22 17:18:46 INFO None 5257007: status FINISHED
2026-07-22 17:18:46 INFO None 5257009: status RUNNING/PENDING
2026-07-22 17:18:46 INFO None 5257011: status RUNNING/PENDING
2026-07-22 17:18:46 INFO None 5257013: status RUNNING/PENDING
2026-07-22 17:18:46 INFO None 5257014: status RUNNING/PENDING
2026-07-22 17:18:46 INFO None 5257016: status FINISHED
2026-07-22 17:18:46 INFO None 5257017: status RUNNING/PENDING
2026-07-22 17:18:46 INFO None 5257020: status RUNNING/PENDING
2026-07-22 17:18:46 INFO None 5257021: status FINISHED
2026-07-22 17:18:46 INFO None 5257022: status RUNNING/PENDING
2026-07-22 17:18:46 INFO None 5257024: status FINISHED
2026-07-22 17:18:46 INFO None 5257025: status RUNNING/PENDING
2026-07-22 17:18:46 INFO Jobs still running: ['5257009', '5257011', '5257013', '5257014', '5257017', '5257020', '5257022', '5257025']. Waiting...
2026-07-22 17:19:03 INFO None 5257003: status FINISHED
2026-07-22 17:19:03 INFO None 5257005: status FINISHED
2026-07-22 17:19:03 INFO None 5257006: status FINISHED
2026-07-22 17:19:03 INFO None 5257007: status FINISHED
2026-07-22 17:19:03 INFO None 5257009: status RUNNING/PENDING
2026-07-22 17:19:03 INFO None 5257011: status FINISHED
2026-07-22 17:19:03 INFO None 5257013: status FINISHED
2026-07-22 17:19:03 INFO None 5257014: status FINISHED
2026-07-22 17:19:03 INFO None 5257016: status FINISHED
2026-07-22 17:19:03 INFO None 5257017: status FINISHED
2026-07-22 17:19:03 INFO None 5257020: status RUNNING/PENDING
2026-07-22 17:19:03 INFO None 5257021: status FINISHED
2026-07-22 17:19:03 INFO None 5257022: status RUNNING/PENDING
2026-07-22 17:19:03 INFO None 5257024: status FINISHED
2026-07-22 17:19:03 INFO None 5257025: status RUNNING/PENDING
2026-07-22 17:19:03 INFO Jobs still running: ['5257009', '5257020', '5257022', '5257025']. Waiting...
2026-07-22 17:19:18 INFO None 5257003: status FINISHED
2026-07-22 17:19:18 INFO None 5257005: status FINISHED
2026-07-22 17:19:18 INFO None 5257006: status FINISHED
2026-07-22 17:19:18 INFO None 5257007: status FINISHED
2026-07-22 17:19:19 INFO None 5257009: status FINISHED
2026-07-22 17:19:19 INFO None 5257011: status FINISHED
2026-07-22 17:19:19 INFO None 5257013: status FINISHED
2026-07-22 17:19:19 INFO None 5257014: status FINISHED
2026-07-22 17:19:19 INFO None 5257016: status FINISHED
2026-07-22 17:19:19 INFO None 5257017: status FINISHED
2026-07-22 17:19:19 INFO None 5257020: status FINISHED
2026-07-22 17:19:19 INFO None 5257021: status FINISHED
2026-07-22 17:19:19 INFO None 5257022: status FINISHED
2026-07-22 17:19:19 INFO None 5257024: status FINISHED
2026-07-22 17:19:21 INFO None 5257025: status RUNNING/PENDING
2026-07-22 17:19:21 INFO Jobs still running: ['5257025']. Waiting...
2026-07-22 17:19:36 INFO None 5257003: status FINISHED
2026-07-22 17:19:36 INFO None 5257005: status FINISHED
2026-07-22 17:19:36 INFO None 5257006: status FINISHED
2026-07-22 17:19:36 INFO None 5257007: status FINISHED
2026-07-22 17:19:36 INFO None 5257009: status FINISHED
2026-07-22 17:19:36 INFO None 5257011: status FINISHED
2026-07-22 17:19:36 INFO None 5257013: status FINISHED
2026-07-22 17:19:36 INFO None 5257014: status FINISHED
2026-07-22 17:19:36 INFO None 5257016: status FINISHED
2026-07-22 17:19:36 INFO None 5257017: status FINISHED
2026-07-22 17:19:36 INFO None 5257020: status FINISHED
2026-07-22 17:19:36 INFO None 5257021: status FINISHED
2026-07-22 17:19:36 INFO None 5257022: status FINISHED
2026-07-22 17:19:36 INFO None 5257024: status FINISHED
2026-07-22 17:19:36 INFO None 5257025: status RUNNING/PENDING
2026-07-22 17:19:36 INFO Jobs still running: ['5257025']. Waiting...
2026-07-22 17:19:51 INFO None 5257003: status FINISHED
2026-07-22 17:19:51 INFO None 5257005: status FINISHED
2026-07-22 17:19:51 INFO None 5257006: status FINISHED
2026-07-22 17:19:51 INFO None 5257007: status FINISHED
2026-07-22 17:19:51 INFO None 5257009: status FINISHED
2026-07-22 17:19:51 INFO None 5257011: status FINISHED
2026-07-22 17:19:51 INFO None 5257013: status FINISHED
2026-07-22 17:19:51 INFO None 5257014: status FINISHED
2026-07-22 17:19:51 INFO None 5257016: status FINISHED
2026-07-22 17:19:51 INFO None 5257017: status FINISHED
2026-07-22 17:19:51 INFO None 5257020: status FINISHED
2026-07-22 17:19:51 INFO None 5257021: status FINISHED
2026-07-22 17:19:51 INFO None 5257022: status FINISHED
2026-07-22 17:19:51 INFO None 5257024: status FINISHED
2026-07-22 17:19:53 INFO None 5257025: status RUNNING/PENDING
2026-07-22 17:19:53 INFO Jobs still running: ['5257025']. Waiting...
2026-07-22 17:20:08 INFO None 5257003: status FINISHED
2026-07-22 17:20:08 INFO None 5257005: status FINISHED
2026-07-22 17:20:08 INFO None 5257006: status FINISHED
2026-07-22 17:20:08 INFO None 5257007: status FINISHED
2026-07-22 17:20:08 INFO None 5257009: status FINISHED
2026-07-22 17:20:08 INFO None 5257011: status FINISHED
2026-07-22 17:20:08 INFO None 5257013: status FINISHED
2026-07-22 17:20:08 INFO None 5257014: status FINISHED
2026-07-22 17:20:09 INFO None 5257016: status FINISHED
2026-07-22 17:20:09 INFO None 5257017: status FINISHED
2026-07-22 17:20:09 INFO None 5257020: status FINISHED
2026-07-22 17:20:09 INFO None 5257021: status FINISHED
2026-07-22 17:20:09 INFO None 5257022: status FINISHED
2026-07-22 17:20:09 INFO None 5257024: status FINISHED
2026-07-22 17:20:09 INFO None 5257025: status FINISHED
2026-07-22 17:20:09 INFO Jobs ['5257003', '5257005', '5257006', '5257007', '5257009', '5257011', '5257013', '5257014', '5257016', '5257017', '5257020', '5257021', '5257022', '5257024', '5257025'] have finished
2026-07-22 17:20:09 INFO Checking restart files were created ...
2026-07-22 17:20:09 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020613_1_ENS1.nc(668832435 bytes)
2026-07-22 17:20:09 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020613_1_ENS2.nc(668832435 bytes)
2026-07-22 17:20:09 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020613_1_ENS3.nc(668832435 bytes)
2026-07-22 17:20:09 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020613_1_ENS4.nc(668832435 bytes)
2026-07-22 17:20:09 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020613_1_ENS5.nc(668832435 bytes)
2026-07-22 17:20:09 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020613_1_ENS6.nc(668832435 bytes)
2026-07-22 17:20:09 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020613_1_ENS7.nc(668832435 bytes)
2026-07-22 17:20:09 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020613_1_ENS8.nc(668832435 bytes)
2026-07-22 17:20:09 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020613_1_ENS9.nc(668832435 bytes)
2026-07-22 17:20:09 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020613_1_ENS10.nc(668832435 bytes)
2026-07-22 17:20:09 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020613_1_ENS11.nc(668832435 bytes)
2026-07-22 17:20:09 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020613_1_ENS12.nc(668832435 bytes)
2026-07-22 17:20:09 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020613_1_ENS13.nc(668832435 bytes)
2026-07-22 17:20:09 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020613_1_ENS14.nc(668832435 bytes)
2026-07-22 17:20:09 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020613_1_ENS15.nc(668832435 bytes)
2026-07-22 17:20:09 INFO  Run_model() completed successfully.
2026-07-22 17:20:09 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 17:20:09 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 14:00:00 days=153072 seconds=50400
2026-07-22 17:20:09 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-22 17:20:09 INFO [TIME] increment current_time 2020-02-06 13:00:00 -> 2020-02-06 14:00:00
2026-07-22 17:20:09 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 17:20:09 INFO ---------->>> Running process_satellite_data()
2026-07-22 17:20:09 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12005.nc
2026-07-22 17:20:09 INFO ---------->>> Running run_obs_converter()
2026-07-22 17:20:09 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_51225_153072.out
2026-07-22 17:20:09 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_51225_153072.out
2026-07-22 17:20:09 INFO ---------->>> Running DART
2026-07-22 17:20:09 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-22 17:20:09 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/chim_ENS1_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/chim_ENS1_2020020614_1_out_toDART.nc
2026-07-22 17:20:09 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/chim_ENS2_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/chim_ENS2_2020020614_1_out_toDART.nc
2026-07-22 17:20:09 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/chim_ENS3_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/chim_ENS3_2020020614_1_out_toDART.nc
2026-07-22 17:20:10 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/chim_ENS4_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/chim_ENS4_2020020614_1_out_toDART.nc
2026-07-22 17:20:10 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/chim_ENS5_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/chim_ENS5_2020020614_1_out_toDART.nc
2026-07-22 17:20:10 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/chim_ENS6_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/chim_ENS6_2020020614_1_out_toDART.nc
2026-07-22 17:20:11 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/chim_ENS7_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/chim_ENS7_2020020614_1_out_toDART.nc
2026-07-22 17:20:11 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/chim_ENS8_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/chim_ENS8_2020020614_1_out_toDART.nc
2026-07-22 17:20:11 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/chim_ENS9_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/chim_ENS9_2020020614_1_out_toDART.nc
2026-07-22 17:20:12 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/chim_ENS10_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/chim_ENS10_2020020614_1_out_toDART.nc
2026-07-22 17:20:12 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/chim_ENS11_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/chim_ENS11_2020020614_1_out_toDART.nc
2026-07-22 17:20:12 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/chim_ENS12_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/chim_ENS12_2020020614_1_out_toDART.nc
2026-07-22 17:20:13 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/chim_ENS13_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/chim_ENS13_2020020614_1_out_toDART.nc
2026-07-22 17:20:13 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/chim_ENS14_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/chim_ENS14_2020020614_1_out_toDART.nc
2026-07-22 17:20:13 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/chim_ENS15_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/chim_ENS15_2020020614_1_out_toDART.nc
2026-07-22 17:20:14 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-22 17:20:14 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-22 17:20:14 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-22 17:20:14 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-22 17:20:14 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-22 17:20:14 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-22 17:20:22 INFO Found: []
2026-07-22 17:20:22 INFO No job id returned by command ./run_filter.bsh
2026-07-22 17:20:22 INFO No monitoring will be performed
2026-07-22 17:20:22 INFO Moving DART output files to analysis and preassim directories for date 2020020614 if present ...
2026-07-22 17:20:22 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 17:20:22 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 17:20:22 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 17:20:22 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 17:20:22 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 17:20:22 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 17:20:22 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 17:20:22 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 17:20:22 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 17:20:22 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 17:20:22 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 17:20:22 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 17:20:22 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 17:20:22 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 17:20:22 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 17:20:22 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 17:20:22 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 17:20:22 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 17:20:22 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 17:20:22 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 17:20:22 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 17:20:22 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 17:20:22 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 17:20:22 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 17:20:22 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 17:20:22 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 17:20:22 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 17:20:22 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 17:20:22 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 17:20:22 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 17:20:22 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 17:20:22 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020614'
2026-07-22 17:20:22 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 17:20:22 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020614'
2026-07-22 17:20:22 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-22 17:20:22 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-22 17:20:22 INFO run_dart() is DONE.
2026-07-22 17:20:22 INFO ---------->>> Running update_pollutant_in_end()
2026-07-22 17:20:23 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:20:23 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS1.nc
2026-07-22 17:20:23 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:20:24 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS1.nc
2026-07-22 17:20:24 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:20:24 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:20:24 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS1.nc
2026-07-22 17:20:25 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:20:25 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS1.nc
2026-07-22 17:20:25 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:20:25 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:20:26 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:20:26 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS2.nc
2026-07-22 17:20:26 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:20:27 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS2.nc
2026-07-22 17:20:27 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:20:27 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:20:27 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS2.nc
2026-07-22 17:20:28 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:20:28 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS2.nc
2026-07-22 17:20:28 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:20:28 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:20:29 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:20:29 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS3.nc
2026-07-22 17:20:29 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:20:30 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS3.nc
2026-07-22 17:20:30 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:20:30 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:20:30 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS3.nc
2026-07-22 17:20:31 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:20:31 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS3.nc
2026-07-22 17:20:31 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:20:31 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:20:32 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:20:32 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS4.nc
2026-07-22 17:20:32 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:20:33 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS4.nc
2026-07-22 17:20:33 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:20:33 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:20:33 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS4.nc
2026-07-22 17:20:34 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:20:35 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS4.nc
2026-07-22 17:20:35 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:20:35 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:20:35 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:20:35 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS5.nc
2026-07-22 17:20:35 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:20:36 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS5.nc
2026-07-22 17:20:36 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:20:37 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:20:37 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS5.nc
2026-07-22 17:20:37 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:20:38 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS5.nc
2026-07-22 17:20:38 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:20:38 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:20:38 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:20:38 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS6.nc
2026-07-22 17:20:39 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:20:39 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS6.nc
2026-07-22 17:20:39 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:20:40 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:20:40 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS6.nc
2026-07-22 17:20:40 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:20:41 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS6.nc
2026-07-22 17:20:41 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:20:41 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:20:41 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:20:41 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS7.nc
2026-07-22 17:20:42 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:20:42 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS7.nc
2026-07-22 17:20:42 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:20:43 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:20:43 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS7.nc
2026-07-22 17:20:43 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:20:44 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS7.nc
2026-07-22 17:20:44 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:20:44 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:20:44 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:20:44 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS8.nc
2026-07-22 17:20:45 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:20:46 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS8.nc
2026-07-22 17:20:46 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:20:46 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:20:46 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS8.nc
2026-07-22 17:20:47 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:20:47 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS8.nc
2026-07-22 17:20:47 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:20:47 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:20:48 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:20:48 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS9.nc
2026-07-22 17:20:48 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:20:49 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS9.nc
2026-07-22 17:20:49 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:20:49 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:20:49 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS9.nc
2026-07-22 17:20:50 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:20:50 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS9.nc
2026-07-22 17:20:50 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:20:50 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:20:50 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:20:50 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS10.nc
2026-07-22 17:20:51 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:20:52 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS10.nc
2026-07-22 17:20:52 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:20:52 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:20:52 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS10.nc
2026-07-22 17:20:53 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:20:54 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS10.nc
2026-07-22 17:20:54 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:20:54 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:20:54 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:20:54 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS11.nc
2026-07-22 17:20:55 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:20:55 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS11.nc
2026-07-22 17:20:55 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:20:56 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:20:56 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS11.nc
2026-07-22 17:20:56 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:20:57 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS11.nc
2026-07-22 17:20:57 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:20:57 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:20:57 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:20:57 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS12.nc
2026-07-22 17:20:58 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:20:58 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS12.nc
2026-07-22 17:20:58 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:20:59 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:20:59 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS12.nc
2026-07-22 17:20:59 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:21:00 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS12.nc
2026-07-22 17:21:00 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:21:00 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:21:00 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:21:00 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS13.nc
2026-07-22 17:21:01 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:21:01 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS13.nc
2026-07-22 17:21:01 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:21:02 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:21:02 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS13.nc
2026-07-22 17:21:02 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:21:03 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS13.nc
2026-07-22 17:21:03 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:21:03 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:21:03 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:21:03 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS14.nc
2026-07-22 17:21:04 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:21:04 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS14.nc
2026-07-22 17:21:04 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:21:05 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:21:05 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS14.nc
2026-07-22 17:21:05 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:21:06 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS14.nc
2026-07-22 17:21:06 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:21:06 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:21:06 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:21:06 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS15.nc
2026-07-22 17:21:07 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:21:07 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS15.nc
2026-07-22 17:21:07 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:21:08 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:21:08 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS15.nc
2026-07-22 17:21:08 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc
2026-07-22 17:21:09 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS15.nc
2026-07-22 17:21:09 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:21:09 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:21:09 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-22 17:21:09 INFO [TIME] step_end current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 17:21:09 INFO [TIME] step_start current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 17:21:09 INFO [TIME] window start=2020-02-06 14:00:00 end=2020-02-07 00:00:00 run_hours=10 has_assimilation=False
2026-07-22 17:21:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:21:10 INFO Hourly dataset computed and listing created
2026-07-22 17:21:26 INFO Hourly dataset computed
2026-07-22 17:21:26 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:21:27 INFO Hourly dataset computed and listing created
2026-07-22 17:21:30 INFO Hourly dataset computed
2026-07-22 17:21:30 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:21:31 INFO Hourly dataset computed and listing created
2026-07-22 17:21:33 INFO Hourly dataset computed
2026-07-22 17:21:33 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:21:34 INFO Hourly dataset computed and listing created
2026-07-22 17:21:37 INFO Hourly dataset computed
2026-07-22 17:21:37 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:21:38 INFO Hourly dataset computed and listing created
2026-07-22 17:21:41 INFO Hourly dataset computed
2026-07-22 17:21:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:21:42 INFO Hourly dataset computed and listing created
2026-07-22 17:21:44 INFO Hourly dataset computed
2026-07-22 17:21:44 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:21:45 INFO Hourly dataset computed and listing created
2026-07-22 17:21:48 INFO Hourly dataset computed
2026-07-22 17:21:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:21:49 INFO Hourly dataset computed and listing created
2026-07-22 17:21:51 INFO Hourly dataset computed
2026-07-22 17:21:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:21:52 INFO Hourly dataset computed and listing created
2026-07-22 17:21:55 INFO Hourly dataset computed
2026-07-22 17:21:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:21:56 INFO Hourly dataset computed and listing created
2026-07-22 17:21:59 INFO Hourly dataset computed
2026-07-22 17:21:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:22:00 INFO Hourly dataset computed and listing created
2026-07-22 17:22:02 INFO Hourly dataset computed
2026-07-22 17:22:02 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:22:04 INFO Hourly dataset computed and listing created
2026-07-22 17:22:06 INFO Hourly dataset computed
2026-07-22 17:22:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:22:07 INFO Hourly dataset computed and listing created
2026-07-22 17:22:10 INFO Hourly dataset computed
2026-07-22 17:22:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:22:11 INFO Hourly dataset computed and listing created
2026-07-22 17:22:13 INFO Hourly dataset computed
2026-07-22 17:22:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:22:14 INFO Hourly dataset computed and listing created
2026-07-22 17:22:17 INFO Hourly dataset computed
2026-07-22 17:22:17 INFO ---------->>> Running CHIMERE model from 2020-02-06 14:00:00 to 2020-02-07 00:00:00
2026-07-22 17:22:17 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:22:17 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1
2026-07-22 17:22:17 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020613_1_ENS1.nc
2026-07-22 17:22:17 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-22 17:22:17 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:22:17 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-22 17:22:17 INFO Queuing job for member 1...
2026-07-22 17:22:17 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:22:17 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-22 17:22:18 INFO Found: ['5257080']
2026-07-22 17:22:23 INFO [TGCC-IRENE] Submitted job with ID:['5257080']
2026-07-22 17:22:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:22:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2
2026-07-22 17:22:23 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020613_1_ENS2.nc
2026-07-22 17:22:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-22 17:22:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:22:23 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-22 17:22:23 INFO Queuing job for member 2...
2026-07-22 17:22:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:22:23 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-22 17:22:24 INFO Found: ['5257086']
2026-07-22 17:22:29 INFO [TGCC-IRENE] Submitted job with ID:['5257086']
2026-07-22 17:22:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:22:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3
2026-07-22 17:22:29 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020613_1_ENS3.nc
2026-07-22 17:22:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-22 17:22:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:22:29 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-22 17:22:29 INFO Queuing job for member 3...
2026-07-22 17:22:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:22:29 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-22 17:22:30 INFO Found: ['5257088']
2026-07-22 17:22:35 INFO [TGCC-IRENE] Submitted job with ID:['5257088']
2026-07-22 17:22:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:22:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4
2026-07-22 17:22:35 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020613_1_ENS4.nc
2026-07-22 17:22:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-22 17:22:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:22:35 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-22 17:22:35 INFO Queuing job for member 4...
2026-07-22 17:22:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:22:35 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-22 17:22:37 INFO Found: ['5257089']
2026-07-22 17:22:42 INFO [TGCC-IRENE] Submitted job with ID:['5257089']
2026-07-22 17:22:42 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:22:42 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5
2026-07-22 17:22:42 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020613_1_ENS5.nc
2026-07-22 17:22:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-22 17:22:42 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:22:42 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-22 17:22:42 INFO Queuing job for member 5...
2026-07-22 17:22:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:22:42 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-22 17:22:45 INFO Found: ['5257091']
2026-07-22 17:22:50 INFO [TGCC-IRENE] Submitted job with ID:['5257091']
2026-07-22 17:22:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:22:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6
2026-07-22 17:22:50 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020613_1_ENS6.nc
2026-07-22 17:22:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-22 17:22:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:22:50 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-22 17:22:50 INFO Queuing job for member 6...
2026-07-22 17:22:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:22:50 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-22 17:22:52 INFO Found: ['5257092']
2026-07-22 17:22:57 INFO [TGCC-IRENE] Submitted job with ID:['5257092']
2026-07-22 17:22:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:22:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7
2026-07-22 17:22:57 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020613_1_ENS7.nc
2026-07-22 17:22:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-22 17:22:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:22:57 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-22 17:22:57 INFO Queuing job for member 7...
2026-07-22 17:22:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:22:57 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-22 17:23:00 INFO Found: ['5257094']
2026-07-22 17:23:05 INFO [TGCC-IRENE] Submitted job with ID:['5257094']
2026-07-22 17:23:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:23:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8
2026-07-22 17:23:05 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020613_1_ENS8.nc
2026-07-22 17:23:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-22 17:23:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:23:05 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-22 17:23:05 INFO Queuing job for member 8...
2026-07-22 17:23:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:23:05 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-22 17:23:07 INFO Found: ['5257096']
2026-07-22 17:23:12 INFO [TGCC-IRENE] Submitted job with ID:['5257096']
2026-07-22 17:23:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:23:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9
2026-07-22 17:23:12 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020613_1_ENS9.nc
2026-07-22 17:23:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-22 17:23:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:23:12 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-22 17:23:12 INFO Queuing job for member 9...
2026-07-22 17:23:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:23:12 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-22 17:23:15 INFO Found: ['5257104']
2026-07-22 17:23:20 INFO [TGCC-IRENE] Submitted job with ID:['5257104']
2026-07-22 17:23:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:23:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10
2026-07-22 17:23:20 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020613_1_ENS10.nc
2026-07-22 17:23:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-22 17:23:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:23:20 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-22 17:23:20 INFO Queuing job for member 10...
2026-07-22 17:23:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:23:20 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-22 17:23:21 INFO Found: ['5257106']
2026-07-22 17:23:26 INFO [TGCC-IRENE] Submitted job with ID:['5257106']
2026-07-22 17:23:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:23:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11
2026-07-22 17:23:26 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020613_1_ENS11.nc
2026-07-22 17:23:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-22 17:23:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:23:26 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-22 17:23:26 INFO Queuing job for member 11...
2026-07-22 17:23:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:23:26 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-22 17:23:27 INFO Found: ['5257108']
2026-07-22 17:23:32 INFO [TGCC-IRENE] Submitted job with ID:['5257108']
2026-07-22 17:23:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:23:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12
2026-07-22 17:23:32 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020613_1_ENS12.nc
2026-07-22 17:23:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-22 17:23:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:23:32 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-22 17:23:32 INFO Queuing job for member 12...
2026-07-22 17:23:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:23:32 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-22 17:23:32 INFO Found: ['5257109']
2026-07-22 17:23:37 INFO [TGCC-IRENE] Submitted job with ID:['5257109']
2026-07-22 17:23:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:23:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13
2026-07-22 17:23:37 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020613_1_ENS13.nc
2026-07-22 17:23:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-22 17:23:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:23:37 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-22 17:23:37 INFO Queuing job for member 13...
2026-07-22 17:23:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:23:37 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-22 17:23:38 INFO Found: ['5257110']
2026-07-22 17:23:43 INFO [TGCC-IRENE] Submitted job with ID:['5257110']
2026-07-22 17:23:43 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:23:43 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14
2026-07-22 17:23:43 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020613_1_ENS14.nc
2026-07-22 17:23:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-22 17:23:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:23:43 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-22 17:23:43 INFO Queuing job for member 14...
2026-07-22 17:23:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:23:43 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-22 17:23:44 INFO Found: ['5257112']
2026-07-22 17:23:49 INFO [TGCC-IRENE] Submitted job with ID:['5257112']
2026-07-22 17:23:49 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:23:49 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15
2026-07-22 17:23:49 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020613_1_ENS15.nc
2026-07-22 17:23:49 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-22 17:23:49 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:23:49 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-22 17:23:49 INFO Queuing job for member 15...
2026-07-22 17:23:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:23:49 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-22 17:23:50 INFO Found: ['5257113']
2026-07-22 17:23:55 INFO [TGCC-IRENE] Submitted job with ID:['5257113']
2026-07-22 17:23:55 INFO Checking job status ...
2026-07-22 17:23:55 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:23:55 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:23:55 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:23:55 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:23:55 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:23:55 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:23:55 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:23:55 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:23:55 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:23:55 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:23:55 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:23:55 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:23:55 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:23:55 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:23:55 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:23:55 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:24:12 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:24:12 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:24:12 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:24:12 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:24:12 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:24:12 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:24:12 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:24:12 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:24:12 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:24:12 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:24:12 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:24:12 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:24:12 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:24:12 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:24:12 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:24:12 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:24:27 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:24:27 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:24:27 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:24:27 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:24:27 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:24:27 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:24:27 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:24:27 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:24:27 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:24:27 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:24:27 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:24:29 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:24:29 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:24:29 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:24:32 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:24:32 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:24:47 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:24:47 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:24:47 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:24:47 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:24:47 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:24:47 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:24:47 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:24:47 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:24:47 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:24:47 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:24:47 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:24:47 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:24:47 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:24:47 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:24:47 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:24:47 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:25:02 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:25:02 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:25:02 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:25:02 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:25:02 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:25:02 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:25:02 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:25:02 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:25:02 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:25:02 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:25:02 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:25:02 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:25:02 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:25:02 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:25:02 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:25:02 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:25:17 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:25:17 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:25:17 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:25:17 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:25:18 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:25:18 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:25:18 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:25:18 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:25:18 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:25:18 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:25:18 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:25:18 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:25:18 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:25:18 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:25:18 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:25:18 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:25:33 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:25:33 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:25:33 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:25:33 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:25:33 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:25:33 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:25:33 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:25:33 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:25:33 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:25:33 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:25:33 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:25:33 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:25:33 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:25:33 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:25:33 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:25:33 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:25:48 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:25:48 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:25:48 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:25:48 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:25:48 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:25:50 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:25:50 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:25:50 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:25:50 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:25:50 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:25:50 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:25:50 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:25:50 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:25:50 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:25:50 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:25:50 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:26:05 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:26:05 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:26:05 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:26:05 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:26:05 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:26:05 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:26:05 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:26:05 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:26:06 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:26:06 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:26:06 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:26:06 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:26:06 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:26:06 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:26:06 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:26:06 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:26:21 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:26:21 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:26:21 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:26:21 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:26:21 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:26:21 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:26:21 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:26:21 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:26:21 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:26:21 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:26:21 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:26:21 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:26:21 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:26:21 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:26:21 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:26:21 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:26:36 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:26:36 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:26:36 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:26:36 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:26:36 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:26:37 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:26:37 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:26:37 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:26:37 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:26:37 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:26:37 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:26:37 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:26:37 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:26:37 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:26:37 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:26:37 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:26:52 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:26:52 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:26:52 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:26:52 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:26:52 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:26:52 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:26:52 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:26:52 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:26:52 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:26:52 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:26:52 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:26:52 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:26:52 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:26:52 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:26:52 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:26:52 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:27:08 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:27:08 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:27:08 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:27:08 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:27:09 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:27:09 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:27:09 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:27:09 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:27:09 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:27:09 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:27:09 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:27:09 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:27:09 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:27:09 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:27:09 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:27:09 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:27:24 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:27:24 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:27:24 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:27:24 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:27:24 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:27:24 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:27:24 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:27:24 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:27:24 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:27:24 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:27:26 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:27:26 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:27:26 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:27:26 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:27:26 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:27:26 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:27:41 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:27:41 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:27:41 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:27:41 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:27:41 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:27:41 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:27:41 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:27:41 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:27:41 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:27:41 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:27:41 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:27:41 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:27:41 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:27:41 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:27:41 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:27:41 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:27:56 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:27:56 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:27:56 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:27:56 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:27:56 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:27:56 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:27:56 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:27:56 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:27:57 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:27:57 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:27:57 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:27:57 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:27:57 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:27:57 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:27:57 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:27:57 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:28:12 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:28:12 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:28:12 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:28:12 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:28:12 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:28:12 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:28:12 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:28:12 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:28:12 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:28:12 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:28:12 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:28:12 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:28:12 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:28:12 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:28:12 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:28:12 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:28:27 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:28:27 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:28:27 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:28:27 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:28:27 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:28:27 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:28:27 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:28:27 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:28:27 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:28:27 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:28:27 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:28:27 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:28:27 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:28:27 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:28:27 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:28:27 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:28:42 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:28:42 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:28:42 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:28:42 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:28:42 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:28:42 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:28:42 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:28:42 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:28:42 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:28:42 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:28:43 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:28:43 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:28:43 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:28:43 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:28:43 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:28:43 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:28:58 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:28:58 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:28:58 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:28:58 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:28:58 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:28:58 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:28:58 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:28:58 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:28:58 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:28:58 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:28:58 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:28:58 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:28:58 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:29:00 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:29:00 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:29:00 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:29:15 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:29:15 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:29:15 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:29:15 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:29:15 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:29:15 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:29:15 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:29:15 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:29:15 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:29:15 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:29:15 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:29:15 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:29:15 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:29:15 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:29:15 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:29:15 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:29:30 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:29:30 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:29:30 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:29:30 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:29:30 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:29:30 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:29:30 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:29:31 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:29:31 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:29:31 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:29:31 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:29:31 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:29:31 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:29:31 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:29:31 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:29:31 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:29:46 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:29:46 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:29:46 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:29:46 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:29:46 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:29:46 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:29:46 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:29:46 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:29:46 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:29:46 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:29:46 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:29:46 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:29:46 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:29:46 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:29:46 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:29:46 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:30:01 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:30:01 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:30:01 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:30:01 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:30:01 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:30:01 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:30:01 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:30:01 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:30:01 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:30:01 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:30:01 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:30:01 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:30:01 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:30:01 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:30:01 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:30:01 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:30:17 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:30:17 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:30:17 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:30:17 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:30:17 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:30:17 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:30:17 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:30:17 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:30:17 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:30:17 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:30:17 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:30:17 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:30:17 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:30:17 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:30:17 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:30:17 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:30:32 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:30:32 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:30:32 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:30:32 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:30:32 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:30:32 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:30:32 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:30:32 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:30:32 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:30:32 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:30:32 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:30:32 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:30:32 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:30:34 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:30:34 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:30:34 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:30:49 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:30:49 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:30:49 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:30:50 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:30:50 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:30:50 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:30:50 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:30:50 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:30:50 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:30:50 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:30:50 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:30:50 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:30:50 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:30:50 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:30:50 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:30:50 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:31:05 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:31:05 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:31:05 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:31:05 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:31:05 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:31:05 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:31:05 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:31:05 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:31:05 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:31:05 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:31:05 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:31:05 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:31:05 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:31:05 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:31:05 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:31:05 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:31:20 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:31:20 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:31:20 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:31:20 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:31:20 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:31:20 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:31:20 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:31:20 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:31:20 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:31:20 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:31:20 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:31:20 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:31:20 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:31:20 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:31:20 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:31:20 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:31:35 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:31:35 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:31:35 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:31:35 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:31:35 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:31:35 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:31:35 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:31:36 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:31:36 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:31:36 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:31:36 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:31:36 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:31:36 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:31:36 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:31:36 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:31:36 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:31:51 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:31:51 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:31:51 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:31:51 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:31:53 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:31:53 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:31:53 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:31:53 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:31:53 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:31:53 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:31:53 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:31:53 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:31:53 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:31:53 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:31:53 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:31:53 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:32:08 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:32:08 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:32:08 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:32:08 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:32:08 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:32:08 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:32:08 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:32:08 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:32:08 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:32:08 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:32:08 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:32:08 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:32:08 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:32:08 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:32:08 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:32:08 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:32:23 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:32:23 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:32:23 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:32:23 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:32:23 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:32:25 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:32:25 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:32:25 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:32:25 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:32:25 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:32:25 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:32:26 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:32:26 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:32:26 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:32:26 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:32:26 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:32:41 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:32:41 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:32:41 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:32:41 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:32:41 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:32:41 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:32:41 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:32:41 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:32:41 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:32:41 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:32:41 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:32:41 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:32:41 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:32:41 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:32:41 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:32:41 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:32:56 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:32:56 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:32:56 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:32:56 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:32:56 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:32:56 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:32:56 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:32:56 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:32:56 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:32:56 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:32:56 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:32:56 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:32:56 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:32:56 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:32:56 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:32:56 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:33:11 INFO None 5257080: status RUNNING/PENDING
2026-07-22 17:33:11 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:33:11 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:33:11 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:33:11 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:33:11 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:33:11 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:33:11 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:33:11 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:33:11 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:33:11 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:33:11 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:33:11 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:33:11 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:33:11 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:33:11 INFO Jobs still running: ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:33:28 INFO None 5257080: status FINISHED
2026-07-22 17:33:28 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:33:28 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:33:28 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:33:28 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:33:28 INFO None 5257092: status RUNNING/PENDING
2026-07-22 17:33:29 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:33:29 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:33:29 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:33:29 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:33:29 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:33:29 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:33:29 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:33:29 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:33:29 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:33:29 INFO Jobs still running: ['5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:33:44 INFO None 5257080: status FINISHED
2026-07-22 17:33:44 INFO None 5257086: status RUNNING/PENDING
2026-07-22 17:33:44 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:33:44 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:33:44 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:33:44 INFO None 5257092: status FINISHED
2026-07-22 17:33:44 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:33:44 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:33:44 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:33:44 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:33:44 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:33:44 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:33:44 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:33:46 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:33:46 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:33:46 INFO Jobs still running: ['5257086', '5257088', '5257089', '5257091', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:34:01 INFO None 5257080: status FINISHED
2026-07-22 17:34:01 INFO None 5257086: status FINISHED
2026-07-22 17:34:01 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:34:01 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:34:01 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:34:01 INFO None 5257092: status FINISHED
2026-07-22 17:34:01 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:34:01 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:34:01 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:34:01 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:34:01 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:34:01 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:34:01 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:34:01 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:34:01 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:34:01 INFO Jobs still running: ['5257088', '5257089', '5257091', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:34:16 INFO None 5257080: status FINISHED
2026-07-22 17:34:16 INFO None 5257086: status FINISHED
2026-07-22 17:34:16 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:34:16 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:34:17 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:34:17 INFO None 5257092: status FINISHED
2026-07-22 17:34:17 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:34:17 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:34:17 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:34:17 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:34:17 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:34:17 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:34:17 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:34:17 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:34:17 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:34:17 INFO Jobs still running: ['5257088', '5257089', '5257091', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:34:32 INFO None 5257080: status FINISHED
2026-07-22 17:34:32 INFO None 5257086: status FINISHED
2026-07-22 17:34:32 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:34:32 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:34:32 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:34:32 INFO None 5257092: status FINISHED
2026-07-22 17:34:32 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:34:32 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:34:32 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:34:32 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:34:32 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:34:32 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:34:32 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:34:32 INFO None 5257112: status RUNNING/PENDING
2026-07-22 17:34:32 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:34:32 INFO Jobs still running: ['5257088', '5257089', '5257091', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113']. Waiting...
2026-07-22 17:34:47 INFO None 5257080: status FINISHED
2026-07-22 17:34:47 INFO None 5257086: status FINISHED
2026-07-22 17:34:47 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:34:47 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:34:47 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:34:47 INFO None 5257092: status FINISHED
2026-07-22 17:34:47 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:34:47 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:34:47 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:34:47 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:34:47 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:34:47 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:34:47 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:34:47 INFO None 5257112: status FINISHED
2026-07-22 17:34:47 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:34:47 INFO Jobs still running: ['5257088', '5257089', '5257091', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257113']. Waiting...
2026-07-22 17:35:04 INFO None 5257080: status FINISHED
2026-07-22 17:35:04 INFO None 5257086: status FINISHED
2026-07-22 17:35:04 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:35:04 INFO None 5257089: status RUNNING/PENDING
2026-07-22 17:35:04 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:35:04 INFO None 5257092: status FINISHED
2026-07-22 17:35:04 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:35:04 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:35:04 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:35:04 INFO None 5257106: status RUNNING/PENDING
2026-07-22 17:35:04 INFO None 5257108: status RUNNING/PENDING
2026-07-22 17:35:04 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:35:04 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:35:04 INFO None 5257112: status FINISHED
2026-07-22 17:35:04 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:35:04 INFO Jobs still running: ['5257088', '5257089', '5257091', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257113']. Waiting...
2026-07-22 17:35:19 INFO None 5257080: status FINISHED
2026-07-22 17:35:19 INFO None 5257086: status FINISHED
2026-07-22 17:35:19 INFO None 5257088: status RUNNING/PENDING
2026-07-22 17:35:19 INFO None 5257089: status FINISHED
2026-07-22 17:35:19 INFO None 5257091: status RUNNING/PENDING
2026-07-22 17:35:19 INFO None 5257092: status FINISHED
2026-07-22 17:35:19 INFO None 5257094: status RUNNING/PENDING
2026-07-22 17:35:19 INFO None 5257096: status RUNNING/PENDING
2026-07-22 17:35:19 INFO None 5257104: status RUNNING/PENDING
2026-07-22 17:35:19 INFO None 5257106: status FINISHED
2026-07-22 17:35:19 INFO None 5257108: status FINISHED
2026-07-22 17:35:19 INFO None 5257109: status RUNNING/PENDING
2026-07-22 17:35:19 INFO None 5257110: status RUNNING/PENDING
2026-07-22 17:35:21 INFO None 5257112: status FINISHED
2026-07-22 17:35:21 INFO None 5257113: status RUNNING/PENDING
2026-07-22 17:35:21 INFO Jobs still running: ['5257088', '5257091', '5257094', '5257096', '5257104', '5257109', '5257110', '5257113']. Waiting...
2026-07-22 17:35:36 INFO None 5257080: status FINISHED
2026-07-22 17:35:36 INFO None 5257086: status FINISHED
2026-07-22 17:35:36 INFO None 5257088: status FINISHED
2026-07-22 17:35:36 INFO None 5257089: status FINISHED
2026-07-22 17:35:36 INFO None 5257091: status FINISHED
2026-07-22 17:35:36 INFO None 5257092: status FINISHED
2026-07-22 17:35:36 INFO None 5257094: status FINISHED
2026-07-22 17:35:36 INFO None 5257096: status FINISHED
2026-07-22 17:35:37 INFO None 5257104: status FINISHED
2026-07-22 17:35:37 INFO None 5257106: status FINISHED
2026-07-22 17:35:37 INFO None 5257108: status FINISHED
2026-07-22 17:35:37 INFO None 5257109: status FINISHED
2026-07-22 17:35:37 INFO None 5257110: status FINISHED
2026-07-22 17:35:37 INFO None 5257112: status FINISHED
2026-07-22 17:35:37 INFO None 5257113: status FINISHED
2026-07-22 17:35:37 INFO Jobs ['5257080', '5257086', '5257088', '5257089', '5257091', '5257092', '5257094', '5257096', '5257104', '5257106', '5257108', '5257109', '5257110', '5257112', '5257113'] have finished
2026-07-22 17:35:37 INFO Checking restart files were created ...
2026-07-22 17:35:37 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020614_10_ENS1.nc(3673513755 bytes)
2026-07-22 17:35:37 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020614_10_ENS2.nc(3673513755 bytes)
2026-07-22 17:35:37 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020614_10_ENS3.nc(3673513755 bytes)
2026-07-22 17:35:37 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020614_10_ENS4.nc(3673513755 bytes)
2026-07-22 17:35:37 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020614_10_ENS5.nc(3673513755 bytes)
2026-07-22 17:35:37 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020614_10_ENS6.nc(3673513755 bytes)
2026-07-22 17:35:37 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020614_10_ENS7.nc(3673513755 bytes)
2026-07-22 17:35:37 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020614_10_ENS8.nc(3673513755 bytes)
2026-07-22 17:35:37 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020614_10_ENS9.nc(3673513755 bytes)
2026-07-22 17:35:37 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020614_10_ENS10.nc(3673513755 bytes)
2026-07-22 17:35:37 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020614_10_ENS11.nc(3673513755 bytes)
2026-07-22 17:35:37 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020614_10_ENS12.nc(3673513755 bytes)
2026-07-22 17:35:37 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020614_10_ENS13.nc(3673513755 bytes)
2026-07-22 17:35:37 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020614_10_ENS14.nc(3673513755 bytes)
2026-07-22 17:35:37 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020614_10_ENS15.nc(3673513755 bytes)
2026-07-22 17:35:37 INFO  Run_model() completed successfully.
2026-07-22 17:35:37 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 14:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 17:35:37 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 00:00:00 days=153073 seconds=0
2026-07-22 17:35:37 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-22 17:35:37 INFO [TIME] increment current_time 2020-02-06 14:00:00 -> 2020-02-07 00:00:00
2026-07-22 17:35:37 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 00:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 17:35:37 INFO ---------->>> Running process_satellite_data()
2026-07-22 17:35:37 INFO [DART] No satellite data found, skipping assimilation
2026-07-22 17:35:37 INFO after_assimilation() skipped
2026-07-22 17:35:37 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-22 17:35:37 INFO [TIME] step_end current_time=2020-02-07 00:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 17:35:37 INFO [TIME] step_start current_time=2020-02-07 00:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 17:35:37 INFO [TIME] window start=2020-02-07 00:00:00 end=2020-02-07 01:00:00 run_hours=1 has_assimilation=False
2026-07-22 17:35:37 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:35:38 INFO Hourly dataset computed and listing created
2026-07-22 17:35:51 INFO Hourly dataset computed
2026-07-22 17:35:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:35:52 INFO Hourly dataset computed and listing created
2026-07-22 17:35:53 INFO Hourly dataset computed
2026-07-22 17:35:53 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:35:54 INFO Hourly dataset computed and listing created
2026-07-22 17:35:54 INFO Hourly dataset computed
2026-07-22 17:35:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:35:55 INFO Hourly dataset computed and listing created
2026-07-22 17:35:55 INFO Hourly dataset computed
2026-07-22 17:35:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:35:56 INFO Hourly dataset computed and listing created
2026-07-22 17:35:57 INFO Hourly dataset computed
2026-07-22 17:35:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:35:58 INFO Hourly dataset computed and listing created
2026-07-22 17:35:58 INFO Hourly dataset computed
2026-07-22 17:35:58 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:35:59 INFO Hourly dataset computed and listing created
2026-07-22 17:35:59 INFO Hourly dataset computed
2026-07-22 17:35:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:36:00 INFO Hourly dataset computed and listing created
2026-07-22 17:36:01 INFO Hourly dataset computed
2026-07-22 17:36:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:36:02 INFO Hourly dataset computed and listing created
2026-07-22 17:36:02 INFO Hourly dataset computed
2026-07-22 17:36:02 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:36:03 INFO Hourly dataset computed and listing created
2026-07-22 17:36:04 INFO Hourly dataset computed
2026-07-22 17:36:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:36:04 INFO Hourly dataset computed and listing created
2026-07-22 17:36:05 INFO Hourly dataset computed
2026-07-22 17:36:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:36:06 INFO Hourly dataset computed and listing created
2026-07-22 17:36:06 INFO Hourly dataset computed
2026-07-22 17:36:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:36:07 INFO Hourly dataset computed and listing created
2026-07-22 17:36:08 INFO Hourly dataset computed
2026-07-22 17:36:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:36:08 INFO Hourly dataset computed and listing created
2026-07-22 17:36:09 INFO Hourly dataset computed
2026-07-22 17:36:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:36:10 INFO Hourly dataset computed and listing created
2026-07-22 17:36:10 INFO Hourly dataset computed
2026-07-22 17:36:10 INFO ---------->>> Running CHIMERE model from 2020-02-07 00:00:00 to 2020-02-07 01:00:00
2026-07-22 17:36:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:36:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1
2026-07-22 17:36:10 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020614_10_ENS1.nc
2026-07-22 17:36:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-22 17:36:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:36:10 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-22 17:36:11 INFO Queuing job for member 1...
2026-07-22 17:36:11 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:36:11 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-22 17:36:11 INFO Found: ['5257213']
2026-07-22 17:36:16 INFO [TGCC-IRENE] Submitted job with ID:['5257213']
2026-07-22 17:36:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:36:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2
2026-07-22 17:36:16 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020614_10_ENS2.nc
2026-07-22 17:36:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-22 17:36:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:36:16 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-22 17:36:16 INFO Queuing job for member 2...
2026-07-22 17:36:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:36:16 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-22 17:36:17 INFO Found: ['5257215']
2026-07-22 17:36:22 INFO [TGCC-IRENE] Submitted job with ID:['5257215']
2026-07-22 17:36:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:36:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3
2026-07-22 17:36:22 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020614_10_ENS3.nc
2026-07-22 17:36:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-22 17:36:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:36:22 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-22 17:36:22 INFO Queuing job for member 3...
2026-07-22 17:36:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:36:22 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-22 17:36:23 INFO Found: ['5257218']
2026-07-22 17:36:28 INFO [TGCC-IRENE] Submitted job with ID:['5257218']
2026-07-22 17:36:28 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:36:28 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4
2026-07-22 17:36:28 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020614_10_ENS4.nc
2026-07-22 17:36:28 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-22 17:36:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:36:28 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-22 17:36:28 INFO Queuing job for member 4...
2026-07-22 17:36:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:36:28 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-22 17:36:30 INFO Found: ['5257220']
2026-07-22 17:36:35 INFO [TGCC-IRENE] Submitted job with ID:['5257220']
2026-07-22 17:36:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:36:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5
2026-07-22 17:36:35 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020614_10_ENS5.nc
2026-07-22 17:36:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-22 17:36:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:36:35 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-22 17:36:35 INFO Queuing job for member 5...
2026-07-22 17:36:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:36:35 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-22 17:36:37 INFO Found: ['5257223']
2026-07-22 17:36:42 INFO [TGCC-IRENE] Submitted job with ID:['5257223']
2026-07-22 17:36:42 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:36:42 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6
2026-07-22 17:36:42 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020614_10_ENS6.nc
2026-07-22 17:36:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-22 17:36:42 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:36:42 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-22 17:36:42 INFO Queuing job for member 6...
2026-07-22 17:36:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:36:42 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-22 17:36:45 INFO Found: ['5257224']
2026-07-22 17:36:50 INFO [TGCC-IRENE] Submitted job with ID:['5257224']
2026-07-22 17:36:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:36:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7
2026-07-22 17:36:50 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020614_10_ENS7.nc
2026-07-22 17:36:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-22 17:36:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:36:50 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-22 17:36:50 INFO Queuing job for member 7...
2026-07-22 17:36:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:36:50 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-22 17:36:52 INFO Found: ['5257225']
2026-07-22 17:36:57 INFO [TGCC-IRENE] Submitted job with ID:['5257225']
2026-07-22 17:36:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:36:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8
2026-07-22 17:36:57 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020614_10_ENS8.nc
2026-07-22 17:36:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-22 17:36:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:36:57 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-22 17:36:57 INFO Queuing job for member 8...
2026-07-22 17:36:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:36:57 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-22 17:37:00 INFO Found: ['5257226']
2026-07-22 17:37:05 INFO [TGCC-IRENE] Submitted job with ID:['5257226']
2026-07-22 17:37:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:37:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9
2026-07-22 17:37:05 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020614_10_ENS9.nc
2026-07-22 17:37:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-22 17:37:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:37:05 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-22 17:37:05 INFO Queuing job for member 9...
2026-07-22 17:37:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:37:05 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-22 17:37:07 INFO Found: ['5257229']
2026-07-22 17:37:12 INFO [TGCC-IRENE] Submitted job with ID:['5257229']
2026-07-22 17:37:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:37:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10
2026-07-22 17:37:12 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020614_10_ENS10.nc
2026-07-22 17:37:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-22 17:37:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:37:12 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-22 17:37:12 INFO Queuing job for member 10...
2026-07-22 17:37:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:37:12 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-22 17:37:15 INFO Found: ['5257231']
2026-07-22 17:37:20 INFO [TGCC-IRENE] Submitted job with ID:['5257231']
2026-07-22 17:37:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:37:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11
2026-07-22 17:37:20 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020614_10_ENS11.nc
2026-07-22 17:37:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-22 17:37:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:37:20 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-22 17:37:20 INFO Queuing job for member 11...
2026-07-22 17:37:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:37:20 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-22 17:37:21 INFO Found: ['5257232']
2026-07-22 17:37:26 INFO [TGCC-IRENE] Submitted job with ID:['5257232']
2026-07-22 17:37:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:37:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12
2026-07-22 17:37:26 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020614_10_ENS12.nc
2026-07-22 17:37:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-22 17:37:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:37:26 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-22 17:37:26 INFO Queuing job for member 12...
2026-07-22 17:37:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:37:26 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-22 17:37:26 INFO Found: ['5257234']
2026-07-22 17:37:31 INFO [TGCC-IRENE] Submitted job with ID:['5257234']
2026-07-22 17:37:31 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:37:31 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13
2026-07-22 17:37:31 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020614_10_ENS13.nc
2026-07-22 17:37:31 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-22 17:37:31 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:37:31 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-22 17:37:31 INFO Queuing job for member 13...
2026-07-22 17:37:31 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:37:31 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-22 17:37:32 INFO Found: ['5257236']
2026-07-22 17:37:37 INFO [TGCC-IRENE] Submitted job with ID:['5257236']
2026-07-22 17:37:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:37:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14
2026-07-22 17:37:37 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020614_10_ENS14.nc
2026-07-22 17:37:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-22 17:37:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:37:37 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-22 17:37:37 INFO Queuing job for member 14...
2026-07-22 17:37:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:37:37 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-22 17:37:38 INFO Found: ['5257237']
2026-07-22 17:37:43 INFO [TGCC-IRENE] Submitted job with ID:['5257237']
2026-07-22 17:37:43 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:37:43 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15
2026-07-22 17:37:43 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020614_10_ENS15.nc
2026-07-22 17:37:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-22 17:37:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:37:43 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-22 17:37:43 INFO Queuing job for member 15...
2026-07-22 17:37:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:37:43 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-22 17:37:44 INFO Found: ['5257238']
2026-07-22 17:37:49 INFO [TGCC-IRENE] Submitted job with ID:['5257238']
2026-07-22 17:37:49 INFO Checking job status ...
2026-07-22 17:37:49 INFO None 5257213: status RUNNING/PENDING
2026-07-22 17:37:49 INFO None 5257215: status RUNNING/PENDING
2026-07-22 17:37:49 INFO None 5257218: status RUNNING/PENDING
2026-07-22 17:37:49 INFO None 5257220: status RUNNING/PENDING
2026-07-22 17:37:49 INFO None 5257223: status RUNNING/PENDING
2026-07-22 17:37:49 INFO None 5257224: status RUNNING/PENDING
2026-07-22 17:37:49 INFO None 5257225: status RUNNING/PENDING
2026-07-22 17:37:49 INFO None 5257226: status RUNNING/PENDING
2026-07-22 17:37:49 INFO None 5257229: status RUNNING/PENDING
2026-07-22 17:37:49 INFO None 5257231: status RUNNING/PENDING
2026-07-22 17:37:49 INFO None 5257232: status RUNNING/PENDING
2026-07-22 17:37:49 INFO None 5257234: status RUNNING/PENDING
2026-07-22 17:37:49 INFO None 5257236: status RUNNING/PENDING
2026-07-22 17:37:49 INFO None 5257237: status RUNNING/PENDING
2026-07-22 17:37:49 INFO None 5257238: status RUNNING/PENDING
2026-07-22 17:37:49 INFO Jobs still running: ['5257213', '5257215', '5257218', '5257220', '5257223', '5257224', '5257225', '5257226', '5257229', '5257231', '5257232', '5257234', '5257236', '5257237', '5257238']. Waiting...
2026-07-22 17:38:04 INFO None 5257213: status RUNNING/PENDING
2026-07-22 17:38:04 INFO None 5257215: status RUNNING/PENDING
2026-07-22 17:38:04 INFO None 5257218: status RUNNING/PENDING
2026-07-22 17:38:04 INFO None 5257220: status RUNNING/PENDING
2026-07-22 17:38:04 INFO None 5257223: status RUNNING/PENDING
2026-07-22 17:38:04 INFO None 5257224: status RUNNING/PENDING
2026-07-22 17:38:04 INFO None 5257225: status RUNNING/PENDING
2026-07-22 17:38:04 INFO None 5257226: status RUNNING/PENDING
2026-07-22 17:38:04 INFO None 5257229: status RUNNING/PENDING
2026-07-22 17:38:04 INFO None 5257231: status RUNNING/PENDING
2026-07-22 17:38:04 INFO None 5257232: status RUNNING/PENDING
2026-07-22 17:38:04 INFO None 5257234: status RUNNING/PENDING
2026-07-22 17:38:04 INFO None 5257236: status RUNNING/PENDING
2026-07-22 17:38:04 INFO None 5257237: status RUNNING/PENDING
2026-07-22 17:38:04 INFO None 5257238: status RUNNING/PENDING
2026-07-22 17:38:04 INFO Jobs still running: ['5257213', '5257215', '5257218', '5257220', '5257223', '5257224', '5257225', '5257226', '5257229', '5257231', '5257232', '5257234', '5257236', '5257237', '5257238']. Waiting...
2026-07-22 17:38:19 INFO None 5257213: status RUNNING/PENDING
2026-07-22 17:38:19 INFO None 5257215: status RUNNING/PENDING
2026-07-22 17:38:19 INFO None 5257218: status RUNNING/PENDING
2026-07-22 17:38:20 INFO None 5257220: status RUNNING/PENDING
2026-07-22 17:38:20 INFO None 5257223: status RUNNING/PENDING
2026-07-22 17:38:20 INFO None 5257224: status RUNNING/PENDING
2026-07-22 17:38:20 INFO None 5257225: status RUNNING/PENDING
2026-07-22 17:38:20 INFO None 5257226: status RUNNING/PENDING
2026-07-22 17:38:20 INFO None 5257229: status RUNNING/PENDING
2026-07-22 17:38:20 INFO None 5257231: status RUNNING/PENDING
2026-07-22 17:38:20 INFO None 5257232: status RUNNING/PENDING
2026-07-22 17:38:20 INFO None 5257234: status RUNNING/PENDING
2026-07-22 17:38:20 INFO None 5257236: status RUNNING/PENDING
2026-07-22 17:38:20 INFO None 5257237: status RUNNING/PENDING
2026-07-22 17:38:20 INFO None 5257238: status RUNNING/PENDING
2026-07-22 17:38:20 INFO Jobs still running: ['5257213', '5257215', '5257218', '5257220', '5257223', '5257224', '5257225', '5257226', '5257229', '5257231', '5257232', '5257234', '5257236', '5257237', '5257238']. Waiting...
2026-07-22 17:38:35 INFO None 5257213: status RUNNING/PENDING
2026-07-22 17:38:35 INFO None 5257215: status RUNNING/PENDING
2026-07-22 17:38:35 INFO None 5257218: status RUNNING/PENDING
2026-07-22 17:38:35 INFO None 5257220: status RUNNING/PENDING
2026-07-22 17:38:35 INFO None 5257223: status RUNNING/PENDING
2026-07-22 17:38:37 INFO None 5257224: status RUNNING/PENDING
2026-07-22 17:38:37 INFO None 5257225: status RUNNING/PENDING
2026-07-22 17:38:37 INFO None 5257226: status RUNNING/PENDING
2026-07-22 17:38:37 INFO None 5257229: status RUNNING/PENDING
2026-07-22 17:38:37 INFO None 5257231: status RUNNING/PENDING
2026-07-22 17:38:37 INFO None 5257232: status RUNNING/PENDING
2026-07-22 17:38:37 INFO None 5257234: status RUNNING/PENDING
2026-07-22 17:38:37 INFO None 5257236: status RUNNING/PENDING
2026-07-22 17:38:37 INFO None 5257237: status RUNNING/PENDING
2026-07-22 17:38:37 INFO None 5257238: status RUNNING/PENDING
2026-07-22 17:38:37 INFO Jobs still running: ['5257213', '5257215', '5257218', '5257220', '5257223', '5257224', '5257225', '5257226', '5257229', '5257231', '5257232', '5257234', '5257236', '5257237', '5257238']. Waiting...
2026-07-22 17:38:52 INFO None 5257213: status FINISHED
2026-07-22 17:38:52 INFO None 5257215: status RUNNING/PENDING
2026-07-22 17:38:52 INFO None 5257218: status RUNNING/PENDING
2026-07-22 17:38:52 INFO None 5257220: status RUNNING/PENDING
2026-07-22 17:38:52 INFO None 5257223: status RUNNING/PENDING
2026-07-22 17:38:52 INFO None 5257224: status RUNNING/PENDING
2026-07-22 17:38:52 INFO None 5257225: status RUNNING/PENDING
2026-07-22 17:38:52 INFO None 5257226: status RUNNING/PENDING
2026-07-22 17:38:52 INFO None 5257229: status RUNNING/PENDING
2026-07-22 17:38:52 INFO None 5257231: status RUNNING/PENDING
2026-07-22 17:38:52 INFO None 5257232: status RUNNING/PENDING
2026-07-22 17:38:52 INFO None 5257234: status RUNNING/PENDING
2026-07-22 17:38:52 INFO None 5257236: status RUNNING/PENDING
2026-07-22 17:38:52 INFO None 5257237: status RUNNING/PENDING
2026-07-22 17:38:52 INFO None 5257238: status RUNNING/PENDING
2026-07-22 17:38:52 INFO Jobs still running: ['5257215', '5257218', '5257220', '5257223', '5257224', '5257225', '5257226', '5257229', '5257231', '5257232', '5257234', '5257236', '5257237', '5257238']. Waiting...
2026-07-22 17:39:07 INFO None 5257213: status FINISHED
2026-07-22 17:39:07 INFO None 5257215: status RUNNING/PENDING
2026-07-22 17:39:07 INFO None 5257218: status RUNNING/PENDING
2026-07-22 17:39:07 INFO None 5257220: status RUNNING/PENDING
2026-07-22 17:39:07 INFO None 5257223: status RUNNING/PENDING
2026-07-22 17:39:07 INFO None 5257224: status RUNNING/PENDING
2026-07-22 17:39:08 INFO None 5257225: status RUNNING/PENDING
2026-07-22 17:39:08 INFO None 5257226: status RUNNING/PENDING
2026-07-22 17:39:08 INFO None 5257229: status RUNNING/PENDING
2026-07-22 17:39:08 INFO None 5257231: status RUNNING/PENDING
2026-07-22 17:39:08 INFO None 5257232: status RUNNING/PENDING
2026-07-22 17:39:08 INFO None 5257234: status RUNNING/PENDING
2026-07-22 17:39:08 INFO None 5257236: status RUNNING/PENDING
2026-07-22 17:39:08 INFO None 5257237: status RUNNING/PENDING
2026-07-22 17:39:08 INFO None 5257238: status RUNNING/PENDING
2026-07-22 17:39:08 INFO Jobs still running: ['5257215', '5257218', '5257220', '5257223', '5257224', '5257225', '5257226', '5257229', '5257231', '5257232', '5257234', '5257236', '5257237', '5257238']. Waiting...
2026-07-22 17:39:23 INFO None 5257213: status FINISHED
2026-07-22 17:39:23 INFO None 5257215: status RUNNING/PENDING
2026-07-22 17:39:23 INFO None 5257218: status RUNNING/PENDING
2026-07-22 17:39:23 INFO None 5257220: status RUNNING/PENDING
2026-07-22 17:39:23 INFO None 5257223: status RUNNING/PENDING
2026-07-22 17:39:23 INFO None 5257224: status RUNNING/PENDING
2026-07-22 17:39:23 INFO None 5257225: status RUNNING/PENDING
2026-07-22 17:39:23 INFO None 5257226: status RUNNING/PENDING
2026-07-22 17:39:23 INFO None 5257229: status RUNNING/PENDING
2026-07-22 17:39:23 INFO None 5257231: status RUNNING/PENDING
2026-07-22 17:39:23 INFO None 5257232: status RUNNING/PENDING
2026-07-22 17:39:23 INFO None 5257234: status RUNNING/PENDING
2026-07-22 17:39:23 INFO None 5257236: status RUNNING/PENDING
2026-07-22 17:39:23 INFO None 5257237: status RUNNING/PENDING
2026-07-22 17:39:23 INFO None 5257238: status RUNNING/PENDING
2026-07-22 17:39:23 INFO Jobs still running: ['5257215', '5257218', '5257220', '5257223', '5257224', '5257225', '5257226', '5257229', '5257231', '5257232', '5257234', '5257236', '5257237', '5257238']. Waiting...
2026-07-22 17:39:38 INFO None 5257213: status FINISHED
2026-07-22 17:39:38 INFO None 5257215: status RUNNING/PENDING
2026-07-22 17:39:38 INFO None 5257218: status RUNNING/PENDING
2026-07-22 17:39:38 INFO None 5257220: status RUNNING/PENDING
2026-07-22 17:39:38 INFO None 5257223: status RUNNING/PENDING
2026-07-22 17:39:38 INFO None 5257224: status RUNNING/PENDING
2026-07-22 17:39:38 INFO None 5257225: status RUNNING/PENDING
2026-07-22 17:39:38 INFO None 5257226: status RUNNING/PENDING
2026-07-22 17:39:38 INFO None 5257229: status RUNNING/PENDING
2026-07-22 17:39:38 INFO None 5257231: status RUNNING/PENDING
2026-07-22 17:39:38 INFO None 5257232: status RUNNING/PENDING
2026-07-22 17:39:38 INFO None 5257234: status RUNNING/PENDING
2026-07-22 17:39:38 INFO None 5257236: status RUNNING/PENDING
2026-07-22 17:39:38 INFO None 5257237: status RUNNING/PENDING
2026-07-22 17:39:38 INFO None 5257238: status RUNNING/PENDING
2026-07-22 17:39:38 INFO Jobs still running: ['5257215', '5257218', '5257220', '5257223', '5257224', '5257225', '5257226', '5257229', '5257231', '5257232', '5257234', '5257236', '5257237', '5257238']. Waiting...
2026-07-22 17:39:53 INFO None 5257213: status FINISHED
2026-07-22 17:39:53 INFO None 5257215: status RUNNING/PENDING
2026-07-22 17:39:53 INFO None 5257218: status RUNNING/PENDING
2026-07-22 17:39:53 INFO None 5257220: status RUNNING/PENDING
2026-07-22 17:39:53 INFO None 5257223: status RUNNING/PENDING
2026-07-22 17:39:53 INFO None 5257224: status RUNNING/PENDING
2026-07-22 17:39:53 INFO None 5257225: status RUNNING/PENDING
2026-07-22 17:39:53 INFO None 5257226: status RUNNING/PENDING
2026-07-22 17:39:53 INFO None 5257229: status RUNNING/PENDING
2026-07-22 17:39:53 INFO None 5257231: status RUNNING/PENDING
2026-07-22 17:39:53 INFO None 5257232: status RUNNING/PENDING
2026-07-22 17:39:56 INFO None 5257234: status RUNNING/PENDING
2026-07-22 17:39:56 INFO None 5257236: status RUNNING/PENDING
2026-07-22 17:39:56 INFO None 5257237: status RUNNING/PENDING
2026-07-22 17:39:56 INFO None 5257238: status RUNNING/PENDING
2026-07-22 17:39:56 INFO Jobs still running: ['5257215', '5257218', '5257220', '5257223', '5257224', '5257225', '5257226', '5257229', '5257231', '5257232', '5257234', '5257236', '5257237', '5257238']. Waiting...
2026-07-22 17:40:11 INFO None 5257213: status FINISHED
2026-07-22 17:40:11 INFO None 5257215: status FINISHED
2026-07-22 17:40:11 INFO None 5257218: status FINISHED
2026-07-22 17:40:11 INFO None 5257220: status RUNNING/PENDING
2026-07-22 17:40:11 INFO None 5257223: status RUNNING/PENDING
2026-07-22 17:40:11 INFO None 5257224: status RUNNING/PENDING
2026-07-22 17:40:11 INFO None 5257225: status RUNNING/PENDING
2026-07-22 17:40:11 INFO None 5257226: status RUNNING/PENDING
2026-07-22 17:40:11 INFO None 5257229: status RUNNING/PENDING
2026-07-22 17:40:11 INFO None 5257231: status RUNNING/PENDING
2026-07-22 17:40:11 INFO None 5257232: status RUNNING/PENDING
2026-07-22 17:40:11 INFO None 5257234: status RUNNING/PENDING
2026-07-22 17:40:11 INFO None 5257236: status RUNNING/PENDING
2026-07-22 17:40:11 INFO None 5257237: status RUNNING/PENDING
2026-07-22 17:40:11 INFO None 5257238: status RUNNING/PENDING
2026-07-22 17:40:11 INFO Jobs still running: ['5257220', '5257223', '5257224', '5257225', '5257226', '5257229', '5257231', '5257232', '5257234', '5257236', '5257237', '5257238']. Waiting...
2026-07-22 17:40:26 INFO None 5257213: status FINISHED
2026-07-22 17:40:26 INFO None 5257215: status FINISHED
2026-07-22 17:40:26 INFO None 5257218: status FINISHED
2026-07-22 17:40:26 INFO None 5257220: status RUNNING/PENDING
2026-07-22 17:40:26 INFO None 5257223: status RUNNING/PENDING
2026-07-22 17:40:26 INFO None 5257224: status RUNNING/PENDING
2026-07-22 17:40:26 INFO None 5257225: status RUNNING/PENDING
2026-07-22 17:40:26 INFO None 5257226: status FINISHED
2026-07-22 17:40:26 INFO None 5257229: status RUNNING/PENDING
2026-07-22 17:40:26 INFO None 5257231: status RUNNING/PENDING
2026-07-22 17:40:26 INFO None 5257232: status RUNNING/PENDING
2026-07-22 17:40:28 INFO None 5257234: status RUNNING/PENDING
2026-07-22 17:40:28 INFO None 5257236: status RUNNING/PENDING
2026-07-22 17:40:28 INFO None 5257237: status RUNNING/PENDING
2026-07-22 17:40:28 INFO None 5257238: status RUNNING/PENDING
2026-07-22 17:40:28 INFO Jobs still running: ['5257220', '5257223', '5257224', '5257225', '5257229', '5257231', '5257232', '5257234', '5257236', '5257237', '5257238']. Waiting...
2026-07-22 17:40:43 INFO None 5257213: status FINISHED
2026-07-22 17:40:43 INFO None 5257215: status FINISHED
2026-07-22 17:40:43 INFO None 5257218: status FINISHED
2026-07-22 17:40:43 INFO None 5257220: status RUNNING/PENDING
2026-07-22 17:40:43 INFO None 5257223: status RUNNING/PENDING
2026-07-22 17:40:43 INFO None 5257224: status RUNNING/PENDING
2026-07-22 17:40:43 INFO None 5257225: status RUNNING/PENDING
2026-07-22 17:40:43 INFO None 5257226: status FINISHED
2026-07-22 17:40:43 INFO None 5257229: status RUNNING/PENDING
2026-07-22 17:40:43 INFO None 5257231: status RUNNING/PENDING
2026-07-22 17:40:43 INFO None 5257232: status RUNNING/PENDING
2026-07-22 17:40:43 INFO None 5257234: status RUNNING/PENDING
2026-07-22 17:40:43 INFO None 5257236: status FINISHED
2026-07-22 17:40:44 INFO None 5257237: status RUNNING/PENDING
2026-07-22 17:40:44 INFO None 5257238: status RUNNING/PENDING
2026-07-22 17:40:44 INFO Jobs still running: ['5257220', '5257223', '5257224', '5257225', '5257229', '5257231', '5257232', '5257234', '5257237', '5257238']. Waiting...
2026-07-22 17:40:59 INFO None 5257213: status FINISHED
2026-07-22 17:40:59 INFO None 5257215: status FINISHED
2026-07-22 17:40:59 INFO None 5257218: status FINISHED
2026-07-22 17:40:59 INFO None 5257220: status RUNNING/PENDING
2026-07-22 17:40:59 INFO None 5257223: status RUNNING/PENDING
2026-07-22 17:40:59 INFO None 5257224: status RUNNING/PENDING
2026-07-22 17:40:59 INFO None 5257225: status RUNNING/PENDING
2026-07-22 17:40:59 INFO None 5257226: status FINISHED
2026-07-22 17:40:59 INFO None 5257229: status RUNNING/PENDING
2026-07-22 17:40:59 INFO None 5257231: status RUNNING/PENDING
2026-07-22 17:40:59 INFO None 5257232: status RUNNING/PENDING
2026-07-22 17:40:59 INFO None 5257234: status RUNNING/PENDING
2026-07-22 17:40:59 INFO None 5257236: status FINISHED
2026-07-22 17:40:59 INFO None 5257237: status RUNNING/PENDING
2026-07-22 17:40:59 INFO None 5257238: status RUNNING/PENDING
2026-07-22 17:40:59 INFO Jobs still running: ['5257220', '5257223', '5257224', '5257225', '5257229', '5257231', '5257232', '5257234', '5257237', '5257238']. Waiting...
2026-07-22 17:41:14 INFO None 5257213: status FINISHED
2026-07-22 17:41:14 INFO None 5257215: status FINISHED
2026-07-22 17:41:14 INFO None 5257218: status FINISHED
2026-07-22 17:41:14 INFO None 5257220: status RUNNING/PENDING
2026-07-22 17:41:14 INFO None 5257223: status RUNNING/PENDING
2026-07-22 17:41:14 INFO None 5257224: status RUNNING/PENDING
2026-07-22 17:41:14 INFO None 5257225: status RUNNING/PENDING
2026-07-22 17:41:14 INFO None 5257226: status FINISHED
2026-07-22 17:41:14 INFO None 5257229: status RUNNING/PENDING
2026-07-22 17:41:14 INFO None 5257231: status RUNNING/PENDING
2026-07-22 17:41:14 INFO None 5257232: status RUNNING/PENDING
2026-07-22 17:41:14 INFO None 5257234: status RUNNING/PENDING
2026-07-22 17:41:14 INFO None 5257236: status FINISHED
2026-07-22 17:41:14 INFO None 5257237: status RUNNING/PENDING
2026-07-22 17:41:14 INFO None 5257238: status RUNNING/PENDING
2026-07-22 17:41:14 INFO Jobs still running: ['5257220', '5257223', '5257224', '5257225', '5257229', '5257231', '5257232', '5257234', '5257237', '5257238']. Waiting...
2026-07-22 17:41:30 INFO None 5257213: status FINISHED
2026-07-22 17:41:30 INFO None 5257215: status FINISHED
2026-07-22 17:41:30 INFO None 5257218: status FINISHED
2026-07-22 17:41:30 INFO None 5257220: status RUNNING/PENDING
2026-07-22 17:41:30 INFO None 5257223: status RUNNING/PENDING
2026-07-22 17:41:30 INFO None 5257224: status RUNNING/PENDING
2026-07-22 17:41:30 INFO None 5257225: status RUNNING/PENDING
2026-07-22 17:41:30 INFO None 5257226: status FINISHED
2026-07-22 17:41:30 INFO None 5257229: status RUNNING/PENDING
2026-07-22 17:41:30 INFO None 5257231: status RUNNING/PENDING
2026-07-22 17:41:30 INFO None 5257232: status RUNNING/PENDING
2026-07-22 17:41:30 INFO None 5257234: status RUNNING/PENDING
2026-07-22 17:41:30 INFO None 5257236: status FINISHED
2026-07-22 17:41:30 INFO None 5257237: status RUNNING/PENDING
2026-07-22 17:41:30 INFO None 5257238: status RUNNING/PENDING
2026-07-22 17:41:30 INFO Jobs still running: ['5257220', '5257223', '5257224', '5257225', '5257229', '5257231', '5257232', '5257234', '5257237', '5257238']. Waiting...
2026-07-22 17:41:45 INFO None 5257213: status FINISHED
2026-07-22 17:41:46 INFO None 5257215: status FINISHED
2026-07-22 17:41:46 INFO None 5257218: status FINISHED
2026-07-22 17:41:46 INFO None 5257220: status FINISHED
2026-07-22 17:41:46 INFO None 5257223: status FINISHED
2026-07-22 17:41:46 INFO None 5257224: status FINISHED
2026-07-22 17:41:46 INFO None 5257225: status RUNNING/PENDING
2026-07-22 17:41:46 INFO None 5257226: status FINISHED
2026-07-22 17:41:46 INFO None 5257229: status FINISHED
2026-07-22 17:41:46 INFO None 5257231: status RUNNING/PENDING
2026-07-22 17:41:46 INFO None 5257232: status FINISHED
2026-07-22 17:41:46 INFO None 5257234: status RUNNING/PENDING
2026-07-22 17:41:46 INFO None 5257236: status FINISHED
2026-07-22 17:41:48 INFO None 5257237: status FINISHED
2026-07-22 17:41:48 INFO None 5257238: status RUNNING/PENDING
2026-07-22 17:41:48 INFO Jobs still running: ['5257225', '5257231', '5257234', '5257238']. Waiting...
2026-07-22 17:42:03 INFO None 5257213: status FINISHED
2026-07-22 17:42:03 INFO None 5257215: status FINISHED
2026-07-22 17:42:03 INFO None 5257218: status FINISHED
2026-07-22 17:42:03 INFO None 5257220: status FINISHED
2026-07-22 17:42:03 INFO None 5257223: status FINISHED
2026-07-22 17:42:03 INFO None 5257224: status FINISHED
2026-07-22 17:42:03 INFO None 5257225: status FINISHED
2026-07-22 17:42:03 INFO None 5257226: status FINISHED
2026-07-22 17:42:03 INFO None 5257229: status FINISHED
2026-07-22 17:42:03 INFO None 5257231: status FINISHED
2026-07-22 17:42:03 INFO None 5257232: status FINISHED
2026-07-22 17:42:03 INFO None 5257234: status FINISHED
2026-07-22 17:42:03 INFO None 5257236: status FINISHED
2026-07-22 17:42:03 INFO None 5257237: status FINISHED
2026-07-22 17:42:03 INFO None 5257238: status FINISHED
2026-07-22 17:42:03 INFO Jobs ['5257213', '5257215', '5257218', '5257220', '5257223', '5257224', '5257225', '5257226', '5257229', '5257231', '5257232', '5257234', '5257236', '5257237', '5257238'] have finished
2026-07-22 17:42:03 INFO Checking restart files were created ...
2026-07-22 17:42:03 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020700_1_ENS1.nc(668832435 bytes)
2026-07-22 17:42:03 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020700_1_ENS2.nc(668832435 bytes)
2026-07-22 17:42:03 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020700_1_ENS3.nc(668832435 bytes)
2026-07-22 17:42:03 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020700_1_ENS4.nc(668832435 bytes)
2026-07-22 17:42:03 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020700_1_ENS5.nc(668832435 bytes)
2026-07-22 17:42:03 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020700_1_ENS6.nc(668832435 bytes)
2026-07-22 17:42:03 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020700_1_ENS7.nc(668832435 bytes)
2026-07-22 17:42:03 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020700_1_ENS8.nc(668832435 bytes)
2026-07-22 17:42:03 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020700_1_ENS9.nc(668832435 bytes)
2026-07-22 17:42:03 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020700_1_ENS10.nc(668832435 bytes)
2026-07-22 17:42:03 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020700_1_ENS11.nc(668832435 bytes)
2026-07-22 17:42:03 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020700_1_ENS12.nc(668832435 bytes)
2026-07-22 17:42:03 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020700_1_ENS13.nc(668832435 bytes)
2026-07-22 17:42:03 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020700_1_ENS14.nc(668832435 bytes)
2026-07-22 17:42:03 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020700_1_ENS15.nc(668832435 bytes)
2026-07-22 17:42:03 INFO  Run_model() completed successfully.
2026-07-22 17:42:03 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 00:00:00 simulated_time=2020-02-07 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 17:42:03 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 01:00:00 days=153073 seconds=3600
2026-07-22 17:42:03 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-22 17:42:03 INFO [TIME] increment current_time 2020-02-07 00:00:00 -> 2020-02-07 01:00:00
2026-07-22 17:42:03 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 01:00:00 simulated_time=2020-02-07 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 17:42:03 INFO ---------->>> Running process_satellite_data()
2026-07-22 17:42:03 INFO [DART] No satellite data found, skipping assimilation
2026-07-22 17:42:03 INFO after_assimilation() skipped
2026-07-22 17:42:03 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-22 17:42:03 INFO [TIME] step_end current_time=2020-02-07 01:00:00 simulated_time=2020-02-07 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 17:42:03 INFO [TIME] step_start current_time=2020-02-07 01:00:00 simulated_time=2020-02-07 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 17:42:03 INFO [TIME] window start=2020-02-07 01:00:00 end=2020-02-07 09:00:00 run_hours=8 has_assimilation=True
2026-07-22 17:42:03 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:42:05 INFO Hourly dataset computed and listing created
2026-07-22 17:42:16 INFO Hourly dataset computed
2026-07-22 17:42:16 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:42:17 INFO Hourly dataset computed and listing created
2026-07-22 17:42:19 INFO Hourly dataset computed
2026-07-22 17:42:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:42:20 INFO Hourly dataset computed and listing created
2026-07-22 17:42:22 INFO Hourly dataset computed
2026-07-22 17:42:22 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:42:23 INFO Hourly dataset computed and listing created
2026-07-22 17:42:25 INFO Hourly dataset computed
2026-07-22 17:42:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:42:26 INFO Hourly dataset computed and listing created
2026-07-22 17:42:28 INFO Hourly dataset computed
2026-07-22 17:42:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:42:29 INFO Hourly dataset computed and listing created
2026-07-22 17:42:31 INFO Hourly dataset computed
2026-07-22 17:42:31 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:42:32 INFO Hourly dataset computed and listing created
2026-07-22 17:42:35 INFO Hourly dataset computed
2026-07-22 17:42:35 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:42:36 INFO Hourly dataset computed and listing created
2026-07-22 17:42:38 INFO Hourly dataset computed
2026-07-22 17:42:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:42:39 INFO Hourly dataset computed and listing created
2026-07-22 17:42:41 INFO Hourly dataset computed
2026-07-22 17:42:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:42:42 INFO Hourly dataset computed and listing created
2026-07-22 17:42:44 INFO Hourly dataset computed
2026-07-22 17:42:44 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:42:45 INFO Hourly dataset computed and listing created
2026-07-22 17:42:47 INFO Hourly dataset computed
2026-07-22 17:42:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:42:49 INFO Hourly dataset computed and listing created
2026-07-22 17:42:51 INFO Hourly dataset computed
2026-07-22 17:42:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:42:52 INFO Hourly dataset computed and listing created
2026-07-22 17:42:54 INFO Hourly dataset computed
2026-07-22 17:42:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:42:55 INFO Hourly dataset computed and listing created
2026-07-22 17:42:58 INFO Hourly dataset computed
2026-07-22 17:42:58 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:42:59 INFO Hourly dataset computed and listing created
2026-07-22 17:43:01 INFO Hourly dataset computed
2026-07-22 17:43:01 INFO ---------->>> Running CHIMERE model from 2020-02-07 01:00:00 to 2020-02-07 09:00:00
2026-07-22 17:43:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:43:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1
2026-07-22 17:43:01 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020700_1_ENS1.nc
2026-07-22 17:43:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-22 17:43:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:43:01 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-22 17:43:01 INFO Queuing job for member 1...
2026-07-22 17:43:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:43:01 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-22 17:43:01 INFO Found: ['5257273']
2026-07-22 17:43:06 INFO [TGCC-IRENE] Submitted job with ID:['5257273']
2026-07-22 17:43:06 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:43:06 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2
2026-07-22 17:43:06 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020700_1_ENS2.nc
2026-07-22 17:43:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-22 17:43:07 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:43:07 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-22 17:43:07 INFO Queuing job for member 2...
2026-07-22 17:43:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:43:07 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-22 17:43:07 INFO Found: ['5257274']
2026-07-22 17:43:12 INFO [TGCC-IRENE] Submitted job with ID:['5257274']
2026-07-22 17:43:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:43:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3
2026-07-22 17:43:12 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020700_1_ENS3.nc
2026-07-22 17:43:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-22 17:43:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:43:12 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-22 17:43:12 INFO Queuing job for member 3...
2026-07-22 17:43:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:43:12 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-22 17:43:15 INFO Found: ['5257275']
2026-07-22 17:43:20 INFO [TGCC-IRENE] Submitted job with ID:['5257275']
2026-07-22 17:43:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:43:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4
2026-07-22 17:43:20 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020700_1_ENS4.nc
2026-07-22 17:43:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-22 17:43:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:43:20 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-22 17:43:20 INFO Queuing job for member 4...
2026-07-22 17:43:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:43:20 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-22 17:43:22 INFO Found: ['5257276']
2026-07-22 17:43:27 INFO [TGCC-IRENE] Submitted job with ID:['5257276']
2026-07-22 17:43:27 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:43:27 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5
2026-07-22 17:43:27 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020700_1_ENS5.nc
2026-07-22 17:43:27 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-22 17:43:27 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:43:27 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-22 17:43:27 INFO Queuing job for member 5...
2026-07-22 17:43:27 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:43:27 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-22 17:43:30 INFO Found: ['5257277']
2026-07-22 17:43:35 INFO [TGCC-IRENE] Submitted job with ID:['5257277']
2026-07-22 17:43:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:43:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6
2026-07-22 17:43:35 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020700_1_ENS6.nc
2026-07-22 17:43:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-22 17:43:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:43:35 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-22 17:43:35 INFO Queuing job for member 6...
2026-07-22 17:43:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:43:35 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-22 17:43:37 INFO Found: ['5257278']
2026-07-22 17:43:42 INFO [TGCC-IRENE] Submitted job with ID:['5257278']
2026-07-22 17:43:42 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:43:42 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7
2026-07-22 17:43:42 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020700_1_ENS7.nc
2026-07-22 17:43:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-22 17:43:42 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:43:42 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-22 17:43:42 INFO Queuing job for member 7...
2026-07-22 17:43:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:43:42 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-22 17:43:45 INFO Found: ['5257279']
2026-07-22 17:43:50 INFO [TGCC-IRENE] Submitted job with ID:['5257279']
2026-07-22 17:43:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:43:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8
2026-07-22 17:43:50 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020700_1_ENS8.nc
2026-07-22 17:43:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-22 17:43:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:43:50 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-22 17:43:50 INFO Queuing job for member 8...
2026-07-22 17:43:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:43:50 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-22 17:43:52 INFO Found: ['5257281']
2026-07-22 17:43:57 INFO [TGCC-IRENE] Submitted job with ID:['5257281']
2026-07-22 17:43:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:43:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9
2026-07-22 17:43:57 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020700_1_ENS9.nc
2026-07-22 17:43:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-22 17:43:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:43:57 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-22 17:43:57 INFO Queuing job for member 9...
2026-07-22 17:43:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:43:57 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-22 17:44:00 INFO Found: ['5257282']
2026-07-22 17:44:05 INFO [TGCC-IRENE] Submitted job with ID:['5257282']
2026-07-22 17:44:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:44:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10
2026-07-22 17:44:05 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020700_1_ENS10.nc
2026-07-22 17:44:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-22 17:44:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:44:05 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-22 17:44:05 INFO Queuing job for member 10...
2026-07-22 17:44:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:44:05 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-22 17:44:05 INFO Found: ['5257284']
2026-07-22 17:44:10 INFO [TGCC-IRENE] Submitted job with ID:['5257284']
2026-07-22 17:44:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:44:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11
2026-07-22 17:44:10 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020700_1_ENS11.nc
2026-07-22 17:44:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-22 17:44:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:44:10 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-22 17:44:10 INFO Queuing job for member 11...
2026-07-22 17:44:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:44:10 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-22 17:44:11 INFO Found: ['5257285']
2026-07-22 17:44:16 INFO [TGCC-IRENE] Submitted job with ID:['5257285']
2026-07-22 17:44:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:44:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12
2026-07-22 17:44:16 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020700_1_ENS12.nc
2026-07-22 17:44:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-22 17:44:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:44:16 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-22 17:44:16 INFO Queuing job for member 12...
2026-07-22 17:44:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:44:16 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-22 17:44:17 INFO Found: ['5257288']
2026-07-22 17:44:22 INFO [TGCC-IRENE] Submitted job with ID:['5257288']
2026-07-22 17:44:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:44:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13
2026-07-22 17:44:22 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020700_1_ENS13.nc
2026-07-22 17:44:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-22 17:44:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:44:22 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-22 17:44:22 INFO Queuing job for member 13...
2026-07-22 17:44:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:44:22 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-22 17:44:23 INFO Found: ['5257289']
2026-07-22 17:44:28 INFO [TGCC-IRENE] Submitted job with ID:['5257289']
2026-07-22 17:44:28 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:44:28 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14
2026-07-22 17:44:28 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020700_1_ENS14.nc
2026-07-22 17:44:28 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-22 17:44:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:44:28 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-22 17:44:28 INFO Queuing job for member 14...
2026-07-22 17:44:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:44:28 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-22 17:44:28 INFO Found: ['5257290']
2026-07-22 17:44:33 INFO [TGCC-IRENE] Submitted job with ID:['5257290']
2026-07-22 17:44:33 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:44:33 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15
2026-07-22 17:44:33 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020700_1_ENS15.nc
2026-07-22 17:44:33 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-22 17:44:33 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:44:33 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-22 17:44:33 INFO Queuing job for member 15...
2026-07-22 17:44:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:44:33 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-22 17:44:34 INFO Found: ['5257291']
2026-07-22 17:44:39 INFO [TGCC-IRENE] Submitted job with ID:['5257291']
2026-07-22 17:44:39 INFO Checking job status ...
2026-07-22 17:44:39 INFO None 5257273: status RUNNING/PENDING
2026-07-22 17:44:39 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:44:39 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:44:39 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:44:39 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:44:39 INFO None 5257278: status RUNNING/PENDING
2026-07-22 17:44:39 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:44:39 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:44:39 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:44:39 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:44:39 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:44:39 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:44:39 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:44:39 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:44:39 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:44:39 INFO Jobs still running: ['5257273', '5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:44:56 INFO None 5257273: status RUNNING/PENDING
2026-07-22 17:44:56 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:44:56 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:44:56 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:44:56 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:44:56 INFO None 5257278: status RUNNING/PENDING
2026-07-22 17:44:56 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:44:56 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:44:56 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:44:56 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:44:56 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:44:56 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:44:56 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:44:56 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:44:56 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:44:56 INFO Jobs still running: ['5257273', '5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:45:11 INFO None 5257273: status RUNNING/PENDING
2026-07-22 17:45:11 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:45:12 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:45:12 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:45:12 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:45:12 INFO None 5257278: status RUNNING/PENDING
2026-07-22 17:45:12 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:45:12 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:45:12 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:45:12 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:45:12 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:45:12 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:45:14 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:45:14 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:45:14 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:45:14 INFO Jobs still running: ['5257273', '5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:45:29 INFO None 5257273: status RUNNING/PENDING
2026-07-22 17:45:29 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:45:29 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:45:29 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:45:29 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:45:29 INFO None 5257278: status RUNNING/PENDING
2026-07-22 17:45:29 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:45:29 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:45:29 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:45:29 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:45:29 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:45:29 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:45:29 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:45:29 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:45:29 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:45:29 INFO Jobs still running: ['5257273', '5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:45:44 INFO None 5257273: status RUNNING/PENDING
2026-07-22 17:45:44 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:45:44 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:45:44 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:45:44 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:45:44 INFO None 5257278: status RUNNING/PENDING
2026-07-22 17:45:44 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:45:44 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:45:44 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:45:44 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:45:44 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:45:44 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:45:44 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:45:44 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:45:44 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:45:44 INFO Jobs still running: ['5257273', '5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:45:59 INFO None 5257273: status RUNNING/PENDING
2026-07-22 17:45:59 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:46:00 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:46:00 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:46:00 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:46:00 INFO None 5257278: status RUNNING/PENDING
2026-07-22 17:46:00 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:46:00 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:46:00 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:46:00 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:46:00 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:46:00 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:46:00 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:46:00 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:46:00 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:46:00 INFO Jobs still running: ['5257273', '5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:46:15 INFO None 5257273: status RUNNING/PENDING
2026-07-22 17:46:15 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:46:15 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:46:15 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:46:15 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:46:15 INFO None 5257278: status RUNNING/PENDING
2026-07-22 17:46:15 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:46:15 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:46:15 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:46:15 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:46:15 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:46:15 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:46:15 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:46:15 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:46:15 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:46:15 INFO Jobs still running: ['5257273', '5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:46:30 INFO None 5257273: status RUNNING/PENDING
2026-07-22 17:46:30 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:46:30 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:46:30 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:46:30 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:46:30 INFO None 5257278: status RUNNING/PENDING
2026-07-22 17:46:30 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:46:30 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:46:30 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:46:31 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:46:31 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:46:31 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:46:31 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:46:31 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:46:31 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:46:31 INFO Jobs still running: ['5257273', '5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:46:46 INFO None 5257273: status RUNNING/PENDING
2026-07-22 17:46:46 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:46:46 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:46:46 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:46:46 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:46:46 INFO None 5257278: status RUNNING/PENDING
2026-07-22 17:46:46 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:46:46 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:46:46 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:46:46 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:46:46 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:46:48 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:46:48 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:46:48 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:46:48 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:46:48 INFO Jobs still running: ['5257273', '5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:47:03 INFO None 5257273: status RUNNING/PENDING
2026-07-22 17:47:03 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:47:03 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:47:03 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:47:03 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:47:03 INFO None 5257278: status RUNNING/PENDING
2026-07-22 17:47:03 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:47:03 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:47:03 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:47:03 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:47:03 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:47:03 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:47:03 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:47:03 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:47:03 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:47:03 INFO Jobs still running: ['5257273', '5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:47:18 INFO None 5257273: status RUNNING/PENDING
2026-07-22 17:47:18 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:47:18 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:47:18 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:47:18 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:47:18 INFO None 5257278: status RUNNING/PENDING
2026-07-22 17:47:18 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:47:18 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:47:18 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:47:21 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:47:21 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:47:21 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:47:21 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:47:21 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:47:21 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:47:21 INFO Jobs still running: ['5257273', '5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:47:36 INFO None 5257273: status RUNNING/PENDING
2026-07-22 17:47:36 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:47:36 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:47:36 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:47:36 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:47:36 INFO None 5257278: status RUNNING/PENDING
2026-07-22 17:47:36 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:47:36 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:47:36 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:47:36 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:47:36 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:47:36 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:47:36 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:47:36 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:47:36 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:47:36 INFO Jobs still running: ['5257273', '5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:47:51 INFO None 5257273: status RUNNING/PENDING
2026-07-22 17:47:51 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:47:51 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:47:51 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:47:51 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:47:51 INFO None 5257278: status RUNNING/PENDING
2026-07-22 17:47:51 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:47:51 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:47:51 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:47:51 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:47:51 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:47:51 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:47:51 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:47:51 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:47:51 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:47:51 INFO Jobs still running: ['5257273', '5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:48:07 INFO None 5257273: status RUNNING/PENDING
2026-07-22 17:48:08 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:48:08 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:48:08 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:48:08 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:48:08 INFO None 5257278: status RUNNING/PENDING
2026-07-22 17:48:08 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:48:08 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:48:08 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:48:08 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:48:08 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:48:08 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:48:08 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:48:08 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:48:08 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:48:08 INFO Jobs still running: ['5257273', '5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:48:23 INFO None 5257273: status RUNNING/PENDING
2026-07-22 17:48:23 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:48:23 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:48:23 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:48:23 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:48:23 INFO None 5257278: status RUNNING/PENDING
2026-07-22 17:48:23 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:48:23 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:48:23 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:48:23 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:48:23 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:48:23 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:48:23 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:48:25 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:48:25 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:48:25 INFO Jobs still running: ['5257273', '5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:48:40 INFO None 5257273: status RUNNING/PENDING
2026-07-22 17:48:40 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:48:40 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:48:40 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:48:40 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:48:40 INFO None 5257278: status RUNNING/PENDING
2026-07-22 17:48:40 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:48:40 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:48:40 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:48:40 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:48:40 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:48:40 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:48:40 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:48:40 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:48:40 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:48:40 INFO Jobs still running: ['5257273', '5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:48:55 INFO None 5257273: status RUNNING/PENDING
2026-07-22 17:48:55 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:48:55 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:48:55 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:48:55 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:48:56 INFO None 5257278: status RUNNING/PENDING
2026-07-22 17:48:56 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:48:56 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:48:56 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:48:56 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:48:56 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:48:56 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:48:56 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:48:56 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:48:56 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:48:56 INFO Jobs still running: ['5257273', '5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:49:12 INFO None 5257273: status RUNNING/PENDING
2026-07-22 17:49:12 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:49:12 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:49:12 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:49:12 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:49:12 INFO None 5257278: status RUNNING/PENDING
2026-07-22 17:49:12 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:49:12 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:49:12 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:49:12 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:49:12 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:49:12 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:49:12 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:49:12 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:49:12 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:49:12 INFO Jobs still running: ['5257273', '5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:49:27 INFO None 5257273: status RUNNING/PENDING
2026-07-22 17:49:27 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:49:27 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:49:27 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:49:27 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:49:27 INFO None 5257278: status RUNNING/PENDING
2026-07-22 17:49:27 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:49:27 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:49:27 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:49:27 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:49:27 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:49:27 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:49:27 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:49:27 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:49:27 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:49:27 INFO Jobs still running: ['5257273', '5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:49:43 INFO None 5257273: status RUNNING/PENDING
2026-07-22 17:49:44 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:49:44 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:49:44 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:49:44 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:49:44 INFO None 5257278: status RUNNING/PENDING
2026-07-22 17:49:44 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:49:44 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:49:44 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:49:44 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:49:44 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:49:44 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:49:44 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:49:44 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:49:44 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:49:44 INFO Jobs still running: ['5257273', '5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:49:59 INFO None 5257273: status RUNNING/PENDING
2026-07-22 17:49:59 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:49:59 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:49:59 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:49:59 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:49:59 INFO None 5257278: status RUNNING/PENDING
2026-07-22 17:49:59 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:49:59 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:49:59 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:49:59 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:49:59 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:49:59 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:49:59 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:49:59 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:50:01 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:50:01 INFO Jobs still running: ['5257273', '5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:50:16 INFO None 5257273: status RUNNING/PENDING
2026-07-22 17:50:16 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:50:16 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:50:16 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:50:16 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:50:16 INFO None 5257278: status RUNNING/PENDING
2026-07-22 17:50:16 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:50:16 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:50:16 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:50:16 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:50:16 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:50:16 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:50:16 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:50:16 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:50:16 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:50:16 INFO Jobs still running: ['5257273', '5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:50:31 INFO None 5257273: status RUNNING/PENDING
2026-07-22 17:50:31 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:50:31 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:50:31 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:50:32 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:50:32 INFO None 5257278: status RUNNING/PENDING
2026-07-22 17:50:32 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:50:32 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:50:32 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:50:32 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:50:32 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:50:32 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:50:32 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:50:34 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:50:34 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:50:34 INFO Jobs still running: ['5257273', '5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:50:49 INFO None 5257273: status RUNNING/PENDING
2026-07-22 17:50:49 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:50:49 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:50:49 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:50:49 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:50:49 INFO None 5257278: status RUNNING/PENDING
2026-07-22 17:50:49 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:50:49 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:50:49 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:50:49 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:50:49 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:50:49 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:50:49 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:50:49 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:50:49 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:50:49 INFO Jobs still running: ['5257273', '5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:51:04 INFO None 5257273: status RUNNING/PENDING
2026-07-22 17:51:04 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:51:04 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:51:04 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:51:04 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:51:04 INFO None 5257278: status RUNNING/PENDING
2026-07-22 17:51:04 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:51:04 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:51:04 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:51:04 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:51:04 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:51:04 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:51:04 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:51:04 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:51:04 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:51:04 INFO Jobs still running: ['5257273', '5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:51:19 INFO None 5257273: status FINISHED
2026-07-22 17:51:19 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:51:19 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:51:19 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:51:19 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:51:19 INFO None 5257278: status RUNNING/PENDING
2026-07-22 17:51:19 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:51:19 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:51:19 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:51:19 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:51:20 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:51:20 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:51:20 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:51:20 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:51:20 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:51:20 INFO Jobs still running: ['5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:51:35 INFO None 5257273: status FINISHED
2026-07-22 17:51:35 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:51:35 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:51:35 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:51:35 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:51:35 INFO None 5257278: status RUNNING/PENDING
2026-07-22 17:51:35 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:51:35 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:51:35 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:51:35 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:51:35 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:51:35 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:51:35 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:51:35 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:51:35 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:51:35 INFO Jobs still running: ['5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:51:50 INFO None 5257273: status FINISHED
2026-07-22 17:51:50 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:51:50 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:51:52 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:51:52 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:51:52 INFO None 5257278: status RUNNING/PENDING
2026-07-22 17:51:52 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:51:52 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:51:52 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:51:52 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:51:52 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:51:52 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:51:52 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:51:52 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:51:52 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:51:52 INFO Jobs still running: ['5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:52:07 INFO None 5257273: status FINISHED
2026-07-22 17:52:07 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:52:07 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:52:07 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:52:07 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:52:07 INFO None 5257278: status RUNNING/PENDING
2026-07-22 17:52:07 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:52:07 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:52:07 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:52:07 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:52:07 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:52:07 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:52:07 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:52:07 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:52:07 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:52:07 INFO Jobs still running: ['5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:52:23 INFO None 5257273: status FINISHED
2026-07-22 17:52:23 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:52:23 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:52:23 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:52:23 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:52:23 INFO None 5257278: status RUNNING/PENDING
2026-07-22 17:52:23 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:52:23 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:52:23 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:52:23 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:52:23 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:52:23 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:52:23 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:52:23 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:52:23 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:52:23 INFO Jobs still running: ['5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:52:38 INFO None 5257273: status FINISHED
2026-07-22 17:52:38 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:52:38 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:52:38 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:52:38 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:52:38 INFO None 5257278: status FINISHED
2026-07-22 17:52:38 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:52:38 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:52:38 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:52:38 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:52:38 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:52:38 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:52:38 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:52:38 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:52:38 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:52:38 INFO Jobs still running: ['5257274', '5257275', '5257276', '5257277', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:52:53 INFO None 5257273: status FINISHED
2026-07-22 17:52:53 INFO None 5257274: status RUNNING/PENDING
2026-07-22 17:52:53 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:52:53 INFO None 5257276: status RUNNING/PENDING
2026-07-22 17:52:53 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:52:53 INFO None 5257278: status FINISHED
2026-07-22 17:52:53 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:52:53 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:52:53 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:52:53 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:52:53 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:52:53 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:52:53 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:52:53 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:52:53 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:52:53 INFO Jobs still running: ['5257274', '5257275', '5257276', '5257277', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:53:08 INFO None 5257273: status FINISHED
2026-07-22 17:53:08 INFO None 5257274: status FINISHED
2026-07-22 17:53:08 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:53:08 INFO None 5257276: status FINISHED
2026-07-22 17:53:08 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:53:09 INFO None 5257278: status FINISHED
2026-07-22 17:53:09 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:53:09 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:53:09 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:53:09 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:53:09 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:53:09 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:53:09 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:53:09 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:53:09 INFO None 5257291: status RUNNING/PENDING
2026-07-22 17:53:09 INFO Jobs still running: ['5257275', '5257277', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291']. Waiting...
2026-07-22 17:53:24 INFO None 5257273: status FINISHED
2026-07-22 17:53:24 INFO None 5257274: status FINISHED
2026-07-22 17:53:24 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:53:24 INFO None 5257276: status FINISHED
2026-07-22 17:53:24 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:53:24 INFO None 5257278: status FINISHED
2026-07-22 17:53:24 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:53:24 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:53:24 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:53:24 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:53:24 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:53:24 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:53:24 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:53:26 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:53:26 INFO None 5257291: status FINISHED
2026-07-22 17:53:26 INFO Jobs still running: ['5257275', '5257277', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290']. Waiting...
2026-07-22 17:53:41 INFO None 5257273: status FINISHED
2026-07-22 17:53:41 INFO None 5257274: status FINISHED
2026-07-22 17:53:41 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:53:41 INFO None 5257276: status FINISHED
2026-07-22 17:53:41 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:53:41 INFO None 5257278: status FINISHED
2026-07-22 17:53:41 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:53:41 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:53:41 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:53:41 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:53:41 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:53:41 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:53:41 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:53:41 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:53:41 INFO None 5257291: status FINISHED
2026-07-22 17:53:41 INFO Jobs still running: ['5257275', '5257277', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290']. Waiting...
2026-07-22 17:53:56 INFO None 5257273: status FINISHED
2026-07-22 17:53:56 INFO None 5257274: status FINISHED
2026-07-22 17:53:56 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:53:56 INFO None 5257276: status FINISHED
2026-07-22 17:53:56 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:53:56 INFO None 5257278: status FINISHED
2026-07-22 17:53:56 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:53:56 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:53:57 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:53:57 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:53:57 INFO None 5257285: status RUNNING/PENDING
2026-07-22 17:53:57 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:53:57 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:53:57 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:53:57 INFO None 5257291: status FINISHED
2026-07-22 17:53:57 INFO Jobs still running: ['5257275', '5257277', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290']. Waiting...
2026-07-22 17:54:12 INFO None 5257273: status FINISHED
2026-07-22 17:54:12 INFO None 5257274: status FINISHED
2026-07-22 17:54:12 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:54:12 INFO None 5257276: status FINISHED
2026-07-22 17:54:12 INFO None 5257277: status RUNNING/PENDING
2026-07-22 17:54:12 INFO None 5257278: status FINISHED
2026-07-22 17:54:12 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:54:12 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:54:12 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:54:12 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:54:12 INFO None 5257285: status FINISHED
2026-07-22 17:54:12 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:54:12 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:54:12 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:54:12 INFO None 5257291: status FINISHED
2026-07-22 17:54:12 INFO Jobs still running: ['5257275', '5257277', '5257279', '5257281', '5257282', '5257284', '5257288', '5257289', '5257290']. Waiting...
2026-07-22 17:54:27 INFO None 5257273: status FINISHED
2026-07-22 17:54:27 INFO None 5257274: status FINISHED
2026-07-22 17:54:27 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:54:27 INFO None 5257276: status FINISHED
2026-07-22 17:54:27 INFO None 5257277: status FINISHED
2026-07-22 17:54:27 INFO None 5257278: status FINISHED
2026-07-22 17:54:27 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:54:27 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:54:27 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:54:27 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:54:27 INFO None 5257285: status FINISHED
2026-07-22 17:54:27 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:54:27 INFO None 5257289: status RUNNING/PENDING
2026-07-22 17:54:27 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:54:27 INFO None 5257291: status FINISHED
2026-07-22 17:54:27 INFO Jobs still running: ['5257275', '5257279', '5257281', '5257282', '5257284', '5257288', '5257289', '5257290']. Waiting...
2026-07-22 17:54:42 INFO None 5257273: status FINISHED
2026-07-22 17:54:42 INFO None 5257274: status FINISHED
2026-07-22 17:54:42 INFO None 5257275: status RUNNING/PENDING
2026-07-22 17:54:42 INFO None 5257276: status FINISHED
2026-07-22 17:54:42 INFO None 5257277: status FINISHED
2026-07-22 17:54:42 INFO None 5257278: status FINISHED
2026-07-22 17:54:42 INFO None 5257279: status RUNNING/PENDING
2026-07-22 17:54:42 INFO None 5257281: status RUNNING/PENDING
2026-07-22 17:54:42 INFO None 5257282: status RUNNING/PENDING
2026-07-22 17:54:42 INFO None 5257284: status RUNNING/PENDING
2026-07-22 17:54:42 INFO None 5257285: status FINISHED
2026-07-22 17:54:42 INFO None 5257288: status RUNNING/PENDING
2026-07-22 17:54:42 INFO None 5257289: status FINISHED
2026-07-22 17:54:42 INFO None 5257290: status RUNNING/PENDING
2026-07-22 17:54:44 INFO None 5257291: status FINISHED
2026-07-22 17:54:44 INFO Jobs still running: ['5257275', '5257279', '5257281', '5257282', '5257284', '5257288', '5257290']. Waiting...
2026-07-22 17:55:00 INFO None 5257273: status FINISHED
2026-07-22 17:55:00 INFO None 5257274: status FINISHED
2026-07-22 17:55:00 INFO None 5257275: status FINISHED
2026-07-22 17:55:00 INFO None 5257276: status FINISHED
2026-07-22 17:55:00 INFO None 5257277: status FINISHED
2026-07-22 17:55:00 INFO None 5257278: status FINISHED
2026-07-22 17:55:00 INFO None 5257279: status FINISHED
2026-07-22 17:55:00 INFO None 5257281: status FINISHED
2026-07-22 17:55:00 INFO None 5257282: status FINISHED
2026-07-22 17:55:00 INFO None 5257284: status FINISHED
2026-07-22 17:55:00 INFO None 5257285: status FINISHED
2026-07-22 17:55:00 INFO None 5257288: status FINISHED
2026-07-22 17:55:00 INFO None 5257289: status FINISHED
2026-07-22 17:55:00 INFO None 5257290: status FINISHED
2026-07-22 17:55:00 INFO None 5257291: status FINISHED
2026-07-22 17:55:00 INFO Jobs ['5257273', '5257274', '5257275', '5257276', '5257277', '5257278', '5257279', '5257281', '5257282', '5257284', '5257285', '5257288', '5257289', '5257290', '5257291'] have finished
2026-07-22 17:55:00 INFO Checking restart files were created ...
2026-07-22 17:55:00 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020701_8_ENS1.nc(3005806795 bytes)
2026-07-22 17:55:00 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020701_8_ENS2.nc(3005806795 bytes)
2026-07-22 17:55:00 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020701_8_ENS3.nc(3005806795 bytes)
2026-07-22 17:55:00 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020701_8_ENS4.nc(3005806795 bytes)
2026-07-22 17:55:00 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020701_8_ENS5.nc(3005806795 bytes)
2026-07-22 17:55:00 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020701_8_ENS6.nc(3005806795 bytes)
2026-07-22 17:55:00 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020701_8_ENS7.nc(3005806795 bytes)
2026-07-22 17:55:00 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020701_8_ENS8.nc(3005806795 bytes)
2026-07-22 17:55:00 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020701_8_ENS9.nc(3005806795 bytes)
2026-07-22 17:55:00 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020701_8_ENS10.nc(3005806795 bytes)
2026-07-22 17:55:00 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020701_8_ENS11.nc(3005806795 bytes)
2026-07-22 17:55:00 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020701_8_ENS12.nc(3005806795 bytes)
2026-07-22 17:55:00 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020701_8_ENS13.nc(3005806795 bytes)
2026-07-22 17:55:00 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020701_8_ENS14.nc(3005806795 bytes)
2026-07-22 17:55:00 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020701_8_ENS15.nc(3005806795 bytes)
2026-07-22 17:55:00 INFO  Run_model() completed successfully.
2026-07-22 17:55:00 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 01:00:00 simulated_time=2020-02-07 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 17:55:00 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 09:00:00 days=153073 seconds=32400
2026-07-22 17:55:00 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-22 17:55:00 INFO [TIME] increment current_time 2020-02-07 01:00:00 -> 2020-02-07 09:00:00
2026-07-22 17:55:00 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 09:00:00 simulated_time=2020-02-07 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 17:55:00 INFO ---------->>> Running process_satellite_data()
2026-07-22 17:55:00 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12016.nc
2026-07-22 17:55:00 INFO ---------->>> Running run_obs_converter()
2026-07-22 17:55:00 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_32090_153073.out
2026-07-22 17:55:00 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_32090_153073.out
2026-07-22 17:55:00 INFO ---------->>> Running DART
2026-07-22 17:55:00 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-22 17:55:00 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/chim_ENS1_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/chim_ENS1_2020020709_1_out_toDART.nc
2026-07-22 17:55:00 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/chim_ENS2_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/chim_ENS2_2020020709_1_out_toDART.nc
2026-07-22 17:55:01 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/chim_ENS3_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/chim_ENS3_2020020709_1_out_toDART.nc
2026-07-22 17:55:01 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/chim_ENS4_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/chim_ENS4_2020020709_1_out_toDART.nc
2026-07-22 17:55:02 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/chim_ENS5_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/chim_ENS5_2020020709_1_out_toDART.nc
2026-07-22 17:55:02 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/chim_ENS6_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/chim_ENS6_2020020709_1_out_toDART.nc
2026-07-22 17:55:02 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/chim_ENS7_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/chim_ENS7_2020020709_1_out_toDART.nc
2026-07-22 17:55:03 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/chim_ENS8_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/chim_ENS8_2020020709_1_out_toDART.nc
2026-07-22 17:55:03 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/chim_ENS9_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/chim_ENS9_2020020709_1_out_toDART.nc
2026-07-22 17:55:04 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/chim_ENS10_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/chim_ENS10_2020020709_1_out_toDART.nc
2026-07-22 17:55:04 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/chim_ENS11_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/chim_ENS11_2020020709_1_out_toDART.nc
2026-07-22 17:55:04 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/chim_ENS12_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/chim_ENS12_2020020709_1_out_toDART.nc
2026-07-22 17:55:05 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/chim_ENS13_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/chim_ENS13_2020020709_1_out_toDART.nc
2026-07-22 17:55:05 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/chim_ENS14_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/chim_ENS14_2020020709_1_out_toDART.nc
2026-07-22 17:55:06 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/chim_ENS15_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/chim_ENS15_2020020709_1_out_toDART.nc
2026-07-22 17:55:06 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-22 17:55:06 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-22 17:55:06 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-22 17:55:06 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-22 17:55:06 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-22 17:55:06 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-22 17:55:13 INFO Found: []
2026-07-22 17:55:13 INFO No job id returned by command ./run_filter.bsh
2026-07-22 17:55:13 INFO No monitoring will be performed
2026-07-22 17:55:13 INFO Moving DART output files to analysis and preassim directories for date 2020020709 if present ...
2026-07-22 17:55:13 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020709'
2026-07-22 17:55:13 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020709'
2026-07-22 17:55:13 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020709'
2026-07-22 17:55:13 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020709'
2026-07-22 17:55:13 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020709'
2026-07-22 17:55:13 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020709'
2026-07-22 17:55:13 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020709'
2026-07-22 17:55:13 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020709'
2026-07-22 17:55:13 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020709'
2026-07-22 17:55:13 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020709'
2026-07-22 17:55:13 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020709'
2026-07-22 17:55:13 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020709'
2026-07-22 17:55:13 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020709'
2026-07-22 17:55:13 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020709'
2026-07-22 17:55:13 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020709'
2026-07-22 17:55:13 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020709'
2026-07-22 17:55:13 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020709'
2026-07-22 17:55:13 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020709'
2026-07-22 17:55:13 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020709'
2026-07-22 17:55:13 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020709'
2026-07-22 17:55:13 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020709'
2026-07-22 17:55:13 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020709'
2026-07-22 17:55:13 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020709'
2026-07-22 17:55:13 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020709'
2026-07-22 17:55:13 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020709'
2026-07-22 17:55:13 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020709'
2026-07-22 17:55:13 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020709'
2026-07-22 17:55:13 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020709'
2026-07-22 17:55:13 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020709'
2026-07-22 17:55:13 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020709'
2026-07-22 17:55:13 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020709'
2026-07-22 17:55:13 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020709'
2026-07-22 17:55:13 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020709'
2026-07-22 17:55:13 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020709'
2026-07-22 17:55:13 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-22 17:55:13 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-22 17:55:13 INFO run_dart() is DONE.
2026-07-22 17:55:13 INFO ---------->>> Running update_pollutant_in_end()
2026-07-22 17:55:14 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens01.nc
2026-07-22 17:55:14 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:55:14 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS1.nc
2026-07-22 17:55:15 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 17:55:16 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS1.nc
2026-07-22 17:55:16 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:55:16 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:55:16 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS1.nc
2026-07-22 17:55:16 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 17:55:17 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS1.nc
2026-07-22 17:55:17 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:55:17 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:55:17 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens02.nc
2026-07-22 17:55:18 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:55:18 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS2.nc
2026-07-22 17:55:18 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 17:55:19 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS2.nc
2026-07-22 17:55:19 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:55:19 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:55:19 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS2.nc
2026-07-22 17:55:20 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 17:55:21 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS2.nc
2026-07-22 17:55:21 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:55:21 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:55:21 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens03.nc
2026-07-22 17:55:22 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:55:22 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS3.nc
2026-07-22 17:55:22 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 17:55:23 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS3.nc
2026-07-22 17:55:23 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:55:23 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:55:23 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS3.nc
2026-07-22 17:55:24 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 17:55:24 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS3.nc
2026-07-22 17:55:24 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:55:24 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:55:25 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens04.nc
2026-07-22 17:55:26 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:55:26 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS4.nc
2026-07-22 17:55:26 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 17:55:27 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS4.nc
2026-07-22 17:55:27 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:55:27 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:55:27 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS4.nc
2026-07-22 17:55:28 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 17:55:28 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS4.nc
2026-07-22 17:55:28 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:55:28 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:55:29 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens05.nc
2026-07-22 17:55:29 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:55:29 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS5.nc
2026-07-22 17:55:30 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 17:55:31 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS5.nc
2026-07-22 17:55:31 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:55:31 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:55:31 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS5.nc
2026-07-22 17:55:32 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 17:55:32 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS5.nc
2026-07-22 17:55:32 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:55:32 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:55:33 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens06.nc
2026-07-22 17:55:33 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:55:33 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS6.nc
2026-07-22 17:55:34 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 17:55:34 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS6.nc
2026-07-22 17:55:34 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:55:35 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:55:35 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS6.nc
2026-07-22 17:55:35 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 17:55:36 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS6.nc
2026-07-22 17:55:36 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:55:36 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:55:37 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens07.nc
2026-07-22 17:55:37 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:55:37 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS7.nc
2026-07-22 17:55:38 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 17:55:38 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS7.nc
2026-07-22 17:55:38 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:55:39 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:55:39 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS7.nc
2026-07-22 17:55:39 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 17:55:40 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS7.nc
2026-07-22 17:55:40 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:55:40 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:55:40 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens08.nc
2026-07-22 17:55:41 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:55:41 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS8.nc
2026-07-22 17:55:41 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 17:55:42 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS8.nc
2026-07-22 17:55:42 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:55:42 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:55:42 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS8.nc
2026-07-22 17:55:43 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 17:55:43 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS8.nc
2026-07-22 17:55:43 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:55:43 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:55:44 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens09.nc
2026-07-22 17:55:44 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:55:44 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS9.nc
2026-07-22 17:55:45 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 17:55:46 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS9.nc
2026-07-22 17:55:46 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:55:46 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:55:46 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS9.nc
2026-07-22 17:55:46 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 17:55:47 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS9.nc
2026-07-22 17:55:47 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:55:47 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:55:48 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens10.nc
2026-07-22 17:55:48 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:55:48 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS10.nc
2026-07-22 17:55:49 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 17:55:50 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS10.nc
2026-07-22 17:55:50 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:55:50 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:55:50 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS10.nc
2026-07-22 17:55:51 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 17:55:52 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS10.nc
2026-07-22 17:55:52 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:55:52 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:55:52 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens11.nc
2026-07-22 17:55:53 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:55:53 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS11.nc
2026-07-22 17:55:53 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 17:55:54 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS11.nc
2026-07-22 17:55:54 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:55:54 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:55:54 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS11.nc
2026-07-22 17:55:55 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 17:55:56 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS11.nc
2026-07-22 17:55:56 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:55:56 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:55:56 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens12.nc
2026-07-22 17:55:57 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:55:57 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS12.nc
2026-07-22 17:55:57 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 17:55:58 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS12.nc
2026-07-22 17:55:58 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:55:58 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:55:58 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS12.nc
2026-07-22 17:55:59 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 17:56:00 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS12.nc
2026-07-22 17:56:00 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:56:00 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:56:01 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens13.nc
2026-07-22 17:56:01 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:56:01 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS13.nc
2026-07-22 17:56:02 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 17:56:02 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS13.nc
2026-07-22 17:56:02 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:56:03 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:56:03 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS13.nc
2026-07-22 17:56:03 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 17:56:04 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS13.nc
2026-07-22 17:56:04 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:56:04 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:56:04 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens14.nc
2026-07-22 17:56:05 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:56:05 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS14.nc
2026-07-22 17:56:06 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 17:56:06 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS14.nc
2026-07-22 17:56:06 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:56:07 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:56:07 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS14.nc
2026-07-22 17:56:08 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 17:56:08 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS14.nc
2026-07-22 17:56:08 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:56:09 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:56:10 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens15.nc
2026-07-22 17:56:10 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:56:10 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS15.nc
2026-07-22 17:56:11 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 17:56:11 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS15.nc
2026-07-22 17:56:11 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:56:11 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 17:56:11 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS15.nc
2026-07-22 17:56:12 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 17:56:13 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS15.nc
2026-07-22 17:56:13 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 17:56:13 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 17:56:13 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-22 17:56:13 INFO [TIME] step_end current_time=2020-02-07 09:00:00 simulated_time=2020-02-07 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 17:56:13 INFO [TIME] step_start current_time=2020-02-07 09:00:00 simulated_time=2020-02-07 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 17:56:13 INFO [TIME] window start=2020-02-07 09:00:00 end=2020-02-07 11:00:00 run_hours=2 has_assimilation=True
2026-07-22 17:56:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:56:14 INFO Hourly dataset computed and listing created
2026-07-22 17:56:20 INFO Hourly dataset computed
2026-07-22 17:56:20 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:56:21 INFO Hourly dataset computed and listing created
2026-07-22 17:56:22 INFO Hourly dataset computed
2026-07-22 17:56:22 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:56:23 INFO Hourly dataset computed and listing created
2026-07-22 17:56:24 INFO Hourly dataset computed
2026-07-22 17:56:24 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:56:25 INFO Hourly dataset computed and listing created
2026-07-22 17:56:26 INFO Hourly dataset computed
2026-07-22 17:56:26 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:56:27 INFO Hourly dataset computed and listing created
2026-07-22 17:56:27 INFO Hourly dataset computed
2026-07-22 17:56:27 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:56:28 INFO Hourly dataset computed and listing created
2026-07-22 17:56:29 INFO Hourly dataset computed
2026-07-22 17:56:29 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:56:30 INFO Hourly dataset computed and listing created
2026-07-22 17:56:31 INFO Hourly dataset computed
2026-07-22 17:56:31 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:56:32 INFO Hourly dataset computed and listing created
2026-07-22 17:56:32 INFO Hourly dataset computed
2026-07-22 17:56:33 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:56:33 INFO Hourly dataset computed and listing created
2026-07-22 17:56:34 INFO Hourly dataset computed
2026-07-22 17:56:34 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:56:35 INFO Hourly dataset computed and listing created
2026-07-22 17:56:36 INFO Hourly dataset computed
2026-07-22 17:56:36 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:56:37 INFO Hourly dataset computed and listing created
2026-07-22 17:56:38 INFO Hourly dataset computed
2026-07-22 17:56:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:56:39 INFO Hourly dataset computed and listing created
2026-07-22 17:56:39 INFO Hourly dataset computed
2026-07-22 17:56:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:56:40 INFO Hourly dataset computed and listing created
2026-07-22 17:56:41 INFO Hourly dataset computed
2026-07-22 17:56:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:56:42 INFO Hourly dataset computed and listing created
2026-07-22 17:56:43 INFO Hourly dataset computed
2026-07-22 17:56:43 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 17:56:44 INFO Hourly dataset computed and listing created
2026-07-22 17:56:45 INFO Hourly dataset computed
2026-07-22 17:56:45 INFO ---------->>> Running CHIMERE model from 2020-02-07 09:00:00 to 2020-02-07 11:00:00
2026-07-22 17:56:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:56:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1
2026-07-22 17:56:45 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020701_8_ENS1.nc
2026-07-22 17:56:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-22 17:56:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:56:45 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-22 17:56:45 INFO Queuing job for member 1...
2026-07-22 17:56:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:56:45 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-22 17:56:46 INFO Found: ['5257338']
2026-07-22 17:56:51 INFO [TGCC-IRENE] Submitted job with ID:['5257338']
2026-07-22 17:56:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:56:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2
2026-07-22 17:56:51 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020701_8_ENS2.nc
2026-07-22 17:56:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-22 17:56:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:56:51 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-22 17:56:51 INFO Queuing job for member 2...
2026-07-22 17:56:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:56:51 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-22 17:56:54 INFO Found: ['5257339']
2026-07-22 17:56:59 INFO [TGCC-IRENE] Submitted job with ID:['5257339']
2026-07-22 17:56:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:56:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3
2026-07-22 17:56:59 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020701_8_ENS3.nc
2026-07-22 17:56:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-22 17:56:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:56:59 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-22 17:56:59 INFO Queuing job for member 3...
2026-07-22 17:56:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:56:59 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-22 17:57:01 INFO Found: ['5257340']
2026-07-22 17:57:06 INFO [TGCC-IRENE] Submitted job with ID:['5257340']
2026-07-22 17:57:06 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:57:06 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4
2026-07-22 17:57:06 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020701_8_ENS4.nc
2026-07-22 17:57:06 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-22 17:57:06 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:57:06 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-22 17:57:06 INFO Queuing job for member 4...
2026-07-22 17:57:06 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:57:06 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-22 17:57:07 INFO Found: ['5257342']
2026-07-22 17:57:12 INFO [TGCC-IRENE] Submitted job with ID:['5257342']
2026-07-22 17:57:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:57:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5
2026-07-22 17:57:12 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020701_8_ENS5.nc
2026-07-22 17:57:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-22 17:57:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:57:12 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-22 17:57:12 INFO Queuing job for member 5...
2026-07-22 17:57:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:57:12 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-22 17:57:13 INFO Found: ['5257343']
2026-07-22 17:57:18 INFO [TGCC-IRENE] Submitted job with ID:['5257343']
2026-07-22 17:57:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:57:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6
2026-07-22 17:57:18 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020701_8_ENS6.nc
2026-07-22 17:57:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-22 17:57:18 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:57:18 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-22 17:57:18 INFO Queuing job for member 6...
2026-07-22 17:57:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:57:18 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-22 17:57:19 INFO Found: ['5257344']
2026-07-22 17:57:24 INFO [TGCC-IRENE] Submitted job with ID:['5257344']
2026-07-22 17:57:24 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:57:24 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7
2026-07-22 17:57:24 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020701_8_ENS7.nc
2026-07-22 17:57:24 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-22 17:57:24 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:57:24 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-22 17:57:24 INFO Queuing job for member 7...
2026-07-22 17:57:24 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:57:24 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-22 17:57:24 INFO Found: ['5257345']
2026-07-22 17:57:29 INFO [TGCC-IRENE] Submitted job with ID:['5257345']
2026-07-22 17:57:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:57:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8
2026-07-22 17:57:29 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020701_8_ENS8.nc
2026-07-22 17:57:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-22 17:57:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:57:29 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-22 17:57:29 INFO Queuing job for member 8...
2026-07-22 17:57:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:57:29 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-22 17:57:30 INFO Found: ['5257346']
2026-07-22 17:57:35 INFO [TGCC-IRENE] Submitted job with ID:['5257346']
2026-07-22 17:57:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:57:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9
2026-07-22 17:57:35 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020701_8_ENS9.nc
2026-07-22 17:57:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-22 17:57:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:57:35 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-22 17:57:35 INFO Queuing job for member 9...
2026-07-22 17:57:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:57:35 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-22 17:57:36 INFO Found: ['5257347']
2026-07-22 17:57:41 INFO [TGCC-IRENE] Submitted job with ID:['5257347']
2026-07-22 17:57:41 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:57:41 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10
2026-07-22 17:57:41 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020701_8_ENS10.nc
2026-07-22 17:57:41 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-22 17:57:41 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:57:41 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-22 17:57:41 INFO Queuing job for member 10...
2026-07-22 17:57:41 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:57:41 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-22 17:57:42 INFO Found: ['5257349']
2026-07-22 17:57:47 INFO [TGCC-IRENE] Submitted job with ID:['5257349']
2026-07-22 17:57:47 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:57:47 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11
2026-07-22 17:57:47 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020701_8_ENS11.nc
2026-07-22 17:57:47 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-22 17:57:47 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:57:47 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-22 17:57:47 INFO Queuing job for member 11...
2026-07-22 17:57:47 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:57:47 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-22 17:57:49 INFO Found: ['5257353']
2026-07-22 17:57:54 INFO [TGCC-IRENE] Submitted job with ID:['5257353']
2026-07-22 17:57:54 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:57:54 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12
2026-07-22 17:57:54 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020701_8_ENS12.nc
2026-07-22 17:57:54 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-22 17:57:54 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:57:54 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-22 17:57:54 INFO Queuing job for member 12...
2026-07-22 17:57:54 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:57:54 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-22 17:57:57 INFO Found: ['5257356']
2026-07-22 17:58:02 INFO [TGCC-IRENE] Submitted job with ID:['5257356']
2026-07-22 17:58:02 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:58:02 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13
2026-07-22 17:58:02 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020701_8_ENS13.nc
2026-07-22 17:58:02 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-22 17:58:02 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:58:02 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-22 17:58:02 INFO Queuing job for member 13...
2026-07-22 17:58:02 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:58:02 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-22 17:58:04 INFO Found: ['5257363']
2026-07-22 17:58:09 INFO [TGCC-IRENE] Submitted job with ID:['5257363']
2026-07-22 17:58:09 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:58:09 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14
2026-07-22 17:58:09 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020701_8_ENS14.nc
2026-07-22 17:58:09 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-22 17:58:09 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:58:09 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-22 17:58:09 INFO Queuing job for member 14...
2026-07-22 17:58:09 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:58:09 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-22 17:58:12 INFO Found: ['5257367']
2026-07-22 17:58:17 INFO [TGCC-IRENE] Submitted job with ID:['5257367']
2026-07-22 17:58:17 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 17:58:17 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15
2026-07-22 17:58:17 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020701_8_ENS15.nc
2026-07-22 17:58:17 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-22 17:58:17 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 17:58:17 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-22 17:58:17 INFO Queuing job for member 15...
2026-07-22 17:58:17 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 17:58:17 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-22 17:58:19 INFO Found: ['5257368']
2026-07-22 17:58:24 INFO [TGCC-IRENE] Submitted job with ID:['5257368']
2026-07-22 17:58:24 INFO Checking job status ...
2026-07-22 17:58:24 INFO None 5257338: status RUNNING/PENDING
2026-07-22 17:58:24 INFO None 5257339: status RUNNING/PENDING
2026-07-22 17:58:24 INFO None 5257340: status RUNNING/PENDING
2026-07-22 17:58:24 INFO None 5257342: status RUNNING/PENDING
2026-07-22 17:58:24 INFO None 5257343: status RUNNING/PENDING
2026-07-22 17:58:24 INFO None 5257344: status RUNNING/PENDING
2026-07-22 17:58:24 INFO None 5257345: status RUNNING/PENDING
2026-07-22 17:58:24 INFO None 5257346: status RUNNING/PENDING
2026-07-22 17:58:24 INFO None 5257347: status RUNNING/PENDING
2026-07-22 17:58:24 INFO None 5257349: status RUNNING/PENDING
2026-07-22 17:58:24 INFO None 5257353: status RUNNING/PENDING
2026-07-22 17:58:24 INFO None 5257356: status RUNNING/PENDING
2026-07-22 17:58:24 INFO None 5257363: status RUNNING/PENDING
2026-07-22 17:58:24 INFO None 5257367: status RUNNING/PENDING
2026-07-22 17:58:24 INFO None 5257368: status RUNNING/PENDING
2026-07-22 17:58:24 INFO Jobs still running: ['5257338', '5257339', '5257340', '5257342', '5257343', '5257344', '5257345', '5257346', '5257347', '5257349', '5257353', '5257356', '5257363', '5257367', '5257368']. Waiting...
2026-07-22 17:58:39 INFO None 5257338: status RUNNING/PENDING
2026-07-22 17:58:39 INFO None 5257339: status RUNNING/PENDING
2026-07-22 17:58:39 INFO None 5257340: status RUNNING/PENDING
2026-07-22 17:58:39 INFO None 5257342: status RUNNING/PENDING
2026-07-22 17:58:39 INFO None 5257343: status RUNNING/PENDING
2026-07-22 17:58:39 INFO None 5257344: status RUNNING/PENDING
2026-07-22 17:58:39 INFO None 5257345: status RUNNING/PENDING
2026-07-22 17:58:39 INFO None 5257346: status RUNNING/PENDING
2026-07-22 17:58:40 INFO None 5257347: status RUNNING/PENDING
2026-07-22 17:58:40 INFO None 5257349: status RUNNING/PENDING
2026-07-22 17:58:40 INFO None 5257353: status RUNNING/PENDING
2026-07-22 17:58:40 INFO None 5257356: status RUNNING/PENDING
2026-07-22 17:58:40 INFO None 5257363: status RUNNING/PENDING
2026-07-22 17:58:40 INFO None 5257367: status RUNNING/PENDING
2026-07-22 17:58:40 INFO None 5257368: status RUNNING/PENDING
2026-07-22 17:58:40 INFO Jobs still running: ['5257338', '5257339', '5257340', '5257342', '5257343', '5257344', '5257345', '5257346', '5257347', '5257349', '5257353', '5257356', '5257363', '5257367', '5257368']. Waiting...
2026-07-22 17:58:55 INFO None 5257338: status RUNNING/PENDING
2026-07-22 17:58:55 INFO None 5257339: status RUNNING/PENDING
2026-07-22 17:58:55 INFO None 5257340: status RUNNING/PENDING
2026-07-22 17:58:55 INFO None 5257342: status RUNNING/PENDING
2026-07-22 17:58:55 INFO None 5257343: status RUNNING/PENDING
2026-07-22 17:58:55 INFO None 5257344: status RUNNING/PENDING
2026-07-22 17:58:55 INFO None 5257345: status RUNNING/PENDING
2026-07-22 17:58:55 INFO None 5257346: status RUNNING/PENDING
2026-07-22 17:58:55 INFO None 5257347: status RUNNING/PENDING
2026-07-22 17:58:55 INFO None 5257349: status RUNNING/PENDING
2026-07-22 17:58:55 INFO None 5257353: status RUNNING/PENDING
2026-07-22 17:58:55 INFO None 5257356: status RUNNING/PENDING
2026-07-22 17:58:55 INFO None 5257363: status RUNNING/PENDING
2026-07-22 17:58:55 INFO None 5257367: status RUNNING/PENDING
2026-07-22 17:58:55 INFO None 5257368: status RUNNING/PENDING
2026-07-22 17:58:55 INFO Jobs still running: ['5257338', '5257339', '5257340', '5257342', '5257343', '5257344', '5257345', '5257346', '5257347', '5257349', '5257353', '5257356', '5257363', '5257367', '5257368']. Waiting...
2026-07-22 17:59:10 INFO None 5257338: status RUNNING/PENDING
2026-07-22 17:59:10 INFO None 5257339: status RUNNING/PENDING
2026-07-22 17:59:10 INFO None 5257340: status RUNNING/PENDING
2026-07-22 17:59:10 INFO None 5257342: status RUNNING/PENDING
2026-07-22 17:59:10 INFO None 5257343: status RUNNING/PENDING
2026-07-22 17:59:10 INFO None 5257344: status RUNNING/PENDING
2026-07-22 17:59:10 INFO None 5257345: status RUNNING/PENDING
2026-07-22 17:59:10 INFO None 5257346: status RUNNING/PENDING
2026-07-22 17:59:10 INFO None 5257347: status RUNNING/PENDING
2026-07-22 17:59:10 INFO None 5257349: status RUNNING/PENDING
2026-07-22 17:59:10 INFO None 5257353: status RUNNING/PENDING
2026-07-22 17:59:10 INFO None 5257356: status RUNNING/PENDING
2026-07-22 17:59:10 INFO None 5257363: status RUNNING/PENDING
2026-07-22 17:59:10 INFO None 5257367: status RUNNING/PENDING
2026-07-22 17:59:10 INFO None 5257368: status RUNNING/PENDING
2026-07-22 17:59:10 INFO Jobs still running: ['5257338', '5257339', '5257340', '5257342', '5257343', '5257344', '5257345', '5257346', '5257347', '5257349', '5257353', '5257356', '5257363', '5257367', '5257368']. Waiting...
2026-07-22 17:59:27 INFO None 5257338: status RUNNING/PENDING
2026-07-22 17:59:27 INFO None 5257339: status RUNNING/PENDING
2026-07-22 17:59:27 INFO None 5257340: status RUNNING/PENDING
2026-07-22 17:59:27 INFO None 5257342: status RUNNING/PENDING
2026-07-22 17:59:27 INFO None 5257343: status RUNNING/PENDING
2026-07-22 17:59:27 INFO None 5257344: status RUNNING/PENDING
2026-07-22 17:59:27 INFO None 5257345: status RUNNING/PENDING
2026-07-22 17:59:27 INFO None 5257346: status RUNNING/PENDING
2026-07-22 17:59:27 INFO None 5257347: status RUNNING/PENDING
2026-07-22 17:59:27 INFO None 5257349: status RUNNING/PENDING
2026-07-22 17:59:27 INFO None 5257353: status RUNNING/PENDING
2026-07-22 17:59:27 INFO None 5257356: status RUNNING/PENDING
2026-07-22 17:59:27 INFO None 5257363: status RUNNING/PENDING
2026-07-22 17:59:27 INFO None 5257367: status RUNNING/PENDING
2026-07-22 17:59:27 INFO None 5257368: status RUNNING/PENDING
2026-07-22 17:59:27 INFO Jobs still running: ['5257338', '5257339', '5257340', '5257342', '5257343', '5257344', '5257345', '5257346', '5257347', '5257349', '5257353', '5257356', '5257363', '5257367', '5257368']. Waiting...
2026-07-22 17:59:42 INFO None 5257338: status RUNNING/PENDING
2026-07-22 17:59:42 INFO None 5257339: status RUNNING/PENDING
2026-07-22 17:59:42 INFO None 5257340: status RUNNING/PENDING
2026-07-22 17:59:42 INFO None 5257342: status RUNNING/PENDING
2026-07-22 17:59:42 INFO None 5257343: status RUNNING/PENDING
2026-07-22 17:59:42 INFO None 5257344: status RUNNING/PENDING
2026-07-22 17:59:42 INFO None 5257345: status RUNNING/PENDING
2026-07-22 17:59:42 INFO None 5257346: status RUNNING/PENDING
2026-07-22 17:59:42 INFO None 5257347: status RUNNING/PENDING
2026-07-22 17:59:42 INFO None 5257349: status RUNNING/PENDING
2026-07-22 17:59:42 INFO None 5257353: status RUNNING/PENDING
2026-07-22 17:59:44 INFO None 5257356: status RUNNING/PENDING
2026-07-22 17:59:44 INFO None 5257363: status RUNNING/PENDING
2026-07-22 17:59:44 INFO None 5257367: status RUNNING/PENDING
2026-07-22 17:59:44 INFO None 5257368: status RUNNING/PENDING
2026-07-22 17:59:44 INFO Jobs still running: ['5257338', '5257339', '5257340', '5257342', '5257343', '5257344', '5257345', '5257346', '5257347', '5257349', '5257353', '5257356', '5257363', '5257367', '5257368']. Waiting...
2026-07-22 17:59:59 INFO None 5257338: status RUNNING/PENDING
2026-07-22 18:00:00 INFO None 5257339: status RUNNING/PENDING
2026-07-22 18:00:00 INFO None 5257340: status RUNNING/PENDING
2026-07-22 18:00:00 INFO None 5257342: status RUNNING/PENDING
2026-07-22 18:00:00 INFO None 5257343: status RUNNING/PENDING
2026-07-22 18:00:00 INFO None 5257344: status RUNNING/PENDING
2026-07-22 18:00:00 INFO None 5257345: status RUNNING/PENDING
2026-07-22 18:00:00 INFO None 5257346: status RUNNING/PENDING
2026-07-22 18:00:00 INFO None 5257347: status RUNNING/PENDING
2026-07-22 18:00:00 INFO None 5257349: status RUNNING/PENDING
2026-07-22 18:00:00 INFO None 5257353: status RUNNING/PENDING
2026-07-22 18:00:00 INFO None 5257356: status RUNNING/PENDING
2026-07-22 18:00:00 INFO None 5257363: status RUNNING/PENDING
2026-07-22 18:00:00 INFO None 5257367: status RUNNING/PENDING
2026-07-22 18:00:00 INFO None 5257368: status RUNNING/PENDING
2026-07-22 18:00:00 INFO Jobs still running: ['5257338', '5257339', '5257340', '5257342', '5257343', '5257344', '5257345', '5257346', '5257347', '5257349', '5257353', '5257356', '5257363', '5257367', '5257368']. Waiting...
2026-07-22 18:00:15 INFO None 5257338: status RUNNING/PENDING
2026-07-22 18:00:15 INFO None 5257339: status RUNNING/PENDING
2026-07-22 18:00:15 INFO None 5257340: status RUNNING/PENDING
2026-07-22 18:00:15 INFO None 5257342: status RUNNING/PENDING
2026-07-22 18:00:15 INFO None 5257343: status RUNNING/PENDING
2026-07-22 18:00:15 INFO None 5257344: status RUNNING/PENDING
2026-07-22 18:00:15 INFO None 5257345: status RUNNING/PENDING
2026-07-22 18:00:15 INFO None 5257346: status RUNNING/PENDING
2026-07-22 18:00:15 INFO None 5257347: status RUNNING/PENDING
2026-07-22 18:00:15 INFO None 5257349: status RUNNING/PENDING
2026-07-22 18:00:15 INFO None 5257353: status RUNNING/PENDING
2026-07-22 18:00:15 INFO None 5257356: status RUNNING/PENDING
2026-07-22 18:00:15 INFO None 5257363: status RUNNING/PENDING
2026-07-22 18:00:15 INFO None 5257367: status RUNNING/PENDING
2026-07-22 18:00:15 INFO None 5257368: status RUNNING/PENDING
2026-07-22 18:00:15 INFO Jobs still running: ['5257338', '5257339', '5257340', '5257342', '5257343', '5257344', '5257345', '5257346', '5257347', '5257349', '5257353', '5257356', '5257363', '5257367', '5257368']. Waiting...
2026-07-22 18:00:30 INFO None 5257338: status RUNNING/PENDING
2026-07-22 18:00:30 INFO None 5257339: status RUNNING/PENDING
2026-07-22 18:00:30 INFO None 5257340: status RUNNING/PENDING
2026-07-22 18:00:30 INFO None 5257342: status RUNNING/PENDING
2026-07-22 18:00:30 INFO None 5257343: status RUNNING/PENDING
2026-07-22 18:00:30 INFO None 5257344: status RUNNING/PENDING
2026-07-22 18:00:30 INFO None 5257345: status RUNNING/PENDING
2026-07-22 18:00:30 INFO None 5257346: status RUNNING/PENDING
2026-07-22 18:00:30 INFO None 5257347: status RUNNING/PENDING
2026-07-22 18:00:30 INFO None 5257349: status RUNNING/PENDING
2026-07-22 18:00:30 INFO None 5257353: status RUNNING/PENDING
2026-07-22 18:00:30 INFO None 5257356: status RUNNING/PENDING
2026-07-22 18:00:30 INFO None 5257363: status RUNNING/PENDING
2026-07-22 18:00:30 INFO None 5257367: status RUNNING/PENDING
2026-07-22 18:00:30 INFO None 5257368: status RUNNING/PENDING
2026-07-22 18:00:30 INFO Jobs still running: ['5257338', '5257339', '5257340', '5257342', '5257343', '5257344', '5257345', '5257346', '5257347', '5257349', '5257353', '5257356', '5257363', '5257367', '5257368']. Waiting...
2026-07-22 18:00:45 INFO None 5257338: status RUNNING/PENDING
2026-07-22 18:00:45 INFO None 5257339: status RUNNING/PENDING
2026-07-22 18:00:45 INFO None 5257340: status RUNNING/PENDING
2026-07-22 18:00:45 INFO None 5257342: status RUNNING/PENDING
2026-07-22 18:00:45 INFO None 5257343: status RUNNING/PENDING
2026-07-22 18:00:46 INFO None 5257344: status RUNNING/PENDING
2026-07-22 18:00:46 INFO None 5257345: status RUNNING/PENDING
2026-07-22 18:00:46 INFO None 5257346: status RUNNING/PENDING
2026-07-22 18:00:46 INFO None 5257347: status RUNNING/PENDING
2026-07-22 18:00:46 INFO None 5257349: status RUNNING/PENDING
2026-07-22 18:00:46 INFO None 5257353: status RUNNING/PENDING
2026-07-22 18:00:46 INFO None 5257356: status RUNNING/PENDING
2026-07-22 18:00:46 INFO None 5257363: status RUNNING/PENDING
2026-07-22 18:00:46 INFO None 5257367: status RUNNING/PENDING
2026-07-22 18:00:46 INFO None 5257368: status RUNNING/PENDING
2026-07-22 18:00:46 INFO Jobs still running: ['5257338', '5257339', '5257340', '5257342', '5257343', '5257344', '5257345', '5257346', '5257347', '5257349', '5257353', '5257356', '5257363', '5257367', '5257368']. Waiting...
2026-07-22 18:01:01 INFO None 5257338: status RUNNING/PENDING
2026-07-22 18:01:01 INFO None 5257339: status RUNNING/PENDING
2026-07-22 18:01:01 INFO None 5257340: status FINISHED
2026-07-22 18:01:01 INFO None 5257342: status FINISHED
2026-07-22 18:01:01 INFO None 5257343: status RUNNING/PENDING
2026-07-22 18:01:01 INFO None 5257344: status RUNNING/PENDING
2026-07-22 18:01:01 INFO None 5257345: status RUNNING/PENDING
2026-07-22 18:01:01 INFO None 5257346: status RUNNING/PENDING
2026-07-22 18:01:01 INFO None 5257347: status RUNNING/PENDING
2026-07-22 18:01:01 INFO None 5257349: status RUNNING/PENDING
2026-07-22 18:01:01 INFO None 5257353: status RUNNING/PENDING
2026-07-22 18:01:01 INFO None 5257356: status RUNNING/PENDING
2026-07-22 18:01:01 INFO None 5257363: status RUNNING/PENDING
2026-07-22 18:01:01 INFO None 5257367: status RUNNING/PENDING
2026-07-22 18:01:01 INFO None 5257368: status RUNNING/PENDING
2026-07-22 18:01:01 INFO Jobs still running: ['5257338', '5257339', '5257343', '5257344', '5257345', '5257346', '5257347', '5257349', '5257353', '5257356', '5257363', '5257367', '5257368']. Waiting...
2026-07-22 18:01:16 INFO None 5257338: status RUNNING/PENDING
2026-07-22 18:01:16 INFO None 5257339: status FINISHED
2026-07-22 18:01:16 INFO None 5257340: status FINISHED
2026-07-22 18:01:16 INFO None 5257342: status FINISHED
2026-07-22 18:01:16 INFO None 5257343: status RUNNING/PENDING
2026-07-22 18:01:16 INFO None 5257344: status RUNNING/PENDING
2026-07-22 18:01:16 INFO None 5257345: status FINISHED
2026-07-22 18:01:16 INFO None 5257346: status RUNNING/PENDING
2026-07-22 18:01:16 INFO None 5257347: status RUNNING/PENDING
2026-07-22 18:01:16 INFO None 5257349: status RUNNING/PENDING
2026-07-22 18:01:16 INFO None 5257353: status RUNNING/PENDING
2026-07-22 18:01:17 INFO None 5257356: status RUNNING/PENDING
2026-07-22 18:01:17 INFO None 5257363: status RUNNING/PENDING
2026-07-22 18:01:19 INFO None 5257367: status RUNNING/PENDING
2026-07-22 18:01:19 INFO None 5257368: status RUNNING/PENDING
2026-07-22 18:01:19 INFO Jobs still running: ['5257338', '5257343', '5257344', '5257346', '5257347', '5257349', '5257353', '5257356', '5257363', '5257367', '5257368']. Waiting...
2026-07-22 18:01:34 INFO None 5257338: status RUNNING/PENDING
2026-07-22 18:01:34 INFO None 5257339: status FINISHED
2026-07-22 18:01:34 INFO None 5257340: status FINISHED
2026-07-22 18:01:34 INFO None 5257342: status FINISHED
2026-07-22 18:01:34 INFO None 5257343: status RUNNING/PENDING
2026-07-22 18:01:34 INFO None 5257344: status RUNNING/PENDING
2026-07-22 18:01:34 INFO None 5257345: status FINISHED
2026-07-22 18:01:34 INFO None 5257346: status RUNNING/PENDING
2026-07-22 18:01:34 INFO None 5257347: status RUNNING/PENDING
2026-07-22 18:01:34 INFO None 5257349: status RUNNING/PENDING
2026-07-22 18:01:34 INFO None 5257353: status RUNNING/PENDING
2026-07-22 18:01:34 INFO None 5257356: status RUNNING/PENDING
2026-07-22 18:01:34 INFO None 5257363: status RUNNING/PENDING
2026-07-22 18:01:34 INFO None 5257367: status RUNNING/PENDING
2026-07-22 18:01:34 INFO None 5257368: status RUNNING/PENDING
2026-07-22 18:01:34 INFO Jobs still running: ['5257338', '5257343', '5257344', '5257346', '5257347', '5257349', '5257353', '5257356', '5257363', '5257367', '5257368']. Waiting...
2026-07-22 18:01:49 INFO None 5257338: status FINISHED
2026-07-22 18:01:49 INFO None 5257339: status FINISHED
2026-07-22 18:01:49 INFO None 5257340: status FINISHED
2026-07-22 18:01:49 INFO None 5257342: status FINISHED
2026-07-22 18:01:49 INFO None 5257343: status RUNNING/PENDING
2026-07-22 18:01:49 INFO None 5257344: status RUNNING/PENDING
2026-07-22 18:01:49 INFO None 5257345: status FINISHED
2026-07-22 18:01:49 INFO None 5257346: status RUNNING/PENDING
2026-07-22 18:01:49 INFO None 5257347: status RUNNING/PENDING
2026-07-22 18:01:49 INFO None 5257349: status RUNNING/PENDING
2026-07-22 18:01:49 INFO None 5257353: status RUNNING/PENDING
2026-07-22 18:01:49 INFO None 5257356: status RUNNING/PENDING
2026-07-22 18:01:49 INFO None 5257363: status RUNNING/PENDING
2026-07-22 18:01:49 INFO None 5257367: status RUNNING/PENDING
2026-07-22 18:01:51 INFO None 5257368: status RUNNING/PENDING
2026-07-22 18:01:51 INFO Jobs still running: ['5257343', '5257344', '5257346', '5257347', '5257349', '5257353', '5257356', '5257363', '5257367', '5257368']. Waiting...
2026-07-22 18:02:06 INFO None 5257338: status FINISHED
2026-07-22 18:02:06 INFO None 5257339: status FINISHED
2026-07-22 18:02:06 INFO None 5257340: status FINISHED
2026-07-22 18:02:06 INFO None 5257342: status FINISHED
2026-07-22 18:02:06 INFO None 5257343: status FINISHED
2026-07-22 18:02:06 INFO None 5257344: status FINISHED
2026-07-22 18:02:06 INFO None 5257345: status FINISHED
2026-07-22 18:02:06 INFO None 5257346: status FINISHED
2026-07-22 18:02:06 INFO None 5257347: status FINISHED
2026-07-22 18:02:06 INFO None 5257349: status FINISHED
2026-07-22 18:02:06 INFO None 5257353: status RUNNING/PENDING
2026-07-22 18:02:06 INFO None 5257356: status RUNNING/PENDING
2026-07-22 18:02:06 INFO None 5257363: status RUNNING/PENDING
2026-07-22 18:02:06 INFO None 5257367: status RUNNING/PENDING
2026-07-22 18:02:06 INFO None 5257368: status RUNNING/PENDING
2026-07-22 18:02:06 INFO Jobs still running: ['5257353', '5257356', '5257363', '5257367', '5257368']. Waiting...
2026-07-22 18:02:22 INFO None 5257338: status FINISHED
2026-07-22 18:02:22 INFO None 5257339: status FINISHED
2026-07-22 18:02:22 INFO None 5257340: status FINISHED
2026-07-22 18:02:22 INFO None 5257342: status FINISHED
2026-07-22 18:02:22 INFO None 5257343: status FINISHED
2026-07-22 18:02:22 INFO None 5257344: status FINISHED
2026-07-22 18:02:22 INFO None 5257345: status FINISHED
2026-07-22 18:02:22 INFO None 5257346: status FINISHED
2026-07-22 18:02:22 INFO None 5257347: status FINISHED
2026-07-22 18:02:22 INFO None 5257349: status FINISHED
2026-07-22 18:02:22 INFO None 5257353: status RUNNING/PENDING
2026-07-22 18:02:22 INFO None 5257356: status RUNNING/PENDING
2026-07-22 18:02:22 INFO None 5257363: status RUNNING/PENDING
2026-07-22 18:02:22 INFO None 5257367: status RUNNING/PENDING
2026-07-22 18:02:22 INFO None 5257368: status RUNNING/PENDING
2026-07-22 18:02:22 INFO Jobs still running: ['5257353', '5257356', '5257363', '5257367', '5257368']. Waiting...
2026-07-22 18:02:37 INFO None 5257338: status FINISHED
2026-07-22 18:02:37 INFO None 5257339: status FINISHED
2026-07-22 18:02:37 INFO None 5257340: status FINISHED
2026-07-22 18:02:37 INFO None 5257342: status FINISHED
2026-07-22 18:02:37 INFO None 5257343: status FINISHED
2026-07-22 18:02:37 INFO None 5257344: status FINISHED
2026-07-22 18:02:37 INFO None 5257345: status FINISHED
2026-07-22 18:02:37 INFO None 5257346: status FINISHED
2026-07-22 18:02:39 INFO None 5257347: status FINISHED
2026-07-22 18:02:39 INFO None 5257349: status FINISHED
2026-07-22 18:02:39 INFO None 5257353: status RUNNING/PENDING
2026-07-22 18:02:39 INFO None 5257356: status RUNNING/PENDING
2026-07-22 18:02:39 INFO None 5257363: status RUNNING/PENDING
2026-07-22 18:02:39 INFO None 5257367: status RUNNING/PENDING
2026-07-22 18:02:39 INFO None 5257368: status RUNNING/PENDING
2026-07-22 18:02:39 INFO Jobs still running: ['5257353', '5257356', '5257363', '5257367', '5257368']. Waiting...
2026-07-22 18:02:54 INFO None 5257338: status FINISHED
2026-07-22 18:02:54 INFO None 5257339: status FINISHED
2026-07-22 18:02:54 INFO None 5257340: status FINISHED
2026-07-22 18:02:54 INFO None 5257342: status FINISHED
2026-07-22 18:02:54 INFO None 5257343: status FINISHED
2026-07-22 18:02:54 INFO None 5257344: status FINISHED
2026-07-22 18:02:54 INFO None 5257345: status FINISHED
2026-07-22 18:02:54 INFO None 5257346: status FINISHED
2026-07-22 18:02:54 INFO None 5257347: status FINISHED
2026-07-22 18:02:54 INFO None 5257349: status FINISHED
2026-07-22 18:02:54 INFO None 5257353: status RUNNING/PENDING
2026-07-22 18:02:54 INFO None 5257356: status RUNNING/PENDING
2026-07-22 18:02:54 INFO None 5257363: status RUNNING/PENDING
2026-07-22 18:02:54 INFO None 5257367: status RUNNING/PENDING
2026-07-22 18:02:54 INFO None 5257368: status RUNNING/PENDING
2026-07-22 18:02:54 INFO Jobs still running: ['5257353', '5257356', '5257363', '5257367', '5257368']. Waiting...
2026-07-22 18:03:09 INFO None 5257338: status FINISHED
2026-07-22 18:03:09 INFO None 5257339: status FINISHED
2026-07-22 18:03:09 INFO None 5257340: status FINISHED
2026-07-22 18:03:09 INFO None 5257342: status FINISHED
2026-07-22 18:03:09 INFO None 5257343: status FINISHED
2026-07-22 18:03:09 INFO None 5257344: status FINISHED
2026-07-22 18:03:09 INFO None 5257345: status FINISHED
2026-07-22 18:03:10 INFO None 5257346: status FINISHED
2026-07-22 18:03:10 INFO None 5257347: status FINISHED
2026-07-22 18:03:10 INFO None 5257349: status FINISHED
2026-07-22 18:03:12 INFO None 5257353: status RUNNING/PENDING
2026-07-22 18:03:12 INFO None 5257356: status RUNNING/PENDING
2026-07-22 18:03:12 INFO None 5257363: status RUNNING/PENDING
2026-07-22 18:03:12 INFO None 5257367: status RUNNING/PENDING
2026-07-22 18:03:12 INFO None 5257368: status RUNNING/PENDING
2026-07-22 18:03:12 INFO Jobs still running: ['5257353', '5257356', '5257363', '5257367', '5257368']. Waiting...
2026-07-22 18:03:27 INFO None 5257338: status FINISHED
2026-07-22 18:03:27 INFO None 5257339: status FINISHED
2026-07-22 18:03:27 INFO None 5257340: status FINISHED
2026-07-22 18:03:27 INFO None 5257342: status FINISHED
2026-07-22 18:03:27 INFO None 5257343: status FINISHED
2026-07-22 18:03:27 INFO None 5257344: status FINISHED
2026-07-22 18:03:27 INFO None 5257345: status FINISHED
2026-07-22 18:03:27 INFO None 5257346: status FINISHED
2026-07-22 18:03:27 INFO None 5257347: status FINISHED
2026-07-22 18:03:27 INFO None 5257349: status FINISHED
2026-07-22 18:03:27 INFO None 5257353: status RUNNING/PENDING
2026-07-22 18:03:27 INFO None 5257356: status RUNNING/PENDING
2026-07-22 18:03:27 INFO None 5257363: status RUNNING/PENDING
2026-07-22 18:03:27 INFO None 5257367: status RUNNING/PENDING
2026-07-22 18:03:27 INFO None 5257368: status FINISHED
2026-07-22 18:03:27 INFO Jobs still running: ['5257353', '5257356', '5257363', '5257367']. Waiting...
2026-07-22 18:03:42 INFO None 5257338: status FINISHED
2026-07-22 18:03:42 INFO None 5257339: status FINISHED
2026-07-22 18:03:42 INFO None 5257340: status FINISHED
2026-07-22 18:03:42 INFO None 5257342: status FINISHED
2026-07-22 18:03:42 INFO None 5257343: status FINISHED
2026-07-22 18:03:42 INFO None 5257344: status FINISHED
2026-07-22 18:03:42 INFO None 5257345: status FINISHED
2026-07-22 18:03:42 INFO None 5257346: status FINISHED
2026-07-22 18:03:42 INFO None 5257347: status FINISHED
2026-07-22 18:03:42 INFO None 5257349: status FINISHED
2026-07-22 18:03:42 INFO None 5257353: status FINISHED
2026-07-22 18:03:42 INFO None 5257356: status FINISHED
2026-07-22 18:03:42 INFO None 5257363: status FINISHED
2026-07-22 18:03:42 INFO None 5257367: status FINISHED
2026-07-22 18:03:42 INFO None 5257368: status FINISHED
2026-07-22 18:03:42 INFO Jobs ['5257338', '5257339', '5257340', '5257342', '5257343', '5257344', '5257345', '5257346', '5257347', '5257349', '5257353', '5257356', '5257363', '5257367', '5257368'] have finished
2026-07-22 18:03:42 INFO Checking restart files were created ...
2026-07-22 18:03:42 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020709_2_ENS1.nc(1002685915 bytes)
2026-07-22 18:03:42 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020709_2_ENS2.nc(1002685915 bytes)
2026-07-22 18:03:42 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020709_2_ENS3.nc(1002685915 bytes)
2026-07-22 18:03:42 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020709_2_ENS4.nc(1002685915 bytes)
2026-07-22 18:03:42 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020709_2_ENS5.nc(1002685915 bytes)
2026-07-22 18:03:42 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020709_2_ENS6.nc(1002685915 bytes)
2026-07-22 18:03:42 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020709_2_ENS7.nc(1002685915 bytes)
2026-07-22 18:03:42 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020709_2_ENS8.nc(1002685915 bytes)
2026-07-22 18:03:42 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020709_2_ENS9.nc(1002685915 bytes)
2026-07-22 18:03:42 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020709_2_ENS10.nc(1002685915 bytes)
2026-07-22 18:03:42 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020709_2_ENS11.nc(1002685915 bytes)
2026-07-22 18:03:42 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020709_2_ENS12.nc(1002685915 bytes)
2026-07-22 18:03:42 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020709_2_ENS13.nc(1002685915 bytes)
2026-07-22 18:03:42 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020709_2_ENS14.nc(1002685915 bytes)
2026-07-22 18:03:42 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020709_2_ENS15.nc(1002685915 bytes)
2026-07-22 18:03:42 INFO  Run_model() completed successfully.
2026-07-22 18:03:42 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 09:00:00 simulated_time=2020-02-07 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 18:03:42 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 11:00:00 days=153073 seconds=39600
2026-07-22 18:03:42 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-22 18:03:42 INFO [TIME] increment current_time 2020-02-07 09:00:00 -> 2020-02-07 11:00:00
2026-07-22 18:03:42 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 11:00:00 simulated_time=2020-02-07 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 18:03:42 INFO ---------->>> Running process_satellite_data()
2026-07-22 18:03:42 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12017.nc
2026-07-22 18:03:42 INFO ---------->>> Running run_obs_converter()
2026-07-22 18:03:42 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_37907_153073.out
2026-07-22 18:03:42 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_37907_153073.out
2026-07-22 18:03:42 INFO ---------->>> Running DART
2026-07-22 18:03:42 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-22 18:03:42 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/chim_ENS1_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/chim_ENS1_2020020711_1_out_toDART.nc
2026-07-22 18:03:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/chim_ENS2_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/chim_ENS2_2020020711_1_out_toDART.nc
2026-07-22 18:03:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/chim_ENS3_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/chim_ENS3_2020020711_1_out_toDART.nc
2026-07-22 18:03:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/chim_ENS4_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/chim_ENS4_2020020711_1_out_toDART.nc
2026-07-22 18:03:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/chim_ENS5_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/chim_ENS5_2020020711_1_out_toDART.nc
2026-07-22 18:03:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/chim_ENS6_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/chim_ENS6_2020020711_1_out_toDART.nc
2026-07-22 18:03:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/chim_ENS7_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/chim_ENS7_2020020711_1_out_toDART.nc
2026-07-22 18:03:45 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/chim_ENS8_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/chim_ENS8_2020020711_1_out_toDART.nc
2026-07-22 18:03:45 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/chim_ENS9_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/chim_ENS9_2020020711_1_out_toDART.nc
2026-07-22 18:03:45 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/chim_ENS10_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/chim_ENS10_2020020711_1_out_toDART.nc
2026-07-22 18:03:46 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/chim_ENS11_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/chim_ENS11_2020020711_1_out_toDART.nc
2026-07-22 18:03:46 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/chim_ENS12_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/chim_ENS12_2020020711_1_out_toDART.nc
2026-07-22 18:03:46 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/chim_ENS13_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/chim_ENS13_2020020711_1_out_toDART.nc
2026-07-22 18:03:47 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/chim_ENS14_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/chim_ENS14_2020020711_1_out_toDART.nc
2026-07-22 18:03:47 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/chim_ENS15_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/chim_ENS15_2020020711_1_out_toDART.nc
2026-07-22 18:03:47 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-22 18:03:47 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-22 18:03:47 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-22 18:03:47 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-22 18:03:47 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-22 18:03:47 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-22 18:03:58 INFO Found: []
2026-07-22 18:03:58 INFO No job id returned by command ./run_filter.bsh
2026-07-22 18:03:58 INFO No monitoring will be performed
2026-07-22 18:03:58 INFO Moving DART output files to analysis and preassim directories for date 2020020711 if present ...
2026-07-22 18:03:58 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020711'
2026-07-22 18:03:58 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020711'
2026-07-22 18:03:58 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020711'
2026-07-22 18:03:58 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020711'
2026-07-22 18:03:58 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020711'
2026-07-22 18:03:59 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020711'
2026-07-22 18:03:59 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020711'
2026-07-22 18:03:59 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020711'
2026-07-22 18:03:59 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020711'
2026-07-22 18:03:59 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020711'
2026-07-22 18:03:59 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020711'
2026-07-22 18:03:59 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020711'
2026-07-22 18:03:59 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020711'
2026-07-22 18:03:59 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020711'
2026-07-22 18:03:59 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020711'
2026-07-22 18:03:59 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020711'
2026-07-22 18:03:59 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020711'
2026-07-22 18:03:59 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020711'
2026-07-22 18:03:59 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020711'
2026-07-22 18:03:59 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020711'
2026-07-22 18:03:59 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020711'
2026-07-22 18:03:59 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020711'
2026-07-22 18:03:59 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020711'
2026-07-22 18:03:59 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020711'
2026-07-22 18:03:59 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020711'
2026-07-22 18:03:59 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020711'
2026-07-22 18:03:59 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020711'
2026-07-22 18:03:59 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020711'
2026-07-22 18:03:59 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020711'
2026-07-22 18:03:59 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020711'
2026-07-22 18:03:59 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020711'
2026-07-22 18:03:59 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020711'
2026-07-22 18:03:59 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020711'
2026-07-22 18:03:59 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020711'
2026-07-22 18:03:59 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-22 18:03:59 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-22 18:03:59 INFO run_dart() is DONE.
2026-07-22 18:03:59 INFO ---------->>> Running update_pollutant_in_end()
2026-07-22 18:03:59 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:03:59 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS1.nc
2026-07-22 18:04:00 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:04:01 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS1.nc
2026-07-22 18:04:01 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:04:01 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:04:01 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS1.nc
2026-07-22 18:04:02 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:04:03 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS1.nc
2026-07-22 18:04:03 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:04:03 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 18:04:03 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:04:03 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS2.nc
2026-07-22 18:04:04 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:04:04 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS2.nc
2026-07-22 18:04:04 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:04:05 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:04:05 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS2.nc
2026-07-22 18:04:05 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:04:06 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS2.nc
2026-07-22 18:04:06 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:04:06 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 18:04:07 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:04:07 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS3.nc
2026-07-22 18:04:07 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:04:08 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS3.nc
2026-07-22 18:04:08 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:04:09 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:04:09 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS3.nc
2026-07-22 18:04:09 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:04:10 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS3.nc
2026-07-22 18:04:10 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:04:10 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 18:04:10 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:04:10 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS4.nc
2026-07-22 18:04:11 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:04:12 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS4.nc
2026-07-22 18:04:12 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:04:12 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:04:12 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS4.nc
2026-07-22 18:04:13 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:04:14 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS4.nc
2026-07-22 18:04:14 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:04:14 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 18:04:14 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:04:14 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS5.nc
2026-07-22 18:04:15 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:04:16 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS5.nc
2026-07-22 18:04:16 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:04:16 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:04:16 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS5.nc
2026-07-22 18:04:17 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:04:17 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS5.nc
2026-07-22 18:04:17 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:04:17 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 18:04:18 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:04:18 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS6.nc
2026-07-22 18:04:18 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:04:19 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS6.nc
2026-07-22 18:04:19 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:04:19 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:04:20 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS6.nc
2026-07-22 18:04:20 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:04:21 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS6.nc
2026-07-22 18:04:21 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:04:21 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 18:04:21 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:04:21 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS7.nc
2026-07-22 18:04:22 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:04:23 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS7.nc
2026-07-22 18:04:23 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:04:23 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:04:23 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS7.nc
2026-07-22 18:04:24 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:04:24 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS7.nc
2026-07-22 18:04:24 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:04:24 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 18:04:25 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:04:25 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS8.nc
2026-07-22 18:04:25 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:04:26 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS8.nc
2026-07-22 18:04:26 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:04:26 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:04:26 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS8.nc
2026-07-22 18:04:27 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:04:28 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS8.nc
2026-07-22 18:04:28 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:04:28 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 18:04:29 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:04:29 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS9.nc
2026-07-22 18:04:29 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:04:30 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS9.nc
2026-07-22 18:04:30 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:04:30 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:04:30 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS9.nc
2026-07-22 18:04:31 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:04:31 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS9.nc
2026-07-22 18:04:31 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:04:31 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 18:04:32 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:04:32 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS10.nc
2026-07-22 18:04:32 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:04:33 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS10.nc
2026-07-22 18:04:33 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:04:34 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:04:34 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS10.nc
2026-07-22 18:04:34 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:04:35 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS10.nc
2026-07-22 18:04:35 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:04:35 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 18:04:35 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:04:35 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS11.nc
2026-07-22 18:04:36 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:04:37 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS11.nc
2026-07-22 18:04:37 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:04:38 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:04:38 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS11.nc
2026-07-22 18:04:38 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:04:39 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS11.nc
2026-07-22 18:04:39 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:04:39 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 18:04:40 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:04:40 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS12.nc
2026-07-22 18:04:40 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:04:41 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS12.nc
2026-07-22 18:04:41 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:04:42 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:04:42 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS12.nc
2026-07-22 18:04:42 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:04:43 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS12.nc
2026-07-22 18:04:43 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:04:43 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 18:04:44 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:04:44 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS13.nc
2026-07-22 18:04:44 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:04:45 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS13.nc
2026-07-22 18:04:45 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:04:45 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:04:45 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS13.nc
2026-07-22 18:04:46 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:04:47 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS13.nc
2026-07-22 18:04:47 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:04:47 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 18:04:47 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:04:47 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS14.nc
2026-07-22 18:04:48 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:04:49 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS14.nc
2026-07-22 18:04:49 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:04:49 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:04:49 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS14.nc
2026-07-22 18:04:50 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:04:50 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS14.nc
2026-07-22 18:04:50 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:04:51 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 18:04:51 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:04:51 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS15.nc
2026-07-22 18:04:52 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:04:52 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS15.nc
2026-07-22 18:04:52 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:04:53 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:04:53 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS15.nc
2026-07-22 18:04:53 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:04:54 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS15.nc
2026-07-22 18:04:54 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:04:54 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 18:04:54 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-22 18:04:54 INFO [TIME] step_end current_time=2020-02-07 11:00:00 simulated_time=2020-02-07 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 18:04:54 INFO [TIME] step_start current_time=2020-02-07 11:00:00 simulated_time=2020-02-07 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 18:04:54 INFO [TIME] window start=2020-02-07 11:00:00 end=2020-02-07 12:00:00 run_hours=1 has_assimilation=True
2026-07-22 18:04:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 18:04:55 INFO Hourly dataset computed and listing created
2026-07-22 18:04:58 INFO Hourly dataset computed
2026-07-22 18:04:58 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 18:04:59 INFO Hourly dataset computed and listing created
2026-07-22 18:04:59 INFO Hourly dataset computed
2026-07-22 18:04:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 18:05:00 INFO Hourly dataset computed and listing created
2026-07-22 18:05:01 INFO Hourly dataset computed
2026-07-22 18:05:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 18:05:02 INFO Hourly dataset computed and listing created
2026-07-22 18:05:02 INFO Hourly dataset computed
2026-07-22 18:05:02 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 18:05:03 INFO Hourly dataset computed and listing created
2026-07-22 18:05:04 INFO Hourly dataset computed
2026-07-22 18:05:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 18:05:05 INFO Hourly dataset computed and listing created
2026-07-22 18:05:05 INFO Hourly dataset computed
2026-07-22 18:05:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 18:05:06 INFO Hourly dataset computed and listing created
2026-07-22 18:05:07 INFO Hourly dataset computed
2026-07-22 18:05:07 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 18:05:08 INFO Hourly dataset computed and listing created
2026-07-22 18:05:08 INFO Hourly dataset computed
2026-07-22 18:05:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 18:05:09 INFO Hourly dataset computed and listing created
2026-07-22 18:05:10 INFO Hourly dataset computed
2026-07-22 18:05:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 18:05:11 INFO Hourly dataset computed and listing created
2026-07-22 18:05:11 INFO Hourly dataset computed
2026-07-22 18:05:11 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 18:05:12 INFO Hourly dataset computed and listing created
2026-07-22 18:05:13 INFO Hourly dataset computed
2026-07-22 18:05:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 18:05:14 INFO Hourly dataset computed and listing created
2026-07-22 18:05:14 INFO Hourly dataset computed
2026-07-22 18:05:14 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 18:05:15 INFO Hourly dataset computed and listing created
2026-07-22 18:05:16 INFO Hourly dataset computed
2026-07-22 18:05:16 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 18:05:17 INFO Hourly dataset computed and listing created
2026-07-22 18:05:17 INFO Hourly dataset computed
2026-07-22 18:05:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 18:05:18 INFO Hourly dataset computed and listing created
2026-07-22 18:05:19 INFO Hourly dataset computed
2026-07-22 18:05:19 INFO ---------->>> Running CHIMERE model from 2020-02-07 11:00:00 to 2020-02-07 12:00:00
2026-07-22 18:05:19 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 18:05:19 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1
2026-07-22 18:05:19 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020709_2_ENS1.nc
2026-07-22 18:05:19 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-22 18:05:19 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 18:05:19 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-22 18:05:20 INFO Queuing job for member 1...
2026-07-22 18:05:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 18:05:20 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-22 18:05:20 INFO Found: ['5257410']
2026-07-22 18:05:25 INFO [TGCC-IRENE] Submitted job with ID:['5257410']
2026-07-22 18:05:25 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 18:05:25 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2
2026-07-22 18:05:25 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020709_2_ENS2.nc
2026-07-22 18:05:25 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-22 18:05:25 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 18:05:25 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-22 18:05:25 INFO Queuing job for member 2...
2026-07-22 18:05:25 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 18:05:25 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-22 18:05:26 INFO Found: ['5257411']
2026-07-22 18:05:31 INFO [TGCC-IRENE] Submitted job with ID:['5257411']
2026-07-22 18:05:31 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 18:05:31 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3
2026-07-22 18:05:31 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020709_2_ENS3.nc
2026-07-22 18:05:31 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-22 18:05:31 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 18:05:31 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-22 18:05:31 INFO Queuing job for member 3...
2026-07-22 18:05:31 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 18:05:31 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-22 18:05:32 INFO Found: ['5257412']
2026-07-22 18:05:37 INFO [TGCC-IRENE] Submitted job with ID:['5257412']
2026-07-22 18:05:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 18:05:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4
2026-07-22 18:05:37 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020709_2_ENS4.nc
2026-07-22 18:05:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-22 18:05:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 18:05:37 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-22 18:05:37 INFO Queuing job for member 4...
2026-07-22 18:05:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 18:05:37 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-22 18:05:37 INFO Found: ['5257413']
2026-07-22 18:05:42 INFO [TGCC-IRENE] Submitted job with ID:['5257413']
2026-07-22 18:05:42 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 18:05:42 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5
2026-07-22 18:05:42 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020709_2_ENS5.nc
2026-07-22 18:05:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-22 18:05:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 18:05:43 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-22 18:05:43 INFO Queuing job for member 5...
2026-07-22 18:05:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 18:05:43 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-22 18:05:43 INFO Found: ['5257414']
2026-07-22 18:05:48 INFO [TGCC-IRENE] Submitted job with ID:['5257414']
2026-07-22 18:05:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 18:05:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6
2026-07-22 18:05:48 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020709_2_ENS6.nc
2026-07-22 18:05:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-22 18:05:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 18:05:48 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-22 18:05:48 INFO Queuing job for member 6...
2026-07-22 18:05:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 18:05:48 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-22 18:05:50 INFO Found: ['5257415']
2026-07-22 18:05:55 INFO [TGCC-IRENE] Submitted job with ID:['5257415']
2026-07-22 18:05:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 18:05:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7
2026-07-22 18:05:55 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020709_2_ENS7.nc
2026-07-22 18:05:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-22 18:05:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 18:05:55 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-22 18:05:55 INFO Queuing job for member 7...
2026-07-22 18:05:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 18:05:55 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-22 18:05:57 INFO Found: ['5257416']
2026-07-22 18:06:02 INFO [TGCC-IRENE] Submitted job with ID:['5257416']
2026-07-22 18:06:02 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 18:06:02 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8
2026-07-22 18:06:02 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020709_2_ENS8.nc
2026-07-22 18:06:03 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-22 18:06:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 18:06:03 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-22 18:06:03 INFO Queuing job for member 8...
2026-07-22 18:06:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 18:06:03 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-22 18:06:05 INFO Found: ['5257419']
2026-07-22 18:06:10 INFO [TGCC-IRENE] Submitted job with ID:['5257419']
2026-07-22 18:06:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 18:06:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9
2026-07-22 18:06:10 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020709_2_ENS9.nc
2026-07-22 18:06:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-22 18:06:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 18:06:10 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-22 18:06:10 INFO Queuing job for member 9...
2026-07-22 18:06:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 18:06:10 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-22 18:06:12 INFO Found: ['5257420']
2026-07-22 18:06:17 INFO [TGCC-IRENE] Submitted job with ID:['5257420']
2026-07-22 18:06:17 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 18:06:17 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10
2026-07-22 18:06:17 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020709_2_ENS10.nc
2026-07-22 18:06:17 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-22 18:06:17 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 18:06:17 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-22 18:06:18 INFO Queuing job for member 10...
2026-07-22 18:06:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 18:06:18 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-22 18:06:20 INFO Found: ['5257421']
2026-07-22 18:06:25 INFO [TGCC-IRENE] Submitted job with ID:['5257421']
2026-07-22 18:06:25 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 18:06:25 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11
2026-07-22 18:06:25 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020709_2_ENS11.nc
2026-07-22 18:06:25 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-22 18:06:25 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 18:06:25 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-22 18:06:25 INFO Queuing job for member 11...
2026-07-22 18:06:25 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 18:06:25 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-22 18:06:28 INFO Found: ['5257423']
2026-07-22 18:06:33 INFO [TGCC-IRENE] Submitted job with ID:['5257423']
2026-07-22 18:06:33 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 18:06:33 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12
2026-07-22 18:06:33 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020709_2_ENS12.nc
2026-07-22 18:06:33 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-22 18:06:33 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 18:06:33 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-22 18:06:33 INFO Queuing job for member 12...
2026-07-22 18:06:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 18:06:33 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-22 18:06:35 INFO Found: ['5257424']
2026-07-22 18:06:40 INFO [TGCC-IRENE] Submitted job with ID:['5257424']
2026-07-22 18:06:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 18:06:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13
2026-07-22 18:06:40 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020709_2_ENS13.nc
2026-07-22 18:06:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-22 18:06:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 18:06:40 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-22 18:06:40 INFO Queuing job for member 13...
2026-07-22 18:06:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 18:06:40 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-22 18:06:41 INFO Found: ['5257426']
2026-07-22 18:06:46 INFO [TGCC-IRENE] Submitted job with ID:['5257426']
2026-07-22 18:06:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 18:06:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14
2026-07-22 18:06:46 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020709_2_ENS14.nc
2026-07-22 18:06:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-22 18:06:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 18:06:46 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-22 18:06:46 INFO Queuing job for member 14...
2026-07-22 18:06:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 18:06:46 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-22 18:06:46 INFO Found: ['5257427']
2026-07-22 18:06:51 INFO [TGCC-IRENE] Submitted job with ID:['5257427']
2026-07-22 18:06:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 18:06:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15
2026-07-22 18:06:51 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020709_2_ENS15.nc
2026-07-22 18:06:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-22 18:06:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 18:06:51 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-22 18:06:51 INFO Queuing job for member 15...
2026-07-22 18:06:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 18:06:51 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-22 18:06:52 INFO Found: ['5257428']
2026-07-22 18:06:57 INFO [TGCC-IRENE] Submitted job with ID:['5257428']
2026-07-22 18:06:57 INFO Checking job status ...
2026-07-22 18:06:57 INFO None 5257410: status RUNNING/PENDING
2026-07-22 18:06:57 INFO None 5257411: status RUNNING/PENDING
2026-07-22 18:06:57 INFO None 5257412: status RUNNING/PENDING
2026-07-22 18:06:57 INFO None 5257413: status RUNNING/PENDING
2026-07-22 18:06:57 INFO None 5257414: status RUNNING/PENDING
2026-07-22 18:06:57 INFO None 5257415: status RUNNING/PENDING
2026-07-22 18:06:57 INFO None 5257416: status RUNNING/PENDING
2026-07-22 18:06:57 INFO None 5257419: status RUNNING/PENDING
2026-07-22 18:06:57 INFO None 5257420: status RUNNING/PENDING
2026-07-22 18:06:57 INFO None 5257421: status RUNNING/PENDING
2026-07-22 18:06:57 INFO None 5257423: status RUNNING/PENDING
2026-07-22 18:06:57 INFO None 5257424: status RUNNING/PENDING
2026-07-22 18:06:57 INFO None 5257426: status RUNNING/PENDING
2026-07-22 18:06:57 INFO None 5257427: status RUNNING/PENDING
2026-07-22 18:06:57 INFO None 5257428: status RUNNING/PENDING
2026-07-22 18:06:57 INFO Jobs still running: ['5257410', '5257411', '5257412', '5257413', '5257414', '5257415', '5257416', '5257419', '5257420', '5257421', '5257423', '5257424', '5257426', '5257427', '5257428']. Waiting...
2026-07-22 18:07:12 INFO None 5257410: status RUNNING/PENDING
2026-07-22 18:07:12 INFO None 5257411: status RUNNING/PENDING
2026-07-22 18:07:12 INFO None 5257412: status RUNNING/PENDING
2026-07-22 18:07:12 INFO None 5257413: status RUNNING/PENDING
2026-07-22 18:07:12 INFO None 5257414: status RUNNING/PENDING
2026-07-22 18:07:12 INFO None 5257415: status RUNNING/PENDING
2026-07-22 18:07:12 INFO None 5257416: status RUNNING/PENDING
2026-07-22 18:07:12 INFO None 5257419: status RUNNING/PENDING
2026-07-22 18:07:13 INFO None 5257420: status RUNNING/PENDING
2026-07-22 18:07:13 INFO None 5257421: status RUNNING/PENDING
2026-07-22 18:07:13 INFO None 5257423: status RUNNING/PENDING
2026-07-22 18:07:13 INFO None 5257424: status RUNNING/PENDING
2026-07-22 18:07:13 INFO None 5257426: status RUNNING/PENDING
2026-07-22 18:07:13 INFO None 5257427: status RUNNING/PENDING
2026-07-22 18:07:13 INFO None 5257428: status RUNNING/PENDING
2026-07-22 18:07:13 INFO Jobs still running: ['5257410', '5257411', '5257412', '5257413', '5257414', '5257415', '5257416', '5257419', '5257420', '5257421', '5257423', '5257424', '5257426', '5257427', '5257428']. Waiting...
2026-07-22 18:07:29 INFO None 5257410: status RUNNING/PENDING
2026-07-22 18:07:29 INFO None 5257411: status RUNNING/PENDING
2026-07-22 18:07:29 INFO None 5257412: status RUNNING/PENDING
2026-07-22 18:07:29 INFO None 5257413: status RUNNING/PENDING
2026-07-22 18:07:29 INFO None 5257414: status RUNNING/PENDING
2026-07-22 18:07:29 INFO None 5257415: status RUNNING/PENDING
2026-07-22 18:07:29 INFO None 5257416: status RUNNING/PENDING
2026-07-22 18:07:29 INFO None 5257419: status RUNNING/PENDING
2026-07-22 18:07:29 INFO None 5257420: status RUNNING/PENDING
2026-07-22 18:07:29 INFO None 5257421: status RUNNING/PENDING
2026-07-22 18:07:29 INFO None 5257423: status RUNNING/PENDING
2026-07-22 18:07:29 INFO None 5257424: status RUNNING/PENDING
2026-07-22 18:07:29 INFO None 5257426: status RUNNING/PENDING
2026-07-22 18:07:29 INFO None 5257427: status RUNNING/PENDING
2026-07-22 18:07:29 INFO None 5257428: status RUNNING/PENDING
2026-07-22 18:07:29 INFO Jobs still running: ['5257410', '5257411', '5257412', '5257413', '5257414', '5257415', '5257416', '5257419', '5257420', '5257421', '5257423', '5257424', '5257426', '5257427', '5257428']. Waiting...
2026-07-22 18:07:44 INFO None 5257410: status RUNNING/PENDING
2026-07-22 18:07:44 INFO None 5257411: status RUNNING/PENDING
2026-07-22 18:07:44 INFO None 5257412: status RUNNING/PENDING
2026-07-22 18:07:45 INFO None 5257413: status RUNNING/PENDING
2026-07-22 18:07:45 INFO None 5257414: status RUNNING/PENDING
2026-07-22 18:07:45 INFO None 5257415: status RUNNING/PENDING
2026-07-22 18:07:45 INFO None 5257416: status RUNNING/PENDING
2026-07-22 18:07:45 INFO None 5257419: status RUNNING/PENDING
2026-07-22 18:07:45 INFO None 5257420: status RUNNING/PENDING
2026-07-22 18:07:45 INFO None 5257421: status RUNNING/PENDING
2026-07-22 18:07:45 INFO None 5257423: status RUNNING/PENDING
2026-07-22 18:07:45 INFO None 5257424: status RUNNING/PENDING
2026-07-22 18:07:45 INFO None 5257426: status RUNNING/PENDING
2026-07-22 18:07:47 INFO None 5257427: status RUNNING/PENDING
2026-07-22 18:07:47 INFO None 5257428: status RUNNING/PENDING
2026-07-22 18:07:47 INFO Jobs still running: ['5257410', '5257411', '5257412', '5257413', '5257414', '5257415', '5257416', '5257419', '5257420', '5257421', '5257423', '5257424', '5257426', '5257427', '5257428']. Waiting...
2026-07-22 18:08:02 INFO None 5257410: status RUNNING/PENDING
2026-07-22 18:08:02 INFO None 5257411: status RUNNING/PENDING
2026-07-22 18:08:02 INFO None 5257412: status RUNNING/PENDING
2026-07-22 18:08:02 INFO None 5257413: status RUNNING/PENDING
2026-07-22 18:08:02 INFO None 5257414: status FINISHED
2026-07-22 18:08:02 INFO None 5257415: status RUNNING/PENDING
2026-07-22 18:08:02 INFO None 5257416: status RUNNING/PENDING
2026-07-22 18:08:02 INFO None 5257419: status RUNNING/PENDING
2026-07-22 18:08:02 INFO None 5257420: status RUNNING/PENDING
2026-07-22 18:08:02 INFO None 5257421: status RUNNING/PENDING
2026-07-22 18:08:02 INFO None 5257423: status RUNNING/PENDING
2026-07-22 18:08:02 INFO None 5257424: status RUNNING/PENDING
2026-07-22 18:08:02 INFO None 5257426: status RUNNING/PENDING
2026-07-22 18:08:02 INFO None 5257427: status RUNNING/PENDING
2026-07-22 18:08:02 INFO None 5257428: status RUNNING/PENDING
2026-07-22 18:08:02 INFO Jobs still running: ['5257410', '5257411', '5257412', '5257413', '5257415', '5257416', '5257419', '5257420', '5257421', '5257423', '5257424', '5257426', '5257427', '5257428']. Waiting...
2026-07-22 18:08:17 INFO None 5257410: status RUNNING/PENDING
2026-07-22 18:08:17 INFO None 5257411: status RUNNING/PENDING
2026-07-22 18:08:17 INFO None 5257412: status RUNNING/PENDING
2026-07-22 18:08:17 INFO None 5257413: status RUNNING/PENDING
2026-07-22 18:08:17 INFO None 5257414: status FINISHED
2026-07-22 18:08:17 INFO None 5257415: status RUNNING/PENDING
2026-07-22 18:08:17 INFO None 5257416: status RUNNING/PENDING
2026-07-22 18:08:17 INFO None 5257419: status RUNNING/PENDING
2026-07-22 18:08:17 INFO None 5257420: status RUNNING/PENDING
2026-07-22 18:08:17 INFO None 5257421: status RUNNING/PENDING
2026-07-22 18:08:17 INFO None 5257423: status RUNNING/PENDING
2026-07-22 18:08:17 INFO None 5257424: status RUNNING/PENDING
2026-07-22 18:08:17 INFO None 5257426: status RUNNING/PENDING
2026-07-22 18:08:17 INFO None 5257427: status RUNNING/PENDING
2026-07-22 18:08:17 INFO None 5257428: status RUNNING/PENDING
2026-07-22 18:08:17 INFO Jobs still running: ['5257410', '5257411', '5257412', '5257413', '5257415', '5257416', '5257419', '5257420', '5257421', '5257423', '5257424', '5257426', '5257427', '5257428']. Waiting...
2026-07-22 18:08:32 INFO None 5257410: status RUNNING/PENDING
2026-07-22 18:08:32 INFO None 5257411: status RUNNING/PENDING
2026-07-22 18:08:32 INFO None 5257412: status RUNNING/PENDING
2026-07-22 18:08:33 INFO None 5257413: status FINISHED
2026-07-22 18:08:33 INFO None 5257414: status FINISHED
2026-07-22 18:08:33 INFO None 5257415: status RUNNING/PENDING
2026-07-22 18:08:33 INFO None 5257416: status RUNNING/PENDING
2026-07-22 18:08:33 INFO None 5257419: status RUNNING/PENDING
2026-07-22 18:08:33 INFO None 5257420: status RUNNING/PENDING
2026-07-22 18:08:33 INFO None 5257421: status RUNNING/PENDING
2026-07-22 18:08:33 INFO None 5257423: status RUNNING/PENDING
2026-07-22 18:08:33 INFO None 5257424: status RUNNING/PENDING
2026-07-22 18:08:33 INFO None 5257426: status RUNNING/PENDING
2026-07-22 18:08:33 INFO None 5257427: status RUNNING/PENDING
2026-07-22 18:08:33 INFO None 5257428: status RUNNING/PENDING
2026-07-22 18:08:33 INFO Jobs still running: ['5257410', '5257411', '5257412', '5257415', '5257416', '5257419', '5257420', '5257421', '5257423', '5257424', '5257426', '5257427', '5257428']. Waiting...
2026-07-22 18:08:48 INFO None 5257410: status RUNNING/PENDING
2026-07-22 18:08:48 INFO None 5257411: status RUNNING/PENDING
2026-07-22 18:08:48 INFO None 5257412: status RUNNING/PENDING
2026-07-22 18:08:48 INFO None 5257413: status FINISHED
2026-07-22 18:08:48 INFO None 5257414: status FINISHED
2026-07-22 18:08:48 INFO None 5257415: status RUNNING/PENDING
2026-07-22 18:08:48 INFO None 5257416: status RUNNING/PENDING
2026-07-22 18:08:48 INFO None 5257419: status RUNNING/PENDING
2026-07-22 18:08:48 INFO None 5257420: status RUNNING/PENDING
2026-07-22 18:08:48 INFO None 5257421: status RUNNING/PENDING
2026-07-22 18:08:48 INFO None 5257423: status RUNNING/PENDING
2026-07-22 18:08:48 INFO None 5257424: status RUNNING/PENDING
2026-07-22 18:08:48 INFO None 5257426: status RUNNING/PENDING
2026-07-22 18:08:48 INFO None 5257427: status RUNNING/PENDING
2026-07-22 18:08:48 INFO None 5257428: status RUNNING/PENDING
2026-07-22 18:08:48 INFO Jobs still running: ['5257410', '5257411', '5257412', '5257415', '5257416', '5257419', '5257420', '5257421', '5257423', '5257424', '5257426', '5257427', '5257428']. Waiting...
2026-07-22 18:09:05 INFO None 5257410: status FINISHED
2026-07-22 18:09:05 INFO None 5257411: status FINISHED
2026-07-22 18:09:05 INFO None 5257412: status FINISHED
2026-07-22 18:09:05 INFO None 5257413: status FINISHED
2026-07-22 18:09:05 INFO None 5257414: status FINISHED
2026-07-22 18:09:05 INFO None 5257415: status RUNNING/PENDING
2026-07-22 18:09:05 INFO None 5257416: status RUNNING/PENDING
2026-07-22 18:09:05 INFO None 5257419: status RUNNING/PENDING
2026-07-22 18:09:05 INFO None 5257420: status RUNNING/PENDING
2026-07-22 18:09:05 INFO None 5257421: status RUNNING/PENDING
2026-07-22 18:09:05 INFO None 5257423: status RUNNING/PENDING
2026-07-22 18:09:05 INFO None 5257424: status RUNNING/PENDING
2026-07-22 18:09:05 INFO None 5257426: status RUNNING/PENDING
2026-07-22 18:09:05 INFO None 5257427: status RUNNING/PENDING
2026-07-22 18:09:05 INFO None 5257428: status RUNNING/PENDING
2026-07-22 18:09:05 INFO Jobs still running: ['5257415', '5257416', '5257419', '5257420', '5257421', '5257423', '5257424', '5257426', '5257427', '5257428']. Waiting...
2026-07-22 18:09:20 INFO None 5257410: status FINISHED
2026-07-22 18:09:20 INFO None 5257411: status FINISHED
2026-07-22 18:09:20 INFO None 5257412: status FINISHED
2026-07-22 18:09:20 INFO None 5257413: status FINISHED
2026-07-22 18:09:20 INFO None 5257414: status FINISHED
2026-07-22 18:09:20 INFO None 5257415: status RUNNING/PENDING
2026-07-22 18:09:20 INFO None 5257416: status RUNNING/PENDING
2026-07-22 18:09:20 INFO None 5257419: status RUNNING/PENDING
2026-07-22 18:09:20 INFO None 5257420: status RUNNING/PENDING
2026-07-22 18:09:20 INFO None 5257421: status RUNNING/PENDING
2026-07-22 18:09:20 INFO None 5257423: status RUNNING/PENDING
2026-07-22 18:09:20 INFO None 5257424: status RUNNING/PENDING
2026-07-22 18:09:20 INFO None 5257426: status RUNNING/PENDING
2026-07-22 18:09:20 INFO None 5257427: status RUNNING/PENDING
2026-07-22 18:09:20 INFO None 5257428: status RUNNING/PENDING
2026-07-22 18:09:20 INFO Jobs still running: ['5257415', '5257416', '5257419', '5257420', '5257421', '5257423', '5257424', '5257426', '5257427', '5257428']. Waiting...
2026-07-22 18:09:35 INFO None 5257410: status FINISHED
2026-07-22 18:09:35 INFO None 5257411: status FINISHED
2026-07-22 18:09:35 INFO None 5257412: status FINISHED
2026-07-22 18:09:35 INFO None 5257413: status FINISHED
2026-07-22 18:09:37 INFO None 5257414: status FINISHED
2026-07-22 18:09:38 INFO None 5257415: status RUNNING/PENDING
2026-07-22 18:09:38 INFO None 5257416: status RUNNING/PENDING
2026-07-22 18:09:38 INFO None 5257419: status RUNNING/PENDING
2026-07-22 18:09:38 INFO None 5257420: status RUNNING/PENDING
2026-07-22 18:09:38 INFO None 5257421: status RUNNING/PENDING
2026-07-22 18:09:38 INFO None 5257423: status RUNNING/PENDING
2026-07-22 18:09:38 INFO None 5257424: status RUNNING/PENDING
2026-07-22 18:09:38 INFO None 5257426: status RUNNING/PENDING
2026-07-22 18:09:38 INFO None 5257427: status RUNNING/PENDING
2026-07-22 18:09:38 INFO None 5257428: status RUNNING/PENDING
2026-07-22 18:09:38 INFO Jobs still running: ['5257415', '5257416', '5257419', '5257420', '5257421', '5257423', '5257424', '5257426', '5257427', '5257428']. Waiting...
2026-07-22 18:09:53 INFO None 5257410: status FINISHED
2026-07-22 18:09:53 INFO None 5257411: status FINISHED
2026-07-22 18:09:53 INFO None 5257412: status FINISHED
2026-07-22 18:09:53 INFO None 5257413: status FINISHED
2026-07-22 18:09:53 INFO None 5257414: status FINISHED
2026-07-22 18:09:53 INFO None 5257415: status RUNNING/PENDING
2026-07-22 18:09:53 INFO None 5257416: status RUNNING/PENDING
2026-07-22 18:09:53 INFO None 5257419: status RUNNING/PENDING
2026-07-22 18:09:53 INFO None 5257420: status RUNNING/PENDING
2026-07-22 18:09:53 INFO None 5257421: status RUNNING/PENDING
2026-07-22 18:09:53 INFO None 5257423: status RUNNING/PENDING
2026-07-22 18:09:53 INFO None 5257424: status RUNNING/PENDING
2026-07-22 18:09:53 INFO None 5257426: status RUNNING/PENDING
2026-07-22 18:09:53 INFO None 5257427: status RUNNING/PENDING
2026-07-22 18:09:53 INFO None 5257428: status RUNNING/PENDING
2026-07-22 18:09:53 INFO Jobs still running: ['5257415', '5257416', '5257419', '5257420', '5257421', '5257423', '5257424', '5257426', '5257427', '5257428']. Waiting...
2026-07-22 18:10:08 INFO None 5257410: status FINISHED
2026-07-22 18:10:08 INFO None 5257411: status FINISHED
2026-07-22 18:10:08 INFO None 5257412: status FINISHED
2026-07-22 18:10:08 INFO None 5257413: status FINISHED
2026-07-22 18:10:08 INFO None 5257414: status FINISHED
2026-07-22 18:10:08 INFO None 5257415: status RUNNING/PENDING
2026-07-22 18:10:08 INFO None 5257416: status FINISHED
2026-07-22 18:10:08 INFO None 5257419: status RUNNING/PENDING
2026-07-22 18:10:08 INFO None 5257420: status RUNNING/PENDING
2026-07-22 18:10:08 INFO None 5257421: status RUNNING/PENDING
2026-07-22 18:10:08 INFO None 5257423: status RUNNING/PENDING
2026-07-22 18:10:08 INFO None 5257424: status RUNNING/PENDING
2026-07-22 18:10:08 INFO None 5257426: status RUNNING/PENDING
2026-07-22 18:10:08 INFO None 5257427: status RUNNING/PENDING
2026-07-22 18:10:08 INFO None 5257428: status RUNNING/PENDING
2026-07-22 18:10:08 INFO Jobs still running: ['5257415', '5257419', '5257420', '5257421', '5257423', '5257424', '5257426', '5257427', '5257428']. Waiting...
2026-07-22 18:10:23 INFO None 5257410: status FINISHED
2026-07-22 18:10:23 INFO None 5257411: status FINISHED
2026-07-22 18:10:23 INFO None 5257412: status FINISHED
2026-07-22 18:10:23 INFO None 5257413: status FINISHED
2026-07-22 18:10:23 INFO None 5257414: status FINISHED
2026-07-22 18:10:23 INFO None 5257415: status RUNNING/PENDING
2026-07-22 18:10:23 INFO None 5257416: status FINISHED
2026-07-22 18:10:23 INFO None 5257419: status RUNNING/PENDING
2026-07-22 18:10:23 INFO None 5257420: status RUNNING/PENDING
2026-07-22 18:10:23 INFO None 5257421: status RUNNING/PENDING
2026-07-22 18:10:23 INFO None 5257423: status RUNNING/PENDING
2026-07-22 18:10:23 INFO None 5257424: status RUNNING/PENDING
2026-07-22 18:10:24 INFO None 5257426: status RUNNING/PENDING
2026-07-22 18:10:24 INFO None 5257427: status RUNNING/PENDING
2026-07-22 18:10:24 INFO None 5257428: status RUNNING/PENDING
2026-07-22 18:10:24 INFO Jobs still running: ['5257415', '5257419', '5257420', '5257421', '5257423', '5257424', '5257426', '5257427', '5257428']. Waiting...
2026-07-22 18:10:40 INFO None 5257410: status FINISHED
2026-07-22 18:10:40 INFO None 5257411: status FINISHED
2026-07-22 18:10:40 INFO None 5257412: status FINISHED
2026-07-22 18:10:41 INFO None 5257413: status FINISHED
2026-07-22 18:10:41 INFO None 5257414: status FINISHED
2026-07-22 18:10:41 INFO None 5257415: status FINISHED
2026-07-22 18:10:41 INFO None 5257416: status FINISHED
2026-07-22 18:10:41 INFO None 5257419: status RUNNING/PENDING
2026-07-22 18:10:41 INFO None 5257420: status RUNNING/PENDING
2026-07-22 18:10:41 INFO None 5257421: status RUNNING/PENDING
2026-07-22 18:10:41 INFO None 5257423: status RUNNING/PENDING
2026-07-22 18:10:41 INFO None 5257424: status RUNNING/PENDING
2026-07-22 18:10:41 INFO None 5257426: status RUNNING/PENDING
2026-07-22 18:10:41 INFO None 5257427: status RUNNING/PENDING
2026-07-22 18:10:41 INFO None 5257428: status RUNNING/PENDING
2026-07-22 18:10:41 INFO Jobs still running: ['5257419', '5257420', '5257421', '5257423', '5257424', '5257426', '5257427', '5257428']. Waiting...
2026-07-22 18:10:56 INFO None 5257410: status FINISHED
2026-07-22 18:10:56 INFO None 5257411: status FINISHED
2026-07-22 18:10:56 INFO None 5257412: status FINISHED
2026-07-22 18:10:56 INFO None 5257413: status FINISHED
2026-07-22 18:10:56 INFO None 5257414: status FINISHED
2026-07-22 18:10:56 INFO None 5257415: status FINISHED
2026-07-22 18:10:56 INFO None 5257416: status FINISHED
2026-07-22 18:10:56 INFO None 5257419: status FINISHED
2026-07-22 18:10:56 INFO None 5257420: status FINISHED
2026-07-22 18:10:56 INFO None 5257421: status FINISHED
2026-07-22 18:10:58 INFO None 5257423: status FINISHED
2026-07-22 18:10:58 INFO None 5257424: status RUNNING/PENDING
2026-07-22 18:10:58 INFO None 5257426: status FINISHED
2026-07-22 18:10:58 INFO None 5257427: status FINISHED
2026-07-22 18:10:58 INFO None 5257428: status RUNNING/PENDING
2026-07-22 18:10:58 INFO Jobs still running: ['5257424', '5257428']. Waiting...
2026-07-22 18:11:13 INFO None 5257410: status FINISHED
2026-07-22 18:11:13 INFO None 5257411: status FINISHED
2026-07-22 18:11:13 INFO None 5257412: status FINISHED
2026-07-22 18:11:13 INFO None 5257413: status FINISHED
2026-07-22 18:11:13 INFO None 5257414: status FINISHED
2026-07-22 18:11:13 INFO None 5257415: status FINISHED
2026-07-22 18:11:13 INFO None 5257416: status FINISHED
2026-07-22 18:11:13 INFO None 5257419: status FINISHED
2026-07-22 18:11:13 INFO None 5257420: status FINISHED
2026-07-22 18:11:13 INFO None 5257421: status FINISHED
2026-07-22 18:11:13 INFO None 5257423: status FINISHED
2026-07-22 18:11:13 INFO None 5257424: status FINISHED
2026-07-22 18:11:13 INFO None 5257426: status FINISHED
2026-07-22 18:11:13 INFO None 5257427: status FINISHED
2026-07-22 18:11:13 INFO None 5257428: status FINISHED
2026-07-22 18:11:13 INFO Jobs ['5257410', '5257411', '5257412', '5257413', '5257414', '5257415', '5257416', '5257419', '5257420', '5257421', '5257423', '5257424', '5257426', '5257427', '5257428'] have finished
2026-07-22 18:11:13 INFO Checking restart files were created ...
2026-07-22 18:11:13 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020711_1_ENS1.nc(668832435 bytes)
2026-07-22 18:11:13 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020711_1_ENS2.nc(668832435 bytes)
2026-07-22 18:11:13 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020711_1_ENS3.nc(668832435 bytes)
2026-07-22 18:11:13 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020711_1_ENS4.nc(668832435 bytes)
2026-07-22 18:11:13 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020711_1_ENS5.nc(668832435 bytes)
2026-07-22 18:11:13 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020711_1_ENS6.nc(668832435 bytes)
2026-07-22 18:11:13 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020711_1_ENS7.nc(668832435 bytes)
2026-07-22 18:11:13 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020711_1_ENS8.nc(668832435 bytes)
2026-07-22 18:11:13 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/end.2020020711_1_ENS9.nc(668832435 bytes)
2026-07-22 18:11:13 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/end.2020020711_1_ENS10.nc(668832435 bytes)
2026-07-22 18:11:13 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/end.2020020711_1_ENS11.nc(668832435 bytes)
2026-07-22 18:11:13 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/end.2020020711_1_ENS12.nc(668832435 bytes)
2026-07-22 18:11:13 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/end.2020020711_1_ENS13.nc(668832435 bytes)
2026-07-22 18:11:13 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/end.2020020711_1_ENS14.nc(668832435 bytes)
2026-07-22 18:11:13 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/end.2020020711_1_ENS15.nc(668832435 bytes)
2026-07-22 18:11:13 INFO  Run_model() completed successfully.
2026-07-22 18:11:13 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 11:00:00 simulated_time=2020-02-07 12:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 18:11:13 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 12:00:00 days=153073 seconds=43200
2026-07-22 18:11:13 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-22 18:11:13 INFO [TIME] increment current_time 2020-02-07 11:00:00 -> 2020-02-07 12:00:00
2026-07-22 18:11:13 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 12:00:00 simulated_time=2020-02-07 12:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 18:11:13 INFO ---------->>> Running process_satellite_data()
2026-07-22 18:11:13 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12018.nc
2026-07-22 18:11:13 INFO ---------->>> Running run_obs_converter()
2026-07-22 18:11:13 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_43993_153073.out
2026-07-22 18:11:13 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_43993_153073.out
2026-07-22 18:11:13 INFO ---------->>> Running DART
2026-07-22 18:11:13 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-22 18:11:13 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/chim_ENS1_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/chim_ENS1_2020020712_1_out_toDART.nc
2026-07-22 18:11:14 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/chim_ENS2_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/chim_ENS2_2020020712_1_out_toDART.nc
2026-07-22 18:11:14 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/chim_ENS3_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/chim_ENS3_2020020712_1_out_toDART.nc
2026-07-22 18:11:15 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/chim_ENS4_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/chim_ENS4_2020020712_1_out_toDART.nc
2026-07-22 18:11:15 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/chim_ENS5_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/chim_ENS5_2020020712_1_out_toDART.nc
2026-07-22 18:11:15 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/chim_ENS6_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/chim_ENS6_2020020712_1_out_toDART.nc
2026-07-22 18:11:16 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/chim_ENS7_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/chim_ENS7_2020020712_1_out_toDART.nc
2026-07-22 18:11:16 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/chim_ENS8_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/chim_ENS8_2020020712_1_out_toDART.nc
2026-07-22 18:11:16 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/chim_ENS9_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/chim_ENS9_2020020712_1_out_toDART.nc
2026-07-22 18:11:17 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/chim_ENS10_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/chim_ENS10_2020020712_1_out_toDART.nc
2026-07-22 18:11:17 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/chim_ENS11_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/chim_ENS11_2020020712_1_out_toDART.nc
2026-07-22 18:11:17 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/chim_ENS12_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/chim_ENS12_2020020712_1_out_toDART.nc
2026-07-22 18:11:18 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/chim_ENS13_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/chim_ENS13_2020020712_1_out_toDART.nc
2026-07-22 18:11:18 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/chim_ENS14_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/chim_ENS14_2020020712_1_out_toDART.nc
2026-07-22 18:11:18 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/chim_ENS15_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/chim_ENS15_2020020712_1_out_toDART.nc
2026-07-22 18:11:19 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-22 18:11:19 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-22 18:11:19 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-22 18:11:19 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-22 18:11:19 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-22 18:11:19 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-22 18:11:33 INFO Found: []
2026-07-22 18:11:33 INFO No job id returned by command ./run_filter.bsh
2026-07-22 18:11:33 INFO No monitoring will be performed
2026-07-22 18:11:33 INFO Moving DART output files to analysis and preassim directories for date 2020020712 if present ...
2026-07-22 18:11:33 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020712'
2026-07-22 18:11:33 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020712'
2026-07-22 18:11:33 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020712'
2026-07-22 18:11:33 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020712'
2026-07-22 18:11:33 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020712'
2026-07-22 18:11:33 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020712'
2026-07-22 18:11:33 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020712'
2026-07-22 18:11:33 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020712'
2026-07-22 18:11:33 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020712'
2026-07-22 18:11:33 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020712'
2026-07-22 18:11:33 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020712'
2026-07-22 18:11:33 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020712'
2026-07-22 18:11:33 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020712'
2026-07-22 18:11:33 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020712'
2026-07-22 18:11:33 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020712'
2026-07-22 18:11:33 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020712'
2026-07-22 18:11:33 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020712'
2026-07-22 18:11:33 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020712'
2026-07-22 18:11:34 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020712'
2026-07-22 18:11:34 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020712'
2026-07-22 18:11:34 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020712'
2026-07-22 18:11:34 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020712'
2026-07-22 18:11:34 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020712'
2026-07-22 18:11:34 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020712'
2026-07-22 18:11:34 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020712'
2026-07-22 18:11:34 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020712'
2026-07-22 18:11:34 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020712'
2026-07-22 18:11:34 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020712'
2026-07-22 18:11:34 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020712'
2026-07-22 18:11:34 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020712'
2026-07-22 18:11:34 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020712'
2026-07-22 18:11:34 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/analysis/2020020712'
2026-07-22 18:11:34 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020712'
2026-07-22 18:11:34 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/preassim/2020020712'
2026-07-22 18:11:34 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-22 18:11:34 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-22 18:11:34 INFO run_dart() is DONE.
2026-07-22 18:11:34 INFO ---------->>> Running update_pollutant_in_end()
2026-07-22 18:11:34 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:11:34 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS1.nc
2026-07-22 18:11:35 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:11:36 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS1.nc
2026-07-22 18:11:36 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:11:36 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:11:36 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS1.nc
2026-07-22 18:11:37 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:11:38 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS1.nc
2026-07-22 18:11:38 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:11:38 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 18:11:38 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:11:38 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS2.nc
2026-07-22 18:11:39 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:11:41 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS2.nc
2026-07-22 18:11:41 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:11:41 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:11:41 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS2.nc
2026-07-22 18:11:42 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:11:42 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS2.nc
2026-07-22 18:11:42 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:11:42 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 18:11:43 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:11:43 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS3.nc
2026-07-22 18:11:44 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:11:45 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS3.nc
2026-07-22 18:11:45 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:11:45 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:11:45 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS3.nc
2026-07-22 18:11:46 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:11:46 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS3.nc
2026-07-22 18:11:46 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:11:46 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 18:11:47 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:11:47 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS4.nc
2026-07-22 18:11:48 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:11:49 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS4.nc
2026-07-22 18:11:49 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:11:49 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:11:49 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS4.nc
2026-07-22 18:11:50 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:11:51 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS4.nc
2026-07-22 18:11:51 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:11:51 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 18:11:51 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:11:51 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS5.nc
2026-07-22 18:11:52 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:11:53 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS5.nc
2026-07-22 18:11:53 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:11:53 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:11:53 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS5.nc
2026-07-22 18:11:53 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:11:54 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS5.nc
2026-07-22 18:11:54 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:11:54 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 18:11:55 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:11:55 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS6.nc
2026-07-22 18:11:55 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:11:56 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS6.nc
2026-07-22 18:11:56 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:11:56 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:11:56 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS6.nc
2026-07-22 18:11:57 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:11:57 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS6.nc
2026-07-22 18:11:57 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:11:57 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 18:11:58 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:11:58 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS7.nc
2026-07-22 18:11:58 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:11:59 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS7.nc
2026-07-22 18:11:59 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:12:00 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:12:00 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS7.nc
2026-07-22 18:12:00 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:12:01 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS7.nc
2026-07-22 18:12:01 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:12:01 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 18:12:02 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:12:02 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS8.nc
2026-07-22 18:12:03 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:12:04 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS8.nc
2026-07-22 18:12:04 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:12:04 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:12:04 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS8.nc
2026-07-22 18:12:05 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:12:06 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS8.nc
2026-07-22 18:12:06 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:12:06 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 18:12:06 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:12:06 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS9.nc
2026-07-22 18:12:07 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:12:08 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS9.nc
2026-07-22 18:12:08 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:12:08 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:12:08 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS9.nc
2026-07-22 18:12:09 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:12:09 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS9.nc
2026-07-22 18:12:09 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:12:09 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 18:12:10 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:12:10 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS10.nc
2026-07-22 18:12:10 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:12:11 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS10.nc
2026-07-22 18:12:11 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:12:12 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:12:12 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS10.nc
2026-07-22 18:12:12 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:12:13 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS10.nc
2026-07-22 18:12:13 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:12:13 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 18:12:14 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:12:14 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS11.nc
2026-07-22 18:12:14 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:12:15 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS11.nc
2026-07-22 18:12:15 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:12:15 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:12:15 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS11.nc
2026-07-22 18:12:16 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:12:17 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS11.nc
2026-07-22 18:12:17 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:12:17 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 18:12:17 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:12:17 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS12.nc
2026-07-22 18:12:18 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:12:19 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS12.nc
2026-07-22 18:12:19 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:12:19 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:12:19 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS12.nc
2026-07-22 18:12:20 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:12:21 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS12.nc
2026-07-22 18:12:21 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:12:21 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 18:12:21 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:12:21 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS13.nc
2026-07-22 18:12:22 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:12:23 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS13.nc
2026-07-22 18:12:23 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:12:23 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:12:23 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS13.nc
2026-07-22 18:12:24 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:12:25 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS13.nc
2026-07-22 18:12:25 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:12:25 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 18:12:25 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:12:25 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS14.nc
2026-07-22 18:12:26 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:12:27 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS14.nc
2026-07-22 18:12:27 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:12:27 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:12:27 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS14.nc
2026-07-22 18:12:28 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:12:28 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS14.nc
2026-07-22 18:12:28 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:12:28 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 18:12:29 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:12:29 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS15.nc
2026-07-22 18:12:29 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:12:30 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS15.nc
2026-07-22 18:12:30 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:12:31 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-22 18:12:31 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS15.nc
2026-07-22 18:12:32 INFO Applying persistent orbit corrections to next day's emission file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc
2026-07-22 18:12:32 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMI_nxday_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS15.nc
2026-07-22 18:12:32 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-22 18:12:32 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-22 18:12:32 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-22 18:12:32 INFO [TIME] step_end current_time=2020-02-07 12:00:00 simulated_time=2020-02-07 12:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 18:12:32 INFO [TIME] step_start current_time=2020-02-07 12:00:00 simulated_time=2020-02-07 12:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-22 18:12:32 INFO [TIME] window start=2020-02-07 12:00:00 end=2020-02-07 14:00:00 run_hours=2 has_assimilation=True
2026-07-22 18:12:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 18:12:34 INFO Hourly dataset computed and listing created
2026-07-22 18:12:39 INFO Hourly dataset computed
2026-07-22 18:12:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 18:12:40 INFO Hourly dataset computed and listing created
2026-07-22 18:12:40 INFO Hourly dataset computed
2026-07-22 18:12:40 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 18:12:42 INFO Hourly dataset computed and listing created
2026-07-22 18:12:42 INFO Hourly dataset computed
2026-07-22 18:12:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 18:12:43 INFO Hourly dataset computed and listing created
2026-07-22 18:12:44 INFO Hourly dataset computed
2026-07-22 18:12:44 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 18:12:45 INFO Hourly dataset computed and listing created
2026-07-22 18:12:46 INFO Hourly dataset computed
2026-07-22 18:12:46 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 18:12:47 INFO Hourly dataset computed and listing created
2026-07-22 18:12:48 INFO Hourly dataset computed
2026-07-22 18:12:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 18:12:49 INFO Hourly dataset computed and listing created
2026-07-22 18:12:50 INFO Hourly dataset computed
2026-07-22 18:12:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 18:12:51 INFO Hourly dataset computed and listing created
2026-07-22 18:12:51 INFO Hourly dataset computed
2026-07-22 18:12:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 18:12:52 INFO Hourly dataset computed and listing created
2026-07-22 18:12:53 INFO Hourly dataset computed
2026-07-22 18:12:53 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 18:12:54 INFO Hourly dataset computed and listing created
2026-07-22 18:12:55 INFO Hourly dataset computed
2026-07-22 18:12:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 18:12:56 INFO Hourly dataset computed and listing created
2026-07-22 18:12:57 INFO Hourly dataset computed
2026-07-22 18:12:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 18:12:58 INFO Hourly dataset computed and listing created
2026-07-22 18:12:59 INFO Hourly dataset computed
2026-07-22 18:12:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 18:13:00 INFO Hourly dataset computed and listing created
2026-07-22 18:13:00 INFO Hourly dataset computed
2026-07-22 18:13:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 18:13:01 INFO Hourly dataset computed and listing created
2026-07-22 18:13:02 INFO Hourly dataset computed
2026-07-22 18:13:02 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-22 18:13:03 INFO Hourly dataset computed and listing created
2026-07-22 18:13:04 INFO Hourly dataset computed
2026-07-22 18:13:04 INFO ---------->>> Running CHIMERE model from 2020-02-07 12:00:00 to 2020-02-07 14:00:00
2026-07-22 18:13:04 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 18:13:04 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1
2026-07-22 18:13:04 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS1/end.2020020711_1_ENS1.nc
2026-07-22 18:13:04 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-22 18:13:04 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 18:13:04 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-22 18:13:04 INFO Queuing job for member 1...
2026-07-22 18:13:04 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 18:13:04 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-22 18:13:05 INFO Found: ['5257451']
2026-07-22 18:13:10 INFO [TGCC-IRENE] Submitted job with ID:['5257451']
2026-07-22 18:13:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 18:13:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2
2026-07-22 18:13:10 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS2/end.2020020711_1_ENS2.nc
2026-07-22 18:13:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-22 18:13:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 18:13:10 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-22 18:13:10 INFO Queuing job for member 2...
2026-07-22 18:13:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 18:13:10 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-22 18:13:11 INFO Found: ['5257452']
2026-07-22 18:13:16 INFO [TGCC-IRENE] Submitted job with ID:['5257452']
2026-07-22 18:13:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 18:13:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3
2026-07-22 18:13:16 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS3/end.2020020711_1_ENS3.nc
2026-07-22 18:13:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-22 18:13:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 18:13:16 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-22 18:13:16 INFO Queuing job for member 3...
2026-07-22 18:13:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 18:13:16 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-22 18:13:16 INFO Found: ['5257453']
2026-07-22 18:13:21 INFO [TGCC-IRENE] Submitted job with ID:['5257453']
2026-07-22 18:13:21 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 18:13:21 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4
2026-07-22 18:13:21 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS4/end.2020020711_1_ENS4.nc
2026-07-22 18:13:21 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-22 18:13:21 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 18:13:21 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-22 18:13:21 INFO Queuing job for member 4...
2026-07-22 18:13:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 18:13:21 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-22 18:13:22 INFO Found: ['5257455']
2026-07-22 18:13:27 INFO [TGCC-IRENE] Submitted job with ID:['5257455']
2026-07-22 18:13:27 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 18:13:27 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5
2026-07-22 18:13:27 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS5/end.2020020711_1_ENS5.nc
2026-07-22 18:13:27 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-22 18:13:27 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 18:13:27 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-22 18:13:27 INFO Queuing job for member 5...
2026-07-22 18:13:27 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 18:13:27 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-22 18:13:28 INFO Found: ['5257456']
2026-07-22 18:13:33 INFO [TGCC-IRENE] Submitted job with ID:['5257456']
2026-07-22 18:13:33 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 18:13:33 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6
2026-07-22 18:13:33 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS6/end.2020020711_1_ENS6.nc
2026-07-22 18:13:33 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-22 18:13:33 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 18:13:33 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-22 18:13:33 INFO Queuing job for member 6...
2026-07-22 18:13:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 18:13:33 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-22 18:13:34 INFO Found: ['5257457']
2026-07-22 18:13:39 INFO [TGCC-IRENE] Submitted job with ID:['5257457']
2026-07-22 18:13:39 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 18:13:39 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7
2026-07-22 18:13:39 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS7/end.2020020711_1_ENS7.nc
2026-07-22 18:13:39 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-22 18:13:39 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 18:13:39 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-22 18:13:39 INFO Queuing job for member 7...
2026-07-22 18:13:39 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 18:13:39 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-22 18:13:41 INFO Found: ['5257458']
2026-07-22 18:13:46 INFO [TGCC-IRENE] Submitted job with ID:['5257458']
2026-07-22 18:13:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-22 18:13:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8
2026-07-22 18:13:46 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMI_nxday_0607_15m_low_v2/ENS8/end.2020020711_1_ENS8.nc
2026-07-22 18:13:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-22 18:13:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-22 18:13:46 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-22 18:13:46 INFO Queuing job for member 8...
2026-07-22 18:13:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-22 18:13:46 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-22 18:13:48 INFO Found: ['5257465']
[2026-07-22T18:13:51.509] error: *** JOB 5255814 ON irene4625 CANCELLED AT 2026-07-22T18:13:51 DUE to SIGNAL Terminated ***
