+ /bin/bash -x /tmp/tmp.Ri3B7v6PrC
+ SCRIPT_PID=859931
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
2026-07-14 09:24:56 INFO 
███    ███ ██ ███    ███ ███████ ███████ ██
████  ████ ██ ████  ████ ██      ██      ██
██ ████ ██ ██ ██ ████ ██ █████   ███████ ██
██  ██  ██ ██ ██  ██  ██ ██           ██ ██
██      ██ ██ ██      ██ ███████ ███████ ██



2026-07-14 09:24:56 INFO [PIPELINE] =======================================
2026-07-14 09:24:56 INFO [PIPELINE] Starting chimere–DART orchestrator
2026-07-14 09:24:56 INFO [PIPELINE] Config file: config/config_irene_IM.yaml
2026-07-14 09:24:56 INFO [PIPELINE] Run dir: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart
2026-07-14 09:24:56 INFO [PIPELINE] Log file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/mimesi_orchestrator_logs/chimere_DART_20260714_092456.log
2026-07-14 09:24:56 INFO [PIPELINE] =======================================
2026-07-14 09:24:57 INFO Running assimilation with model_type=ModelType.CHIMERE
2026-07-14 09:24:57 INFO Using scheduler=Scheduler.SLURM, queue=rome
2026-07-14 09:24:57 INFO [STEP] ---- TIME LOOP START ----
2026-07-14 09:24:57 INFO [TIME] step_start current_time=2020-02-06 00:00:00 simulated_time=None dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 09:24:57 INFO [TIME] window start=2020-02-06 00:00:00 end=2020-02-06 01:00:00 run_hours=1 has_assimilation=False
2026-07-14 09:24:57 INFO Creating directories and links for ENS1 to run chimere's parallel part
2026-07-14 09:24:57 INFO Copying EMIS ...
2026-07-14 09:24:57 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens01.nc
2026-07-14 09:24:57 INFO Linking first END ...
2026-07-14 09:24:57 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020500_24_ENS1.nc
2026-07-14 09:24:57 INFO >> Checking links...
2026-07-14 09:24:57 INFO >> All links are good for ENS1  ...
2026-07-14 09:24:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:25:05 INFO Hourly dataset computed and listing created
2026-07-14 09:25:17 INFO Hourly dataset computed
2026-07-14 09:25:17 INFO Creating directories and links for ENS2 to run chimere's parallel part
2026-07-14 09:25:17 INFO Copying EMIS ...
2026-07-14 09:25:17 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens02.nc
2026-07-14 09:25:17 INFO Linking first END ...
2026-07-14 09:25:17 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020500_24_ENS2.nc
2026-07-14 09:25:17 INFO >> Checking links...
2026-07-14 09:25:18 INFO >> All links are good for ENS2  ...
2026-07-14 09:25:18 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:25:20 INFO Hourly dataset computed and listing created
2026-07-14 09:25:33 INFO Hourly dataset computed
2026-07-14 09:25:33 INFO Creating directories and links for ENS3 to run chimere's parallel part
2026-07-14 09:25:33 INFO Copying EMIS ...
2026-07-14 09:25:34 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens03.nc
2026-07-14 09:25:34 INFO Linking first END ...
2026-07-14 09:25:34 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020500_24_ENS3.nc
2026-07-14 09:25:34 INFO >> Checking links...
2026-07-14 09:25:34 INFO >> All links are good for ENS3  ...
2026-07-14 09:25:34 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:25:36 INFO Hourly dataset computed and listing created
2026-07-14 09:25:40 INFO Hourly dataset computed
2026-07-14 09:25:40 INFO Creating directories and links for ENS4 to run chimere's parallel part
2026-07-14 09:25:40 INFO Copying EMIS ...
2026-07-14 09:25:41 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens04.nc
2026-07-14 09:25:41 INFO Linking first END ...
2026-07-14 09:25:41 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020500_24_ENS4.nc
2026-07-14 09:25:41 INFO >> Checking links...
2026-07-14 09:25:41 INFO >> All links are good for ENS4  ...
2026-07-14 09:25:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:25:42 INFO Hourly dataset computed and listing created
2026-07-14 09:25:44 INFO Hourly dataset computed
2026-07-14 09:25:44 INFO Creating directories and links for ENS5 to run chimere's parallel part
2026-07-14 09:25:44 INFO Copying EMIS ...
2026-07-14 09:25:45 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens05.nc
2026-07-14 09:25:45 INFO Linking first END ...
2026-07-14 09:25:45 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020500_24_ENS5.nc
2026-07-14 09:25:45 INFO >> Checking links...
2026-07-14 09:25:45 INFO >> All links are good for ENS5  ...
2026-07-14 09:25:45 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:25:46 INFO Hourly dataset computed and listing created
2026-07-14 09:25:49 INFO Hourly dataset computed
2026-07-14 09:25:49 INFO Creating directories and links for ENS6 to run chimere's parallel part
2026-07-14 09:25:49 INFO Copying EMIS ...
2026-07-14 09:25:50 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens06.nc
2026-07-14 09:25:50 INFO Linking first END ...
2026-07-14 09:25:50 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020500_24_ENS6.nc
2026-07-14 09:25:50 INFO >> Checking links...
2026-07-14 09:25:50 INFO >> All links are good for ENS6  ...
2026-07-14 09:25:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:25:51 INFO Hourly dataset computed and listing created
2026-07-14 09:25:54 INFO Hourly dataset computed
2026-07-14 09:25:54 INFO Creating directories and links for ENS7 to run chimere's parallel part
2026-07-14 09:25:54 INFO Copying EMIS ...
2026-07-14 09:25:55 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens07.nc
2026-07-14 09:25:55 INFO Linking first END ...
2026-07-14 09:25:55 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020500_24_ENS7.nc
2026-07-14 09:25:55 INFO >> Checking links...
2026-07-14 09:25:55 INFO >> All links are good for ENS7  ...
2026-07-14 09:25:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:25:56 INFO Hourly dataset computed and listing created
2026-07-14 09:25:57 INFO Hourly dataset computed
2026-07-14 09:25:57 INFO Creating directories and links for ENS8 to run chimere's parallel part
2026-07-14 09:25:57 INFO Copying EMIS ...
2026-07-14 09:25:57 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens08.nc
2026-07-14 09:25:57 INFO Linking first END ...
2026-07-14 09:25:57 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020500_24_ENS8.nc
2026-07-14 09:25:57 INFO >> Checking links...
2026-07-14 09:25:58 INFO >> All links are good for ENS8  ...
2026-07-14 09:25:58 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:25:59 INFO Hourly dataset computed and listing created
2026-07-14 09:28:11 INFO Hourly dataset computed
2026-07-14 09:28:11 INFO Creating directories and links for ENS9 to run chimere's parallel part
2026-07-14 09:28:11 INFO Copying EMIS ...
2026-07-14 09:28:12 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens09.nc
2026-07-14 09:28:12 INFO Linking first END ...
2026-07-14 09:28:12 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020500_24_ENS9.nc
2026-07-14 09:28:12 INFO >> Checking links...
2026-07-14 09:28:12 INFO >> All links are good for ENS9  ...
2026-07-14 09:28:12 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:28:14 INFO Hourly dataset computed and listing created
2026-07-14 09:28:17 INFO Hourly dataset computed
2026-07-14 09:28:17 INFO Creating directories and links for ENS10 to run chimere's parallel part
2026-07-14 09:28:17 INFO Copying EMIS ...
2026-07-14 09:28:18 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens10.nc
2026-07-14 09:28:18 INFO Linking first END ...
2026-07-14 09:28:18 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020500_24_ENS10.nc
2026-07-14 09:28:18 INFO >> Checking links...
2026-07-14 09:28:18 INFO >> All links are good for ENS10  ...
2026-07-14 09:28:18 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:28:19 INFO Hourly dataset computed and listing created
2026-07-14 09:28:20 INFO Hourly dataset computed
2026-07-14 09:28:20 INFO Creating directories and links for ENS11 to run chimere's parallel part
2026-07-14 09:28:20 INFO Copying EMIS ...
2026-07-14 09:28:20 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens11.nc
2026-07-14 09:28:20 INFO Linking first END ...
2026-07-14 09:28:20 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020500_24_ENS11.nc
2026-07-14 09:28:20 INFO >> Checking links...
2026-07-14 09:28:21 INFO >> All links are good for ENS11  ...
2026-07-14 09:28:21 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:28:21 INFO Hourly dataset computed and listing created
2026-07-14 09:28:22 INFO Hourly dataset computed
2026-07-14 09:28:22 INFO Creating directories and links for ENS12 to run chimere's parallel part
2026-07-14 09:28:22 INFO Copying EMIS ...
2026-07-14 09:28:23 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens12.nc
2026-07-14 09:28:23 INFO Linking first END ...
2026-07-14 09:28:23 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020500_24_ENS12.nc
2026-07-14 09:28:23 INFO >> Checking links...
2026-07-14 09:28:23 INFO >> All links are good for ENS12  ...
2026-07-14 09:28:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:28:24 INFO Hourly dataset computed and listing created
2026-07-14 09:28:25 INFO Hourly dataset computed
2026-07-14 09:28:25 INFO Creating directories and links for ENS13 to run chimere's parallel part
2026-07-14 09:28:25 INFO Copying EMIS ...
2026-07-14 09:28:25 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens13.nc
2026-07-14 09:28:25 INFO Linking first END ...
2026-07-14 09:28:25 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020500_24_ENS13.nc
2026-07-14 09:28:25 INFO >> Checking links...
2026-07-14 09:28:26 INFO >> All links are good for ENS13  ...
2026-07-14 09:28:26 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:28:27 INFO Hourly dataset computed and listing created
2026-07-14 09:28:27 INFO Hourly dataset computed
2026-07-14 09:28:27 INFO Creating directories and links for ENS14 to run chimere's parallel part
2026-07-14 09:28:27 INFO Copying EMIS ...
2026-07-14 09:28:28 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens14.nc
2026-07-14 09:28:28 INFO Linking first END ...
2026-07-14 09:28:28 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020500_24_ENS14.nc
2026-07-14 09:28:28 INFO >> Checking links...
2026-07-14 09:28:28 INFO >> All links are good for ENS14  ...
2026-07-14 09:28:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:28:29 INFO Hourly dataset computed and listing created
2026-07-14 09:28:30 INFO Hourly dataset computed
2026-07-14 09:28:30 INFO Creating directories and links for ENS15 to run chimere's parallel part
2026-07-14 09:28:30 INFO Copying EMIS ...
2026-07-14 09:28:30 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens15.nc
2026-07-14 09:28:30 INFO Linking first END ...
2026-07-14 09:28:30 INFO Symlink already exists and is correct: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020500_24_ENS15.nc
2026-07-14 09:28:30 INFO >> Checking links...
2026-07-14 09:28:31 INFO >> All links are good for ENS15  ...
2026-07-14 09:28:31 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:28:32 INFO Hourly dataset computed and listing created
2026-07-14 09:28:32 INFO Hourly dataset computed
2026-07-14 09:28:32 INFO ---------->>> Running CHIMERE model from 2020-02-06 00:00:00 to 2020-02-06 01:00:00
2026-07-14 09:28:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:28:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1
2026-07-14 09:28:32 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020500_24_ENS1.nc
2026-07-14 09:28:33 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-14 09:28:33 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:28:33 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-14 09:28:33 INFO Queuing job for member 1...
2026-07-14 09:28:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:28:33 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-14 09:28:34 INFO Found: ['5154530']
2026-07-14 09:28:39 INFO [TGCC-IRENE] Submitted job with ID:['5154530']
2026-07-14 09:28:39 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:28:39 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2
2026-07-14 09:28:39 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020500_24_ENS2.nc
2026-07-14 09:28:39 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-14 09:28:39 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:28:39 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-14 09:28:39 INFO Queuing job for member 2...
2026-07-14 09:28:39 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:28:39 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-14 09:28:40 INFO Found: ['5154531']
2026-07-14 09:28:45 INFO [TGCC-IRENE] Submitted job with ID:['5154531']
2026-07-14 09:28:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:28:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3
2026-07-14 09:28:45 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020500_24_ENS3.nc
2026-07-14 09:28:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-14 09:28:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:28:45 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-14 09:28:45 INFO Queuing job for member 3...
2026-07-14 09:28:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:28:45 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-14 09:28:45 INFO Found: ['5154532']
2026-07-14 09:28:50 INFO [TGCC-IRENE] Submitted job with ID:['5154532']
2026-07-14 09:28:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:28:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4
2026-07-14 09:28:50 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020500_24_ENS4.nc
2026-07-14 09:28:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-14 09:28:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:28:51 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-14 09:28:51 INFO Queuing job for member 4...
2026-07-14 09:28:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:28:51 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-14 09:28:51 INFO Found: ['5154533']
2026-07-14 09:28:56 INFO [TGCC-IRENE] Submitted job with ID:['5154533']
2026-07-14 09:28:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:28:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5
2026-07-14 09:28:57 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020500_24_ENS5.nc
2026-07-14 09:28:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-14 09:28:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:28:57 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-14 09:28:57 INFO Queuing job for member 5...
2026-07-14 09:28:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:28:57 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-14 09:28:57 INFO Found: ['5154534']
2026-07-14 09:29:02 INFO [TGCC-IRENE] Submitted job with ID:['5154534']
2026-07-14 09:29:02 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:29:02 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6
2026-07-14 09:29:02 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020500_24_ENS6.nc
2026-07-14 09:29:02 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-14 09:29:02 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:29:02 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-14 09:29:02 INFO Queuing job for member 6...
2026-07-14 09:29:02 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:29:02 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-14 09:29:03 INFO Found: ['5154536']
2026-07-14 09:29:08 INFO [TGCC-IRENE] Submitted job with ID:['5154536']
2026-07-14 09:29:08 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:29:08 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7
2026-07-14 09:29:08 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020500_24_ENS7.nc
2026-07-14 09:29:08 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-14 09:29:08 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:29:08 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-14 09:29:08 INFO Queuing job for member 7...
2026-07-14 09:29:08 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:29:08 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-14 09:29:09 INFO Found: ['5154537']
2026-07-14 09:29:14 INFO [TGCC-IRENE] Submitted job with ID:['5154537']
2026-07-14 09:29:14 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:29:14 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8
2026-07-14 09:29:14 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020500_24_ENS8.nc
2026-07-14 09:29:14 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-14 09:29:14 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:29:14 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-14 09:29:14 INFO Queuing job for member 8...
2026-07-14 09:29:14 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:29:14 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-14 09:29:15 INFO Found: ['5154538']
2026-07-14 09:29:20 INFO [TGCC-IRENE] Submitted job with ID:['5154538']
2026-07-14 09:29:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:29:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9
2026-07-14 09:29:20 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020500_24_ENS9.nc
2026-07-14 09:29:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-14 09:29:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:29:20 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-14 09:29:20 INFO Queuing job for member 9...
2026-07-14 09:29:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:29:20 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-14 09:29:21 INFO Found: ['5154539']
2026-07-14 09:29:26 INFO [TGCC-IRENE] Submitted job with ID:['5154539']
2026-07-14 09:29:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:29:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10
2026-07-14 09:29:26 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020500_24_ENS10.nc
2026-07-14 09:29:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-14 09:29:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:29:26 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-14 09:29:26 INFO Queuing job for member 10...
2026-07-14 09:29:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:29:26 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-14 09:29:27 INFO Found: ['5154540']
2026-07-14 09:29:32 INFO [TGCC-IRENE] Submitted job with ID:['5154540']
2026-07-14 09:29:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:29:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11
2026-07-14 09:29:32 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020500_24_ENS11.nc
2026-07-14 09:29:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-14 09:29:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:29:32 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-14 09:29:32 INFO Queuing job for member 11...
2026-07-14 09:29:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:29:32 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-14 09:29:33 INFO Found: ['5154541']
2026-07-14 09:29:38 INFO [TGCC-IRENE] Submitted job with ID:['5154541']
2026-07-14 09:29:38 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:29:38 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12
2026-07-14 09:29:38 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020500_24_ENS12.nc
2026-07-14 09:29:38 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-14 09:29:38 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:29:38 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-14 09:29:38 INFO Queuing job for member 12...
2026-07-14 09:29:38 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:29:38 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-14 09:29:38 INFO Found: ['5154542']
2026-07-14 09:29:43 INFO [TGCC-IRENE] Submitted job with ID:['5154542']
2026-07-14 09:29:43 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:29:43 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13
2026-07-14 09:29:43 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020500_24_ENS13.nc
2026-07-14 09:29:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-14 09:29:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:29:43 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-14 09:29:43 INFO Queuing job for member 13...
2026-07-14 09:29:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:29:43 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-14 09:29:44 INFO Found: ['5154543']
2026-07-14 09:29:49 INFO [TGCC-IRENE] Submitted job with ID:['5154543']
2026-07-14 09:29:49 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:29:49 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14
2026-07-14 09:29:49 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020500_24_ENS14.nc
2026-07-14 09:29:49 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-14 09:29:49 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:29:49 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-14 09:29:49 INFO Queuing job for member 14...
2026-07-14 09:29:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:29:49 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-14 09:29:50 INFO Found: ['5154544']
2026-07-14 09:29:55 INFO [TGCC-IRENE] Submitted job with ID:['5154544']
2026-07-14 09:29:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:29:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15
2026-07-14 09:29:55 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020500_24_ENS15.nc
2026-07-14 09:29:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-14 09:29:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:29:55 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-14 09:29:55 INFO Queuing job for member 15...
2026-07-14 09:29:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:29:55 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-14 09:29:56 INFO Found: ['5154545']
2026-07-14 09:30:01 INFO [TGCC-IRENE] Submitted job with ID:['5154545']
2026-07-14 09:30:01 INFO Checking job status ...
2026-07-14 09:30:01 INFO None 5154530: status RUNNING/PENDING
2026-07-14 09:30:01 INFO None 5154531: status RUNNING/PENDING
2026-07-14 09:30:01 INFO None 5154532: status RUNNING/PENDING
2026-07-14 09:30:01 INFO None 5154533: status RUNNING/PENDING
2026-07-14 09:30:01 INFO None 5154534: status RUNNING/PENDING
2026-07-14 09:30:01 INFO None 5154536: status RUNNING/PENDING
2026-07-14 09:30:01 INFO None 5154537: status RUNNING/PENDING
2026-07-14 09:30:01 INFO None 5154538: status RUNNING/PENDING
2026-07-14 09:30:01 INFO None 5154539: status RUNNING/PENDING
2026-07-14 09:30:01 INFO None 5154540: status RUNNING/PENDING
2026-07-14 09:30:01 INFO None 5154541: status RUNNING/PENDING
2026-07-14 09:30:01 INFO None 5154542: status RUNNING/PENDING
2026-07-14 09:30:01 INFO None 5154543: status RUNNING/PENDING
2026-07-14 09:30:01 INFO None 5154544: status RUNNING/PENDING
2026-07-14 09:30:01 INFO None 5154545: status RUNNING/PENDING
2026-07-14 09:30:01 INFO Jobs still running: ['5154530', '5154531', '5154532', '5154533', '5154534', '5154536', '5154537', '5154538', '5154539', '5154540', '5154541', '5154542', '5154543', '5154544', '5154545']. Waiting...
2026-07-14 09:30:16 INFO None 5154530: status RUNNING/PENDING
2026-07-14 09:30:16 INFO None 5154531: status RUNNING/PENDING
2026-07-14 09:30:16 INFO None 5154532: status RUNNING/PENDING
2026-07-14 09:30:16 INFO None 5154533: status RUNNING/PENDING
2026-07-14 09:30:17 INFO None 5154534: status RUNNING/PENDING
2026-07-14 09:30:17 INFO None 5154536: status RUNNING/PENDING
2026-07-14 09:30:17 INFO None 5154537: status RUNNING/PENDING
2026-07-14 09:30:17 INFO None 5154538: status RUNNING/PENDING
2026-07-14 09:30:17 INFO None 5154539: status RUNNING/PENDING
2026-07-14 09:30:17 INFO None 5154540: status RUNNING/PENDING
2026-07-14 09:30:17 INFO None 5154541: status RUNNING/PENDING
2026-07-14 09:30:17 INFO None 5154542: status RUNNING/PENDING
2026-07-14 09:30:17 INFO None 5154543: status RUNNING/PENDING
2026-07-14 09:30:17 INFO None 5154544: status RUNNING/PENDING
2026-07-14 09:30:17 INFO None 5154545: status RUNNING/PENDING
2026-07-14 09:30:17 INFO Jobs still running: ['5154530', '5154531', '5154532', '5154533', '5154534', '5154536', '5154537', '5154538', '5154539', '5154540', '5154541', '5154542', '5154543', '5154544', '5154545']. Waiting...
2026-07-14 09:30:32 INFO None 5154530: status RUNNING/PENDING
2026-07-14 09:30:32 INFO None 5154531: status RUNNING/PENDING
2026-07-14 09:30:32 INFO None 5154532: status RUNNING/PENDING
2026-07-14 09:30:32 INFO None 5154533: status RUNNING/PENDING
2026-07-14 09:30:32 INFO None 5154534: status RUNNING/PENDING
2026-07-14 09:30:32 INFO None 5154536: status RUNNING/PENDING
2026-07-14 09:30:32 INFO None 5154537: status RUNNING/PENDING
2026-07-14 09:30:32 INFO None 5154538: status RUNNING/PENDING
2026-07-14 09:30:32 INFO None 5154539: status RUNNING/PENDING
2026-07-14 09:30:32 INFO None 5154540: status RUNNING/PENDING
2026-07-14 09:30:32 INFO None 5154541: status RUNNING/PENDING
2026-07-14 09:30:32 INFO None 5154542: status RUNNING/PENDING
2026-07-14 09:30:32 INFO None 5154543: status RUNNING/PENDING
2026-07-14 09:30:32 INFO None 5154544: status RUNNING/PENDING
2026-07-14 09:30:32 INFO None 5154545: status RUNNING/PENDING
2026-07-14 09:30:32 INFO Jobs still running: ['5154530', '5154531', '5154532', '5154533', '5154534', '5154536', '5154537', '5154538', '5154539', '5154540', '5154541', '5154542', '5154543', '5154544', '5154545']. Waiting...
2026-07-14 09:30:47 INFO None 5154530: status RUNNING/PENDING
2026-07-14 09:30:47 INFO None 5154531: status RUNNING/PENDING
2026-07-14 09:30:47 INFO None 5154532: status RUNNING/PENDING
2026-07-14 09:30:47 INFO None 5154533: status RUNNING/PENDING
2026-07-14 09:30:47 INFO None 5154534: status RUNNING/PENDING
2026-07-14 09:30:47 INFO None 5154536: status RUNNING/PENDING
2026-07-14 09:30:47 INFO None 5154537: status RUNNING/PENDING
2026-07-14 09:30:47 INFO None 5154538: status RUNNING/PENDING
2026-07-14 09:30:47 INFO None 5154539: status RUNNING/PENDING
2026-07-14 09:30:47 INFO None 5154540: status RUNNING/PENDING
2026-07-14 09:30:47 INFO None 5154541: status RUNNING/PENDING
2026-07-14 09:30:47 INFO None 5154542: status RUNNING/PENDING
2026-07-14 09:30:47 INFO None 5154543: status RUNNING/PENDING
2026-07-14 09:30:47 INFO None 5154544: status RUNNING/PENDING
2026-07-14 09:30:47 INFO None 5154545: status RUNNING/PENDING
2026-07-14 09:30:47 INFO Jobs still running: ['5154530', '5154531', '5154532', '5154533', '5154534', '5154536', '5154537', '5154538', '5154539', '5154540', '5154541', '5154542', '5154543', '5154544', '5154545']. Waiting...
2026-07-14 09:31:02 INFO None 5154530: status RUNNING/PENDING
2026-07-14 09:31:02 INFO None 5154531: status RUNNING/PENDING
2026-07-14 09:31:02 INFO None 5154532: status RUNNING/PENDING
2026-07-14 09:31:02 INFO None 5154533: status RUNNING/PENDING
2026-07-14 09:31:03 INFO None 5154534: status RUNNING/PENDING
2026-07-14 09:31:03 INFO None 5154536: status RUNNING/PENDING
2026-07-14 09:31:03 INFO None 5154537: status RUNNING/PENDING
2026-07-14 09:31:03 INFO None 5154538: status RUNNING/PENDING
2026-07-14 09:31:03 INFO None 5154539: status RUNNING/PENDING
2026-07-14 09:31:03 INFO None 5154540: status RUNNING/PENDING
2026-07-14 09:31:03 INFO None 5154541: status RUNNING/PENDING
2026-07-14 09:31:03 INFO None 5154542: status RUNNING/PENDING
2026-07-14 09:31:03 INFO None 5154543: status RUNNING/PENDING
2026-07-14 09:31:03 INFO None 5154544: status RUNNING/PENDING
2026-07-14 09:31:03 INFO None 5154545: status RUNNING/PENDING
2026-07-14 09:31:03 INFO Jobs still running: ['5154530', '5154531', '5154532', '5154533', '5154534', '5154536', '5154537', '5154538', '5154539', '5154540', '5154541', '5154542', '5154543', '5154544', '5154545']. Waiting...
2026-07-14 09:31:18 INFO None 5154530: status FINISHED
2026-07-14 09:31:18 INFO None 5154531: status FINISHED
2026-07-14 09:31:18 INFO None 5154532: status FINISHED
2026-07-14 09:31:18 INFO None 5154533: status FINISHED
2026-07-14 09:31:18 INFO None 5154534: status FINISHED
2026-07-14 09:31:18 INFO None 5154536: status RUNNING/PENDING
2026-07-14 09:31:18 INFO None 5154537: status RUNNING/PENDING
2026-07-14 09:31:18 INFO None 5154538: status RUNNING/PENDING
2026-07-14 09:31:18 INFO None 5154539: status RUNNING/PENDING
2026-07-14 09:31:18 INFO None 5154540: status RUNNING/PENDING
2026-07-14 09:31:18 INFO None 5154541: status RUNNING/PENDING
2026-07-14 09:31:18 INFO None 5154542: status RUNNING/PENDING
2026-07-14 09:31:18 INFO None 5154543: status RUNNING/PENDING
2026-07-14 09:31:18 INFO None 5154544: status RUNNING/PENDING
2026-07-14 09:31:18 INFO None 5154545: status RUNNING/PENDING
2026-07-14 09:31:18 INFO Jobs still running: ['5154536', '5154537', '5154538', '5154539', '5154540', '5154541', '5154542', '5154543', '5154544', '5154545']. Waiting...
2026-07-14 09:31:33 INFO None 5154530: status FINISHED
2026-07-14 09:31:33 INFO None 5154531: status FINISHED
2026-07-14 09:31:33 INFO None 5154532: status FINISHED
2026-07-14 09:31:33 INFO None 5154533: status FINISHED
2026-07-14 09:31:33 INFO None 5154534: status FINISHED
2026-07-14 09:31:33 INFO None 5154536: status FINISHED
2026-07-14 09:31:33 INFO None 5154537: status RUNNING/PENDING
2026-07-14 09:31:33 INFO None 5154538: status RUNNING/PENDING
2026-07-14 09:31:33 INFO None 5154539: status RUNNING/PENDING
2026-07-14 09:31:33 INFO None 5154540: status RUNNING/PENDING
2026-07-14 09:31:33 INFO None 5154541: status RUNNING/PENDING
2026-07-14 09:31:33 INFO None 5154542: status RUNNING/PENDING
2026-07-14 09:31:33 INFO None 5154543: status RUNNING/PENDING
2026-07-14 09:31:33 INFO None 5154544: status RUNNING/PENDING
2026-07-14 09:31:33 INFO None 5154545: status RUNNING/PENDING
2026-07-14 09:31:33 INFO Jobs still running: ['5154537', '5154538', '5154539', '5154540', '5154541', '5154542', '5154543', '5154544', '5154545']. Waiting...
2026-07-14 09:31:48 INFO None 5154530: status FINISHED
2026-07-14 09:31:48 INFO None 5154531: status FINISHED
2026-07-14 09:31:48 INFO None 5154532: status FINISHED
2026-07-14 09:31:48 INFO None 5154533: status FINISHED
2026-07-14 09:31:48 INFO None 5154534: status FINISHED
2026-07-14 09:31:48 INFO None 5154536: status FINISHED
2026-07-14 09:31:49 INFO None 5154537: status RUNNING/PENDING
2026-07-14 09:31:49 INFO None 5154538: status RUNNING/PENDING
2026-07-14 09:31:49 INFO None 5154539: status RUNNING/PENDING
2026-07-14 09:31:49 INFO None 5154540: status RUNNING/PENDING
2026-07-14 09:31:49 INFO None 5154541: status RUNNING/PENDING
2026-07-14 09:31:49 INFO None 5154542: status RUNNING/PENDING
2026-07-14 09:31:49 INFO None 5154543: status RUNNING/PENDING
2026-07-14 09:31:49 INFO None 5154544: status FINISHED
2026-07-14 09:31:49 INFO None 5154545: status FINISHED
2026-07-14 09:31:49 INFO Jobs still running: ['5154537', '5154538', '5154539', '5154540', '5154541', '5154542', '5154543']. Waiting...
2026-07-14 09:32:04 INFO None 5154530: status FINISHED
2026-07-14 09:33:31 INFO None 5154531: status FINISHED
2026-07-14 09:33:31 INFO None 5154532: status FINISHED
2026-07-14 09:33:31 INFO None 5154533: status FINISHED
2026-07-14 09:33:31 INFO None 5154534: status FINISHED
2026-07-14 09:33:31 INFO None 5154536: status FINISHED
2026-07-14 09:33:31 INFO None 5154537: status FINISHED
2026-07-14 09:33:31 INFO None 5154538: status FINISHED
2026-07-14 09:33:31 INFO None 5154539: status FINISHED
2026-07-14 09:33:31 INFO None 5154540: status FINISHED
2026-07-14 09:33:31 INFO None 5154541: status FINISHED
2026-07-14 09:33:31 INFO None 5154542: status FINISHED
2026-07-14 09:33:31 INFO None 5154543: status FINISHED
2026-07-14 09:33:31 INFO None 5154544: status FINISHED
2026-07-14 09:33:31 INFO None 5154545: status FINISHED
2026-07-14 09:33:31 INFO Jobs ['5154530', '5154531', '5154532', '5154533', '5154534', '5154536', '5154537', '5154538', '5154539', '5154540', '5154541', '5154542', '5154543', '5154544', '5154545'] have finished
2026-07-14 09:33:31 INFO Checking restart files were created ...
2026-07-14 09:33:31 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020600_1_ENS1.nc(668832435 bytes)
2026-07-14 09:33:31 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020600_1_ENS2.nc(668832435 bytes)
2026-07-14 09:33:31 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020600_1_ENS3.nc(668832435 bytes)
2026-07-14 09:33:31 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020600_1_ENS4.nc(668832435 bytes)
2026-07-14 09:33:31 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020600_1_ENS5.nc(668832435 bytes)
2026-07-14 09:33:31 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020600_1_ENS6.nc(668832435 bytes)
2026-07-14 09:33:31 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020600_1_ENS7.nc(668832435 bytes)
2026-07-14 09:33:31 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020600_1_ENS8.nc(668832435 bytes)
2026-07-14 09:33:31 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020600_1_ENS9.nc(668832435 bytes)
2026-07-14 09:33:31 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020600_1_ENS10.nc(668832435 bytes)
2026-07-14 09:33:31 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020600_1_ENS11.nc(668832435 bytes)
2026-07-14 09:33:31 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020600_1_ENS12.nc(668832435 bytes)
2026-07-14 09:33:31 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020600_1_ENS13.nc(668832435 bytes)
2026-07-14 09:33:31 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020600_1_ENS14.nc(668832435 bytes)
2026-07-14 09:33:31 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020600_1_ENS15.nc(668832435 bytes)
2026-07-14 09:33:31 INFO  Run_model() completed successfully.
2026-07-14 09:33:31 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 00:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 09:33:31 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 01:00:00 days=153072 seconds=3600
2026-07-14 09:33:31 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-14 09:33:31 INFO [TIME] increment current_time 2020-02-06 00:00:00 -> 2020-02-06 01:00:00
2026-07-14 09:33:31 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 09:33:31 INFO ---------->>> Running process_satellite_data()
2026-07-14 09:33:31 INFO [DART] No satellite data found, skipping assimilation
2026-07-14 09:33:31 INFO after_assimilation() skipped
2026-07-14 09:33:31 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-14 09:33:31 INFO [TIME] step_end current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 09:33:31 INFO [TIME] step_start current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 09:33:31 INFO [TIME] window start=2020-02-06 01:00:00 end=2020-02-06 09:00:00 run_hours=8 has_assimilation=True
2026-07-14 09:33:31 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:33:33 INFO Hourly dataset computed and listing created
2026-07-14 09:33:45 INFO Hourly dataset computed
2026-07-14 09:33:45 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:33:46 INFO Hourly dataset computed and listing created
2026-07-14 09:33:56 INFO Hourly dataset computed
2026-07-14 09:33:56 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:33:57 INFO Hourly dataset computed and listing created
2026-07-14 09:34:07 INFO Hourly dataset computed
2026-07-14 09:34:07 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:34:08 INFO Hourly dataset computed and listing created
2026-07-14 09:34:11 INFO Hourly dataset computed
2026-07-14 09:34:11 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:34:12 INFO Hourly dataset computed and listing created
2026-07-14 09:34:15 INFO Hourly dataset computed
2026-07-14 09:34:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:34:16 INFO Hourly dataset computed and listing created
2026-07-14 09:34:18 INFO Hourly dataset computed
2026-07-14 09:34:18 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:34:20 INFO Hourly dataset computed and listing created
2026-07-14 09:34:22 INFO Hourly dataset computed
2026-07-14 09:34:22 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:34:23 INFO Hourly dataset computed and listing created
2026-07-14 09:34:26 INFO Hourly dataset computed
2026-07-14 09:34:26 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:34:27 INFO Hourly dataset computed and listing created
2026-07-14 09:34:30 INFO Hourly dataset computed
2026-07-14 09:34:30 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:34:31 INFO Hourly dataset computed and listing created
2026-07-14 09:34:34 INFO Hourly dataset computed
2026-07-14 09:34:34 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:34:35 INFO Hourly dataset computed and listing created
2026-07-14 09:34:37 INFO Hourly dataset computed
2026-07-14 09:34:37 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:34:39 INFO Hourly dataset computed and listing created
2026-07-14 09:34:41 INFO Hourly dataset computed
2026-07-14 09:34:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:34:43 INFO Hourly dataset computed and listing created
2026-07-14 09:34:45 INFO Hourly dataset computed
2026-07-14 09:34:45 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:34:47 INFO Hourly dataset computed and listing created
2026-07-14 09:34:49 INFO Hourly dataset computed
2026-07-14 09:34:49 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:34:50 INFO Hourly dataset computed and listing created
2026-07-14 09:34:53 INFO Hourly dataset computed
2026-07-14 09:34:53 INFO ---------->>> Running CHIMERE model from 2020-02-06 01:00:00 to 2020-02-06 09:00:00
2026-07-14 09:34:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:34:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1
2026-07-14 09:34:53 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020600_1_ENS1.nc
2026-07-14 09:34:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-14 09:34:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:34:53 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-14 09:34:53 INFO Queuing job for member 1...
2026-07-14 09:34:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:34:53 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-14 09:34:54 INFO Found: ['5154558']
2026-07-14 09:34:59 INFO [TGCC-IRENE] Submitted job with ID:['5154558']
2026-07-14 09:34:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:34:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2
2026-07-14 09:34:59 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020600_1_ENS2.nc
2026-07-14 09:34:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-14 09:34:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:34:59 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-14 09:34:59 INFO Queuing job for member 2...
2026-07-14 09:34:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:34:59 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-14 09:35:00 INFO Found: ['5154559']
2026-07-14 09:35:05 INFO [TGCC-IRENE] Submitted job with ID:['5154559']
2026-07-14 09:35:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:35:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3
2026-07-14 09:35:05 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020600_1_ENS3.nc
2026-07-14 09:35:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-14 09:35:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:35:05 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-14 09:35:05 INFO Queuing job for member 3...
2026-07-14 09:35:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:35:05 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-14 09:35:06 INFO Found: ['5154561']
2026-07-14 09:35:11 INFO [TGCC-IRENE] Submitted job with ID:['5154561']
2026-07-14 09:35:11 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:35:11 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4
2026-07-14 09:35:11 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020600_1_ENS4.nc
2026-07-14 09:35:11 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-14 09:35:11 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:35:11 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-14 09:35:11 INFO Queuing job for member 4...
2026-07-14 09:35:11 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:35:11 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-14 09:35:11 INFO Found: ['5154562']
2026-07-14 09:35:16 INFO [TGCC-IRENE] Submitted job with ID:['5154562']
2026-07-14 09:35:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:35:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5
2026-07-14 09:35:16 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020600_1_ENS5.nc
2026-07-14 09:35:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-14 09:35:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:35:17 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-14 09:35:17 INFO Queuing job for member 5...
2026-07-14 09:35:17 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:35:17 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-14 09:35:17 INFO Found: ['5154563']
2026-07-14 09:35:22 INFO [TGCC-IRENE] Submitted job with ID:['5154563']
2026-07-14 09:35:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:35:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6
2026-07-14 09:35:22 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020600_1_ENS6.nc
2026-07-14 09:35:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-14 09:35:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:35:22 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-14 09:35:22 INFO Queuing job for member 6...
2026-07-14 09:35:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:35:22 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-14 09:35:23 INFO Found: ['5154564']
2026-07-14 09:35:28 INFO [TGCC-IRENE] Submitted job with ID:['5154564']
2026-07-14 09:35:28 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:35:28 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7
2026-07-14 09:35:28 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020600_1_ENS7.nc
2026-07-14 09:35:28 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-14 09:35:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:35:28 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-14 09:35:28 INFO Queuing job for member 7...
2026-07-14 09:35:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:35:28 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-14 09:35:29 INFO Found: ['5154565']
2026-07-14 09:35:34 INFO [TGCC-IRENE] Submitted job with ID:['5154565']
2026-07-14 09:35:34 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:35:34 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8
2026-07-14 09:35:34 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020600_1_ENS8.nc
2026-07-14 09:35:34 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-14 09:35:34 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:35:34 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-14 09:35:34 INFO Queuing job for member 8...
2026-07-14 09:35:34 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:35:34 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-14 09:35:35 INFO Found: ['5154566']
2026-07-14 09:35:40 INFO [TGCC-IRENE] Submitted job with ID:['5154566']
2026-07-14 09:35:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:35:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9
2026-07-14 09:35:40 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020600_1_ENS9.nc
2026-07-14 09:35:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-14 09:35:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:35:40 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-14 09:35:40 INFO Queuing job for member 9...
2026-07-14 09:35:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:35:40 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-14 09:35:40 INFO Found: ['5154567']
2026-07-14 09:35:45 INFO [TGCC-IRENE] Submitted job with ID:['5154567']
2026-07-14 09:35:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:35:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10
2026-07-14 09:35:46 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020600_1_ENS10.nc
2026-07-14 09:35:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-14 09:35:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:35:46 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-14 09:35:46 INFO Queuing job for member 10...
2026-07-14 09:35:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:35:46 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-14 09:35:46 INFO Found: ['5154568']
2026-07-14 09:35:51 INFO [TGCC-IRENE] Submitted job with ID:['5154568']
2026-07-14 09:35:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:35:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11
2026-07-14 09:35:51 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020600_1_ENS11.nc
2026-07-14 09:35:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-14 09:35:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:35:51 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-14 09:35:51 INFO Queuing job for member 11...
2026-07-14 09:35:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:35:51 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-14 09:35:52 INFO Found: ['5154570']
2026-07-14 09:35:57 INFO [TGCC-IRENE] Submitted job with ID:['5154570']
2026-07-14 09:35:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:35:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12
2026-07-14 09:35:57 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020600_1_ENS12.nc
2026-07-14 09:35:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-14 09:35:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:35:57 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-14 09:35:57 INFO Queuing job for member 12...
2026-07-14 09:35:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:35:57 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-14 09:35:58 INFO Found: ['5154572']
2026-07-14 09:36:03 INFO [TGCC-IRENE] Submitted job with ID:['5154572']
2026-07-14 09:36:03 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:36:03 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13
2026-07-14 09:36:03 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020600_1_ENS13.nc
2026-07-14 09:38:11 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-14 09:38:11 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:38:11 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-14 09:39:23 INFO Queuing job for member 13...
2026-07-14 09:39:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:39:23 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-14 09:39:24 INFO Found: ['5154584']
2026-07-14 09:39:29 INFO [TGCC-IRENE] Submitted job with ID:['5154584']
2026-07-14 09:39:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:39:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14
2026-07-14 09:39:29 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020600_1_ENS14.nc
2026-07-14 09:39:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-14 09:39:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:39:29 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-14 09:39:29 INFO Queuing job for member 14...
2026-07-14 09:39:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:39:29 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-14 09:39:30 INFO Found: ['5154585']
2026-07-14 09:39:35 INFO [TGCC-IRENE] Submitted job with ID:['5154585']
2026-07-14 09:39:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:39:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15
2026-07-14 09:39:35 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020600_1_ENS15.nc
2026-07-14 09:39:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-14 09:39:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:39:35 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-14 09:39:35 INFO Queuing job for member 15...
2026-07-14 09:39:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:39:35 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-14 09:39:36 INFO Found: ['5154586']
2026-07-14 09:39:41 INFO [TGCC-IRENE] Submitted job with ID:['5154586']
2026-07-14 09:39:41 INFO Checking job status ...
2026-07-14 09:39:41 INFO None 5154558: status RUNNING/PENDING
2026-07-14 09:39:41 INFO None 5154559: status RUNNING/PENDING
2026-07-14 09:39:41 INFO None 5154561: status RUNNING/PENDING
2026-07-14 09:39:41 INFO None 5154562: status RUNNING/PENDING
2026-07-14 09:39:41 INFO None 5154563: status RUNNING/PENDING
2026-07-14 09:39:41 INFO None 5154564: status RUNNING/PENDING
2026-07-14 09:39:41 INFO None 5154565: status RUNNING/PENDING
2026-07-14 09:39:41 INFO None 5154566: status RUNNING/PENDING
2026-07-14 09:39:41 INFO None 5154567: status RUNNING/PENDING
2026-07-14 09:39:41 INFO None 5154568: status RUNNING/PENDING
2026-07-14 09:39:41 INFO None 5154570: status RUNNING/PENDING
2026-07-14 09:39:41 INFO None 5154572: status RUNNING/PENDING
2026-07-14 09:39:41 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:39:41 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:39:41 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:39:41 INFO Jobs still running: ['5154558', '5154559', '5154561', '5154562', '5154563', '5154564', '5154565', '5154566', '5154567', '5154568', '5154570', '5154572', '5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:39:56 INFO None 5154558: status RUNNING/PENDING
2026-07-14 09:39:56 INFO None 5154559: status RUNNING/PENDING
2026-07-14 09:39:56 INFO None 5154561: status RUNNING/PENDING
2026-07-14 09:39:56 INFO None 5154562: status RUNNING/PENDING
2026-07-14 09:39:56 INFO None 5154563: status RUNNING/PENDING
2026-07-14 09:39:56 INFO None 5154564: status RUNNING/PENDING
2026-07-14 09:39:56 INFO None 5154565: status RUNNING/PENDING
2026-07-14 09:39:56 INFO None 5154566: status RUNNING/PENDING
2026-07-14 09:39:56 INFO None 5154567: status RUNNING/PENDING
2026-07-14 09:39:56 INFO None 5154568: status RUNNING/PENDING
2026-07-14 09:39:56 INFO None 5154570: status RUNNING/PENDING
2026-07-14 09:39:56 INFO None 5154572: status RUNNING/PENDING
2026-07-14 09:39:56 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:39:56 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:39:56 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:39:56 INFO Jobs still running: ['5154558', '5154559', '5154561', '5154562', '5154563', '5154564', '5154565', '5154566', '5154567', '5154568', '5154570', '5154572', '5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:40:12 INFO None 5154558: status RUNNING/PENDING
2026-07-14 09:40:12 INFO None 5154559: status RUNNING/PENDING
2026-07-14 09:40:12 INFO None 5154561: status RUNNING/PENDING
2026-07-14 09:40:12 INFO None 5154562: status RUNNING/PENDING
2026-07-14 09:40:12 INFO None 5154563: status RUNNING/PENDING
2026-07-14 09:40:12 INFO None 5154564: status RUNNING/PENDING
2026-07-14 09:40:12 INFO None 5154565: status RUNNING/PENDING
2026-07-14 09:40:12 INFO None 5154566: status RUNNING/PENDING
2026-07-14 09:40:12 INFO None 5154567: status RUNNING/PENDING
2026-07-14 09:40:12 INFO None 5154568: status RUNNING/PENDING
2026-07-14 09:40:12 INFO None 5154570: status RUNNING/PENDING
2026-07-14 09:40:12 INFO None 5154572: status RUNNING/PENDING
2026-07-14 09:40:12 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:40:12 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:40:12 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:40:12 INFO Jobs still running: ['5154558', '5154559', '5154561', '5154562', '5154563', '5154564', '5154565', '5154566', '5154567', '5154568', '5154570', '5154572', '5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:40:27 INFO None 5154558: status RUNNING/PENDING
2026-07-14 09:40:27 INFO None 5154559: status RUNNING/PENDING
2026-07-14 09:40:27 INFO None 5154561: status RUNNING/PENDING
2026-07-14 09:40:27 INFO None 5154562: status RUNNING/PENDING
2026-07-14 09:40:27 INFO None 5154563: status RUNNING/PENDING
2026-07-14 09:40:27 INFO None 5154564: status RUNNING/PENDING
2026-07-14 09:40:27 INFO None 5154565: status RUNNING/PENDING
2026-07-14 09:40:27 INFO None 5154566: status RUNNING/PENDING
2026-07-14 09:40:27 INFO None 5154567: status RUNNING/PENDING
2026-07-14 09:40:27 INFO None 5154568: status RUNNING/PENDING
2026-07-14 09:40:27 INFO None 5154570: status RUNNING/PENDING
2026-07-14 09:40:27 INFO None 5154572: status RUNNING/PENDING
2026-07-14 09:40:27 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:40:27 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:40:27 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:40:27 INFO Jobs still running: ['5154558', '5154559', '5154561', '5154562', '5154563', '5154564', '5154565', '5154566', '5154567', '5154568', '5154570', '5154572', '5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:40:42 INFO None 5154558: status RUNNING/PENDING
2026-07-14 09:40:42 INFO None 5154559: status RUNNING/PENDING
2026-07-14 09:40:42 INFO None 5154561: status RUNNING/PENDING
2026-07-14 09:40:42 INFO None 5154562: status RUNNING/PENDING
2026-07-14 09:40:42 INFO None 5154563: status RUNNING/PENDING
2026-07-14 09:40:42 INFO None 5154564: status RUNNING/PENDING
2026-07-14 09:40:42 INFO None 5154565: status RUNNING/PENDING
2026-07-14 09:40:42 INFO None 5154566: status RUNNING/PENDING
2026-07-14 09:40:42 INFO None 5154567: status RUNNING/PENDING
2026-07-14 09:40:42 INFO None 5154568: status RUNNING/PENDING
2026-07-14 09:40:42 INFO None 5154570: status RUNNING/PENDING
2026-07-14 09:40:42 INFO None 5154572: status RUNNING/PENDING
2026-07-14 09:40:42 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:40:42 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:40:42 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:40:42 INFO Jobs still running: ['5154558', '5154559', '5154561', '5154562', '5154563', '5154564', '5154565', '5154566', '5154567', '5154568', '5154570', '5154572', '5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:40:57 INFO None 5154558: status RUNNING/PENDING
2026-07-14 09:40:58 INFO None 5154559: status RUNNING/PENDING
2026-07-14 09:40:58 INFO None 5154561: status RUNNING/PENDING
2026-07-14 09:40:58 INFO None 5154562: status RUNNING/PENDING
2026-07-14 09:40:58 INFO None 5154563: status RUNNING/PENDING
2026-07-14 09:40:58 INFO None 5154564: status RUNNING/PENDING
2026-07-14 09:40:58 INFO None 5154565: status RUNNING/PENDING
2026-07-14 09:40:58 INFO None 5154566: status RUNNING/PENDING
2026-07-14 09:40:58 INFO None 5154567: status RUNNING/PENDING
2026-07-14 09:40:58 INFO None 5154568: status RUNNING/PENDING
2026-07-14 09:40:58 INFO None 5154570: status RUNNING/PENDING
2026-07-14 09:40:58 INFO None 5154572: status RUNNING/PENDING
2026-07-14 09:40:58 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:40:58 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:40:58 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:40:58 INFO Jobs still running: ['5154558', '5154559', '5154561', '5154562', '5154563', '5154564', '5154565', '5154566', '5154567', '5154568', '5154570', '5154572', '5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:41:13 INFO None 5154558: status RUNNING/PENDING
2026-07-14 09:41:13 INFO None 5154559: status RUNNING/PENDING
2026-07-14 09:41:13 INFO None 5154561: status RUNNING/PENDING
2026-07-14 09:41:13 INFO None 5154562: status RUNNING/PENDING
2026-07-14 09:41:13 INFO None 5154563: status RUNNING/PENDING
2026-07-14 09:41:13 INFO None 5154564: status RUNNING/PENDING
2026-07-14 09:41:13 INFO None 5154565: status RUNNING/PENDING
2026-07-14 09:41:13 INFO None 5154566: status RUNNING/PENDING
2026-07-14 09:41:13 INFO None 5154567: status RUNNING/PENDING
2026-07-14 09:41:13 INFO None 5154568: status RUNNING/PENDING
2026-07-14 09:41:13 INFO None 5154570: status RUNNING/PENDING
2026-07-14 09:41:13 INFO None 5154572: status RUNNING/PENDING
2026-07-14 09:41:13 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:41:13 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:41:13 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:41:13 INFO Jobs still running: ['5154558', '5154559', '5154561', '5154562', '5154563', '5154564', '5154565', '5154566', '5154567', '5154568', '5154570', '5154572', '5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:41:28 INFO None 5154558: status RUNNING/PENDING
2026-07-14 09:43:14 INFO None 5154559: status RUNNING/PENDING
2026-07-14 09:43:14 INFO None 5154561: status RUNNING/PENDING
2026-07-14 09:43:14 INFO None 5154562: status RUNNING/PENDING
2026-07-14 09:43:14 INFO None 5154563: status RUNNING/PENDING
2026-07-14 09:43:14 INFO None 5154564: status RUNNING/PENDING
2026-07-14 09:43:14 INFO None 5154565: status RUNNING/PENDING
2026-07-14 09:43:14 INFO None 5154566: status RUNNING/PENDING
2026-07-14 09:43:14 INFO None 5154567: status RUNNING/PENDING
2026-07-14 09:43:14 INFO None 5154568: status RUNNING/PENDING
2026-07-14 09:43:14 INFO None 5154570: status RUNNING/PENDING
2026-07-14 09:43:14 INFO None 5154572: status RUNNING/PENDING
2026-07-14 09:43:14 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:43:14 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:43:14 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:43:14 INFO Jobs still running: ['5154558', '5154559', '5154561', '5154562', '5154563', '5154564', '5154565', '5154566', '5154567', '5154568', '5154570', '5154572', '5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:43:29 INFO None 5154558: status RUNNING/PENDING
2026-07-14 09:43:29 INFO None 5154559: status RUNNING/PENDING
2026-07-14 09:43:29 INFO None 5154561: status RUNNING/PENDING
2026-07-14 09:43:29 INFO None 5154562: status RUNNING/PENDING
2026-07-14 09:43:29 INFO None 5154563: status RUNNING/PENDING
2026-07-14 09:43:29 INFO None 5154564: status RUNNING/PENDING
2026-07-14 09:43:30 INFO None 5154565: status RUNNING/PENDING
2026-07-14 09:43:30 INFO None 5154566: status RUNNING/PENDING
2026-07-14 09:43:30 INFO None 5154567: status RUNNING/PENDING
2026-07-14 09:43:30 INFO None 5154568: status RUNNING/PENDING
2026-07-14 09:43:30 INFO None 5154570: status RUNNING/PENDING
2026-07-14 09:43:30 INFO None 5154572: status RUNNING/PENDING
2026-07-14 09:43:30 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:43:30 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:43:30 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:43:30 INFO Jobs still running: ['5154558', '5154559', '5154561', '5154562', '5154563', '5154564', '5154565', '5154566', '5154567', '5154568', '5154570', '5154572', '5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:43:45 INFO None 5154558: status RUNNING/PENDING
2026-07-14 09:43:45 INFO None 5154559: status RUNNING/PENDING
2026-07-14 09:43:45 INFO None 5154561: status RUNNING/PENDING
2026-07-14 09:43:45 INFO None 5154562: status RUNNING/PENDING
2026-07-14 09:43:45 INFO None 5154563: status RUNNING/PENDING
2026-07-14 09:43:45 INFO None 5154564: status RUNNING/PENDING
2026-07-14 09:43:45 INFO None 5154565: status RUNNING/PENDING
2026-07-14 09:43:45 INFO None 5154566: status RUNNING/PENDING
2026-07-14 09:43:45 INFO None 5154567: status RUNNING/PENDING
2026-07-14 09:43:45 INFO None 5154568: status RUNNING/PENDING
2026-07-14 09:43:45 INFO None 5154570: status RUNNING/PENDING
2026-07-14 09:43:45 INFO None 5154572: status RUNNING/PENDING
2026-07-14 09:43:45 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:43:45 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:43:45 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:43:45 INFO Jobs still running: ['5154558', '5154559', '5154561', '5154562', '5154563', '5154564', '5154565', '5154566', '5154567', '5154568', '5154570', '5154572', '5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:44:00 INFO None 5154558: status RUNNING/PENDING
2026-07-14 09:44:00 INFO None 5154559: status RUNNING/PENDING
2026-07-14 09:44:00 INFO None 5154561: status RUNNING/PENDING
2026-07-14 09:44:00 INFO None 5154562: status RUNNING/PENDING
2026-07-14 09:44:00 INFO None 5154563: status RUNNING/PENDING
2026-07-14 09:44:00 INFO None 5154564: status RUNNING/PENDING
2026-07-14 09:44:00 INFO None 5154565: status RUNNING/PENDING
2026-07-14 09:44:00 INFO None 5154566: status RUNNING/PENDING
2026-07-14 09:44:00 INFO None 5154567: status RUNNING/PENDING
2026-07-14 09:44:00 INFO None 5154568: status RUNNING/PENDING
2026-07-14 09:44:00 INFO None 5154570: status RUNNING/PENDING
2026-07-14 09:44:00 INFO None 5154572: status RUNNING/PENDING
2026-07-14 09:44:00 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:44:00 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:44:00 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:44:00 INFO Jobs still running: ['5154558', '5154559', '5154561', '5154562', '5154563', '5154564', '5154565', '5154566', '5154567', '5154568', '5154570', '5154572', '5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:44:15 INFO None 5154558: status RUNNING/PENDING
2026-07-14 09:44:15 INFO None 5154559: status RUNNING/PENDING
2026-07-14 09:44:15 INFO None 5154561: status RUNNING/PENDING
2026-07-14 09:44:15 INFO None 5154562: status RUNNING/PENDING
2026-07-14 09:44:15 INFO None 5154563: status RUNNING/PENDING
2026-07-14 09:44:15 INFO None 5154564: status RUNNING/PENDING
2026-07-14 09:44:15 INFO None 5154565: status RUNNING/PENDING
2026-07-14 09:44:15 INFO None 5154566: status RUNNING/PENDING
2026-07-14 09:44:15 INFO None 5154567: status RUNNING/PENDING
2026-07-14 09:44:15 INFO None 5154568: status RUNNING/PENDING
2026-07-14 09:44:16 INFO None 5154570: status RUNNING/PENDING
2026-07-14 09:44:16 INFO None 5154572: status RUNNING/PENDING
2026-07-14 09:44:16 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:44:16 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:44:16 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:44:16 INFO Jobs still running: ['5154558', '5154559', '5154561', '5154562', '5154563', '5154564', '5154565', '5154566', '5154567', '5154568', '5154570', '5154572', '5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:44:31 INFO None 5154558: status RUNNING/PENDING
2026-07-14 09:44:31 INFO None 5154559: status RUNNING/PENDING
2026-07-14 09:44:31 INFO None 5154561: status RUNNING/PENDING
2026-07-14 09:44:31 INFO None 5154562: status RUNNING/PENDING
2026-07-14 09:44:31 INFO None 5154563: status RUNNING/PENDING
2026-07-14 09:44:31 INFO None 5154564: status RUNNING/PENDING
2026-07-14 09:44:31 INFO None 5154565: status RUNNING/PENDING
2026-07-14 09:44:31 INFO None 5154566: status RUNNING/PENDING
2026-07-14 09:44:31 INFO None 5154567: status RUNNING/PENDING
2026-07-14 09:44:31 INFO None 5154568: status RUNNING/PENDING
2026-07-14 09:44:31 INFO None 5154570: status RUNNING/PENDING
2026-07-14 09:44:31 INFO None 5154572: status RUNNING/PENDING
2026-07-14 09:44:31 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:44:31 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:44:31 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:44:31 INFO Jobs still running: ['5154558', '5154559', '5154561', '5154562', '5154563', '5154564', '5154565', '5154566', '5154567', '5154568', '5154570', '5154572', '5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:44:46 INFO None 5154558: status RUNNING/PENDING
2026-07-14 09:44:46 INFO None 5154559: status RUNNING/PENDING
2026-07-14 09:44:46 INFO None 5154561: status RUNNING/PENDING
2026-07-14 09:44:46 INFO None 5154562: status RUNNING/PENDING
2026-07-14 09:44:46 INFO None 5154563: status RUNNING/PENDING
2026-07-14 09:44:46 INFO None 5154564: status RUNNING/PENDING
2026-07-14 09:44:46 INFO None 5154565: status RUNNING/PENDING
2026-07-14 09:44:46 INFO None 5154566: status RUNNING/PENDING
2026-07-14 09:44:46 INFO None 5154567: status RUNNING/PENDING
2026-07-14 09:44:46 INFO None 5154568: status RUNNING/PENDING
2026-07-14 09:44:46 INFO None 5154570: status RUNNING/PENDING
2026-07-14 09:44:46 INFO None 5154572: status RUNNING/PENDING
2026-07-14 09:44:46 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:44:46 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:44:46 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:44:46 INFO Jobs still running: ['5154558', '5154559', '5154561', '5154562', '5154563', '5154564', '5154565', '5154566', '5154567', '5154568', '5154570', '5154572', '5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:45:01 INFO None 5154558: status RUNNING/PENDING
2026-07-14 09:45:01 INFO None 5154559: status RUNNING/PENDING
2026-07-14 09:45:01 INFO None 5154561: status RUNNING/PENDING
2026-07-14 09:45:01 INFO None 5154562: status RUNNING/PENDING
2026-07-14 09:45:01 INFO None 5154563: status RUNNING/PENDING
2026-07-14 09:45:01 INFO None 5154564: status RUNNING/PENDING
2026-07-14 09:45:01 INFO None 5154565: status RUNNING/PENDING
2026-07-14 09:45:01 INFO None 5154566: status RUNNING/PENDING
2026-07-14 09:45:01 INFO None 5154567: status RUNNING/PENDING
2026-07-14 09:45:01 INFO None 5154568: status RUNNING/PENDING
2026-07-14 09:45:01 INFO None 5154570: status RUNNING/PENDING
2026-07-14 09:45:01 INFO None 5154572: status RUNNING/PENDING
2026-07-14 09:45:01 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:45:01 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:45:02 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:45:02 INFO Jobs still running: ['5154558', '5154559', '5154561', '5154562', '5154563', '5154564', '5154565', '5154566', '5154567', '5154568', '5154570', '5154572', '5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:45:17 INFO None 5154558: status RUNNING/PENDING
2026-07-14 09:45:17 INFO None 5154559: status RUNNING/PENDING
2026-07-14 09:45:17 INFO None 5154561: status RUNNING/PENDING
2026-07-14 09:45:17 INFO None 5154562: status RUNNING/PENDING
2026-07-14 09:45:17 INFO None 5154563: status RUNNING/PENDING
2026-07-14 09:45:17 INFO None 5154564: status RUNNING/PENDING
2026-07-14 09:45:17 INFO None 5154565: status RUNNING/PENDING
2026-07-14 09:45:17 INFO None 5154566: status RUNNING/PENDING
2026-07-14 09:45:17 INFO None 5154567: status RUNNING/PENDING
2026-07-14 09:45:17 INFO None 5154568: status RUNNING/PENDING
2026-07-14 09:45:17 INFO None 5154570: status RUNNING/PENDING
2026-07-14 09:45:17 INFO None 5154572: status RUNNING/PENDING
2026-07-14 09:45:17 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:45:17 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:45:17 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:45:17 INFO Jobs still running: ['5154558', '5154559', '5154561', '5154562', '5154563', '5154564', '5154565', '5154566', '5154567', '5154568', '5154570', '5154572', '5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:45:32 INFO None 5154558: status RUNNING/PENDING
2026-07-14 09:45:32 INFO None 5154559: status RUNNING/PENDING
2026-07-14 09:45:32 INFO None 5154561: status RUNNING/PENDING
2026-07-14 09:45:32 INFO None 5154562: status RUNNING/PENDING
2026-07-14 09:45:32 INFO None 5154563: status RUNNING/PENDING
2026-07-14 09:45:32 INFO None 5154564: status RUNNING/PENDING
2026-07-14 09:45:32 INFO None 5154565: status RUNNING/PENDING
2026-07-14 09:45:32 INFO None 5154566: status RUNNING/PENDING
2026-07-14 09:45:32 INFO None 5154567: status RUNNING/PENDING
2026-07-14 09:45:32 INFO None 5154568: status RUNNING/PENDING
2026-07-14 09:45:32 INFO None 5154570: status RUNNING/PENDING
2026-07-14 09:45:32 INFO None 5154572: status RUNNING/PENDING
2026-07-14 09:45:32 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:45:32 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:45:32 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:45:32 INFO Jobs still running: ['5154558', '5154559', '5154561', '5154562', '5154563', '5154564', '5154565', '5154566', '5154567', '5154568', '5154570', '5154572', '5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:45:47 INFO None 5154558: status RUNNING/PENDING
2026-07-14 09:45:47 INFO None 5154559: status RUNNING/PENDING
2026-07-14 09:45:47 INFO None 5154561: status RUNNING/PENDING
2026-07-14 09:45:47 INFO None 5154562: status RUNNING/PENDING
2026-07-14 09:45:47 INFO None 5154563: status RUNNING/PENDING
2026-07-14 09:45:47 INFO None 5154564: status RUNNING/PENDING
2026-07-14 09:45:47 INFO None 5154565: status RUNNING/PENDING
2026-07-14 09:45:47 INFO None 5154566: status RUNNING/PENDING
2026-07-14 09:45:47 INFO None 5154567: status RUNNING/PENDING
2026-07-14 09:45:47 INFO None 5154568: status RUNNING/PENDING
2026-07-14 09:45:47 INFO None 5154570: status RUNNING/PENDING
2026-07-14 09:45:47 INFO None 5154572: status RUNNING/PENDING
2026-07-14 09:45:47 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:45:47 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:45:47 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:45:47 INFO Jobs still running: ['5154558', '5154559', '5154561', '5154562', '5154563', '5154564', '5154565', '5154566', '5154567', '5154568', '5154570', '5154572', '5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:46:02 INFO None 5154558: status FINISHED
2026-07-14 09:46:02 INFO None 5154559: status FINISHED
2026-07-14 09:46:02 INFO None 5154561: status RUNNING/PENDING
2026-07-14 09:46:03 INFO None 5154562: status RUNNING/PENDING
2026-07-14 09:46:03 INFO None 5154563: status RUNNING/PENDING
2026-07-14 09:46:03 INFO None 5154564: status RUNNING/PENDING
2026-07-14 09:46:03 INFO None 5154565: status RUNNING/PENDING
2026-07-14 09:46:03 INFO None 5154566: status RUNNING/PENDING
2026-07-14 09:46:03 INFO None 5154567: status RUNNING/PENDING
2026-07-14 09:46:03 INFO None 5154568: status RUNNING/PENDING
2026-07-14 09:46:03 INFO None 5154570: status RUNNING/PENDING
2026-07-14 09:46:03 INFO None 5154572: status RUNNING/PENDING
2026-07-14 09:46:03 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:46:03 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:46:03 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:46:03 INFO Jobs still running: ['5154561', '5154562', '5154563', '5154564', '5154565', '5154566', '5154567', '5154568', '5154570', '5154572', '5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:46:18 INFO None 5154558: status FINISHED
2026-07-14 09:48:14 INFO None 5154559: status FINISHED
2026-07-14 09:48:14 INFO None 5154561: status FINISHED
2026-07-14 09:48:14 INFO None 5154562: status FINISHED
2026-07-14 09:48:14 INFO None 5154563: status FINISHED
2026-07-14 09:48:14 INFO None 5154564: status FINISHED
2026-07-14 09:48:14 INFO None 5154565: status FINISHED
2026-07-14 09:48:14 INFO None 5154566: status FINISHED
2026-07-14 09:48:14 INFO None 5154567: status FINISHED
2026-07-14 09:48:14 INFO None 5154568: status FINISHED
2026-07-14 09:48:14 INFO None 5154570: status FINISHED
2026-07-14 09:48:14 INFO None 5154572: status FINISHED
2026-07-14 09:48:14 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:48:14 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:48:14 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:48:14 INFO Jobs still running: ['5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:48:29 INFO None 5154558: status FINISHED
2026-07-14 09:48:29 INFO None 5154559: status FINISHED
2026-07-14 09:48:29 INFO None 5154561: status FINISHED
2026-07-14 09:48:29 INFO None 5154562: status FINISHED
2026-07-14 09:48:29 INFO None 5154563: status FINISHED
2026-07-14 09:48:29 INFO None 5154564: status FINISHED
2026-07-14 09:48:29 INFO None 5154565: status FINISHED
2026-07-14 09:48:29 INFO None 5154566: status FINISHED
2026-07-14 09:48:29 INFO None 5154567: status FINISHED
2026-07-14 09:48:29 INFO None 5154568: status FINISHED
2026-07-14 09:48:29 INFO None 5154570: status FINISHED
2026-07-14 09:48:29 INFO None 5154572: status FINISHED
2026-07-14 09:48:29 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:48:29 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:48:29 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:48:29 INFO Jobs still running: ['5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:48:44 INFO None 5154558: status FINISHED
2026-07-14 09:48:44 INFO None 5154559: status FINISHED
2026-07-14 09:48:45 INFO None 5154561: status FINISHED
2026-07-14 09:48:45 INFO None 5154562: status FINISHED
2026-07-14 09:48:45 INFO None 5154563: status FINISHED
2026-07-14 09:48:45 INFO None 5154564: status FINISHED
2026-07-14 09:48:45 INFO None 5154565: status FINISHED
2026-07-14 09:48:45 INFO None 5154566: status FINISHED
2026-07-14 09:48:45 INFO None 5154567: status FINISHED
2026-07-14 09:48:45 INFO None 5154568: status FINISHED
2026-07-14 09:48:45 INFO None 5154570: status FINISHED
2026-07-14 09:48:45 INFO None 5154572: status FINISHED
2026-07-14 09:48:45 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:48:45 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:48:45 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:48:45 INFO Jobs still running: ['5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:49:00 INFO None 5154558: status FINISHED
2026-07-14 09:49:00 INFO None 5154559: status FINISHED
2026-07-14 09:49:00 INFO None 5154561: status FINISHED
2026-07-14 09:49:00 INFO None 5154562: status FINISHED
2026-07-14 09:49:00 INFO None 5154563: status FINISHED
2026-07-14 09:49:00 INFO None 5154564: status FINISHED
2026-07-14 09:49:00 INFO None 5154565: status FINISHED
2026-07-14 09:49:00 INFO None 5154566: status FINISHED
2026-07-14 09:49:00 INFO None 5154567: status FINISHED
2026-07-14 09:49:00 INFO None 5154568: status FINISHED
2026-07-14 09:49:00 INFO None 5154570: status FINISHED
2026-07-14 09:49:00 INFO None 5154572: status FINISHED
2026-07-14 09:49:00 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:49:00 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:49:00 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:49:00 INFO Jobs still running: ['5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:49:15 INFO None 5154558: status FINISHED
2026-07-14 09:49:15 INFO None 5154559: status FINISHED
2026-07-14 09:49:15 INFO None 5154561: status FINISHED
2026-07-14 09:49:15 INFO None 5154562: status FINISHED
2026-07-14 09:49:15 INFO None 5154563: status FINISHED
2026-07-14 09:49:15 INFO None 5154564: status FINISHED
2026-07-14 09:49:15 INFO None 5154565: status FINISHED
2026-07-14 09:49:15 INFO None 5154566: status FINISHED
2026-07-14 09:49:15 INFO None 5154567: status FINISHED
2026-07-14 09:49:15 INFO None 5154568: status FINISHED
2026-07-14 09:49:15 INFO None 5154570: status FINISHED
2026-07-14 09:49:15 INFO None 5154572: status FINISHED
2026-07-14 09:49:15 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:49:15 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:49:15 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:49:15 INFO Jobs still running: ['5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:49:30 INFO None 5154558: status FINISHED
2026-07-14 09:49:30 INFO None 5154559: status FINISHED
2026-07-14 09:49:30 INFO None 5154561: status FINISHED
2026-07-14 09:49:30 INFO None 5154562: status FINISHED
2026-07-14 09:49:30 INFO None 5154563: status FINISHED
2026-07-14 09:49:30 INFO None 5154564: status FINISHED
2026-07-14 09:49:30 INFO None 5154565: status FINISHED
2026-07-14 09:49:30 INFO None 5154566: status FINISHED
2026-07-14 09:49:30 INFO None 5154567: status FINISHED
2026-07-14 09:49:31 INFO None 5154568: status FINISHED
2026-07-14 09:49:31 INFO None 5154570: status FINISHED
2026-07-14 09:49:31 INFO None 5154572: status FINISHED
2026-07-14 09:49:31 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:49:31 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:49:31 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:49:31 INFO Jobs still running: ['5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:49:46 INFO None 5154558: status FINISHED
2026-07-14 09:49:46 INFO None 5154559: status FINISHED
2026-07-14 09:49:46 INFO None 5154561: status FINISHED
2026-07-14 09:49:46 INFO None 5154562: status FINISHED
2026-07-14 09:49:46 INFO None 5154563: status FINISHED
2026-07-14 09:49:46 INFO None 5154564: status FINISHED
2026-07-14 09:49:46 INFO None 5154565: status FINISHED
2026-07-14 09:49:46 INFO None 5154566: status FINISHED
2026-07-14 09:49:46 INFO None 5154567: status FINISHED
2026-07-14 09:49:46 INFO None 5154568: status FINISHED
2026-07-14 09:49:46 INFO None 5154570: status FINISHED
2026-07-14 09:49:46 INFO None 5154572: status FINISHED
2026-07-14 09:49:46 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:49:46 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:49:46 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:49:46 INFO Jobs still running: ['5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:50:01 INFO None 5154558: status FINISHED
2026-07-14 09:50:01 INFO None 5154559: status FINISHED
2026-07-14 09:50:01 INFO None 5154561: status FINISHED
2026-07-14 09:50:01 INFO None 5154562: status FINISHED
2026-07-14 09:50:01 INFO None 5154563: status FINISHED
2026-07-14 09:50:01 INFO None 5154564: status FINISHED
2026-07-14 09:50:01 INFO None 5154565: status FINISHED
2026-07-14 09:50:01 INFO None 5154566: status FINISHED
2026-07-14 09:50:01 INFO None 5154567: status FINISHED
2026-07-14 09:50:01 INFO None 5154568: status FINISHED
2026-07-14 09:50:01 INFO None 5154570: status FINISHED
2026-07-14 09:50:01 INFO None 5154572: status FINISHED
2026-07-14 09:50:01 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:50:01 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:50:01 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:50:01 INFO Jobs still running: ['5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:50:16 INFO None 5154558: status FINISHED
2026-07-14 09:50:16 INFO None 5154559: status FINISHED
2026-07-14 09:50:16 INFO None 5154561: status FINISHED
2026-07-14 09:50:16 INFO None 5154562: status FINISHED
2026-07-14 09:50:16 INFO None 5154563: status FINISHED
2026-07-14 09:50:16 INFO None 5154564: status FINISHED
2026-07-14 09:50:16 INFO None 5154565: status FINISHED
2026-07-14 09:50:16 INFO None 5154566: status FINISHED
2026-07-14 09:50:16 INFO None 5154567: status FINISHED
2026-07-14 09:50:16 INFO None 5154568: status FINISHED
2026-07-14 09:50:16 INFO None 5154570: status FINISHED
2026-07-14 09:50:16 INFO None 5154572: status FINISHED
2026-07-14 09:50:16 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:50:16 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:50:16 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:50:16 INFO Jobs still running: ['5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:50:32 INFO None 5154558: status FINISHED
2026-07-14 09:50:32 INFO None 5154559: status FINISHED
2026-07-14 09:50:32 INFO None 5154561: status FINISHED
2026-07-14 09:50:32 INFO None 5154562: status FINISHED
2026-07-14 09:50:32 INFO None 5154563: status FINISHED
2026-07-14 09:50:32 INFO None 5154564: status FINISHED
2026-07-14 09:50:32 INFO None 5154565: status FINISHED
2026-07-14 09:50:32 INFO None 5154566: status FINISHED
2026-07-14 09:50:32 INFO None 5154567: status FINISHED
2026-07-14 09:50:32 INFO None 5154568: status FINISHED
2026-07-14 09:50:32 INFO None 5154570: status FINISHED
2026-07-14 09:50:32 INFO None 5154572: status FINISHED
2026-07-14 09:50:32 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:50:32 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:50:32 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:50:32 INFO Jobs still running: ['5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:50:47 INFO None 5154558: status FINISHED
2026-07-14 09:50:47 INFO None 5154559: status FINISHED
2026-07-14 09:50:47 INFO None 5154561: status FINISHED
2026-07-14 09:50:47 INFO None 5154562: status FINISHED
2026-07-14 09:50:47 INFO None 5154563: status FINISHED
2026-07-14 09:50:47 INFO None 5154564: status FINISHED
2026-07-14 09:50:47 INFO None 5154565: status FINISHED
2026-07-14 09:50:47 INFO None 5154566: status FINISHED
2026-07-14 09:50:47 INFO None 5154567: status FINISHED
2026-07-14 09:50:47 INFO None 5154568: status FINISHED
2026-07-14 09:50:47 INFO None 5154570: status FINISHED
2026-07-14 09:50:47 INFO None 5154572: status FINISHED
2026-07-14 09:50:47 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:50:47 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:50:47 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:50:47 INFO Jobs still running: ['5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:51:02 INFO None 5154558: status FINISHED
2026-07-14 09:51:02 INFO None 5154559: status FINISHED
2026-07-14 09:51:02 INFO None 5154561: status FINISHED
2026-07-14 09:51:02 INFO None 5154562: status FINISHED
2026-07-14 09:51:02 INFO None 5154563: status FINISHED
2026-07-14 09:51:02 INFO None 5154564: status FINISHED
2026-07-14 09:51:02 INFO None 5154565: status FINISHED
2026-07-14 09:51:02 INFO None 5154566: status FINISHED
2026-07-14 09:51:02 INFO None 5154567: status FINISHED
2026-07-14 09:51:02 INFO None 5154568: status FINISHED
2026-07-14 09:51:02 INFO None 5154570: status FINISHED
2026-07-14 09:51:02 INFO None 5154572: status FINISHED
2026-07-14 09:51:02 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:51:02 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:51:02 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:51:02 INFO Jobs still running: ['5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:51:17 INFO None 5154558: status FINISHED
2026-07-14 09:51:17 INFO None 5154559: status FINISHED
2026-07-14 09:51:17 INFO None 5154561: status FINISHED
2026-07-14 09:51:17 INFO None 5154562: status FINISHED
2026-07-14 09:51:17 INFO None 5154563: status FINISHED
2026-07-14 09:51:17 INFO None 5154564: status FINISHED
2026-07-14 09:51:17 INFO None 5154565: status FINISHED
2026-07-14 09:51:18 INFO None 5154566: status FINISHED
2026-07-14 09:51:18 INFO None 5154567: status FINISHED
2026-07-14 09:51:18 INFO None 5154568: status FINISHED
2026-07-14 09:51:18 INFO None 5154570: status FINISHED
2026-07-14 09:51:18 INFO None 5154572: status FINISHED
2026-07-14 09:51:18 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:51:18 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:51:18 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:51:18 INFO Jobs still running: ['5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:51:33 INFO None 5154558: status FINISHED
2026-07-14 09:53:31 INFO None 5154559: status FINISHED
2026-07-14 09:53:31 INFO None 5154561: status FINISHED
2026-07-14 09:53:31 INFO None 5154562: status FINISHED
2026-07-14 09:53:31 INFO None 5154563: status FINISHED
2026-07-14 09:53:31 INFO None 5154564: status FINISHED
2026-07-14 09:53:31 INFO None 5154565: status FINISHED
2026-07-14 09:53:31 INFO None 5154566: status FINISHED
2026-07-14 09:53:31 INFO None 5154567: status FINISHED
2026-07-14 09:53:31 INFO None 5154568: status FINISHED
2026-07-14 09:53:31 INFO None 5154570: status FINISHED
2026-07-14 09:53:31 INFO None 5154572: status FINISHED
2026-07-14 09:53:31 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:53:31 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:53:31 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:53:31 INFO Jobs still running: ['5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:53:46 INFO None 5154558: status FINISHED
2026-07-14 09:53:46 INFO None 5154559: status FINISHED
2026-07-14 09:53:46 INFO None 5154561: status FINISHED
2026-07-14 09:53:46 INFO None 5154562: status FINISHED
2026-07-14 09:53:46 INFO None 5154563: status FINISHED
2026-07-14 09:53:46 INFO None 5154564: status FINISHED
2026-07-14 09:53:46 INFO None 5154565: status FINISHED
2026-07-14 09:53:46 INFO None 5154566: status FINISHED
2026-07-14 09:53:46 INFO None 5154567: status FINISHED
2026-07-14 09:53:46 INFO None 5154568: status FINISHED
2026-07-14 09:53:46 INFO None 5154570: status FINISHED
2026-07-14 09:53:46 INFO None 5154572: status FINISHED
2026-07-14 09:53:46 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:53:46 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:53:46 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:53:46 INFO Jobs still running: ['5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:54:01 INFO None 5154558: status FINISHED
2026-07-14 09:54:02 INFO None 5154559: status FINISHED
2026-07-14 09:54:02 INFO None 5154561: status FINISHED
2026-07-14 09:54:02 INFO None 5154562: status FINISHED
2026-07-14 09:54:02 INFO None 5154563: status FINISHED
2026-07-14 09:54:02 INFO None 5154564: status FINISHED
2026-07-14 09:54:02 INFO None 5154565: status FINISHED
2026-07-14 09:54:02 INFO None 5154566: status FINISHED
2026-07-14 09:54:02 INFO None 5154567: status FINISHED
2026-07-14 09:54:02 INFO None 5154568: status FINISHED
2026-07-14 09:54:02 INFO None 5154570: status FINISHED
2026-07-14 09:54:02 INFO None 5154572: status FINISHED
2026-07-14 09:54:02 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:54:02 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:54:02 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:54:02 INFO Jobs still running: ['5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:54:17 INFO None 5154558: status FINISHED
2026-07-14 09:54:17 INFO None 5154559: status FINISHED
2026-07-14 09:54:17 INFO None 5154561: status FINISHED
2026-07-14 09:54:17 INFO None 5154562: status FINISHED
2026-07-14 09:54:17 INFO None 5154563: status FINISHED
2026-07-14 09:54:17 INFO None 5154564: status FINISHED
2026-07-14 09:54:17 INFO None 5154565: status FINISHED
2026-07-14 09:54:17 INFO None 5154566: status FINISHED
2026-07-14 09:54:17 INFO None 5154567: status FINISHED
2026-07-14 09:54:17 INFO None 5154568: status FINISHED
2026-07-14 09:54:17 INFO None 5154570: status FINISHED
2026-07-14 09:54:17 INFO None 5154572: status FINISHED
2026-07-14 09:54:17 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:54:17 INFO None 5154585: status RUNNING/PENDING
2026-07-14 09:54:17 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:54:17 INFO Jobs still running: ['5154584', '5154585', '5154586']. Waiting...
2026-07-14 09:54:32 INFO None 5154558: status FINISHED
2026-07-14 09:54:32 INFO None 5154559: status FINISHED
2026-07-14 09:54:32 INFO None 5154561: status FINISHED
2026-07-14 09:54:32 INFO None 5154562: status FINISHED
2026-07-14 09:54:32 INFO None 5154563: status FINISHED
2026-07-14 09:54:32 INFO None 5154564: status FINISHED
2026-07-14 09:54:32 INFO None 5154565: status FINISHED
2026-07-14 09:54:32 INFO None 5154566: status FINISHED
2026-07-14 09:54:32 INFO None 5154567: status FINISHED
2026-07-14 09:54:32 INFO None 5154568: status FINISHED
2026-07-14 09:54:32 INFO None 5154570: status FINISHED
2026-07-14 09:54:32 INFO None 5154572: status FINISHED
2026-07-14 09:54:32 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:54:32 INFO None 5154585: status FINISHED
2026-07-14 09:54:32 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:54:32 INFO Jobs still running: ['5154584', '5154586']. Waiting...
2026-07-14 09:54:47 INFO None 5154558: status FINISHED
2026-07-14 09:54:47 INFO None 5154559: status FINISHED
2026-07-14 09:54:47 INFO None 5154561: status FINISHED
2026-07-14 09:54:47 INFO None 5154562: status FINISHED
2026-07-14 09:54:47 INFO None 5154563: status FINISHED
2026-07-14 09:54:47 INFO None 5154564: status FINISHED
2026-07-14 09:54:47 INFO None 5154565: status FINISHED
2026-07-14 09:54:47 INFO None 5154566: status FINISHED
2026-07-14 09:54:48 INFO None 5154567: status FINISHED
2026-07-14 09:54:48 INFO None 5154568: status FINISHED
2026-07-14 09:54:48 INFO None 5154570: status FINISHED
2026-07-14 09:54:48 INFO None 5154572: status FINISHED
2026-07-14 09:54:48 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:54:48 INFO None 5154585: status FINISHED
2026-07-14 09:54:48 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:54:48 INFO Jobs still running: ['5154584', '5154586']. Waiting...
2026-07-14 09:55:03 INFO None 5154558: status FINISHED
2026-07-14 09:55:03 INFO None 5154559: status FINISHED
2026-07-14 09:55:03 INFO None 5154561: status FINISHED
2026-07-14 09:55:03 INFO None 5154562: status FINISHED
2026-07-14 09:55:03 INFO None 5154563: status FINISHED
2026-07-14 09:55:03 INFO None 5154564: status FINISHED
2026-07-14 09:55:03 INFO None 5154565: status FINISHED
2026-07-14 09:55:03 INFO None 5154566: status FINISHED
2026-07-14 09:55:03 INFO None 5154567: status FINISHED
2026-07-14 09:55:03 INFO None 5154568: status FINISHED
2026-07-14 09:55:03 INFO None 5154570: status FINISHED
2026-07-14 09:55:03 INFO None 5154572: status FINISHED
2026-07-14 09:55:03 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:55:03 INFO None 5154585: status FINISHED
2026-07-14 09:55:03 INFO None 5154586: status RUNNING/PENDING
2026-07-14 09:55:03 INFO Jobs still running: ['5154584', '5154586']. Waiting...
2026-07-14 09:55:18 INFO None 5154558: status FINISHED
2026-07-14 09:55:18 INFO None 5154559: status FINISHED
2026-07-14 09:55:18 INFO None 5154561: status FINISHED
2026-07-14 09:55:18 INFO None 5154562: status FINISHED
2026-07-14 09:55:18 INFO None 5154563: status FINISHED
2026-07-14 09:55:18 INFO None 5154564: status FINISHED
2026-07-14 09:55:18 INFO None 5154565: status FINISHED
2026-07-14 09:55:18 INFO None 5154566: status FINISHED
2026-07-14 09:55:18 INFO None 5154567: status FINISHED
2026-07-14 09:55:18 INFO None 5154568: status FINISHED
2026-07-14 09:55:18 INFO None 5154570: status FINISHED
2026-07-14 09:55:18 INFO None 5154572: status FINISHED
2026-07-14 09:55:18 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:55:18 INFO None 5154585: status FINISHED
2026-07-14 09:55:18 INFO None 5154586: status FINISHED
2026-07-14 09:55:18 INFO Jobs still running: ['5154584']. Waiting...
2026-07-14 09:55:33 INFO None 5154558: status FINISHED
2026-07-14 09:55:33 INFO None 5154559: status FINISHED
2026-07-14 09:55:33 INFO None 5154561: status FINISHED
2026-07-14 09:55:33 INFO None 5154562: status FINISHED
2026-07-14 09:55:33 INFO None 5154563: status FINISHED
2026-07-14 09:55:33 INFO None 5154564: status FINISHED
2026-07-14 09:55:33 INFO None 5154565: status FINISHED
2026-07-14 09:55:33 INFO None 5154566: status FINISHED
2026-07-14 09:55:33 INFO None 5154567: status FINISHED
2026-07-14 09:55:33 INFO None 5154568: status FINISHED
2026-07-14 09:55:33 INFO None 5154570: status FINISHED
2026-07-14 09:55:33 INFO None 5154572: status FINISHED
2026-07-14 09:55:33 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:55:33 INFO None 5154585: status FINISHED
2026-07-14 09:55:33 INFO None 5154586: status FINISHED
2026-07-14 09:55:33 INFO Jobs still running: ['5154584']. Waiting...
2026-07-14 09:55:48 INFO None 5154558: status FINISHED
2026-07-14 09:55:48 INFO None 5154559: status FINISHED
2026-07-14 09:55:48 INFO None 5154561: status FINISHED
2026-07-14 09:55:49 INFO None 5154562: status FINISHED
2026-07-14 09:55:49 INFO None 5154563: status FINISHED
2026-07-14 09:55:49 INFO None 5154564: status FINISHED
2026-07-14 09:55:49 INFO None 5154565: status FINISHED
2026-07-14 09:55:49 INFO None 5154566: status FINISHED
2026-07-14 09:55:49 INFO None 5154567: status FINISHED
2026-07-14 09:55:49 INFO None 5154568: status FINISHED
2026-07-14 09:55:49 INFO None 5154570: status FINISHED
2026-07-14 09:55:49 INFO None 5154572: status FINISHED
2026-07-14 09:55:49 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:55:49 INFO None 5154585: status FINISHED
2026-07-14 09:55:49 INFO None 5154586: status FINISHED
2026-07-14 09:55:49 INFO Jobs still running: ['5154584']. Waiting...
2026-07-14 09:56:04 INFO None 5154558: status FINISHED
2026-07-14 09:56:04 INFO None 5154559: status FINISHED
2026-07-14 09:56:04 INFO None 5154561: status FINISHED
2026-07-14 09:56:04 INFO None 5154562: status FINISHED
2026-07-14 09:56:04 INFO None 5154563: status FINISHED
2026-07-14 09:56:04 INFO None 5154564: status FINISHED
2026-07-14 09:56:04 INFO None 5154565: status FINISHED
2026-07-14 09:56:04 INFO None 5154566: status FINISHED
2026-07-14 09:56:04 INFO None 5154567: status FINISHED
2026-07-14 09:56:04 INFO None 5154568: status FINISHED
2026-07-14 09:56:04 INFO None 5154570: status FINISHED
2026-07-14 09:56:04 INFO None 5154572: status FINISHED
2026-07-14 09:56:04 INFO None 5154584: status RUNNING/PENDING
2026-07-14 09:56:04 INFO None 5154585: status FINISHED
2026-07-14 09:56:04 INFO None 5154586: status FINISHED
2026-07-14 09:56:04 INFO Jobs still running: ['5154584']. Waiting...
2026-07-14 09:56:19 INFO None 5154558: status FINISHED
2026-07-14 09:58:17 INFO None 5154559: status FINISHED (not in squeue)
2026-07-14 09:58:17 INFO None 5154561: status FINISHED (not in squeue)
2026-07-14 09:58:17 INFO None 5154562: status FINISHED (not in squeue)
2026-07-14 09:58:17 INFO None 5154563: status FINISHED (not in squeue)
2026-07-14 09:58:17 INFO None 5154564: status FINISHED (not in squeue)
2026-07-14 09:58:17 INFO None 5154565: status FINISHED (not in squeue)
2026-07-14 09:58:17 INFO None 5154566: status FINISHED (not in squeue)
2026-07-14 09:58:17 INFO None 5154567: status FINISHED (not in squeue)
2026-07-14 09:58:17 INFO None 5154568: status FINISHED (not in squeue)
2026-07-14 09:58:17 INFO None 5154570: status FINISHED (not in squeue)
2026-07-14 09:58:17 INFO None 5154572: status FINISHED (not in squeue)
2026-07-14 09:58:17 INFO None 5154584: status FINISHED
2026-07-14 09:58:17 INFO None 5154585: status FINISHED
2026-07-14 09:58:17 INFO None 5154586: status FINISHED
2026-07-14 09:58:17 INFO Jobs ['5154558', '5154559', '5154561', '5154562', '5154563', '5154564', '5154565', '5154566', '5154567', '5154568', '5154570', '5154572', '5154584', '5154585', '5154586'] have finished
2026-07-14 09:58:17 INFO Checking restart files were created ...
2026-07-14 09:58:17 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020601_8_ENS1.nc(3005806795 bytes)
2026-07-14 09:58:17 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020601_8_ENS2.nc(3005806795 bytes)
2026-07-14 09:58:17 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020601_8_ENS3.nc(3005806795 bytes)
2026-07-14 09:58:17 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020601_8_ENS4.nc(3005806795 bytes)
2026-07-14 09:58:17 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020601_8_ENS5.nc(3005806795 bytes)
2026-07-14 09:58:17 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020601_8_ENS6.nc(3005806795 bytes)
2026-07-14 09:58:17 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020601_8_ENS7.nc(3005806795 bytes)
2026-07-14 09:58:17 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020601_8_ENS8.nc(3005806795 bytes)
2026-07-14 09:58:17 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020601_8_ENS9.nc(3005806795 bytes)
2026-07-14 09:58:17 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020601_8_ENS10.nc(3005806795 bytes)
2026-07-14 09:58:17 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020601_8_ENS11.nc(3005806795 bytes)
2026-07-14 09:58:17 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020601_8_ENS12.nc(3005806795 bytes)
2026-07-14 09:58:17 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020601_8_ENS13.nc(3005806795 bytes)
2026-07-14 09:58:17 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020601_8_ENS14.nc(3005806795 bytes)
2026-07-14 09:58:17 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020601_8_ENS15.nc(3005806795 bytes)
2026-07-14 09:58:17 INFO  Run_model() completed successfully.
2026-07-14 09:58:17 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 09:58:17 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 09:00:00 days=153072 seconds=32400
2026-07-14 09:58:17 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-14 09:58:17 INFO [TIME] increment current_time 2020-02-06 01:00:00 -> 2020-02-06 09:00:00
2026-07-14 09:58:17 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 09:58:17 INFO ---------->>> Running process_satellite_data()
2026-07-14 09:58:17 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12002.nc
2026-07-14 09:58:17 INFO ---------->>> Running run_obs_converter()
2026-07-14 09:58:17 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_33067_153072.out
2026-07-14 09:58:17 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_33067_153072.out
2026-07-14 09:58:17 INFO ---------->>> Running DART
2026-07-14 09:58:17 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-14 09:58:17 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/chim_ENS1_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/chim_ENS1_2020020609_1_out_toDART.nc
2026-07-14 09:58:17 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/chim_ENS2_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/chim_ENS2_2020020609_1_out_toDART.nc
2026-07-14 09:58:18 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/chim_ENS3_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/chim_ENS3_2020020609_1_out_toDART.nc
2026-07-14 09:58:18 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/chim_ENS4_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/chim_ENS4_2020020609_1_out_toDART.nc
2026-07-14 09:58:19 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/chim_ENS5_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/chim_ENS5_2020020609_1_out_toDART.nc
2026-07-14 09:58:19 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/chim_ENS6_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/chim_ENS6_2020020609_1_out_toDART.nc
2026-07-14 09:58:19 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/chim_ENS7_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/chim_ENS7_2020020609_1_out_toDART.nc
2026-07-14 09:58:20 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/chim_ENS8_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/chim_ENS8_2020020609_1_out_toDART.nc
2026-07-14 09:58:20 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/chim_ENS9_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/chim_ENS9_2020020609_1_out_toDART.nc
2026-07-14 09:58:20 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/chim_ENS10_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/chim_ENS10_2020020609_1_out_toDART.nc
2026-07-14 09:58:21 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/chim_ENS11_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/chim_ENS11_2020020609_1_out_toDART.nc
2026-07-14 09:58:21 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/chim_ENS12_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/chim_ENS12_2020020609_1_out_toDART.nc
2026-07-14 09:58:22 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/chim_ENS13_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/chim_ENS13_2020020609_1_out_toDART.nc
2026-07-14 09:58:22 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/chim_ENS14_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/chim_ENS14_2020020609_1_out_toDART.nc
2026-07-14 09:58:22 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/chim_ENS15_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/chim_ENS15_2020020609_1_out_toDART.nc
2026-07-14 09:58:23 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-14 09:58:23 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-14 09:58:23 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-14 09:58:23 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-14 09:58:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-14 09:58:23 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-14 09:58:31 INFO Found: []
2026-07-14 09:58:31 INFO No job id returned by command ./run_filter.bsh
2026-07-14 09:58:31 INFO No monitoring will be performed
2026-07-14 09:58:31 INFO Moving DART output files to analysis and preassim directories for date 2020020609 if present ...
2026-07-14 09:58:31 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-14 09:58:31 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-14 09:58:31 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-14 09:58:31 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-14 09:58:31 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-14 09:58:31 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-14 09:58:31 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-14 09:58:31 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-14 09:58:31 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-14 09:58:31 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-14 09:58:31 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-14 09:58:31 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-14 09:58:31 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-14 09:58:31 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-14 09:58:31 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-14 09:58:31 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-14 09:58:31 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-14 09:58:31 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-14 09:58:31 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-14 09:58:31 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-14 09:58:31 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-14 09:58:31 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-14 09:58:31 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-14 09:58:31 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-14 09:58:31 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-14 09:58:31 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-14 09:58:31 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-14 09:58:31 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-14 09:58:31 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-14 09:58:31 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-14 09:58:31 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-14 09:58:31 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-14 09:58:31 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-14 09:58:32 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-14 09:58:32 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-14 09:58:34 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-14 09:58:34 INFO run_dart() is DONE.
2026-07-14 09:58:34 INFO ---------->>> Running update_pollutant_in_end()
2026-07-14 09:58:34 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 09:58:34 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS1_2020020609.nc
2026-07-14 09:58:34 INFO Emission correction applied with pixel-based damping.
2026-07-14 09:58:34 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS1_2020020609.nc
2026-07-14 09:58:35 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 09:58:35 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 09:58:35 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS1_2020020609.nc
2026-07-14 09:58:35 INFO Emission correction applied with pixel-based damping.
2026-07-14 09:58:35 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS1_2020020609.nc
2026-07-14 09:58:36 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 09:58:36 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 09:58:36 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 09:58:36 INFO No previous orbit memory found.
2026-07-14 09:58:36 INFO Emission correction applied with pixel-based damping.
2026-07-14 09:58:36 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS2_2020020609.nc
2026-07-14 09:58:37 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 09:58:37 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 09:58:37 INFO No previous orbit memory found.
2026-07-14 09:58:37 INFO Emission correction applied with pixel-based damping.
2026-07-14 09:58:37 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS2_2020020609.nc
2026-07-14 09:58:38 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 09:58:38 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 09:58:38 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 09:58:38 INFO No previous orbit memory found.
2026-07-14 09:58:38 INFO Emission correction applied with pixel-based damping.
2026-07-14 09:58:38 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS3_2020020609.nc
2026-07-14 09:58:38 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 09:58:39 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 09:58:39 INFO No previous orbit memory found.
2026-07-14 09:58:39 INFO Emission correction applied with pixel-based damping.
2026-07-14 09:58:39 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS3_2020020609.nc
2026-07-14 09:58:39 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 09:58:39 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 09:58:39 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 09:58:39 INFO No previous orbit memory found.
2026-07-14 09:58:39 INFO Emission correction applied with pixel-based damping.
2026-07-14 09:58:39 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS4_2020020609.nc
2026-07-14 09:58:40 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 09:58:40 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 09:58:40 INFO No previous orbit memory found.
2026-07-14 09:58:40 INFO Emission correction applied with pixel-based damping.
2026-07-14 09:58:40 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS4_2020020609.nc
2026-07-14 09:58:41 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 09:58:41 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 09:58:41 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 09:58:41 INFO No previous orbit memory found.
2026-07-14 09:58:41 INFO Emission correction applied with pixel-based damping.
2026-07-14 09:58:41 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS5_2020020609.nc
2026-07-14 09:58:42 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 09:58:42 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 09:58:42 INFO No previous orbit memory found.
2026-07-14 09:58:42 INFO Emission correction applied with pixel-based damping.
2026-07-14 09:58:42 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS5_2020020609.nc
2026-07-14 09:58:43 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 09:58:43 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 09:58:43 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 09:58:43 INFO No previous orbit memory found.
2026-07-14 09:58:43 INFO Emission correction applied with pixel-based damping.
2026-07-14 09:58:43 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS6_2020020609.nc
2026-07-14 09:58:44 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 09:58:44 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 09:58:44 INFO No previous orbit memory found.
2026-07-14 09:58:44 INFO Emission correction applied with pixel-based damping.
2026-07-14 09:58:44 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS6_2020020609.nc
2026-07-14 09:58:45 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 09:58:45 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 09:58:45 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 09:58:45 INFO No previous orbit memory found.
2026-07-14 09:58:45 INFO Emission correction applied with pixel-based damping.
2026-07-14 09:58:45 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS7_2020020609.nc
2026-07-14 09:58:46 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 09:58:46 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 09:58:46 INFO No previous orbit memory found.
2026-07-14 09:58:46 INFO Emission correction applied with pixel-based damping.
2026-07-14 09:58:46 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS7_2020020609.nc
2026-07-14 09:58:47 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 09:58:47 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 09:58:47 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 09:58:47 INFO No previous orbit memory found.
2026-07-14 09:58:47 INFO Emission correction applied with pixel-based damping.
2026-07-14 09:58:47 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS8_2020020609.nc
2026-07-14 09:58:48 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 09:58:48 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 09:58:48 INFO No previous orbit memory found.
2026-07-14 09:58:48 INFO Emission correction applied with pixel-based damping.
2026-07-14 09:58:48 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS8_2020020609.nc
2026-07-14 09:58:48 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 09:58:48 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 09:58:49 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 09:58:49 INFO No previous orbit memory found.
2026-07-14 09:58:49 INFO Emission correction applied with pixel-based damping.
2026-07-14 09:58:49 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS9_2020020609.nc
2026-07-14 09:58:49 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 09:58:50 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 09:58:50 INFO No previous orbit memory found.
2026-07-14 09:58:50 INFO Emission correction applied with pixel-based damping.
2026-07-14 09:58:50 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS9_2020020609.nc
2026-07-14 09:58:50 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 09:58:50 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 09:58:50 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 09:58:50 INFO No previous orbit memory found.
2026-07-14 09:58:50 INFO Emission correction applied with pixel-based damping.
2026-07-14 09:58:51 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS10_2020020609.nc
2026-07-14 09:58:51 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 09:58:51 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 09:58:51 INFO No previous orbit memory found.
2026-07-14 09:58:51 INFO Emission correction applied with pixel-based damping.
2026-07-14 09:58:51 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS10_2020020609.nc
2026-07-14 09:58:52 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 09:58:52 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 09:58:52 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 09:58:52 INFO No previous orbit memory found.
2026-07-14 09:58:52 INFO Emission correction applied with pixel-based damping.
2026-07-14 09:58:52 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS11_2020020609.nc
2026-07-14 09:58:53 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 09:58:53 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 09:58:53 INFO No previous orbit memory found.
2026-07-14 09:58:53 INFO Emission correction applied with pixel-based damping.
2026-07-14 09:58:53 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS11_2020020609.nc
2026-07-14 09:58:54 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 09:58:54 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 09:58:54 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 09:58:54 INFO No previous orbit memory found.
2026-07-14 09:58:54 INFO Emission correction applied with pixel-based damping.
2026-07-14 09:58:54 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS12_2020020609.nc
2026-07-14 09:58:55 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 09:58:55 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 09:58:55 INFO No previous orbit memory found.
2026-07-14 09:58:55 INFO Emission correction applied with pixel-based damping.
2026-07-14 09:58:55 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS12_2020020609.nc
2026-07-14 09:58:56 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 09:58:56 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 09:58:56 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 09:58:56 INFO No previous orbit memory found.
2026-07-14 09:58:56 INFO Emission correction applied with pixel-based damping.
2026-07-14 09:58:56 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS13_2020020609.nc
2026-07-14 09:58:57 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 09:58:57 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 09:58:57 INFO No previous orbit memory found.
2026-07-14 09:58:57 INFO Emission correction applied with pixel-based damping.
2026-07-14 09:58:57 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS13_2020020609.nc
2026-07-14 09:58:58 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 09:58:58 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 09:58:58 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 09:58:58 INFO No previous orbit memory found.
2026-07-14 09:58:58 INFO Emission correction applied with pixel-based damping.
2026-07-14 09:58:58 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS14_2020020609.nc
2026-07-14 09:58:59 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 09:58:59 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 09:58:59 INFO No previous orbit memory found.
2026-07-14 09:58:59 INFO Emission correction applied with pixel-based damping.
2026-07-14 09:58:59 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS14_2020020609.nc
2026-07-14 09:58:59 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 09:58:59 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 09:59:00 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 09:59:00 INFO No previous orbit memory found.
2026-07-14 09:59:00 INFO Emission correction applied with pixel-based damping.
2026-07-14 09:59:00 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS15_2020020609.nc
2026-07-14 09:59:00 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 09:59:01 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 09:59:01 INFO No previous orbit memory found.
2026-07-14 09:59:01 INFO Emission correction applied with pixel-based damping.
2026-07-14 09:59:01 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS15_2020020609.nc
2026-07-14 09:59:01 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 09:59:01 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 09:59:01 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-14 09:59:01 INFO [TIME] step_end current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 09:59:01 INFO [TIME] step_start current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 09:59:01 INFO [TIME] window start=2020-02-06 09:00:00 end=2020-02-06 11:00:00 run_hours=2 has_assimilation=True
2026-07-14 09:59:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:59:03 INFO Hourly dataset computed and listing created
2026-07-14 09:59:06 INFO Hourly dataset computed
2026-07-14 09:59:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:59:07 INFO Hourly dataset computed and listing created
2026-07-14 09:59:08 INFO Hourly dataset computed
2026-07-14 09:59:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:59:09 INFO Hourly dataset computed and listing created
2026-07-14 09:59:10 INFO Hourly dataset computed
2026-07-14 09:59:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:59:11 INFO Hourly dataset computed and listing created
2026-07-14 09:59:11 INFO Hourly dataset computed
2026-07-14 09:59:11 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:59:12 INFO Hourly dataset computed and listing created
2026-07-14 09:59:13 INFO Hourly dataset computed
2026-07-14 09:59:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:59:14 INFO Hourly dataset computed and listing created
2026-07-14 09:59:15 INFO Hourly dataset computed
2026-07-14 09:59:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:59:16 INFO Hourly dataset computed and listing created
2026-07-14 09:59:16 INFO Hourly dataset computed
2026-07-14 09:59:16 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:59:17 INFO Hourly dataset computed and listing created
2026-07-14 09:59:19 INFO Hourly dataset computed
2026-07-14 09:59:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:59:20 INFO Hourly dataset computed and listing created
2026-07-14 09:59:21 INFO Hourly dataset computed
2026-07-14 09:59:21 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:59:22 INFO Hourly dataset computed and listing created
2026-07-14 09:59:22 INFO Hourly dataset computed
2026-07-14 09:59:22 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:59:23 INFO Hourly dataset computed and listing created
2026-07-14 09:59:24 INFO Hourly dataset computed
2026-07-14 09:59:24 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:59:25 INFO Hourly dataset computed and listing created
2026-07-14 09:59:26 INFO Hourly dataset computed
2026-07-14 09:59:26 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:59:27 INFO Hourly dataset computed and listing created
2026-07-14 09:59:27 INFO Hourly dataset computed
2026-07-14 09:59:27 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:59:28 INFO Hourly dataset computed and listing created
2026-07-14 09:59:29 INFO Hourly dataset computed
2026-07-14 09:59:29 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 09:59:30 INFO Hourly dataset computed and listing created
2026-07-14 09:59:31 INFO Hourly dataset computed
2026-07-14 09:59:31 INFO ---------->>> Running CHIMERE model from 2020-02-06 09:00:00 to 2020-02-06 11:00:00
2026-07-14 09:59:31 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:59:31 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1
2026-07-14 09:59:31 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020601_8_ENS1.nc
2026-07-14 09:59:31 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-14 09:59:31 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:59:31 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-14 09:59:31 INFO Queuing job for member 1...
2026-07-14 09:59:31 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:59:31 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-14 09:59:32 INFO Found: ['5154705']
2026-07-14 09:59:37 INFO [TGCC-IRENE] Submitted job with ID:['5154705']
2026-07-14 09:59:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:59:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2
2026-07-14 09:59:37 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020601_8_ENS2.nc
2026-07-14 09:59:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-14 09:59:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:59:37 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-14 09:59:37 INFO Queuing job for member 2...
2026-07-14 09:59:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:59:37 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-14 09:59:37 INFO Found: ['5154707']
2026-07-14 09:59:42 INFO [TGCC-IRENE] Submitted job with ID:['5154707']
2026-07-14 09:59:42 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:59:42 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3
2026-07-14 09:59:42 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020601_8_ENS3.nc
2026-07-14 09:59:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-14 09:59:42 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:59:42 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-14 09:59:42 INFO Queuing job for member 3...
2026-07-14 09:59:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:59:42 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-14 09:59:43 INFO Found: ['5154708']
2026-07-14 09:59:48 INFO [TGCC-IRENE] Submitted job with ID:['5154708']
2026-07-14 09:59:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:59:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4
2026-07-14 09:59:48 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020601_8_ENS4.nc
2026-07-14 09:59:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-14 09:59:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:59:48 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-14 09:59:48 INFO Queuing job for member 4...
2026-07-14 09:59:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:59:48 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-14 09:59:49 INFO Found: ['5154711']
2026-07-14 09:59:54 INFO [TGCC-IRENE] Submitted job with ID:['5154711']
2026-07-14 09:59:54 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 09:59:54 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5
2026-07-14 09:59:54 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020601_8_ENS5.nc
2026-07-14 09:59:54 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-14 09:59:54 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 09:59:54 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-14 09:59:54 INFO Queuing job for member 5...
2026-07-14 09:59:54 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 09:59:54 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-14 09:59:55 INFO Found: ['5154712']
2026-07-14 10:00:00 INFO [TGCC-IRENE] Submitted job with ID:['5154712']
2026-07-14 10:00:00 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:00:00 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6
2026-07-14 10:00:00 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020601_8_ENS6.nc
2026-07-14 10:00:00 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-14 10:00:00 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:00:00 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-14 10:00:00 INFO Queuing job for member 6...
2026-07-14 10:00:00 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:00:00 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-14 10:00:00 INFO Found: ['5154713']
2026-07-14 10:00:05 INFO [TGCC-IRENE] Submitted job with ID:['5154713']
2026-07-14 10:00:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:00:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7
2026-07-14 10:00:05 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020601_8_ENS7.nc
2026-07-14 10:00:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-14 10:00:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:00:05 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-14 10:00:05 INFO Queuing job for member 7...
2026-07-14 10:00:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:00:05 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-14 10:00:06 INFO Found: ['5154715']
2026-07-14 10:00:11 INFO [TGCC-IRENE] Submitted job with ID:['5154715']
2026-07-14 10:00:11 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:00:11 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8
2026-07-14 10:00:11 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020601_8_ENS8.nc
2026-07-14 10:00:11 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-14 10:00:11 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:00:11 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-14 10:00:11 INFO Queuing job for member 8...
2026-07-14 10:00:11 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:00:11 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-14 10:00:12 INFO Found: ['5154716']
2026-07-14 10:00:17 INFO [TGCC-IRENE] Submitted job with ID:['5154716']
2026-07-14 10:00:17 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:00:17 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9
2026-07-14 10:00:17 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020601_8_ENS9.nc
2026-07-14 10:00:17 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-14 10:00:17 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:00:17 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-14 10:00:17 INFO Queuing job for member 9...
2026-07-14 10:00:17 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:00:17 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-14 10:00:18 INFO Found: ['5154717']
2026-07-14 10:00:23 INFO [TGCC-IRENE] Submitted job with ID:['5154717']
2026-07-14 10:00:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:00:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10
2026-07-14 10:00:23 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020601_8_ENS10.nc
2026-07-14 10:00:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-14 10:00:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:00:23 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-14 10:00:23 INFO Queuing job for member 10...
2026-07-14 10:00:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:00:23 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-14 10:00:24 INFO Found: ['5154718']
2026-07-14 10:00:29 INFO [TGCC-IRENE] Submitted job with ID:['5154718']
2026-07-14 10:00:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:00:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11
2026-07-14 10:00:29 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020601_8_ENS11.nc
2026-07-14 10:00:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-14 10:00:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:00:29 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-14 10:00:29 INFO Queuing job for member 11...
2026-07-14 10:00:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:00:29 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-14 10:00:29 INFO Found: ['5154720']
2026-07-14 10:00:34 INFO [TGCC-IRENE] Submitted job with ID:['5154720']
2026-07-14 10:00:34 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:00:34 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12
2026-07-14 10:00:34 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020601_8_ENS12.nc
2026-07-14 10:00:34 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-14 10:00:34 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:00:34 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-14 10:00:34 INFO Queuing job for member 12...
2026-07-14 10:00:34 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:00:34 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-14 10:00:35 INFO Found: ['5154722']
2026-07-14 10:00:40 INFO [TGCC-IRENE] Submitted job with ID:['5154722']
2026-07-14 10:00:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:00:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13
2026-07-14 10:00:40 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020601_8_ENS13.nc
2026-07-14 10:00:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-14 10:00:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:00:40 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-14 10:00:40 INFO Queuing job for member 13...
2026-07-14 10:00:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:00:40 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-14 10:00:41 INFO Found: ['5154723']
2026-07-14 10:00:46 INFO [TGCC-IRENE] Submitted job with ID:['5154723']
2026-07-14 10:00:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:00:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14
2026-07-14 10:00:46 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020601_8_ENS14.nc
2026-07-14 10:00:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-14 10:00:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:00:46 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-14 10:00:46 INFO Queuing job for member 14...
2026-07-14 10:00:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:00:46 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-14 10:00:47 INFO Found: ['5154724']
2026-07-14 10:00:52 INFO [TGCC-IRENE] Submitted job with ID:['5154724']
2026-07-14 10:00:52 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:00:52 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15
2026-07-14 10:00:52 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020601_8_ENS15.nc
2026-07-14 10:00:52 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-14 10:00:52 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:00:52 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-14 10:00:52 INFO Queuing job for member 15...
2026-07-14 10:00:52 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:00:52 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-14 10:00:52 INFO Found: ['5154727']
2026-07-14 10:00:57 INFO [TGCC-IRENE] Submitted job with ID:['5154727']
2026-07-14 10:00:57 INFO Checking job status ...
2026-07-14 10:00:57 INFO None 5154705: status RUNNING/PENDING
2026-07-14 10:00:57 INFO None 5154707: status RUNNING/PENDING
2026-07-14 10:00:57 INFO None 5154708: status RUNNING/PENDING
2026-07-14 10:00:57 INFO None 5154711: status RUNNING/PENDING
2026-07-14 10:00:57 INFO None 5154712: status RUNNING/PENDING
2026-07-14 10:00:57 INFO None 5154713: status RUNNING/PENDING
2026-07-14 10:00:57 INFO None 5154715: status RUNNING/PENDING
2026-07-14 10:00:58 INFO None 5154716: status RUNNING/PENDING
2026-07-14 10:00:58 INFO None 5154717: status RUNNING/PENDING
2026-07-14 10:00:58 INFO None 5154718: status RUNNING/PENDING
2026-07-14 10:00:58 INFO None 5154720: status RUNNING/PENDING
2026-07-14 10:00:58 INFO None 5154722: status RUNNING/PENDING
2026-07-14 10:00:58 INFO None 5154723: status RUNNING/PENDING
2026-07-14 10:00:58 INFO None 5154724: status RUNNING/PENDING
2026-07-14 10:00:58 INFO None 5154727: status RUNNING/PENDING
2026-07-14 10:00:58 INFO Jobs still running: ['5154705', '5154707', '5154708', '5154711', '5154712', '5154713', '5154715', '5154716', '5154717', '5154718', '5154720', '5154722', '5154723', '5154724', '5154727']. Waiting...
2026-07-14 10:01:13 INFO None 5154705: status RUNNING/PENDING
2026-07-14 10:01:13 INFO None 5154707: status RUNNING/PENDING
2026-07-14 10:01:13 INFO None 5154708: status RUNNING/PENDING
2026-07-14 10:01:13 INFO None 5154711: status RUNNING/PENDING
2026-07-14 10:01:13 INFO None 5154712: status RUNNING/PENDING
2026-07-14 10:01:13 INFO None 5154713: status RUNNING/PENDING
2026-07-14 10:01:13 INFO None 5154715: status RUNNING/PENDING
2026-07-14 10:01:13 INFO None 5154716: status RUNNING/PENDING
2026-07-14 10:01:13 INFO None 5154717: status RUNNING/PENDING
2026-07-14 10:01:13 INFO None 5154718: status RUNNING/PENDING
2026-07-14 10:01:13 INFO None 5154720: status RUNNING/PENDING
2026-07-14 10:01:13 INFO None 5154722: status RUNNING/PENDING
2026-07-14 10:01:13 INFO None 5154723: status RUNNING/PENDING
2026-07-14 10:01:13 INFO None 5154724: status RUNNING/PENDING
2026-07-14 10:01:13 INFO None 5154727: status RUNNING/PENDING
2026-07-14 10:01:13 INFO Jobs still running: ['5154705', '5154707', '5154708', '5154711', '5154712', '5154713', '5154715', '5154716', '5154717', '5154718', '5154720', '5154722', '5154723', '5154724', '5154727']. Waiting...
2026-07-14 10:01:28 INFO None 5154705: status RUNNING/PENDING
2026-07-14 10:01:28 INFO None 5154707: status RUNNING/PENDING
2026-07-14 10:01:28 INFO None 5154708: status RUNNING/PENDING
2026-07-14 10:01:28 INFO None 5154711: status RUNNING/PENDING
2026-07-14 10:01:28 INFO None 5154712: status RUNNING/PENDING
2026-07-14 10:01:28 INFO None 5154713: status RUNNING/PENDING
2026-07-14 10:01:28 INFO None 5154715: status RUNNING/PENDING
2026-07-14 10:01:28 INFO None 5154716: status RUNNING/PENDING
2026-07-14 10:01:28 INFO None 5154717: status RUNNING/PENDING
2026-07-14 10:01:28 INFO None 5154718: status RUNNING/PENDING
2026-07-14 10:01:28 INFO None 5154720: status RUNNING/PENDING
2026-07-14 10:01:28 INFO None 5154722: status RUNNING/PENDING
2026-07-14 10:01:28 INFO None 5154723: status RUNNING/PENDING
2026-07-14 10:01:28 INFO None 5154724: status RUNNING/PENDING
2026-07-14 10:01:28 INFO None 5154727: status RUNNING/PENDING
2026-07-14 10:01:28 INFO Jobs still running: ['5154705', '5154707', '5154708', '5154711', '5154712', '5154713', '5154715', '5154716', '5154717', '5154718', '5154720', '5154722', '5154723', '5154724', '5154727']. Waiting...
2026-07-14 10:01:43 INFO None 5154705: status RUNNING/PENDING
2026-07-14 10:03:10 INFO None 5154707: status RUNNING/PENDING
2026-07-14 10:03:10 INFO None 5154708: status RUNNING/PENDING
2026-07-14 10:03:10 INFO None 5154711: status FINISHED
2026-07-14 10:03:10 INFO None 5154712: status FINISHED
2026-07-14 10:03:10 INFO None 5154713: status RUNNING/PENDING
2026-07-14 10:03:10 INFO None 5154715: status RUNNING/PENDING
2026-07-14 10:03:10 INFO None 5154716: status RUNNING/PENDING
2026-07-14 10:03:10 INFO None 5154717: status RUNNING/PENDING
2026-07-14 10:03:11 INFO None 5154718: status RUNNING/PENDING
2026-07-14 10:03:11 INFO None 5154720: status RUNNING/PENDING
2026-07-14 10:03:11 INFO None 5154722: status RUNNING/PENDING
2026-07-14 10:03:11 INFO None 5154723: status RUNNING/PENDING
2026-07-14 10:03:11 INFO None 5154724: status RUNNING/PENDING
2026-07-14 10:03:11 INFO None 5154727: status RUNNING/PENDING
2026-07-14 10:03:11 INFO Jobs still running: ['5154705', '5154707', '5154708', '5154713', '5154715', '5154716', '5154717', '5154718', '5154720', '5154722', '5154723', '5154724', '5154727']. Waiting...
2026-07-14 10:03:26 INFO None 5154705: status RUNNING/PENDING
2026-07-14 10:03:26 INFO None 5154707: status RUNNING/PENDING
2026-07-14 10:03:26 INFO None 5154708: status RUNNING/PENDING
2026-07-14 10:03:26 INFO None 5154711: status FINISHED
2026-07-14 10:03:26 INFO None 5154712: status FINISHED
2026-07-14 10:03:26 INFO None 5154713: status RUNNING/PENDING
2026-07-14 10:03:26 INFO None 5154715: status RUNNING/PENDING
2026-07-14 10:03:26 INFO None 5154716: status RUNNING/PENDING
2026-07-14 10:03:26 INFO None 5154717: status RUNNING/PENDING
2026-07-14 10:03:26 INFO None 5154718: status RUNNING/PENDING
2026-07-14 10:03:26 INFO None 5154720: status RUNNING/PENDING
2026-07-14 10:03:26 INFO None 5154722: status RUNNING/PENDING
2026-07-14 10:03:26 INFO None 5154723: status RUNNING/PENDING
2026-07-14 10:03:26 INFO None 5154724: status RUNNING/PENDING
2026-07-14 10:03:26 INFO None 5154727: status RUNNING/PENDING
2026-07-14 10:03:26 INFO Jobs still running: ['5154705', '5154707', '5154708', '5154713', '5154715', '5154716', '5154717', '5154718', '5154720', '5154722', '5154723', '5154724', '5154727']. Waiting...
2026-07-14 10:03:41 INFO None 5154705: status RUNNING/PENDING
2026-07-14 10:03:41 INFO None 5154707: status RUNNING/PENDING
2026-07-14 10:03:41 INFO None 5154708: status RUNNING/PENDING
2026-07-14 10:03:41 INFO None 5154711: status FINISHED
2026-07-14 10:03:41 INFO None 5154712: status FINISHED
2026-07-14 10:03:41 INFO None 5154713: status RUNNING/PENDING
2026-07-14 10:03:41 INFO None 5154715: status RUNNING/PENDING
2026-07-14 10:03:41 INFO None 5154716: status RUNNING/PENDING
2026-07-14 10:03:41 INFO None 5154717: status RUNNING/PENDING
2026-07-14 10:03:41 INFO None 5154718: status RUNNING/PENDING
2026-07-14 10:03:41 INFO None 5154720: status RUNNING/PENDING
2026-07-14 10:03:41 INFO None 5154722: status RUNNING/PENDING
2026-07-14 10:03:41 INFO None 5154723: status RUNNING/PENDING
2026-07-14 10:03:41 INFO None 5154724: status RUNNING/PENDING
2026-07-14 10:03:41 INFO None 5154727: status RUNNING/PENDING
2026-07-14 10:03:41 INFO Jobs still running: ['5154705', '5154707', '5154708', '5154713', '5154715', '5154716', '5154717', '5154718', '5154720', '5154722', '5154723', '5154724', '5154727']. Waiting...
2026-07-14 10:03:56 INFO None 5154705: status RUNNING/PENDING
2026-07-14 10:03:56 INFO None 5154707: status RUNNING/PENDING
2026-07-14 10:03:56 INFO None 5154708: status RUNNING/PENDING
2026-07-14 10:03:56 INFO None 5154711: status FINISHED
2026-07-14 10:03:56 INFO None 5154712: status FINISHED
2026-07-14 10:03:56 INFO None 5154713: status RUNNING/PENDING
2026-07-14 10:03:56 INFO None 5154715: status RUNNING/PENDING
2026-07-14 10:03:56 INFO None 5154716: status RUNNING/PENDING
2026-07-14 10:03:56 INFO None 5154717: status RUNNING/PENDING
2026-07-14 10:03:56 INFO None 5154718: status RUNNING/PENDING
2026-07-14 10:03:56 INFO None 5154720: status RUNNING/PENDING
2026-07-14 10:03:56 INFO None 5154722: status RUNNING/PENDING
2026-07-14 10:03:56 INFO None 5154723: status RUNNING/PENDING
2026-07-14 10:03:56 INFO None 5154724: status RUNNING/PENDING
2026-07-14 10:03:56 INFO None 5154727: status RUNNING/PENDING
2026-07-14 10:03:56 INFO Jobs still running: ['5154705', '5154707', '5154708', '5154713', '5154715', '5154716', '5154717', '5154718', '5154720', '5154722', '5154723', '5154724', '5154727']. Waiting...
2026-07-14 10:04:12 INFO None 5154705: status RUNNING/PENDING
2026-07-14 10:04:12 INFO None 5154707: status RUNNING/PENDING
2026-07-14 10:04:12 INFO None 5154708: status RUNNING/PENDING
2026-07-14 10:04:12 INFO None 5154711: status FINISHED
2026-07-14 10:04:12 INFO None 5154712: status FINISHED
2026-07-14 10:04:12 INFO None 5154713: status RUNNING/PENDING
2026-07-14 10:04:12 INFO None 5154715: status RUNNING/PENDING
2026-07-14 10:04:12 INFO None 5154716: status RUNNING/PENDING
2026-07-14 10:04:12 INFO None 5154717: status RUNNING/PENDING
2026-07-14 10:04:12 INFO None 5154718: status RUNNING/PENDING
2026-07-14 10:04:12 INFO None 5154720: status RUNNING/PENDING
2026-07-14 10:04:12 INFO None 5154722: status RUNNING/PENDING
2026-07-14 10:04:12 INFO None 5154723: status RUNNING/PENDING
2026-07-14 10:04:12 INFO None 5154724: status RUNNING/PENDING
2026-07-14 10:04:12 INFO None 5154727: status RUNNING/PENDING
2026-07-14 10:04:12 INFO Jobs still running: ['5154705', '5154707', '5154708', '5154713', '5154715', '5154716', '5154717', '5154718', '5154720', '5154722', '5154723', '5154724', '5154727']. Waiting...
2026-07-14 10:04:27 INFO None 5154705: status RUNNING/PENDING
2026-07-14 10:04:27 INFO None 5154707: status RUNNING/PENDING
2026-07-14 10:04:27 INFO None 5154708: status RUNNING/PENDING
2026-07-14 10:04:27 INFO None 5154711: status FINISHED
2026-07-14 10:04:27 INFO None 5154712: status FINISHED
2026-07-14 10:04:27 INFO None 5154713: status RUNNING/PENDING
2026-07-14 10:04:27 INFO None 5154715: status RUNNING/PENDING
2026-07-14 10:04:27 INFO None 5154716: status RUNNING/PENDING
2026-07-14 10:04:27 INFO None 5154717: status RUNNING/PENDING
2026-07-14 10:04:27 INFO None 5154718: status RUNNING/PENDING
2026-07-14 10:04:27 INFO None 5154720: status RUNNING/PENDING
2026-07-14 10:04:27 INFO None 5154722: status RUNNING/PENDING
2026-07-14 10:04:27 INFO None 5154723: status RUNNING/PENDING
2026-07-14 10:04:27 INFO None 5154724: status RUNNING/PENDING
2026-07-14 10:04:27 INFO None 5154727: status RUNNING/PENDING
2026-07-14 10:04:27 INFO Jobs still running: ['5154705', '5154707', '5154708', '5154713', '5154715', '5154716', '5154717', '5154718', '5154720', '5154722', '5154723', '5154724', '5154727']. Waiting...
2026-07-14 10:04:42 INFO None 5154705: status RUNNING/PENDING
2026-07-14 10:04:42 INFO None 5154707: status RUNNING/PENDING
2026-07-14 10:04:42 INFO None 5154708: status FINISHED
2026-07-14 10:04:42 INFO None 5154711: status FINISHED
2026-07-14 10:04:42 INFO None 5154712: status FINISHED
2026-07-14 10:04:42 INFO None 5154713: status RUNNING/PENDING
2026-07-14 10:04:42 INFO None 5154715: status RUNNING/PENDING
2026-07-14 10:04:42 INFO None 5154716: status RUNNING/PENDING
2026-07-14 10:04:42 INFO None 5154717: status RUNNING/PENDING
2026-07-14 10:04:42 INFO None 5154718: status RUNNING/PENDING
2026-07-14 10:04:42 INFO None 5154720: status RUNNING/PENDING
2026-07-14 10:04:42 INFO None 5154722: status RUNNING/PENDING
2026-07-14 10:04:42 INFO None 5154723: status RUNNING/PENDING
2026-07-14 10:04:42 INFO None 5154724: status RUNNING/PENDING
2026-07-14 10:04:42 INFO None 5154727: status RUNNING/PENDING
2026-07-14 10:04:42 INFO Jobs still running: ['5154705', '5154707', '5154713', '5154715', '5154716', '5154717', '5154718', '5154720', '5154722', '5154723', '5154724', '5154727']. Waiting...
2026-07-14 10:04:57 INFO None 5154705: status RUNNING/PENDING
2026-07-14 10:04:57 INFO None 5154707: status RUNNING/PENDING
2026-07-14 10:04:57 INFO None 5154708: status FINISHED
2026-07-14 10:04:57 INFO None 5154711: status FINISHED
2026-07-14 10:04:57 INFO None 5154712: status FINISHED
2026-07-14 10:04:57 INFO None 5154713: status RUNNING/PENDING
2026-07-14 10:04:58 INFO None 5154715: status RUNNING/PENDING
2026-07-14 10:04:58 INFO None 5154716: status RUNNING/PENDING
2026-07-14 10:04:58 INFO None 5154717: status RUNNING/PENDING
2026-07-14 10:04:58 INFO None 5154718: status RUNNING/PENDING
2026-07-14 10:04:58 INFO None 5154720: status RUNNING/PENDING
2026-07-14 10:04:58 INFO None 5154722: status RUNNING/PENDING
2026-07-14 10:04:58 INFO None 5154723: status RUNNING/PENDING
2026-07-14 10:04:58 INFO None 5154724: status RUNNING/PENDING
2026-07-14 10:04:58 INFO None 5154727: status RUNNING/PENDING
2026-07-14 10:04:58 INFO Jobs still running: ['5154705', '5154707', '5154713', '5154715', '5154716', '5154717', '5154718', '5154720', '5154722', '5154723', '5154724', '5154727']. Waiting...
2026-07-14 10:05:13 INFO None 5154705: status RUNNING/PENDING
2026-07-14 10:05:13 INFO None 5154707: status RUNNING/PENDING
2026-07-14 10:05:13 INFO None 5154708: status FINISHED
2026-07-14 10:05:13 INFO None 5154711: status FINISHED
2026-07-14 10:05:13 INFO None 5154712: status FINISHED
2026-07-14 10:05:13 INFO None 5154713: status RUNNING/PENDING
2026-07-14 10:05:13 INFO None 5154715: status RUNNING/PENDING
2026-07-14 10:05:13 INFO None 5154716: status RUNNING/PENDING
2026-07-14 10:05:13 INFO None 5154717: status RUNNING/PENDING
2026-07-14 10:05:13 INFO None 5154718: status RUNNING/PENDING
2026-07-14 10:05:13 INFO None 5154720: status RUNNING/PENDING
2026-07-14 10:05:13 INFO None 5154722: status RUNNING/PENDING
2026-07-14 10:05:13 INFO None 5154723: status RUNNING/PENDING
2026-07-14 10:05:13 INFO None 5154724: status RUNNING/PENDING
2026-07-14 10:05:13 INFO None 5154727: status RUNNING/PENDING
2026-07-14 10:05:13 INFO Jobs still running: ['5154705', '5154707', '5154713', '5154715', '5154716', '5154717', '5154718', '5154720', '5154722', '5154723', '5154724', '5154727']. Waiting...
2026-07-14 10:05:28 INFO None 5154705: status RUNNING/PENDING
2026-07-14 10:05:28 INFO None 5154707: status RUNNING/PENDING
2026-07-14 10:05:28 INFO None 5154708: status FINISHED
2026-07-14 10:05:28 INFO None 5154711: status FINISHED
2026-07-14 10:05:28 INFO None 5154712: status FINISHED
2026-07-14 10:05:28 INFO None 5154713: status RUNNING/PENDING
2026-07-14 10:05:28 INFO None 5154715: status RUNNING/PENDING
2026-07-14 10:05:28 INFO None 5154716: status RUNNING/PENDING
2026-07-14 10:05:28 INFO None 5154717: status RUNNING/PENDING
2026-07-14 10:05:28 INFO None 5154718: status RUNNING/PENDING
2026-07-14 10:05:28 INFO None 5154720: status RUNNING/PENDING
2026-07-14 10:05:28 INFO None 5154722: status RUNNING/PENDING
2026-07-14 10:05:28 INFO None 5154723: status RUNNING/PENDING
2026-07-14 10:05:28 INFO None 5154724: status RUNNING/PENDING
2026-07-14 10:05:28 INFO None 5154727: status RUNNING/PENDING
2026-07-14 10:05:28 INFO Jobs still running: ['5154705', '5154707', '5154713', '5154715', '5154716', '5154717', '5154718', '5154720', '5154722', '5154723', '5154724', '5154727']. Waiting...
2026-07-14 10:05:43 INFO None 5154705: status RUNNING/PENDING
2026-07-14 10:05:43 INFO None 5154707: status RUNNING/PENDING
2026-07-14 10:05:43 INFO None 5154708: status FINISHED
2026-07-14 10:05:43 INFO None 5154711: status FINISHED
2026-07-14 10:05:43 INFO None 5154712: status FINISHED
2026-07-14 10:05:43 INFO None 5154713: status RUNNING/PENDING
2026-07-14 10:05:43 INFO None 5154715: status RUNNING/PENDING
2026-07-14 10:05:43 INFO None 5154716: status FINISHED
2026-07-14 10:05:43 INFO None 5154717: status FINISHED
2026-07-14 10:05:43 INFO None 5154718: status FINISHED
2026-07-14 10:05:43 INFO None 5154720: status FINISHED
2026-07-14 10:05:43 INFO None 5154722: status FINISHED
2026-07-14 10:05:43 INFO None 5154723: status FINISHED
2026-07-14 10:05:44 INFO None 5154724: status RUNNING/PENDING
2026-07-14 10:05:44 INFO None 5154727: status RUNNING/PENDING
2026-07-14 10:05:44 INFO Jobs still running: ['5154705', '5154707', '5154713', '5154715', '5154724', '5154727']. Waiting...
2026-07-14 10:05:59 INFO None 5154705: status RUNNING/PENDING
2026-07-14 10:05:59 INFO None 5154707: status RUNNING/PENDING
2026-07-14 10:05:59 INFO None 5154708: status FINISHED
2026-07-14 10:05:59 INFO None 5154711: status FINISHED
2026-07-14 10:05:59 INFO None 5154712: status FINISHED
2026-07-14 10:05:59 INFO None 5154713: status RUNNING/PENDING
2026-07-14 10:05:59 INFO None 5154715: status RUNNING/PENDING
2026-07-14 10:05:59 INFO None 5154716: status FINISHED
2026-07-14 10:05:59 INFO None 5154717: status FINISHED
2026-07-14 10:05:59 INFO None 5154718: status FINISHED
2026-07-14 10:05:59 INFO None 5154720: status FINISHED
2026-07-14 10:05:59 INFO None 5154722: status FINISHED
2026-07-14 10:05:59 INFO None 5154723: status FINISHED
2026-07-14 10:05:59 INFO None 5154724: status RUNNING/PENDING
2026-07-14 10:05:59 INFO None 5154727: status RUNNING/PENDING
2026-07-14 10:05:59 INFO Jobs still running: ['5154705', '5154707', '5154713', '5154715', '5154724', '5154727']. Waiting...
2026-07-14 10:06:14 INFO None 5154705: status RUNNING/PENDING
2026-07-14 10:08:12 INFO None 5154707: status FINISHED
2026-07-14 10:08:12 INFO None 5154708: status FINISHED
2026-07-14 10:08:12 INFO None 5154711: status FINISHED
2026-07-14 10:08:12 INFO None 5154712: status FINISHED
2026-07-14 10:08:12 INFO None 5154713: status FINISHED
2026-07-14 10:08:12 INFO None 5154715: status FINISHED
2026-07-14 10:08:12 INFO None 5154716: status FINISHED
2026-07-14 10:08:12 INFO None 5154717: status FINISHED
2026-07-14 10:08:12 INFO None 5154718: status FINISHED
2026-07-14 10:08:12 INFO None 5154720: status FINISHED
2026-07-14 10:08:12 INFO None 5154722: status FINISHED
2026-07-14 10:08:12 INFO None 5154723: status FINISHED
2026-07-14 10:08:12 INFO None 5154724: status FINISHED
2026-07-14 10:08:12 INFO None 5154727: status FINISHED
2026-07-14 10:08:12 INFO Jobs still running: ['5154705']. Waiting...
2026-07-14 10:08:27 INFO None 5154705: status FINISHED
2026-07-14 10:08:27 INFO None 5154707: status FINISHED
2026-07-14 10:08:27 INFO None 5154708: status FINISHED
2026-07-14 10:08:27 INFO None 5154711: status FINISHED
2026-07-14 10:08:27 INFO None 5154712: status FINISHED
2026-07-14 10:08:27 INFO None 5154713: status FINISHED
2026-07-14 10:08:27 INFO None 5154715: status FINISHED
2026-07-14 10:08:27 INFO None 5154716: status FINISHED
2026-07-14 10:08:27 INFO None 5154717: status FINISHED
2026-07-14 10:08:27 INFO None 5154718: status FINISHED
2026-07-14 10:08:27 INFO None 5154720: status FINISHED
2026-07-14 10:08:27 INFO None 5154722: status FINISHED
2026-07-14 10:08:27 INFO None 5154723: status FINISHED
2026-07-14 10:08:27 INFO None 5154724: status FINISHED
2026-07-14 10:08:27 INFO None 5154727: status FINISHED
2026-07-14 10:08:27 INFO Jobs ['5154705', '5154707', '5154708', '5154711', '5154712', '5154713', '5154715', '5154716', '5154717', '5154718', '5154720', '5154722', '5154723', '5154724', '5154727'] have finished
2026-07-14 10:08:27 INFO Checking restart files were created ...
2026-07-14 10:08:27 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc(1002685915 bytes)
2026-07-14 10:08:27 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc(1002685915 bytes)
2026-07-14 10:08:27 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc(1002685915 bytes)
2026-07-14 10:08:27 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc(1002685915 bytes)
2026-07-14 10:08:27 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc(1002685915 bytes)
2026-07-14 10:08:27 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc(1002685915 bytes)
2026-07-14 10:08:27 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc(1002685915 bytes)
2026-07-14 10:08:27 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc(1002685915 bytes)
2026-07-14 10:08:27 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc(1002685915 bytes)
2026-07-14 10:08:27 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc(1002685915 bytes)
2026-07-14 10:08:27 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc(1002685915 bytes)
2026-07-14 10:08:27 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc(1002685915 bytes)
2026-07-14 10:08:27 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc(1002685915 bytes)
2026-07-14 10:08:27 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc(1002685915 bytes)
2026-07-14 10:08:27 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc(1002685915 bytes)
2026-07-14 10:08:27 INFO  Run_model() completed successfully.
2026-07-14 10:08:27 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 10:08:27 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 11:00:00 days=153072 seconds=39600
2026-07-14 10:08:27 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-14 10:08:27 INFO [TIME] increment current_time 2020-02-06 09:00:00 -> 2020-02-06 11:00:00
2026-07-14 10:08:27 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 10:08:27 INFO ---------->>> Running process_satellite_data()
2026-07-14 10:08:28 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12003.nc
2026-07-14 10:08:28 INFO ---------->>> Running run_obs_converter()
2026-07-14 10:08:28 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_39045_153072.out
2026-07-14 10:08:28 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_39045_153072.out
2026-07-14 10:08:28 INFO ---------->>> Running DART
2026-07-14 10:08:28 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-14 10:08:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/chim_ENS1_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/chim_ENS1_2020020611_1_out_toDART.nc
2026-07-14 10:08:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/chim_ENS2_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/chim_ENS2_2020020611_1_out_toDART.nc
2026-07-14 10:08:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/chim_ENS3_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/chim_ENS3_2020020611_1_out_toDART.nc
2026-07-14 10:08:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/chim_ENS4_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/chim_ENS4_2020020611_1_out_toDART.nc
2026-07-14 10:08:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/chim_ENS5_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/chim_ENS5_2020020611_1_out_toDART.nc
2026-07-14 10:08:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/chim_ENS6_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/chim_ENS6_2020020611_1_out_toDART.nc
2026-07-14 10:08:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/chim_ENS7_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/chim_ENS7_2020020611_1_out_toDART.nc
2026-07-14 10:08:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/chim_ENS8_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/chim_ENS8_2020020611_1_out_toDART.nc
2026-07-14 10:08:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/chim_ENS9_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/chim_ENS9_2020020611_1_out_toDART.nc
2026-07-14 10:08:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/chim_ENS10_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/chim_ENS10_2020020611_1_out_toDART.nc
2026-07-14 10:08:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/chim_ENS11_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/chim_ENS11_2020020611_1_out_toDART.nc
2026-07-14 10:08:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/chim_ENS12_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/chim_ENS12_2020020611_1_out_toDART.nc
2026-07-14 10:08:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/chim_ENS13_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/chim_ENS13_2020020611_1_out_toDART.nc
2026-07-14 10:08:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/chim_ENS14_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/chim_ENS14_2020020611_1_out_toDART.nc
2026-07-14 10:08:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/chim_ENS15_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/chim_ENS15_2020020611_1_out_toDART.nc
2026-07-14 10:08:32 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-14 10:08:32 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-14 10:08:32 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-14 10:08:32 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-14 10:08:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-14 10:08:32 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-14 10:08:44 INFO Found: []
2026-07-14 10:08:44 INFO No job id returned by command ./run_filter.bsh
2026-07-14 10:08:44 INFO No monitoring will be performed
2026-07-14 10:08:44 INFO Moving DART output files to analysis and preassim directories for date 2020020611 if present ...
2026-07-14 10:08:44 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-14 10:08:44 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-14 10:08:44 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-14 10:08:44 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-14 10:08:44 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-14 10:08:44 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-14 10:08:44 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-14 10:08:44 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-14 10:08:44 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-14 10:08:44 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-14 10:08:44 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-14 10:08:44 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-14 10:08:44 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-14 10:08:44 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-14 10:08:44 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-14 10:08:44 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-14 10:08:44 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-14 10:08:44 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-14 10:08:44 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-14 10:08:44 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-14 10:08:44 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-14 10:08:44 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-14 10:08:44 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-14 10:08:44 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-14 10:08:44 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-14 10:08:44 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-14 10:08:44 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-14 10:08:44 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-14 10:08:44 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-14 10:08:44 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-14 10:08:44 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-14 10:08:45 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-14 10:08:45 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-14 10:08:45 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-14 10:08:45 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-14 10:08:45 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-14 10:08:45 INFO run_dart() is DONE.
2026-07-14 10:08:45 INFO ---------->>> Running update_pollutant_in_end()
2026-07-14 10:08:45 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:08:45 INFO No previous orbit memory found.
2026-07-14 10:08:45 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:08:45 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS1_2020020611.nc
2026-07-14 10:08:46 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:08:46 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:08:46 INFO No previous orbit memory found.
2026-07-14 10:08:46 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:08:46 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS1_2020020611.nc
2026-07-14 10:08:47 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:08:47 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:08:47 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:08:47 INFO No previous orbit memory found.
2026-07-14 10:08:47 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:08:47 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS2_2020020611.nc
2026-07-14 10:08:48 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:08:48 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:08:48 INFO No previous orbit memory found.
2026-07-14 10:08:48 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:08:48 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS2_2020020611.nc
2026-07-14 10:08:48 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:08:49 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:08:49 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:08:49 INFO No previous orbit memory found.
2026-07-14 10:08:49 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:08:49 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS3_2020020611.nc
2026-07-14 10:08:50 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:08:50 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:08:50 INFO No previous orbit memory found.
2026-07-14 10:08:50 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:08:50 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS3_2020020611.nc
2026-07-14 10:08:51 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:08:51 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:08:51 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:08:51 INFO No previous orbit memory found.
2026-07-14 10:08:51 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:08:51 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS4_2020020611.nc
2026-07-14 10:08:51 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:08:52 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:08:52 INFO No previous orbit memory found.
2026-07-14 10:08:52 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:08:52 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS4_2020020611.nc
2026-07-14 10:08:52 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:08:52 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:08:53 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:08:53 INFO No previous orbit memory found.
2026-07-14 10:08:53 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:08:53 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS5_2020020611.nc
2026-07-14 10:08:53 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:08:54 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:08:54 INFO No previous orbit memory found.
2026-07-14 10:08:54 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:08:54 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS5_2020020611.nc
2026-07-14 10:08:54 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:08:54 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:08:54 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:08:54 INFO No previous orbit memory found.
2026-07-14 10:08:54 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:08:54 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS6_2020020611.nc
2026-07-14 10:08:55 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:08:55 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:08:55 INFO No previous orbit memory found.
2026-07-14 10:08:55 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:08:55 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS6_2020020611.nc
2026-07-14 10:08:56 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:08:56 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:08:56 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:08:56 INFO No previous orbit memory found.
2026-07-14 10:08:56 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:08:56 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS7_2020020611.nc
2026-07-14 10:08:57 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:08:57 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:08:57 INFO No previous orbit memory found.
2026-07-14 10:08:57 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:08:57 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS7_2020020611.nc
2026-07-14 10:08:58 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:08:58 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:08:58 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:08:58 INFO No previous orbit memory found.
2026-07-14 10:08:58 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:08:58 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS8_2020020611.nc
2026-07-14 10:08:59 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:08:59 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:08:59 INFO No previous orbit memory found.
2026-07-14 10:08:59 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:08:59 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS8_2020020611.nc
2026-07-14 10:09:00 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:09:00 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:09:00 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:09:00 INFO No previous orbit memory found.
2026-07-14 10:09:00 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:09:00 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS9_2020020611.nc
2026-07-14 10:09:00 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:09:01 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:09:01 INFO No previous orbit memory found.
2026-07-14 10:09:01 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:09:01 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS9_2020020611.nc
2026-07-14 10:09:01 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:09:01 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:09:02 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:09:02 INFO No previous orbit memory found.
2026-07-14 10:09:02 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:09:02 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS10_2020020611.nc
2026-07-14 10:09:02 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:09:02 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:09:02 INFO No previous orbit memory found.
2026-07-14 10:09:02 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:09:02 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS10_2020020611.nc
2026-07-14 10:09:03 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:09:03 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:09:03 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:09:03 INFO No previous orbit memory found.
2026-07-14 10:09:03 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:09:03 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS11_2020020611.nc
2026-07-14 10:09:04 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:09:04 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:09:04 INFO No previous orbit memory found.
2026-07-14 10:09:04 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:09:04 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS11_2020020611.nc
2026-07-14 10:09:05 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:09:05 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:09:05 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:09:05 INFO No previous orbit memory found.
2026-07-14 10:09:05 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:09:05 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS12_2020020611.nc
2026-07-14 10:09:06 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:09:06 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:09:06 INFO No previous orbit memory found.
2026-07-14 10:09:06 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:09:06 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS12_2020020611.nc
2026-07-14 10:09:06 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:09:06 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:09:07 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:09:07 INFO No previous orbit memory found.
2026-07-14 10:09:07 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:09:07 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS13_2020020611.nc
2026-07-14 10:09:07 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:09:08 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:09:08 INFO No previous orbit memory found.
2026-07-14 10:09:08 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:09:08 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS13_2020020611.nc
2026-07-14 10:09:08 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:09:08 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:09:09 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:09:09 INFO No previous orbit memory found.
2026-07-14 10:09:09 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:09:09 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS14_2020020611.nc
2026-07-14 10:09:09 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:09:10 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:09:10 INFO No previous orbit memory found.
2026-07-14 10:09:10 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:09:10 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS14_2020020611.nc
2026-07-14 10:09:10 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:09:10 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:09:10 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:09:10 INFO No previous orbit memory found.
2026-07-14 10:09:10 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:09:10 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS15_2020020611.nc
2026-07-14 10:09:11 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:09:11 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:09:11 INFO No previous orbit memory found.
2026-07-14 10:09:11 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:09:11 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS15_2020020611.nc
2026-07-14 10:09:12 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:09:12 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:09:12 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-14 10:09:12 INFO [TIME] step_end current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 10:09:12 INFO [TIME] step_start current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 10:09:12 INFO [TIME] window start=2020-02-06 11:00:00 end=2020-02-06 13:00:00 run_hours=2 has_assimilation=True
2026-07-14 10:09:12 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:09:13 INFO Hourly dataset computed and listing created
2026-07-14 10:09:18 INFO Hourly dataset computed
2026-07-14 10:09:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:09:19 INFO Hourly dataset computed and listing created
2026-07-14 10:09:20 INFO Hourly dataset computed
2026-07-14 10:09:20 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:09:21 INFO Hourly dataset computed and listing created
2026-07-14 10:09:22 INFO Hourly dataset computed
2026-07-14 10:09:22 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:09:23 INFO Hourly dataset computed and listing created
2026-07-14 10:09:23 INFO Hourly dataset computed
2026-07-14 10:09:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:09:24 INFO Hourly dataset computed and listing created
2026-07-14 10:09:25 INFO Hourly dataset computed
2026-07-14 10:09:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:09:26 INFO Hourly dataset computed and listing created
2026-07-14 10:09:27 INFO Hourly dataset computed
2026-07-14 10:09:27 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:09:28 INFO Hourly dataset computed and listing created
2026-07-14 10:09:28 INFO Hourly dataset computed
2026-07-14 10:09:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:09:29 INFO Hourly dataset computed and listing created
2026-07-14 10:09:30 INFO Hourly dataset computed
2026-07-14 10:09:30 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:09:31 INFO Hourly dataset computed and listing created
2026-07-14 10:09:32 INFO Hourly dataset computed
2026-07-14 10:09:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:09:33 INFO Hourly dataset computed and listing created
2026-07-14 10:09:33 INFO Hourly dataset computed
2026-07-14 10:09:33 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:09:34 INFO Hourly dataset computed and listing created
2026-07-14 10:09:35 INFO Hourly dataset computed
2026-07-14 10:09:35 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:09:36 INFO Hourly dataset computed and listing created
2026-07-14 10:09:37 INFO Hourly dataset computed
2026-07-14 10:09:37 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:09:38 INFO Hourly dataset computed and listing created
2026-07-14 10:09:38 INFO Hourly dataset computed
2026-07-14 10:09:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:09:39 INFO Hourly dataset computed and listing created
2026-07-14 10:09:40 INFO Hourly dataset computed
2026-07-14 10:09:40 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:09:41 INFO Hourly dataset computed and listing created
2026-07-14 10:09:41 INFO Hourly dataset computed
2026-07-14 10:09:42 INFO ---------->>> Running CHIMERE model from 2020-02-06 11:00:00 to 2020-02-06 13:00:00
2026-07-14 10:09:42 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:09:42 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1
2026-07-14 10:09:42 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc
2026-07-14 10:09:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-14 10:09:42 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:09:42 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-14 10:09:42 INFO Queuing job for member 1...
2026-07-14 10:09:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:09:42 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-14 10:09:42 INFO Found: ['5154754']
2026-07-14 10:09:47 INFO [TGCC-IRENE] Submitted job with ID:['5154754']
2026-07-14 10:09:47 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:09:47 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2
2026-07-14 10:09:47 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc
2026-07-14 10:09:47 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-14 10:09:47 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:09:47 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-14 10:09:47 INFO Queuing job for member 2...
2026-07-14 10:09:47 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:09:47 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-14 10:09:48 INFO Found: ['5154756']
2026-07-14 10:09:53 INFO [TGCC-IRENE] Submitted job with ID:['5154756']
2026-07-14 10:09:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:09:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3
2026-07-14 10:09:53 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc
2026-07-14 10:09:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-14 10:09:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:09:53 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-14 10:09:53 INFO Queuing job for member 3...
2026-07-14 10:09:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:09:53 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-14 10:09:54 INFO Found: ['5154760']
2026-07-14 10:09:59 INFO [TGCC-IRENE] Submitted job with ID:['5154760']
2026-07-14 10:09:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:09:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4
2026-07-14 10:09:59 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc
2026-07-14 10:09:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-14 10:09:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:09:59 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-14 10:09:59 INFO Queuing job for member 4...
2026-07-14 10:09:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:09:59 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-14 10:10:00 INFO Found: ['5154762']
2026-07-14 10:10:05 INFO [TGCC-IRENE] Submitted job with ID:['5154762']
2026-07-14 10:10:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:10:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5
2026-07-14 10:10:05 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc
2026-07-14 10:10:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-14 10:10:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:10:05 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-14 10:10:05 INFO Queuing job for member 5...
2026-07-14 10:10:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:10:05 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-14 10:10:05 INFO Found: ['5154764']
2026-07-14 10:10:10 INFO [TGCC-IRENE] Submitted job with ID:['5154764']
2026-07-14 10:10:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:10:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6
2026-07-14 10:10:10 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc
2026-07-14 10:10:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-14 10:10:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:10:10 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-14 10:10:10 INFO Queuing job for member 6...
2026-07-14 10:10:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:10:10 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-14 10:10:12 INFO Found: ['5154765']
2026-07-14 10:10:17 INFO [TGCC-IRENE] Submitted job with ID:['5154765']
2026-07-14 10:10:17 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:10:17 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7
2026-07-14 10:10:17 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc
2026-07-14 10:10:17 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-14 10:10:17 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:10:17 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-14 10:10:17 INFO Queuing job for member 7...
2026-07-14 10:10:17 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:10:17 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-14 10:10:18 INFO Found: ['5154766']
2026-07-14 10:10:23 INFO [TGCC-IRENE] Submitted job with ID:['5154766']
2026-07-14 10:10:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:10:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8
2026-07-14 10:10:23 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc
2026-07-14 10:10:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-14 10:10:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:10:23 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-14 10:10:23 INFO Queuing job for member 8...
2026-07-14 10:10:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:10:23 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-14 10:10:24 INFO Found: ['5154767']
2026-07-14 10:10:29 INFO [TGCC-IRENE] Submitted job with ID:['5154767']
2026-07-14 10:10:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:10:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9
2026-07-14 10:10:29 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc
2026-07-14 10:10:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-14 10:10:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:10:29 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-14 10:10:29 INFO Queuing job for member 9...
2026-07-14 10:10:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:10:29 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-14 10:10:30 INFO Found: ['5154768']
2026-07-14 10:10:35 INFO [TGCC-IRENE] Submitted job with ID:['5154768']
2026-07-14 10:10:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:10:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10
2026-07-14 10:10:35 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc
2026-07-14 10:10:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-14 10:10:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:10:35 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-14 10:10:35 INFO Queuing job for member 10...
2026-07-14 10:10:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:10:35 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-14 10:10:36 INFO Found: ['5154769']
2026-07-14 10:10:41 INFO [TGCC-IRENE] Submitted job with ID:['5154769']
2026-07-14 10:10:41 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:10:41 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11
2026-07-14 10:10:41 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc
2026-07-14 10:10:41 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-14 10:10:41 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:10:41 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-14 10:10:41 INFO Queuing job for member 11...
2026-07-14 10:10:41 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:10:41 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-14 10:10:41 INFO Found: ['5154770']
2026-07-14 10:10:46 INFO [TGCC-IRENE] Submitted job with ID:['5154770']
2026-07-14 10:10:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:10:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12
2026-07-14 10:10:46 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc
2026-07-14 10:10:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-14 10:10:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:10:46 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-14 10:10:46 INFO Queuing job for member 12...
2026-07-14 10:10:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:10:46 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-14 10:10:47 INFO Found: ['5154771']
2026-07-14 10:10:52 INFO [TGCC-IRENE] Submitted job with ID:['5154771']
2026-07-14 10:10:52 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:10:52 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13
2026-07-14 10:10:52 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc
2026-07-14 10:10:52 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-14 10:10:52 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:10:52 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-14 10:10:52 INFO Queuing job for member 13...
2026-07-14 10:10:52 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:10:52 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-14 10:10:53 INFO Found: ['5154772']
2026-07-14 10:10:58 INFO [TGCC-IRENE] Submitted job with ID:['5154772']
2026-07-14 10:10:58 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:10:58 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14
2026-07-14 10:10:58 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc
2026-07-14 10:10:58 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-14 10:10:58 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:10:58 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-14 10:10:58 INFO Queuing job for member 14...
2026-07-14 10:10:58 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:10:58 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-14 10:15:11 INFO Found: ['5154786']
2026-07-14 10:15:16 INFO [TGCC-IRENE] Submitted job with ID:['5154786']
2026-07-14 10:15:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:15:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15
2026-07-14 10:15:16 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc
2026-07-14 10:15:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-14 10:15:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:15:16 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-14 10:15:16 INFO Queuing job for member 15...
2026-07-14 10:15:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:15:16 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-14 10:15:17 INFO Found: ['5154787']
2026-07-14 10:15:22 INFO [TGCC-IRENE] Submitted job with ID:['5154787']
2026-07-14 10:15:22 INFO Checking job status ...
2026-07-14 10:15:22 INFO None 5154754: status FINISHED
2026-07-14 10:15:22 INFO None 5154756: status FINISHED
2026-07-14 10:15:22 INFO None 5154760: status RUNNING/PENDING
2026-07-14 10:15:22 INFO None 5154762: status RUNNING/PENDING
2026-07-14 10:15:22 INFO None 5154764: status RUNNING/PENDING
2026-07-14 10:15:22 INFO None 5154765: status RUNNING/PENDING
2026-07-14 10:15:22 INFO None 5154766: status RUNNING/PENDING
2026-07-14 10:15:22 INFO None 5154767: status RUNNING/PENDING
2026-07-14 10:15:22 INFO None 5154768: status FINISHED
2026-07-14 10:15:22 INFO None 5154769: status RUNNING/PENDING
2026-07-14 10:15:22 INFO None 5154770: status RUNNING/PENDING
2026-07-14 10:15:22 INFO None 5154771: status RUNNING/PENDING
2026-07-14 10:15:22 INFO None 5154772: status RUNNING/PENDING
2026-07-14 10:15:22 INFO None 5154786: status RUNNING/PENDING
2026-07-14 10:15:22 INFO None 5154787: status RUNNING/PENDING
2026-07-14 10:15:22 INFO Jobs still running: ['5154760', '5154762', '5154764', '5154765', '5154766', '5154767', '5154769', '5154770', '5154771', '5154772', '5154786', '5154787']. Waiting...
2026-07-14 10:15:38 INFO None 5154754: status FINISHED
2026-07-14 10:15:38 INFO None 5154756: status FINISHED
2026-07-14 10:15:38 INFO None 5154760: status FINISHED
2026-07-14 10:15:38 INFO None 5154762: status FINISHED
2026-07-14 10:15:38 INFO None 5154764: status RUNNING/PENDING
2026-07-14 10:15:38 INFO None 5154765: status RUNNING/PENDING
2026-07-14 10:15:38 INFO None 5154766: status RUNNING/PENDING
2026-07-14 10:15:38 INFO None 5154767: status RUNNING/PENDING
2026-07-14 10:15:38 INFO None 5154768: status FINISHED
2026-07-14 10:15:38 INFO None 5154769: status RUNNING/PENDING
2026-07-14 10:15:38 INFO None 5154770: status RUNNING/PENDING
2026-07-14 10:15:38 INFO None 5154771: status RUNNING/PENDING
2026-07-14 10:15:38 INFO None 5154772: status RUNNING/PENDING
2026-07-14 10:15:38 INFO None 5154786: status RUNNING/PENDING
2026-07-14 10:15:38 INFO None 5154787: status RUNNING/PENDING
2026-07-14 10:15:38 INFO Jobs still running: ['5154764', '5154765', '5154766', '5154767', '5154769', '5154770', '5154771', '5154772', '5154786', '5154787']. Waiting...
2026-07-14 10:15:53 INFO None 5154754: status FINISHED
2026-07-14 10:15:53 INFO None 5154756: status FINISHED
2026-07-14 10:15:53 INFO None 5154760: status FINISHED
2026-07-14 10:15:53 INFO None 5154762: status FINISHED
2026-07-14 10:15:53 INFO None 5154764: status RUNNING/PENDING
2026-07-14 10:15:53 INFO None 5154765: status RUNNING/PENDING
2026-07-14 10:15:53 INFO None 5154766: status RUNNING/PENDING
2026-07-14 10:15:53 INFO None 5154767: status RUNNING/PENDING
2026-07-14 10:15:53 INFO None 5154768: status FINISHED
2026-07-14 10:15:53 INFO None 5154769: status RUNNING/PENDING
2026-07-14 10:15:53 INFO None 5154770: status FINISHED
2026-07-14 10:15:53 INFO None 5154771: status FINISHED
2026-07-14 10:15:53 INFO None 5154772: status RUNNING/PENDING
2026-07-14 10:15:53 INFO None 5154786: status RUNNING/PENDING
2026-07-14 10:15:53 INFO None 5154787: status RUNNING/PENDING
2026-07-14 10:15:53 INFO Jobs still running: ['5154764', '5154765', '5154766', '5154767', '5154769', '5154772', '5154786', '5154787']. Waiting...
2026-07-14 10:16:08 INFO None 5154754: status FINISHED
2026-07-14 10:16:08 INFO None 5154756: status FINISHED
2026-07-14 10:16:08 INFO None 5154760: status FINISHED
2026-07-14 10:16:08 INFO None 5154762: status FINISHED
2026-07-14 10:16:08 INFO None 5154764: status RUNNING/PENDING
2026-07-14 10:16:08 INFO None 5154765: status RUNNING/PENDING
2026-07-14 10:16:08 INFO None 5154766: status RUNNING/PENDING
2026-07-14 10:16:08 INFO None 5154767: status RUNNING/PENDING
2026-07-14 10:16:08 INFO None 5154768: status FINISHED
2026-07-14 10:16:08 INFO None 5154769: status FINISHED
2026-07-14 10:16:08 INFO None 5154770: status FINISHED
2026-07-14 10:16:08 INFO None 5154771: status FINISHED
2026-07-14 10:16:08 INFO None 5154772: status RUNNING/PENDING
2026-07-14 10:16:08 INFO None 5154786: status RUNNING/PENDING
2026-07-14 10:16:09 INFO None 5154787: status RUNNING/PENDING
2026-07-14 10:16:09 INFO Jobs still running: ['5154764', '5154765', '5154766', '5154767', '5154772', '5154786', '5154787']. Waiting...
2026-07-14 10:16:24 INFO None 5154754: status FINISHED
2026-07-14 10:16:24 INFO None 5154756: status FINISHED
2026-07-14 10:16:24 INFO None 5154760: status FINISHED
2026-07-14 10:16:24 INFO None 5154762: status FINISHED
2026-07-14 10:16:24 INFO None 5154764: status FINISHED
2026-07-14 10:16:24 INFO None 5154765: status FINISHED
2026-07-14 10:16:24 INFO None 5154766: status FINISHED
2026-07-14 10:16:24 INFO None 5154767: status RUNNING/PENDING
2026-07-14 10:16:24 INFO None 5154768: status FINISHED
2026-07-14 10:16:24 INFO None 5154769: status FINISHED
2026-07-14 10:16:24 INFO None 5154770: status FINISHED
2026-07-14 10:16:24 INFO None 5154771: status FINISHED
2026-07-14 10:16:24 INFO None 5154772: status FINISHED
2026-07-14 10:16:24 INFO None 5154786: status RUNNING/PENDING
2026-07-14 10:16:24 INFO None 5154787: status RUNNING/PENDING
2026-07-14 10:16:24 INFO Jobs still running: ['5154767', '5154786', '5154787']. Waiting...
2026-07-14 10:16:39 INFO None 5154754: status FINISHED
2026-07-14 10:16:39 INFO None 5154756: status FINISHED
2026-07-14 10:16:39 INFO None 5154760: status FINISHED
2026-07-14 10:16:39 INFO None 5154762: status FINISHED
2026-07-14 10:16:39 INFO None 5154764: status FINISHED
2026-07-14 10:16:39 INFO None 5154765: status FINISHED
2026-07-14 10:16:39 INFO None 5154766: status FINISHED
2026-07-14 10:16:39 INFO None 5154767: status RUNNING/PENDING
2026-07-14 10:16:39 INFO None 5154768: status FINISHED
2026-07-14 10:16:39 INFO None 5154769: status FINISHED
2026-07-14 10:16:39 INFO None 5154770: status FINISHED
2026-07-14 10:16:39 INFO None 5154771: status FINISHED
2026-07-14 10:16:39 INFO None 5154772: status FINISHED
2026-07-14 10:16:39 INFO None 5154786: status RUNNING/PENDING
2026-07-14 10:16:39 INFO None 5154787: status RUNNING/PENDING
2026-07-14 10:16:39 INFO Jobs still running: ['5154767', '5154786', '5154787']. Waiting...
2026-07-14 10:16:54 INFO None 5154754: status FINISHED
2026-07-14 10:16:54 INFO None 5154756: status FINISHED
2026-07-14 10:16:54 INFO None 5154760: status FINISHED
2026-07-14 10:16:54 INFO None 5154762: status FINISHED
2026-07-14 10:16:54 INFO None 5154764: status FINISHED
2026-07-14 10:16:54 INFO None 5154765: status FINISHED
2026-07-14 10:16:54 INFO None 5154766: status FINISHED
2026-07-14 10:16:54 INFO None 5154767: status FINISHED
2026-07-14 10:16:54 INFO None 5154768: status FINISHED
2026-07-14 10:16:54 INFO None 5154769: status FINISHED
2026-07-14 10:16:54 INFO None 5154770: status FINISHED
2026-07-14 10:16:54 INFO None 5154771: status FINISHED
2026-07-14 10:16:54 INFO None 5154772: status FINISHED
2026-07-14 10:16:54 INFO None 5154786: status RUNNING/PENDING
2026-07-14 10:16:54 INFO None 5154787: status RUNNING/PENDING
2026-07-14 10:16:54 INFO Jobs still running: ['5154786', '5154787']. Waiting...
2026-07-14 10:17:09 INFO None 5154754: status FINISHED
2026-07-14 10:17:09 INFO None 5154756: status FINISHED
2026-07-14 10:17:09 INFO None 5154760: status FINISHED
2026-07-14 10:17:09 INFO None 5154762: status FINISHED
2026-07-14 10:17:09 INFO None 5154764: status FINISHED
2026-07-14 10:17:09 INFO None 5154765: status FINISHED
2026-07-14 10:17:09 INFO None 5154766: status FINISHED
2026-07-14 10:17:09 INFO None 5154767: status FINISHED
2026-07-14 10:17:10 INFO None 5154768: status FINISHED
2026-07-14 10:17:10 INFO None 5154769: status FINISHED
2026-07-14 10:17:10 INFO None 5154770: status FINISHED
2026-07-14 10:17:10 INFO None 5154771: status FINISHED
2026-07-14 10:17:10 INFO None 5154772: status FINISHED
2026-07-14 10:17:10 INFO None 5154786: status RUNNING/PENDING
2026-07-14 10:17:10 INFO None 5154787: status RUNNING/PENDING
2026-07-14 10:17:10 INFO Jobs still running: ['5154786', '5154787']. Waiting...
2026-07-14 10:17:25 INFO None 5154754: status FINISHED
2026-07-14 10:17:25 INFO None 5154756: status FINISHED
2026-07-14 10:17:25 INFO None 5154760: status FINISHED
2026-07-14 10:17:25 INFO None 5154762: status FINISHED
2026-07-14 10:17:25 INFO None 5154764: status FINISHED
2026-07-14 10:17:25 INFO None 5154765: status FINISHED
2026-07-14 10:17:25 INFO None 5154766: status FINISHED
2026-07-14 10:17:25 INFO None 5154767: status FINISHED
2026-07-14 10:17:25 INFO None 5154768: status FINISHED
2026-07-14 10:17:25 INFO None 5154769: status FINISHED
2026-07-14 10:17:25 INFO None 5154770: status FINISHED
2026-07-14 10:17:25 INFO None 5154771: status FINISHED
2026-07-14 10:17:25 INFO None 5154772: status FINISHED
2026-07-14 10:17:25 INFO None 5154786: status RUNNING/PENDING
2026-07-14 10:17:25 INFO None 5154787: status RUNNING/PENDING
2026-07-14 10:17:25 INFO Jobs still running: ['5154786', '5154787']. Waiting...
2026-07-14 10:17:40 INFO None 5154754: status FINISHED
2026-07-14 10:17:40 INFO None 5154756: status FINISHED
2026-07-14 10:17:40 INFO None 5154760: status FINISHED
2026-07-14 10:17:40 INFO None 5154762: status FINISHED
2026-07-14 10:17:40 INFO None 5154764: status FINISHED
2026-07-14 10:17:40 INFO None 5154765: status FINISHED
2026-07-14 10:17:40 INFO None 5154766: status FINISHED
2026-07-14 10:17:40 INFO None 5154767: status FINISHED
2026-07-14 10:17:40 INFO None 5154768: status FINISHED
2026-07-14 10:17:40 INFO None 5154769: status FINISHED
2026-07-14 10:17:40 INFO None 5154770: status FINISHED
2026-07-14 10:17:40 INFO None 5154771: status FINISHED
2026-07-14 10:17:40 INFO None 5154772: status FINISHED
2026-07-14 10:17:40 INFO None 5154786: status RUNNING/PENDING
2026-07-14 10:17:40 INFO None 5154787: status RUNNING/PENDING
2026-07-14 10:17:40 INFO Jobs still running: ['5154786', '5154787']. Waiting...
2026-07-14 10:17:55 INFO None 5154754: status FINISHED
2026-07-14 10:17:55 INFO None 5154756: status FINISHED
2026-07-14 10:17:55 INFO None 5154760: status FINISHED
2026-07-14 10:17:55 INFO None 5154762: status FINISHED
2026-07-14 10:17:55 INFO None 5154764: status FINISHED
2026-07-14 10:17:55 INFO None 5154765: status FINISHED
2026-07-14 10:17:55 INFO None 5154766: status FINISHED
2026-07-14 10:17:55 INFO None 5154767: status FINISHED
2026-07-14 10:17:55 INFO None 5154768: status FINISHED
2026-07-14 10:17:55 INFO None 5154769: status FINISHED
2026-07-14 10:17:55 INFO None 5154770: status FINISHED
2026-07-14 10:17:55 INFO None 5154771: status FINISHED
2026-07-14 10:17:55 INFO None 5154772: status FINISHED
2026-07-14 10:17:55 INFO None 5154786: status RUNNING/PENDING
2026-07-14 10:17:55 INFO None 5154787: status RUNNING/PENDING
2026-07-14 10:17:55 INFO Jobs still running: ['5154786', '5154787']. Waiting...
2026-07-14 10:18:10 INFO None 5154754: status FINISHED
2026-07-14 10:18:11 INFO None 5154756: status FINISHED
2026-07-14 10:18:11 INFO None 5154760: status FINISHED
2026-07-14 10:18:11 INFO None 5154762: status FINISHED
2026-07-14 10:18:11 INFO None 5154764: status FINISHED
2026-07-14 10:18:11 INFO None 5154765: status FINISHED
2026-07-14 10:18:11 INFO None 5154766: status FINISHED
2026-07-14 10:18:11 INFO None 5154767: status FINISHED
2026-07-14 10:18:11 INFO None 5154768: status FINISHED
2026-07-14 10:18:11 INFO None 5154769: status FINISHED
2026-07-14 10:18:11 INFO None 5154770: status FINISHED
2026-07-14 10:18:11 INFO None 5154771: status FINISHED
2026-07-14 10:18:11 INFO None 5154772: status FINISHED
2026-07-14 10:18:11 INFO None 5154786: status RUNNING/PENDING
2026-07-14 10:18:11 INFO None 5154787: status RUNNING/PENDING
2026-07-14 10:18:11 INFO Jobs still running: ['5154786', '5154787']. Waiting...
2026-07-14 10:18:26 INFO None 5154754: status FINISHED
2026-07-14 10:18:26 INFO None 5154756: status FINISHED
2026-07-14 10:18:26 INFO None 5154760: status FINISHED
2026-07-14 10:18:26 INFO None 5154762: status FINISHED
2026-07-14 10:18:26 INFO None 5154764: status FINISHED
2026-07-14 10:18:26 INFO None 5154765: status FINISHED
2026-07-14 10:18:26 INFO None 5154766: status FINISHED
2026-07-14 10:18:26 INFO None 5154767: status FINISHED
2026-07-14 10:18:26 INFO None 5154768: status FINISHED
2026-07-14 10:18:26 INFO None 5154769: status FINISHED
2026-07-14 10:18:26 INFO None 5154770: status FINISHED
2026-07-14 10:18:26 INFO None 5154771: status FINISHED
2026-07-14 10:18:26 INFO None 5154772: status FINISHED
2026-07-14 10:18:26 INFO None 5154786: status RUNNING/PENDING
2026-07-14 10:18:27 INFO None 5154787: status RUNNING/PENDING
2026-07-14 10:18:27 INFO Jobs still running: ['5154786', '5154787']. Waiting...
2026-07-14 10:18:42 INFO None 5154754: status FINISHED
2026-07-14 10:18:42 INFO None 5154756: status FINISHED
2026-07-14 10:18:42 INFO None 5154760: status FINISHED
2026-07-14 10:18:42 INFO None 5154762: status FINISHED
2026-07-14 10:18:42 INFO None 5154764: status FINISHED
2026-07-14 10:18:42 INFO None 5154765: status FINISHED
2026-07-14 10:18:42 INFO None 5154766: status FINISHED
2026-07-14 10:18:42 INFO None 5154767: status FINISHED
2026-07-14 10:18:42 INFO None 5154768: status FINISHED
2026-07-14 10:18:42 INFO None 5154769: status FINISHED
2026-07-14 10:18:42 INFO None 5154770: status FINISHED
2026-07-14 10:18:42 INFO None 5154771: status FINISHED
2026-07-14 10:18:42 INFO None 5154772: status FINISHED
2026-07-14 10:18:42 INFO None 5154786: status RUNNING/PENDING
2026-07-14 10:18:42 INFO None 5154787: status RUNNING/PENDING
2026-07-14 10:18:42 INFO Jobs still running: ['5154786', '5154787']. Waiting...
2026-07-14 10:18:57 INFO None 5154754: status FINISHED
2026-07-14 10:18:57 INFO None 5154756: status FINISHED
2026-07-14 10:18:57 INFO None 5154760: status FINISHED
2026-07-14 10:18:57 INFO None 5154762: status FINISHED
2026-07-14 10:18:57 INFO None 5154764: status FINISHED
2026-07-14 10:18:57 INFO None 5154765: status FINISHED
2026-07-14 10:18:57 INFO None 5154766: status FINISHED
2026-07-14 10:18:57 INFO None 5154767: status FINISHED
2026-07-14 10:18:57 INFO None 5154768: status FINISHED
2026-07-14 10:18:57 INFO None 5154769: status FINISHED
2026-07-14 10:18:57 INFO None 5154770: status FINISHED
2026-07-14 10:18:57 INFO None 5154771: status FINISHED
2026-07-14 10:18:57 INFO None 5154772: status FINISHED
2026-07-14 10:18:57 INFO None 5154786: status RUNNING/PENDING
2026-07-14 10:18:57 INFO None 5154787: status RUNNING/PENDING
2026-07-14 10:18:57 INFO Jobs still running: ['5154786', '5154787']. Waiting...
2026-07-14 10:19:12 INFO None 5154754: status FINISHED
2026-07-14 10:19:12 INFO None 5154756: status FINISHED
2026-07-14 10:19:12 INFO None 5154760: status FINISHED
2026-07-14 10:19:12 INFO None 5154762: status FINISHED
2026-07-14 10:19:12 INFO None 5154764: status FINISHED
2026-07-14 10:19:12 INFO None 5154765: status FINISHED
2026-07-14 10:19:12 INFO None 5154766: status FINISHED
2026-07-14 10:19:12 INFO None 5154767: status FINISHED
2026-07-14 10:19:12 INFO None 5154768: status FINISHED
2026-07-14 10:19:12 INFO None 5154769: status FINISHED
2026-07-14 10:19:12 INFO None 5154770: status FINISHED
2026-07-14 10:19:12 INFO None 5154771: status FINISHED
2026-07-14 10:19:12 INFO None 5154772: status FINISHED
2026-07-14 10:19:12 INFO None 5154786: status RUNNING/PENDING
2026-07-14 10:19:12 INFO None 5154787: status RUNNING/PENDING
2026-07-14 10:19:12 INFO Jobs still running: ['5154786', '5154787']. Waiting...
2026-07-14 10:19:27 INFO None 5154754: status FINISHED
2026-07-14 10:19:27 INFO None 5154756: status FINISHED
2026-07-14 10:19:27 INFO None 5154760: status FINISHED
2026-07-14 10:19:28 INFO None 5154762: status FINISHED
2026-07-14 10:19:28 INFO None 5154764: status FINISHED
2026-07-14 10:19:28 INFO None 5154765: status FINISHED
2026-07-14 10:19:28 INFO None 5154766: status FINISHED
2026-07-14 10:19:28 INFO None 5154767: status FINISHED
2026-07-14 10:19:28 INFO None 5154768: status FINISHED
2026-07-14 10:19:28 INFO None 5154769: status FINISHED
2026-07-14 10:19:28 INFO None 5154770: status FINISHED
2026-07-14 10:19:28 INFO None 5154771: status FINISHED
2026-07-14 10:19:28 INFO None 5154772: status FINISHED
2026-07-14 10:19:28 INFO None 5154786: status RUNNING/PENDING
2026-07-14 10:19:28 INFO None 5154787: status RUNNING/PENDING
2026-07-14 10:19:28 INFO Jobs still running: ['5154786', '5154787']. Waiting...
2026-07-14 10:19:43 INFO None 5154754: status FINISHED
2026-07-14 10:19:43 INFO None 5154756: status FINISHED
2026-07-14 10:19:43 INFO None 5154760: status FINISHED
2026-07-14 10:19:43 INFO None 5154762: status FINISHED
2026-07-14 10:19:43 INFO None 5154764: status FINISHED
2026-07-14 10:19:43 INFO None 5154765: status FINISHED
2026-07-14 10:19:43 INFO None 5154766: status FINISHED
2026-07-14 10:19:43 INFO None 5154767: status FINISHED
2026-07-14 10:19:43 INFO None 5154768: status FINISHED
2026-07-14 10:19:43 INFO None 5154769: status FINISHED
2026-07-14 10:19:43 INFO None 5154770: status FINISHED
2026-07-14 10:19:43 INFO None 5154771: status FINISHED
2026-07-14 10:19:43 INFO None 5154772: status FINISHED
2026-07-14 10:19:43 INFO None 5154786: status RUNNING/PENDING
2026-07-14 10:19:43 INFO None 5154787: status RUNNING/PENDING
2026-07-14 10:19:43 INFO Jobs still running: ['5154786', '5154787']. Waiting...
2026-07-14 10:19:58 INFO None 5154754: status FINISHED
2026-07-14 10:19:58 INFO None 5154756: status FINISHED
2026-07-14 10:19:58 INFO None 5154760: status FINISHED
2026-07-14 10:19:58 INFO None 5154762: status FINISHED
2026-07-14 10:19:58 INFO None 5154764: status FINISHED
2026-07-14 10:19:58 INFO None 5154765: status FINISHED
2026-07-14 10:19:58 INFO None 5154766: status FINISHED
2026-07-14 10:19:58 INFO None 5154767: status FINISHED
2026-07-14 10:19:58 INFO None 5154768: status FINISHED
2026-07-14 10:19:58 INFO None 5154769: status FINISHED
2026-07-14 10:19:58 INFO None 5154770: status FINISHED
2026-07-14 10:19:58 INFO None 5154771: status FINISHED
2026-07-14 10:19:58 INFO None 5154772: status FINISHED
2026-07-14 10:19:58 INFO None 5154786: status RUNNING/PENDING
2026-07-14 10:19:58 INFO None 5154787: status RUNNING/PENDING
2026-07-14 10:19:58 INFO Jobs still running: ['5154786', '5154787']. Waiting...
2026-07-14 10:20:13 INFO None 5154754: status FINISHED
2026-07-14 10:20:13 INFO None 5154756: status FINISHED
2026-07-14 10:20:13 INFO None 5154760: status FINISHED
2026-07-14 10:20:13 INFO None 5154762: status FINISHED
2026-07-14 10:20:13 INFO None 5154764: status FINISHED
2026-07-14 10:20:13 INFO None 5154765: status FINISHED
2026-07-14 10:20:13 INFO None 5154766: status FINISHED
2026-07-14 10:20:13 INFO None 5154767: status FINISHED
2026-07-14 10:20:13 INFO None 5154768: status FINISHED
2026-07-14 10:20:13 INFO None 5154769: status FINISHED
2026-07-14 10:20:13 INFO None 5154770: status FINISHED
2026-07-14 10:20:13 INFO None 5154771: status FINISHED
2026-07-14 10:20:14 INFO None 5154772: status FINISHED
2026-07-14 10:20:14 INFO None 5154786: status RUNNING/PENDING
2026-07-14 10:20:14 INFO None 5154787: status RUNNING/PENDING
2026-07-14 10:20:14 INFO Jobs still running: ['5154786', '5154787']. Waiting...
2026-07-14 10:20:29 INFO None 5154754: status FINISHED
2026-07-14 10:20:29 INFO None 5154756: status FINISHED
2026-07-14 10:20:29 INFO None 5154760: status FINISHED
2026-07-14 10:20:29 INFO None 5154762: status FINISHED
2026-07-14 10:20:29 INFO None 5154764: status FINISHED
2026-07-14 10:20:29 INFO None 5154765: status FINISHED
2026-07-14 10:20:29 INFO None 5154766: status FINISHED
2026-07-14 10:20:29 INFO None 5154767: status FINISHED
2026-07-14 10:20:29 INFO None 5154768: status FINISHED
2026-07-14 10:20:29 INFO None 5154769: status FINISHED
2026-07-14 10:20:29 INFO None 5154770: status FINISHED
2026-07-14 10:20:29 INFO None 5154771: status FINISHED
2026-07-14 10:20:29 INFO None 5154772: status FINISHED
2026-07-14 10:20:29 INFO None 5154786: status RUNNING/PENDING
2026-07-14 10:20:29 INFO None 5154787: status RUNNING/PENDING
2026-07-14 10:20:29 INFO Jobs still running: ['5154786', '5154787']. Waiting...
2026-07-14 10:20:44 INFO None 5154754: status FINISHED
2026-07-14 10:20:44 INFO None 5154756: status FINISHED
2026-07-14 10:20:44 INFO None 5154760: status FINISHED
2026-07-14 10:20:44 INFO None 5154762: status FINISHED
2026-07-14 10:20:44 INFO None 5154764: status FINISHED
2026-07-14 10:20:44 INFO None 5154765: status FINISHED
2026-07-14 10:20:44 INFO None 5154766: status FINISHED
2026-07-14 10:20:44 INFO None 5154767: status FINISHED
2026-07-14 10:20:44 INFO None 5154768: status FINISHED
2026-07-14 10:20:44 INFO None 5154769: status FINISHED
2026-07-14 10:20:44 INFO None 5154770: status FINISHED
2026-07-14 10:20:44 INFO None 5154771: status FINISHED
2026-07-14 10:20:44 INFO None 5154772: status FINISHED
2026-07-14 10:20:44 INFO None 5154786: status FINISHED
2026-07-14 10:20:44 INFO None 5154787: status FINISHED
2026-07-14 10:20:44 INFO Jobs ['5154754', '5154756', '5154760', '5154762', '5154764', '5154765', '5154766', '5154767', '5154768', '5154769', '5154770', '5154771', '5154772', '5154786', '5154787'] have finished
2026-07-14 10:20:44 INFO Checking restart files were created ...
2026-07-14 10:20:44 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020611_2_ENS1.nc(1002685915 bytes)
2026-07-14 10:20:44 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020611_2_ENS2.nc(1002685915 bytes)
2026-07-14 10:20:44 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020611_2_ENS3.nc(1002685915 bytes)
2026-07-14 10:20:44 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020611_2_ENS4.nc(1002685915 bytes)
2026-07-14 10:20:44 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020611_2_ENS5.nc(1002685915 bytes)
2026-07-14 10:20:44 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020611_2_ENS6.nc(1002685915 bytes)
2026-07-14 10:20:44 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020611_2_ENS7.nc(1002685915 bytes)
2026-07-14 10:20:44 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020611_2_ENS8.nc(1002685915 bytes)
2026-07-14 10:20:44 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020611_2_ENS9.nc(1002685915 bytes)
2026-07-14 10:20:44 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020611_2_ENS10.nc(1002685915 bytes)
2026-07-14 10:20:44 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020611_2_ENS11.nc(1002685915 bytes)
2026-07-14 10:20:44 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020611_2_ENS12.nc(1002685915 bytes)
2026-07-14 10:20:44 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020611_2_ENS13.nc(1002685915 bytes)
2026-07-14 10:20:44 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020611_2_ENS14.nc(1002685915 bytes)
2026-07-14 10:20:44 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020611_2_ENS15.nc(1002685915 bytes)
2026-07-14 10:20:44 INFO  Run_model() completed successfully.
2026-07-14 10:20:44 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 10:20:44 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 13:00:00 days=153072 seconds=46800
2026-07-14 10:20:44 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-14 10:20:44 INFO [TIME] increment current_time 2020-02-06 11:00:00 -> 2020-02-06 13:00:00
2026-07-14 10:20:44 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 10:20:44 INFO ---------->>> Running process_satellite_data()
2026-07-14 10:20:44 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12004.nc
2026-07-14 10:20:44 INFO ---------->>> Running run_obs_converter()
2026-07-14 10:20:44 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_45135_153072.out
2026-07-14 10:20:44 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_45135_153072.out
2026-07-14 10:20:44 INFO ---------->>> Running DART
2026-07-14 10:20:44 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-14 10:20:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/chim_ENS1_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/chim_ENS1_2020020613_1_out_toDART.nc
2026-07-14 10:20:45 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/chim_ENS2_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/chim_ENS2_2020020613_1_out_toDART.nc
2026-07-14 10:20:45 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/chim_ENS3_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/chim_ENS3_2020020613_1_out_toDART.nc
2026-07-14 10:20:45 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/chim_ENS4_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/chim_ENS4_2020020613_1_out_toDART.nc
2026-07-14 10:20:46 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/chim_ENS5_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/chim_ENS5_2020020613_1_out_toDART.nc
2026-07-14 10:20:46 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/chim_ENS6_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/chim_ENS6_2020020613_1_out_toDART.nc
2026-07-14 10:20:46 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/chim_ENS7_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/chim_ENS7_2020020613_1_out_toDART.nc
2026-07-14 10:20:47 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/chim_ENS8_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/chim_ENS8_2020020613_1_out_toDART.nc
2026-07-14 10:20:47 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/chim_ENS9_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/chim_ENS9_2020020613_1_out_toDART.nc
2026-07-14 10:20:47 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/chim_ENS10_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/chim_ENS10_2020020613_1_out_toDART.nc
2026-07-14 10:20:48 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/chim_ENS11_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/chim_ENS11_2020020613_1_out_toDART.nc
2026-07-14 10:20:48 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/chim_ENS12_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/chim_ENS12_2020020613_1_out_toDART.nc
2026-07-14 10:20:48 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/chim_ENS13_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/chim_ENS13_2020020613_1_out_toDART.nc
2026-07-14 10:20:49 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/chim_ENS14_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/chim_ENS14_2020020613_1_out_toDART.nc
2026-07-14 10:20:49 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/chim_ENS15_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/chim_ENS15_2020020613_1_out_toDART.nc
2026-07-14 10:20:49 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-14 10:20:49 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-14 10:20:49 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-14 10:20:49 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-14 10:20:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-14 10:20:49 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-14 10:23:38 INFO Found: []
2026-07-14 10:23:38 INFO No job id returned by command ./run_filter.bsh
2026-07-14 10:23:38 INFO No monitoring will be performed
2026-07-14 10:23:38 INFO Moving DART output files to analysis and preassim directories for date 2020020613 if present ...
2026-07-14 10:23:38 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-14 10:23:38 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-14 10:23:38 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-14 10:23:38 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-14 10:23:38 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-14 10:23:38 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-14 10:23:38 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-14 10:23:38 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-14 10:23:38 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-14 10:23:39 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-14 10:23:39 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-14 10:23:39 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-14 10:23:39 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-14 10:23:39 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-14 10:23:39 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-14 10:23:39 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-14 10:23:39 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-14 10:23:39 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-14 10:23:39 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-14 10:23:39 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-14 10:23:39 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-14 10:23:39 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-14 10:23:39 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-14 10:23:39 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-14 10:23:39 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-14 10:23:39 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-14 10:23:39 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-14 10:23:39 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-14 10:23:39 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-14 10:23:39 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-14 10:23:39 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-14 10:23:39 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-14 10:23:39 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-14 10:23:39 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-14 10:23:39 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-14 10:23:39 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-14 10:23:39 INFO run_dart() is DONE.
2026-07-14 10:23:39 INFO ---------->>> Running update_pollutant_in_end()
2026-07-14 10:23:39 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:23:39 INFO No previous orbit memory found.
2026-07-14 10:23:39 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:23:39 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS1_2020020613.nc
2026-07-14 10:23:40 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:23:40 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:23:40 INFO No previous orbit memory found.
2026-07-14 10:23:40 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:23:40 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS1_2020020613.nc
2026-07-14 10:23:41 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:23:41 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:23:41 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:23:41 INFO No previous orbit memory found.
2026-07-14 10:23:41 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:23:41 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS2_2020020613.nc
2026-07-14 10:23:42 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:23:42 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:23:42 INFO No previous orbit memory found.
2026-07-14 10:23:42 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:23:42 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS2_2020020613.nc
2026-07-14 10:23:42 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:23:42 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:23:43 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:23:43 INFO No previous orbit memory found.
2026-07-14 10:23:43 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:23:43 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS3_2020020613.nc
2026-07-14 10:23:43 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:23:44 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:23:44 INFO No previous orbit memory found.
2026-07-14 10:23:44 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:23:44 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS3_2020020613.nc
2026-07-14 10:23:44 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:23:44 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:23:44 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:23:44 INFO No previous orbit memory found.
2026-07-14 10:23:44 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:23:44 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS4_2020020613.nc
2026-07-14 10:23:45 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:23:45 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:23:45 INFO No previous orbit memory found.
2026-07-14 10:23:45 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:23:45 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS4_2020020613.nc
2026-07-14 10:23:46 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:23:46 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:23:46 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:23:46 INFO No previous orbit memory found.
2026-07-14 10:23:46 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:23:46 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS5_2020020613.nc
2026-07-14 10:23:47 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:23:47 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:23:47 INFO No previous orbit memory found.
2026-07-14 10:23:47 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:23:47 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS5_2020020613.nc
2026-07-14 10:23:47 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:23:47 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:23:48 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:23:48 INFO No previous orbit memory found.
2026-07-14 10:23:48 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:23:48 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS6_2020020613.nc
2026-07-14 10:23:48 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:23:49 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:23:49 INFO No previous orbit memory found.
2026-07-14 10:23:49 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:23:49 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS6_2020020613.nc
2026-07-14 10:23:49 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:23:49 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:23:50 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:23:50 INFO No previous orbit memory found.
2026-07-14 10:23:50 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:23:50 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS7_2020020613.nc
2026-07-14 10:23:50 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:23:51 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:23:51 INFO No previous orbit memory found.
2026-07-14 10:23:51 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:23:51 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS7_2020020613.nc
2026-07-14 10:23:51 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:23:51 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:23:51 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:23:51 INFO No previous orbit memory found.
2026-07-14 10:23:51 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:23:51 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS8_2020020613.nc
2026-07-14 10:23:52 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:23:52 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:23:52 INFO No previous orbit memory found.
2026-07-14 10:23:52 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:23:52 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS8_2020020613.nc
2026-07-14 10:23:53 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:23:53 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:23:53 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:23:53 INFO No previous orbit memory found.
2026-07-14 10:23:53 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:23:53 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS9_2020020613.nc
2026-07-14 10:23:54 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:23:54 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:23:54 INFO No previous orbit memory found.
2026-07-14 10:23:54 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:23:54 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS9_2020020613.nc
2026-07-14 10:23:55 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:23:55 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:23:55 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:23:55 INFO No previous orbit memory found.
2026-07-14 10:23:55 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:23:55 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS10_2020020613.nc
2026-07-14 10:23:55 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:23:56 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:23:56 INFO No previous orbit memory found.
2026-07-14 10:23:56 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:23:56 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS10_2020020613.nc
2026-07-14 10:23:56 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:23:56 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:23:57 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:23:57 INFO No previous orbit memory found.
2026-07-14 10:23:57 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:23:57 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS11_2020020613.nc
2026-07-14 10:23:57 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:23:57 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:23:57 INFO No previous orbit memory found.
2026-07-14 10:23:57 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:23:57 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS11_2020020613.nc
2026-07-14 10:23:58 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:23:58 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:23:58 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:23:58 INFO No previous orbit memory found.
2026-07-14 10:23:58 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:23:58 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS12_2020020613.nc
2026-07-14 10:23:59 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:23:59 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:23:59 INFO No previous orbit memory found.
2026-07-14 10:23:59 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:23:59 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS12_2020020613.nc
2026-07-14 10:24:00 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:24:00 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:24:00 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:24:00 INFO No previous orbit memory found.
2026-07-14 10:24:00 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:24:00 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS13_2020020613.nc
2026-07-14 10:24:01 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:24:01 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:24:01 INFO No previous orbit memory found.
2026-07-14 10:24:01 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:24:01 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS13_2020020613.nc
2026-07-14 10:24:01 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:24:01 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:24:02 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:24:02 INFO No previous orbit memory found.
2026-07-14 10:24:02 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:24:02 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS14_2020020613.nc
2026-07-14 10:24:02 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:24:03 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:24:03 INFO No previous orbit memory found.
2026-07-14 10:24:03 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:24:03 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS14_2020020613.nc
2026-07-14 10:24:03 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:24:03 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:24:04 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:24:04 INFO No previous orbit memory found.
2026-07-14 10:24:04 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:24:04 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS15_2020020613.nc
2026-07-14 10:24:04 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:24:04 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:24:04 INFO No previous orbit memory found.
2026-07-14 10:24:04 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:24:04 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS15_2020020613.nc
2026-07-14 10:24:05 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:24:05 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:24:05 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-14 10:24:05 INFO [TIME] step_end current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 10:24:05 INFO [TIME] step_start current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 10:24:05 INFO [TIME] window start=2020-02-06 13:00:00 end=2020-02-06 14:00:00 run_hours=1 has_assimilation=True
2026-07-14 10:24:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:24:06 INFO Hourly dataset computed and listing created
2026-07-14 10:24:08 INFO Hourly dataset computed
2026-07-14 10:24:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:24:09 INFO Hourly dataset computed and listing created
2026-07-14 10:24:10 INFO Hourly dataset computed
2026-07-14 10:24:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:24:10 INFO Hourly dataset computed and listing created
2026-07-14 10:24:11 INFO Hourly dataset computed
2026-07-14 10:24:11 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:24:12 INFO Hourly dataset computed and listing created
2026-07-14 10:24:12 INFO Hourly dataset computed
2026-07-14 10:24:12 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:24:13 INFO Hourly dataset computed and listing created
2026-07-14 10:24:14 INFO Hourly dataset computed
2026-07-14 10:24:14 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:24:15 INFO Hourly dataset computed and listing created
2026-07-14 10:24:15 INFO Hourly dataset computed
2026-07-14 10:24:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:24:16 INFO Hourly dataset computed and listing created
2026-07-14 10:24:17 INFO Hourly dataset computed
2026-07-14 10:24:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:24:18 INFO Hourly dataset computed and listing created
2026-07-14 10:24:18 INFO Hourly dataset computed
2026-07-14 10:24:18 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:24:19 INFO Hourly dataset computed and listing created
2026-07-14 10:24:19 INFO Hourly dataset computed
2026-07-14 10:24:20 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:24:20 INFO Hourly dataset computed and listing created
2026-07-14 10:24:21 INFO Hourly dataset computed
2026-07-14 10:24:21 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:24:22 INFO Hourly dataset computed and listing created
2026-07-14 10:24:22 INFO Hourly dataset computed
2026-07-14 10:24:22 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:24:23 INFO Hourly dataset computed and listing created
2026-07-14 10:24:24 INFO Hourly dataset computed
2026-07-14 10:24:24 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:24:25 INFO Hourly dataset computed and listing created
2026-07-14 10:24:25 INFO Hourly dataset computed
2026-07-14 10:24:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:24:26 INFO Hourly dataset computed and listing created
2026-07-14 10:24:27 INFO Hourly dataset computed
2026-07-14 10:24:27 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:24:28 INFO Hourly dataset computed and listing created
2026-07-14 10:24:28 INFO Hourly dataset computed
2026-07-14 10:24:28 INFO ---------->>> Running CHIMERE model from 2020-02-06 13:00:00 to 2020-02-06 14:00:00
2026-07-14 10:24:28 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:24:28 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1
2026-07-14 10:24:28 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020611_2_ENS1.nc
2026-07-14 10:24:28 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-14 10:24:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:24:28 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-14 10:24:28 INFO Queuing job for member 1...
2026-07-14 10:24:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:24:28 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-14 10:24:29 INFO Found: ['5154819']
2026-07-14 10:24:34 INFO [TGCC-IRENE] Submitted job with ID:['5154819']
2026-07-14 10:24:34 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:24:34 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2
2026-07-14 10:24:34 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020611_2_ENS2.nc
2026-07-14 10:24:34 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-14 10:24:34 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:24:34 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-14 10:24:34 INFO Queuing job for member 2...
2026-07-14 10:24:34 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:24:34 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-14 10:24:35 INFO Found: ['5154820']
2026-07-14 10:24:40 INFO [TGCC-IRENE] Submitted job with ID:['5154820']
2026-07-14 10:24:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:24:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3
2026-07-14 10:24:40 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020611_2_ENS3.nc
2026-07-14 10:24:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-14 10:24:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:24:40 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-14 10:24:40 INFO Queuing job for member 3...
2026-07-14 10:24:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:24:40 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-14 10:24:40 INFO Found: ['5154821']
2026-07-14 10:24:45 INFO [TGCC-IRENE] Submitted job with ID:['5154821']
2026-07-14 10:24:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:24:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4
2026-07-14 10:24:45 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020611_2_ENS4.nc
2026-07-14 10:24:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-14 10:24:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:24:45 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-14 10:24:45 INFO Queuing job for member 4...
2026-07-14 10:24:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:24:45 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-14 10:24:46 INFO Found: ['5154822']
2026-07-14 10:24:51 INFO [TGCC-IRENE] Submitted job with ID:['5154822']
2026-07-14 10:24:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:24:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5
2026-07-14 10:24:51 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020611_2_ENS5.nc
2026-07-14 10:24:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-14 10:24:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:24:51 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-14 10:24:51 INFO Queuing job for member 5...
2026-07-14 10:24:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:24:51 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-14 10:24:52 INFO Found: ['5154824']
2026-07-14 10:24:57 INFO [TGCC-IRENE] Submitted job with ID:['5154824']
2026-07-14 10:24:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:24:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6
2026-07-14 10:24:57 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020611_2_ENS6.nc
2026-07-14 10:24:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-14 10:24:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:24:57 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-14 10:24:57 INFO Queuing job for member 6...
2026-07-14 10:24:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:24:57 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-14 10:24:58 INFO Found: ['5154827']
2026-07-14 10:25:03 INFO [TGCC-IRENE] Submitted job with ID:['5154827']
2026-07-14 10:25:03 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:25:03 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7
2026-07-14 10:25:03 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020611_2_ENS7.nc
2026-07-14 10:25:03 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-14 10:25:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:25:03 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-14 10:25:03 INFO Queuing job for member 7...
2026-07-14 10:25:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:25:03 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-14 10:25:04 INFO Found: ['5154830']
2026-07-14 10:25:09 INFO [TGCC-IRENE] Submitted job with ID:['5154830']
2026-07-14 10:25:09 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:25:09 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8
2026-07-14 10:25:09 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020611_2_ENS8.nc
2026-07-14 10:25:09 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-14 10:25:09 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:25:09 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-14 10:25:09 INFO Queuing job for member 8...
2026-07-14 10:25:09 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:25:09 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-14 10:25:10 INFO Found: ['5154831']
2026-07-14 10:25:15 INFO [TGCC-IRENE] Submitted job with ID:['5154831']
2026-07-14 10:25:15 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:25:15 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9
2026-07-14 10:25:15 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020611_2_ENS9.nc
2026-07-14 10:25:15 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-14 10:25:15 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:25:15 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-14 10:25:15 INFO Queuing job for member 9...
2026-07-14 10:25:15 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:25:15 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-14 10:25:15 INFO Found: ['5154832']
2026-07-14 10:25:20 INFO [TGCC-IRENE] Submitted job with ID:['5154832']
2026-07-14 10:25:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:25:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10
2026-07-14 10:25:20 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020611_2_ENS10.nc
2026-07-14 10:25:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-14 10:25:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:25:20 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-14 10:25:20 INFO Queuing job for member 10...
2026-07-14 10:25:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:25:20 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-14 10:25:22 INFO Found: ['5154833']
2026-07-14 10:25:27 INFO [TGCC-IRENE] Submitted job with ID:['5154833']
2026-07-14 10:25:27 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:25:27 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11
2026-07-14 10:25:27 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020611_2_ENS11.nc
2026-07-14 10:25:27 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-14 10:25:27 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:25:27 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-14 10:25:27 INFO Queuing job for member 11...
2026-07-14 10:25:27 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:25:27 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-14 10:25:27 INFO Found: ['5154834']
2026-07-14 10:25:32 INFO [TGCC-IRENE] Submitted job with ID:['5154834']
2026-07-14 10:25:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:25:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12
2026-07-14 10:25:32 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020611_2_ENS12.nc
2026-07-14 10:25:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-14 10:25:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:25:32 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-14 10:25:32 INFO Queuing job for member 12...
2026-07-14 10:25:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:25:32 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-14 10:25:33 INFO Found: ['5154835']
2026-07-14 10:25:38 INFO [TGCC-IRENE] Submitted job with ID:['5154835']
2026-07-14 10:25:38 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:25:38 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13
2026-07-14 10:25:38 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020611_2_ENS13.nc
2026-07-14 10:25:38 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-14 10:25:38 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:25:38 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-14 10:25:38 INFO Queuing job for member 13...
2026-07-14 10:25:38 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:25:38 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-14 10:25:39 INFO Found: ['5154836']
2026-07-14 10:25:44 INFO [TGCC-IRENE] Submitted job with ID:['5154836']
2026-07-14 10:25:44 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:25:44 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14
2026-07-14 10:25:44 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020611_2_ENS14.nc
2026-07-14 10:25:44 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-14 10:25:44 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:25:44 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-14 10:25:44 INFO Queuing job for member 14...
2026-07-14 10:25:44 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:25:44 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-14 10:25:45 INFO Found: ['5154837']
2026-07-14 10:25:50 INFO [TGCC-IRENE] Submitted job with ID:['5154837']
2026-07-14 10:25:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:25:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15
2026-07-14 10:25:50 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020611_2_ENS15.nc
2026-07-14 10:25:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-14 10:25:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:25:50 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-14 10:25:50 INFO Queuing job for member 15...
2026-07-14 10:25:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:25:50 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-14 10:25:50 INFO Found: ['5154838']
2026-07-14 10:25:55 INFO [TGCC-IRENE] Submitted job with ID:['5154838']
2026-07-14 10:25:55 INFO Checking job status ...
2026-07-14 10:25:55 INFO None 5154819: status RUNNING/PENDING
2026-07-14 10:25:55 INFO None 5154820: status RUNNING/PENDING
2026-07-14 10:25:55 INFO None 5154821: status RUNNING/PENDING
2026-07-14 10:25:55 INFO None 5154822: status RUNNING/PENDING
2026-07-14 10:25:55 INFO None 5154824: status RUNNING/PENDING
2026-07-14 10:25:55 INFO None 5154827: status RUNNING/PENDING
2026-07-14 10:25:56 INFO None 5154830: status RUNNING/PENDING
2026-07-14 10:25:56 INFO None 5154831: status RUNNING/PENDING
2026-07-14 10:25:56 INFO None 5154832: status RUNNING/PENDING
2026-07-14 10:25:56 INFO None 5154833: status RUNNING/PENDING
2026-07-14 10:25:56 INFO None 5154834: status RUNNING/PENDING
2026-07-14 10:25:56 INFO None 5154835: status RUNNING/PENDING
2026-07-14 10:25:56 INFO None 5154836: status RUNNING/PENDING
2026-07-14 10:25:56 INFO None 5154837: status RUNNING/PENDING
2026-07-14 10:25:56 INFO None 5154838: status RUNNING/PENDING
2026-07-14 10:25:56 INFO Jobs still running: ['5154819', '5154820', '5154821', '5154822', '5154824', '5154827', '5154830', '5154831', '5154832', '5154833', '5154834', '5154835', '5154836', '5154837', '5154838']. Waiting...
2026-07-14 10:26:11 INFO None 5154819: status RUNNING/PENDING
2026-07-14 10:26:11 INFO None 5154820: status RUNNING/PENDING
2026-07-14 10:26:11 INFO None 5154821: status RUNNING/PENDING
2026-07-14 10:26:11 INFO None 5154822: status RUNNING/PENDING
2026-07-14 10:26:11 INFO None 5154824: status RUNNING/PENDING
2026-07-14 10:26:11 INFO None 5154827: status RUNNING/PENDING
2026-07-14 10:26:11 INFO None 5154830: status RUNNING/PENDING
2026-07-14 10:26:11 INFO None 5154831: status RUNNING/PENDING
2026-07-14 10:26:11 INFO None 5154832: status RUNNING/PENDING
2026-07-14 10:26:11 INFO None 5154833: status RUNNING/PENDING
2026-07-14 10:26:11 INFO None 5154834: status RUNNING/PENDING
2026-07-14 10:26:11 INFO None 5154835: status RUNNING/PENDING
2026-07-14 10:26:11 INFO None 5154836: status RUNNING/PENDING
2026-07-14 10:26:11 INFO None 5154837: status RUNNING/PENDING
2026-07-14 10:26:11 INFO None 5154838: status RUNNING/PENDING
2026-07-14 10:26:11 INFO Jobs still running: ['5154819', '5154820', '5154821', '5154822', '5154824', '5154827', '5154830', '5154831', '5154832', '5154833', '5154834', '5154835', '5154836', '5154837', '5154838']. Waiting...
2026-07-14 10:26:26 INFO None 5154819: status RUNNING/PENDING
2026-07-14 10:26:26 INFO None 5154820: status RUNNING/PENDING
2026-07-14 10:26:26 INFO None 5154821: status RUNNING/PENDING
2026-07-14 10:26:26 INFO None 5154822: status RUNNING/PENDING
2026-07-14 10:26:26 INFO None 5154824: status RUNNING/PENDING
2026-07-14 10:26:26 INFO None 5154827: status RUNNING/PENDING
2026-07-14 10:26:26 INFO None 5154830: status RUNNING/PENDING
2026-07-14 10:26:26 INFO None 5154831: status RUNNING/PENDING
2026-07-14 10:26:26 INFO None 5154832: status RUNNING/PENDING
2026-07-14 10:26:26 INFO None 5154833: status RUNNING/PENDING
2026-07-14 10:26:26 INFO None 5154834: status RUNNING/PENDING
2026-07-14 10:26:26 INFO None 5154835: status RUNNING/PENDING
2026-07-14 10:26:26 INFO None 5154836: status RUNNING/PENDING
2026-07-14 10:26:26 INFO None 5154837: status RUNNING/PENDING
2026-07-14 10:26:26 INFO None 5154838: status RUNNING/PENDING
2026-07-14 10:26:26 INFO Jobs still running: ['5154819', '5154820', '5154821', '5154822', '5154824', '5154827', '5154830', '5154831', '5154832', '5154833', '5154834', '5154835', '5154836', '5154837', '5154838']. Waiting...
2026-07-14 10:26:41 INFO None 5154819: status RUNNING/PENDING
2026-07-14 10:28:58 INFO None 5154820: status FINISHED
2026-07-14 10:28:58 INFO None 5154821: status FINISHED
2026-07-14 10:28:58 INFO None 5154822: status FINISHED
2026-07-14 10:28:58 INFO None 5154824: status FINISHED
2026-07-14 10:28:58 INFO None 5154827: status FINISHED
2026-07-14 10:28:58 INFO None 5154830: status FINISHED
2026-07-14 10:28:58 INFO None 5154831: status FINISHED
2026-07-14 10:28:58 INFO None 5154832: status RUNNING/PENDING
2026-07-14 10:28:58 INFO None 5154833: status RUNNING/PENDING
2026-07-14 10:28:58 INFO None 5154834: status RUNNING/PENDING
2026-07-14 10:28:58 INFO None 5154835: status RUNNING/PENDING
2026-07-14 10:28:58 INFO None 5154836: status RUNNING/PENDING
2026-07-14 10:28:58 INFO None 5154837: status RUNNING/PENDING
2026-07-14 10:28:58 INFO None 5154838: status RUNNING/PENDING
2026-07-14 10:28:58 INFO Jobs still running: ['5154819', '5154832', '5154833', '5154834', '5154835', '5154836', '5154837', '5154838']. Waiting...
2026-07-14 10:29:13 INFO None 5154819: status FINISHED
2026-07-14 10:29:13 INFO None 5154820: status FINISHED
2026-07-14 10:29:13 INFO None 5154821: status FINISHED
2026-07-14 10:29:13 INFO None 5154822: status FINISHED
2026-07-14 10:29:14 INFO None 5154824: status FINISHED
2026-07-14 10:29:14 INFO None 5154827: status FINISHED
2026-07-14 10:29:14 INFO None 5154830: status FINISHED
2026-07-14 10:29:14 INFO None 5154831: status FINISHED
2026-07-14 10:29:14 INFO None 5154832: status RUNNING/PENDING
2026-07-14 10:29:14 INFO None 5154833: status RUNNING/PENDING
2026-07-14 10:29:14 INFO None 5154834: status RUNNING/PENDING
2026-07-14 10:29:14 INFO None 5154835: status RUNNING/PENDING
2026-07-14 10:29:14 INFO None 5154836: status RUNNING/PENDING
2026-07-14 10:29:14 INFO None 5154837: status RUNNING/PENDING
2026-07-14 10:29:14 INFO None 5154838: status RUNNING/PENDING
2026-07-14 10:29:14 INFO Jobs still running: ['5154832', '5154833', '5154834', '5154835', '5154836', '5154837', '5154838']. Waiting...
2026-07-14 10:29:29 INFO None 5154819: status FINISHED
2026-07-14 10:29:29 INFO None 5154820: status FINISHED
2026-07-14 10:29:29 INFO None 5154821: status FINISHED
2026-07-14 10:29:29 INFO None 5154822: status FINISHED
2026-07-14 10:29:29 INFO None 5154824: status FINISHED
2026-07-14 10:29:29 INFO None 5154827: status FINISHED
2026-07-14 10:29:29 INFO None 5154830: status FINISHED
2026-07-14 10:29:29 INFO None 5154831: status FINISHED
2026-07-14 10:29:29 INFO None 5154832: status FINISHED
2026-07-14 10:29:29 INFO None 5154833: status FINISHED
2026-07-14 10:29:29 INFO None 5154834: status RUNNING/PENDING
2026-07-14 10:29:29 INFO None 5154835: status RUNNING/PENDING
2026-07-14 10:29:29 INFO None 5154836: status FINISHED
2026-07-14 10:29:29 INFO None 5154837: status RUNNING/PENDING
2026-07-14 10:29:29 INFO None 5154838: status RUNNING/PENDING
2026-07-14 10:29:29 INFO Jobs still running: ['5154834', '5154835', '5154837', '5154838']. Waiting...
2026-07-14 10:29:44 INFO None 5154819: status FINISHED
2026-07-14 10:29:44 INFO None 5154820: status FINISHED
2026-07-14 10:29:44 INFO None 5154821: status FINISHED
2026-07-14 10:29:44 INFO None 5154822: status FINISHED
2026-07-14 10:29:44 INFO None 5154824: status FINISHED
2026-07-14 10:29:44 INFO None 5154827: status FINISHED
2026-07-14 10:29:44 INFO None 5154830: status FINISHED
2026-07-14 10:29:44 INFO None 5154831: status FINISHED
2026-07-14 10:29:44 INFO None 5154832: status FINISHED
2026-07-14 10:29:44 INFO None 5154833: status FINISHED
2026-07-14 10:29:44 INFO None 5154834: status RUNNING/PENDING
2026-07-14 10:29:44 INFO None 5154835: status RUNNING/PENDING
2026-07-14 10:29:44 INFO None 5154836: status FINISHED
2026-07-14 10:29:44 INFO None 5154837: status RUNNING/PENDING
2026-07-14 10:29:44 INFO None 5154838: status RUNNING/PENDING
2026-07-14 10:29:44 INFO Jobs still running: ['5154834', '5154835', '5154837', '5154838']. Waiting...
2026-07-14 10:29:59 INFO None 5154819: status FINISHED
2026-07-14 10:29:59 INFO None 5154820: status FINISHED
2026-07-14 10:29:59 INFO None 5154821: status FINISHED
2026-07-14 10:29:59 INFO None 5154822: status FINISHED
2026-07-14 10:30:00 INFO None 5154824: status FINISHED
2026-07-14 10:30:00 INFO None 5154827: status FINISHED
2026-07-14 10:30:00 INFO None 5154830: status FINISHED
2026-07-14 10:30:00 INFO None 5154831: status FINISHED
2026-07-14 10:30:00 INFO None 5154832: status FINISHED
2026-07-14 10:30:00 INFO None 5154833: status FINISHED
2026-07-14 10:30:00 INFO None 5154834: status RUNNING/PENDING
2026-07-14 10:30:00 INFO None 5154835: status RUNNING/PENDING
2026-07-14 10:30:00 INFO None 5154836: status FINISHED
2026-07-14 10:30:00 INFO None 5154837: status RUNNING/PENDING
2026-07-14 10:30:00 INFO None 5154838: status RUNNING/PENDING
2026-07-14 10:30:00 INFO Jobs still running: ['5154834', '5154835', '5154837', '5154838']. Waiting...
2026-07-14 10:30:15 INFO None 5154819: status FINISHED
2026-07-14 10:30:15 INFO None 5154820: status FINISHED
2026-07-14 10:30:15 INFO None 5154821: status FINISHED
2026-07-14 10:30:15 INFO None 5154822: status FINISHED
2026-07-14 10:30:15 INFO None 5154824: status FINISHED
2026-07-14 10:30:15 INFO None 5154827: status FINISHED
2026-07-14 10:30:15 INFO None 5154830: status FINISHED
2026-07-14 10:30:15 INFO None 5154831: status FINISHED
2026-07-14 10:30:15 INFO None 5154832: status FINISHED
2026-07-14 10:30:15 INFO None 5154833: status FINISHED
2026-07-14 10:30:15 INFO None 5154834: status RUNNING/PENDING
2026-07-14 10:30:15 INFO None 5154835: status RUNNING/PENDING
2026-07-14 10:30:15 INFO None 5154836: status FINISHED
2026-07-14 10:30:15 INFO None 5154837: status RUNNING/PENDING
2026-07-14 10:30:15 INFO None 5154838: status RUNNING/PENDING
2026-07-14 10:30:15 INFO Jobs still running: ['5154834', '5154835', '5154837', '5154838']. Waiting...
2026-07-14 10:30:30 INFO None 5154819: status FINISHED
2026-07-14 10:30:30 INFO None 5154820: status FINISHED
2026-07-14 10:30:30 INFO None 5154821: status FINISHED
2026-07-14 10:30:30 INFO None 5154822: status FINISHED
2026-07-14 10:30:30 INFO None 5154824: status FINISHED
2026-07-14 10:30:30 INFO None 5154827: status FINISHED
2026-07-14 10:30:30 INFO None 5154830: status FINISHED
2026-07-14 10:30:30 INFO None 5154831: status FINISHED
2026-07-14 10:30:30 INFO None 5154832: status FINISHED
2026-07-14 10:30:30 INFO None 5154833: status FINISHED
2026-07-14 10:30:30 INFO None 5154834: status RUNNING/PENDING
2026-07-14 10:30:30 INFO None 5154835: status RUNNING/PENDING
2026-07-14 10:30:30 INFO None 5154836: status FINISHED
2026-07-14 10:30:30 INFO None 5154837: status RUNNING/PENDING
2026-07-14 10:30:30 INFO None 5154838: status RUNNING/PENDING
2026-07-14 10:30:30 INFO Jobs still running: ['5154834', '5154835', '5154837', '5154838']. Waiting...
2026-07-14 10:30:45 INFO None 5154819: status FINISHED
2026-07-14 10:30:45 INFO None 5154820: status FINISHED
2026-07-14 10:30:45 INFO None 5154821: status FINISHED
2026-07-14 10:30:45 INFO None 5154822: status FINISHED
2026-07-14 10:30:45 INFO None 5154824: status FINISHED
2026-07-14 10:30:45 INFO None 5154827: status FINISHED
2026-07-14 10:30:45 INFO None 5154830: status FINISHED
2026-07-14 10:30:46 INFO None 5154831: status FINISHED
2026-07-14 10:30:46 INFO None 5154832: status FINISHED
2026-07-14 10:30:46 INFO None 5154833: status FINISHED
2026-07-14 10:30:46 INFO None 5154834: status FINISHED
2026-07-14 10:30:46 INFO None 5154835: status FINISHED
2026-07-14 10:30:46 INFO None 5154836: status FINISHED
2026-07-14 10:30:46 INFO None 5154837: status RUNNING/PENDING
2026-07-14 10:30:46 INFO None 5154838: status RUNNING/PENDING
2026-07-14 10:30:46 INFO Jobs still running: ['5154837', '5154838']. Waiting...
2026-07-14 10:31:01 INFO None 5154819: status FINISHED
2026-07-14 10:31:01 INFO None 5154820: status FINISHED
2026-07-14 10:31:01 INFO None 5154821: status FINISHED
2026-07-14 10:31:01 INFO None 5154822: status FINISHED
2026-07-14 10:31:01 INFO None 5154824: status FINISHED
2026-07-14 10:31:01 INFO None 5154827: status FINISHED
2026-07-14 10:31:01 INFO None 5154830: status FINISHED
2026-07-14 10:31:01 INFO None 5154831: status FINISHED
2026-07-14 10:31:01 INFO None 5154832: status FINISHED
2026-07-14 10:31:01 INFO None 5154833: status FINISHED
2026-07-14 10:31:01 INFO None 5154834: status FINISHED
2026-07-14 10:31:01 INFO None 5154835: status FINISHED
2026-07-14 10:31:01 INFO None 5154836: status FINISHED
2026-07-14 10:31:01 INFO None 5154837: status FINISHED
2026-07-14 10:31:01 INFO None 5154838: status FINISHED
2026-07-14 10:31:01 INFO Jobs ['5154819', '5154820', '5154821', '5154822', '5154824', '5154827', '5154830', '5154831', '5154832', '5154833', '5154834', '5154835', '5154836', '5154837', '5154838'] have finished
2026-07-14 10:31:01 INFO Checking restart files were created ...
2026-07-14 10:33:39 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020613_1_ENS1.nc(668832435 bytes)
2026-07-14 10:33:39 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020613_1_ENS2.nc(668832435 bytes)
2026-07-14 10:33:39 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020613_1_ENS3.nc(668832435 bytes)
2026-07-14 10:33:39 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020613_1_ENS4.nc(668832435 bytes)
2026-07-14 10:33:39 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020613_1_ENS5.nc(668832435 bytes)
2026-07-14 10:33:39 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020613_1_ENS6.nc(668832435 bytes)
2026-07-14 10:33:39 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020613_1_ENS7.nc(668832435 bytes)
2026-07-14 10:33:39 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020613_1_ENS8.nc(668832435 bytes)
2026-07-14 10:33:39 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020613_1_ENS9.nc(668832435 bytes)
2026-07-14 10:33:39 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020613_1_ENS10.nc(668832435 bytes)
2026-07-14 10:33:39 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020613_1_ENS11.nc(668832435 bytes)
2026-07-14 10:33:39 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020613_1_ENS12.nc(668832435 bytes)
2026-07-14 10:33:39 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020613_1_ENS13.nc(668832435 bytes)
2026-07-14 10:33:39 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020613_1_ENS14.nc(668832435 bytes)
2026-07-14 10:33:39 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020613_1_ENS15.nc(668832435 bytes)
2026-07-14 10:33:39 INFO  Run_model() completed successfully.
2026-07-14 10:33:39 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 10:33:39 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 14:00:00 days=153072 seconds=50400
2026-07-14 10:33:39 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-14 10:33:39 INFO [TIME] increment current_time 2020-02-06 13:00:00 -> 2020-02-06 14:00:00
2026-07-14 10:33:39 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 10:33:39 INFO ---------->>> Running process_satellite_data()
2026-07-14 10:33:39 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12005.nc
2026-07-14 10:33:39 INFO ---------->>> Running run_obs_converter()
2026-07-14 10:33:39 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_51225_153072.out
2026-07-14 10:33:39 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_51225_153072.out
2026-07-14 10:33:39 INFO ---------->>> Running DART
2026-07-14 10:33:39 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-14 10:33:39 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/chim_ENS1_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/chim_ENS1_2020020614_1_out_toDART.nc
2026-07-14 10:33:40 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/chim_ENS2_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/chim_ENS2_2020020614_1_out_toDART.nc
2026-07-14 10:33:40 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/chim_ENS3_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/chim_ENS3_2020020614_1_out_toDART.nc
2026-07-14 10:33:40 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/chim_ENS4_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/chim_ENS4_2020020614_1_out_toDART.nc
2026-07-14 10:33:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/chim_ENS5_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/chim_ENS5_2020020614_1_out_toDART.nc
2026-07-14 10:33:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/chim_ENS6_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/chim_ENS6_2020020614_1_out_toDART.nc
2026-07-14 10:33:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/chim_ENS7_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/chim_ENS7_2020020614_1_out_toDART.nc
2026-07-14 10:33:42 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/chim_ENS8_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/chim_ENS8_2020020614_1_out_toDART.nc
2026-07-14 10:33:42 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/chim_ENS9_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/chim_ENS9_2020020614_1_out_toDART.nc
2026-07-14 10:33:42 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/chim_ENS10_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/chim_ENS10_2020020614_1_out_toDART.nc
2026-07-14 10:33:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/chim_ENS11_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/chim_ENS11_2020020614_1_out_toDART.nc
2026-07-14 10:33:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/chim_ENS12_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/chim_ENS12_2020020614_1_out_toDART.nc
2026-07-14 10:33:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/chim_ENS13_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/chim_ENS13_2020020614_1_out_toDART.nc
2026-07-14 10:33:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/chim_ENS14_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/chim_ENS14_2020020614_1_out_toDART.nc
2026-07-14 10:33:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/chim_ENS15_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/chim_ENS15_2020020614_1_out_toDART.nc
2026-07-14 10:33:44 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-14 10:33:44 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-14 10:33:44 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-14 10:33:44 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-14 10:33:44 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-14 10:33:44 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-14 10:33:52 INFO Found: []
2026-07-14 10:33:52 INFO No job id returned by command ./run_filter.bsh
2026-07-14 10:33:52 INFO No monitoring will be performed
2026-07-14 10:33:52 INFO Moving DART output files to analysis and preassim directories for date 2020020614 if present ...
2026-07-14 10:33:52 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-14 10:33:52 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-14 10:33:52 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-14 10:33:52 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-14 10:33:52 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-14 10:33:52 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-14 10:33:52 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-14 10:33:52 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-14 10:33:52 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-14 10:33:52 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-14 10:33:52 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-14 10:33:52 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-14 10:33:52 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-14 10:33:52 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-14 10:33:52 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-14 10:33:52 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-14 10:33:52 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-14 10:33:52 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-14 10:33:52 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-14 10:33:52 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-14 10:33:52 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-14 10:33:52 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-14 10:33:52 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-14 10:33:52 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-14 10:33:52 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-14 10:33:52 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-14 10:33:52 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-14 10:33:52 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-14 10:33:52 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-14 10:33:52 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-14 10:33:52 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-14 10:33:52 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-14 10:33:52 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-14 10:33:52 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-14 10:33:52 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-14 10:33:52 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-14 10:33:52 INFO run_dart() is DONE.
2026-07-14 10:33:52 INFO ---------->>> Running update_pollutant_in_end()
2026-07-14 10:33:53 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:33:53 INFO No previous orbit memory found.
2026-07-14 10:33:53 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:33:53 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS1_2020020614.nc
2026-07-14 10:33:53 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:33:54 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:33:54 INFO No previous orbit memory found.
2026-07-14 10:33:54 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:33:54 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS1_2020020614.nc
2026-07-14 10:33:54 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:33:54 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:33:54 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:33:54 INFO No previous orbit memory found.
2026-07-14 10:33:54 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:33:54 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS2_2020020614.nc
2026-07-14 10:33:55 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:33:55 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:33:55 INFO No previous orbit memory found.
2026-07-14 10:33:55 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:33:55 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS2_2020020614.nc
2026-07-14 10:33:56 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:33:56 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:33:56 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:33:56 INFO No previous orbit memory found.
2026-07-14 10:33:56 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:33:56 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS3_2020020614.nc
2026-07-14 10:33:57 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:33:57 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:33:57 INFO No previous orbit memory found.
2026-07-14 10:33:57 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:33:57 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS3_2020020614.nc
2026-07-14 10:33:57 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:33:57 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:33:58 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:33:58 INFO No previous orbit memory found.
2026-07-14 10:33:58 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:33:58 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS4_2020020614.nc
2026-07-14 10:33:58 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:33:59 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:33:59 INFO No previous orbit memory found.
2026-07-14 10:33:59 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:33:59 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS4_2020020614.nc
2026-07-14 10:33:59 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:33:59 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:33:59 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:33:59 INFO No previous orbit memory found.
2026-07-14 10:33:59 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:33:59 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS5_2020020614.nc
2026-07-14 10:34:00 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:34:00 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:34:00 INFO No previous orbit memory found.
2026-07-14 10:34:00 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:34:00 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS5_2020020614.nc
2026-07-14 10:34:01 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:34:01 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:34:01 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:34:01 INFO No previous orbit memory found.
2026-07-14 10:34:01 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:34:01 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS6_2020020614.nc
2026-07-14 10:34:02 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:34:02 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:34:02 INFO No previous orbit memory found.
2026-07-14 10:34:02 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:34:02 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS6_2020020614.nc
2026-07-14 10:34:02 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:34:02 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:34:03 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:34:03 INFO No previous orbit memory found.
2026-07-14 10:34:03 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:34:03 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS7_2020020614.nc
2026-07-14 10:34:03 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:34:04 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:34:04 INFO No previous orbit memory found.
2026-07-14 10:34:04 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:34:04 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS7_2020020614.nc
2026-07-14 10:34:04 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:34:04 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:34:05 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:34:05 INFO No previous orbit memory found.
2026-07-14 10:34:05 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:34:05 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS8_2020020614.nc
2026-07-14 10:34:05 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:34:06 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:34:06 INFO No previous orbit memory found.
2026-07-14 10:34:06 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:34:06 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS8_2020020614.nc
2026-07-14 10:34:06 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:34:06 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:34:06 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:34:06 INFO No previous orbit memory found.
2026-07-14 10:34:06 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:34:06 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS9_2020020614.nc
2026-07-14 10:34:07 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:34:07 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:34:07 INFO No previous orbit memory found.
2026-07-14 10:34:07 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:34:07 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS9_2020020614.nc
2026-07-14 10:34:08 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:34:08 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:34:08 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:34:08 INFO No previous orbit memory found.
2026-07-14 10:34:08 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:34:08 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS10_2020020614.nc
2026-07-14 10:34:09 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:34:09 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:34:09 INFO No previous orbit memory found.
2026-07-14 10:34:09 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:34:09 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS10_2020020614.nc
2026-07-14 10:34:10 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:34:10 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:34:10 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:34:10 INFO No previous orbit memory found.
2026-07-14 10:34:10 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:34:10 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS11_2020020614.nc
2026-07-14 10:34:11 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:34:11 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:34:11 INFO No previous orbit memory found.
2026-07-14 10:34:11 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:34:11 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS11_2020020614.nc
2026-07-14 10:34:11 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:34:11 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:34:12 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:34:12 INFO No previous orbit memory found.
2026-07-14 10:34:12 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:34:12 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS12_2020020614.nc
2026-07-14 10:34:12 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:34:12 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:34:12 INFO No previous orbit memory found.
2026-07-14 10:34:12 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:34:12 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS12_2020020614.nc
2026-07-14 10:34:13 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:34:13 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:34:13 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:34:13 INFO No previous orbit memory found.
2026-07-14 10:34:13 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:34:13 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS13_2020020614.nc
2026-07-14 10:34:14 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:34:14 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:34:14 INFO No previous orbit memory found.
2026-07-14 10:34:14 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:34:14 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS13_2020020614.nc
2026-07-14 10:34:15 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:34:15 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:34:15 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:34:15 INFO No previous orbit memory found.
2026-07-14 10:34:15 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:34:15 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS14_2020020614.nc
2026-07-14 10:34:16 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:34:16 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:34:16 INFO No previous orbit memory found.
2026-07-14 10:34:16 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:34:16 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS14_2020020614.nc
2026-07-14 10:34:17 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:34:17 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:34:17 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:34:17 INFO No previous orbit memory found.
2026-07-14 10:34:17 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:34:17 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS15_2020020614.nc
2026-07-14 10:34:18 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:34:18 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 10:34:18 INFO No previous orbit memory found.
2026-07-14 10:34:18 INFO Emission correction applied with pixel-based damping.
2026-07-14 10:34:18 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS15_2020020614.nc
2026-07-14 10:34:18 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 10:34:18 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 10:34:18 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-14 10:34:18 INFO [TIME] step_end current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 10:34:18 INFO [TIME] step_start current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 10:34:18 INFO [TIME] window start=2020-02-06 14:00:00 end=2020-02-07 00:00:00 run_hours=10 has_assimilation=False
2026-07-14 10:34:18 INFO Copying EMIS of next day ...
2026-07-14 10:34:19 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens01.nc
2026-07-14 10:34:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:34:20 INFO Hourly dataset computed and listing created
2026-07-14 10:39:53 INFO Hourly dataset computed
2026-07-14 10:39:53 INFO Copying EMIS of next day ...
2026-07-14 10:39:54 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens02.nc
2026-07-14 10:39:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:39:55 INFO Hourly dataset computed and listing created
2026-07-14 10:40:21 INFO Hourly dataset computed
2026-07-14 10:40:21 INFO Copying EMIS of next day ...
2026-07-14 10:40:22 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens03.nc
2026-07-14 10:40:22 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:40:24 INFO Hourly dataset computed and listing created
2026-07-14 10:40:45 INFO Hourly dataset computed
2026-07-14 10:40:45 INFO Copying EMIS of next day ...
2026-07-14 10:40:45 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens04.nc
2026-07-14 10:40:45 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:40:47 INFO Hourly dataset computed and listing created
2026-07-14 10:43:12 INFO Hourly dataset computed
2026-07-14 10:43:12 INFO Copying EMIS of next day ...
2026-07-14 10:43:12 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens05.nc
2026-07-14 10:43:12 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:43:14 INFO Hourly dataset computed and listing created
2026-07-14 10:43:31 INFO Hourly dataset computed
2026-07-14 10:43:31 INFO Copying EMIS of next day ...
2026-07-14 10:43:31 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens06.nc
2026-07-14 10:43:31 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:43:33 INFO Hourly dataset computed and listing created
2026-07-14 10:43:48 INFO Hourly dataset computed
2026-07-14 10:43:48 INFO Copying EMIS of next day ...
2026-07-14 10:43:48 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens07.nc
2026-07-14 10:43:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:43:50 INFO Hourly dataset computed and listing created
2026-07-14 10:44:05 INFO Hourly dataset computed
2026-07-14 10:44:05 INFO Copying EMIS of next day ...
2026-07-14 10:44:06 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens08.nc
2026-07-14 10:44:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:44:08 INFO Hourly dataset computed and listing created
2026-07-14 10:44:23 INFO Hourly dataset computed
2026-07-14 10:44:23 INFO Copying EMIS of next day ...
2026-07-14 10:44:24 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens09.nc
2026-07-14 10:44:24 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:44:25 INFO Hourly dataset computed and listing created
2026-07-14 10:44:41 INFO Hourly dataset computed
2026-07-14 10:44:41 INFO Copying EMIS of next day ...
2026-07-14 10:44:42 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens10.nc
2026-07-14 10:44:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:44:43 INFO Hourly dataset computed and listing created
2026-07-14 10:44:57 INFO Hourly dataset computed
2026-07-14 10:44:57 INFO Copying EMIS of next day ...
2026-07-14 10:44:57 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens11.nc
2026-07-14 10:44:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:44:59 INFO Hourly dataset computed and listing created
2026-07-14 10:45:16 INFO Hourly dataset computed
2026-07-14 10:45:16 INFO Copying EMIS of next day ...
2026-07-14 10:45:16 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens12.nc
2026-07-14 10:45:16 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:45:19 INFO Hourly dataset computed and listing created
2026-07-14 10:45:43 INFO Hourly dataset computed
2026-07-14 10:45:43 INFO Copying EMIS of next day ...
2026-07-14 10:45:43 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens13.nc
2026-07-14 10:45:43 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:45:45 INFO Hourly dataset computed and listing created
2026-07-14 10:48:09 INFO Hourly dataset computed
2026-07-14 10:48:09 INFO Copying EMIS of next day ...
2026-07-14 10:48:10 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens14.nc
2026-07-14 10:48:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:48:11 INFO Hourly dataset computed and listing created
2026-07-14 10:48:29 INFO Hourly dataset computed
2026-07-14 10:48:29 INFO Copying EMIS of next day ...
2026-07-14 10:48:30 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens15.nc
2026-07-14 10:48:30 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 10:48:31 INFO Hourly dataset computed and listing created
2026-07-14 10:48:48 INFO Hourly dataset computed
2026-07-14 10:48:48 INFO ---------->>> Running CHIMERE model from 2020-02-06 14:00:00 to 2020-02-07 00:00:00
2026-07-14 10:48:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:48:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1
2026-07-14 10:48:48 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020613_1_ENS1.nc
2026-07-14 10:48:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-14 10:48:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:48:48 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-14 10:48:48 INFO Queuing job for member 1...
2026-07-14 10:48:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:48:48 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-14 10:48:48 INFO Found: ['5154903']
2026-07-14 10:48:53 INFO [TGCC-IRENE] Submitted job with ID:['5154903']
2026-07-14 10:48:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:48:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2
2026-07-14 10:48:53 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020613_1_ENS2.nc
2026-07-14 10:48:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-14 10:48:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:48:53 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-14 10:48:53 INFO Queuing job for member 2...
2026-07-14 10:48:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:48:53 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-14 10:48:54 INFO Found: ['5154905']
2026-07-14 10:48:59 INFO [TGCC-IRENE] Submitted job with ID:['5154905']
2026-07-14 10:48:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:48:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3
2026-07-14 10:48:59 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020613_1_ENS3.nc
2026-07-14 10:48:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-14 10:48:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:48:59 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-14 10:48:59 INFO Queuing job for member 3...
2026-07-14 10:48:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:48:59 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-14 10:49:00 INFO Found: ['5154906']
2026-07-14 10:49:05 INFO [TGCC-IRENE] Submitted job with ID:['5154906']
2026-07-14 10:49:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:49:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4
2026-07-14 10:49:05 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020613_1_ENS4.nc
2026-07-14 10:49:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-14 10:49:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:49:05 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-14 10:49:05 INFO Queuing job for member 4...
2026-07-14 10:49:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:49:05 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-14 10:49:06 INFO Found: ['5154923']
2026-07-14 10:49:11 INFO [TGCC-IRENE] Submitted job with ID:['5154923']
2026-07-14 10:49:11 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:49:11 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5
2026-07-14 10:49:11 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020613_1_ENS5.nc
2026-07-14 10:49:11 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-14 10:49:11 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:49:11 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-14 10:49:11 INFO Queuing job for member 5...
2026-07-14 10:49:11 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:49:11 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-14 10:49:11 INFO Found: ['5154928']
2026-07-14 10:49:16 INFO [TGCC-IRENE] Submitted job with ID:['5154928']
2026-07-14 10:49:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:49:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6
2026-07-14 10:49:16 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020613_1_ENS6.nc
2026-07-14 10:49:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-14 10:49:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:49:16 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-14 10:49:16 INFO Queuing job for member 6...
2026-07-14 10:49:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:49:16 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-14 10:49:17 INFO Found: ['5154929']
2026-07-14 10:49:22 INFO [TGCC-IRENE] Submitted job with ID:['5154929']
2026-07-14 10:49:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:49:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7
2026-07-14 10:49:22 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020613_1_ENS7.nc
2026-07-14 10:49:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-14 10:49:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:49:22 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-14 10:49:22 INFO Queuing job for member 7...
2026-07-14 10:49:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:49:22 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-14 10:49:23 INFO Found: ['5154930']
2026-07-14 10:49:28 INFO [TGCC-IRENE] Submitted job with ID:['5154930']
2026-07-14 10:49:28 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:49:28 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8
2026-07-14 10:49:28 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020613_1_ENS8.nc
2026-07-14 10:49:28 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-14 10:49:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:49:28 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-14 10:49:28 INFO Queuing job for member 8...
2026-07-14 10:49:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:49:28 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-14 10:49:29 INFO Found: ['5154931']
2026-07-14 10:49:34 INFO [TGCC-IRENE] Submitted job with ID:['5154931']
2026-07-14 10:49:34 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:49:34 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9
2026-07-14 10:49:34 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020613_1_ENS9.nc
2026-07-14 10:49:34 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-14 10:49:34 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:49:34 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-14 10:49:34 INFO Queuing job for member 9...
2026-07-14 10:49:34 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:49:34 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-14 10:49:34 INFO Found: ['5154932']
2026-07-14 10:49:39 INFO [TGCC-IRENE] Submitted job with ID:['5154932']
2026-07-14 10:49:39 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:49:39 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10
2026-07-14 10:49:39 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020613_1_ENS10.nc
2026-07-14 10:49:39 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-14 10:49:39 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:49:39 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-14 10:49:39 INFO Queuing job for member 10...
2026-07-14 10:49:39 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:49:39 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-14 10:49:40 INFO Found: ['5154933']
2026-07-14 10:49:45 INFO [TGCC-IRENE] Submitted job with ID:['5154933']
2026-07-14 10:49:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:49:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11
2026-07-14 10:49:45 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020613_1_ENS11.nc
2026-07-14 10:49:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-14 10:49:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:49:45 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-14 10:49:45 INFO Queuing job for member 11...
2026-07-14 10:49:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:49:45 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-14 10:49:46 INFO Found: ['5154934']
2026-07-14 10:49:51 INFO [TGCC-IRENE] Submitted job with ID:['5154934']
2026-07-14 10:49:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:49:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12
2026-07-14 10:49:51 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020613_1_ENS12.nc
2026-07-14 10:49:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-14 10:49:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:49:51 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-14 10:49:51 INFO Queuing job for member 12...
2026-07-14 10:49:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:49:51 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-14 10:49:52 INFO Found: ['5154949']
2026-07-14 10:49:57 INFO [TGCC-IRENE] Submitted job with ID:['5154949']
2026-07-14 10:49:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:49:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13
2026-07-14 10:49:57 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020613_1_ENS13.nc
2026-07-14 10:49:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-14 10:49:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:49:57 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-14 10:49:57 INFO Queuing job for member 13...
2026-07-14 10:49:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:49:57 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-14 10:49:57 INFO Found: ['5154970']
2026-07-14 10:50:02 INFO [TGCC-IRENE] Submitted job with ID:['5154970']
2026-07-14 10:50:02 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:50:02 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14
2026-07-14 10:50:02 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020613_1_ENS14.nc
2026-07-14 10:50:02 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-14 10:50:02 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:50:02 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-14 10:50:02 INFO Queuing job for member 14...
2026-07-14 10:50:02 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:50:02 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-14 10:50:03 INFO Found: ['5154989']
2026-07-14 10:50:08 INFO [TGCC-IRENE] Submitted job with ID:['5154989']
2026-07-14 10:50:08 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 10:50:08 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15
2026-07-14 10:50:08 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020613_1_ENS15.nc
2026-07-14 10:50:08 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-14 10:50:08 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 10:50:08 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-14 10:50:08 INFO Queuing job for member 15...
2026-07-14 10:50:08 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 10:50:08 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-14 10:50:09 INFO Found: ['5154994']
2026-07-14 10:50:14 INFO [TGCC-IRENE] Submitted job with ID:['5154994']
2026-07-14 10:50:14 INFO Checking job status ...
2026-07-14 10:50:14 INFO None 5154903: status RUNNING/PENDING
2026-07-14 10:50:14 INFO None 5154905: status RUNNING/PENDING
2026-07-14 10:50:14 INFO None 5154906: status RUNNING/PENDING
2026-07-14 10:50:14 INFO None 5154923: status RUNNING/PENDING
2026-07-14 10:50:14 INFO None 5154928: status RUNNING/PENDING
2026-07-14 10:50:14 INFO None 5154929: status RUNNING/PENDING
2026-07-14 10:50:14 INFO None 5154930: status RUNNING/PENDING
2026-07-14 10:50:14 INFO None 5154931: status RUNNING/PENDING
2026-07-14 10:50:14 INFO None 5154932: status RUNNING/PENDING
2026-07-14 10:50:14 INFO None 5154933: status RUNNING/PENDING
2026-07-14 10:50:14 INFO None 5154934: status RUNNING/PENDING
2026-07-14 10:50:14 INFO None 5154949: status RUNNING/PENDING
2026-07-14 10:50:14 INFO None 5154970: status RUNNING/PENDING
2026-07-14 10:50:14 INFO None 5154989: status RUNNING/PENDING
2026-07-14 10:50:14 INFO None 5154994: status RUNNING/PENDING
2026-07-14 10:50:14 INFO Jobs still running: ['5154903', '5154905', '5154906', '5154923', '5154928', '5154929', '5154930', '5154931', '5154932', '5154933', '5154934', '5154949', '5154970', '5154989', '5154994']. Waiting...
2026-07-14 10:50:29 INFO None 5154903: status RUNNING/PENDING
2026-07-14 10:50:29 INFO None 5154905: status RUNNING/PENDING
2026-07-14 10:50:29 INFO None 5154906: status RUNNING/PENDING
2026-07-14 10:50:29 INFO None 5154923: status RUNNING/PENDING
2026-07-14 10:50:29 INFO None 5154928: status RUNNING/PENDING
2026-07-14 10:50:29 INFO None 5154929: status RUNNING/PENDING
2026-07-14 10:50:29 INFO None 5154930: status RUNNING/PENDING
2026-07-14 10:50:29 INFO None 5154931: status RUNNING/PENDING
2026-07-14 10:50:29 INFO None 5154932: status RUNNING/PENDING
2026-07-14 10:50:29 INFO None 5154933: status RUNNING/PENDING
2026-07-14 10:50:29 INFO None 5154934: status RUNNING/PENDING
2026-07-14 10:50:29 INFO None 5154949: status RUNNING/PENDING
2026-07-14 10:50:29 INFO None 5154970: status RUNNING/PENDING
2026-07-14 10:50:29 INFO None 5154989: status RUNNING/PENDING
2026-07-14 10:50:29 INFO None 5154994: status RUNNING/PENDING
2026-07-14 10:50:29 INFO Jobs still running: ['5154903', '5154905', '5154906', '5154923', '5154928', '5154929', '5154930', '5154931', '5154932', '5154933', '5154934', '5154949', '5154970', '5154989', '5154994']. Waiting...
2026-07-14 10:50:46 INFO None 5154903: status RUNNING/PENDING
2026-07-14 10:50:46 INFO None 5154905: status RUNNING/PENDING
2026-07-14 10:50:46 INFO None 5154906: status RUNNING/PENDING
2026-07-14 10:50:46 INFO None 5154923: status RUNNING/PENDING
2026-07-14 10:50:46 INFO None 5154928: status RUNNING/PENDING
2026-07-14 10:50:46 INFO None 5154929: status RUNNING/PENDING
2026-07-14 10:50:46 INFO None 5154930: status RUNNING/PENDING
2026-07-14 10:50:46 INFO None 5154931: status RUNNING/PENDING
2026-07-14 10:50:46 INFO None 5154932: status RUNNING/PENDING
2026-07-14 10:50:46 INFO None 5154933: status RUNNING/PENDING
2026-07-14 10:50:46 INFO None 5154934: status RUNNING/PENDING
2026-07-14 10:50:46 INFO None 5154949: status RUNNING/PENDING
2026-07-14 10:50:46 INFO None 5154970: status RUNNING/PENDING
2026-07-14 10:50:46 INFO None 5154989: status RUNNING/PENDING
2026-07-14 10:50:46 INFO None 5154994: status RUNNING/PENDING
2026-07-14 10:50:46 INFO Jobs still running: ['5154903', '5154905', '5154906', '5154923', '5154928', '5154929', '5154930', '5154931', '5154932', '5154933', '5154934', '5154949', '5154970', '5154989', '5154994']. Waiting...
2026-07-14 10:51:01 INFO None 5154903: status RUNNING/PENDING
2026-07-14 10:51:01 INFO None 5154905: status RUNNING/PENDING
2026-07-14 10:51:01 INFO None 5154906: status RUNNING/PENDING
2026-07-14 10:51:01 INFO None 5154923: status RUNNING/PENDING
2026-07-14 10:51:01 INFO None 5154928: status RUNNING/PENDING
2026-07-14 10:51:01 INFO None 5154929: status RUNNING/PENDING
2026-07-14 10:51:01 INFO None 5154930: status RUNNING/PENDING
2026-07-14 10:51:01 INFO None 5154931: status RUNNING/PENDING
2026-07-14 10:51:01 INFO None 5154932: status RUNNING/PENDING
2026-07-14 10:51:01 INFO None 5154933: status RUNNING/PENDING
2026-07-14 10:51:01 INFO None 5154934: status RUNNING/PENDING
2026-07-14 10:51:01 INFO None 5154949: status RUNNING/PENDING
2026-07-14 10:51:01 INFO None 5154970: status RUNNING/PENDING
2026-07-14 10:51:01 INFO None 5154989: status RUNNING/PENDING
2026-07-14 10:51:01 INFO None 5154994: status RUNNING/PENDING
2026-07-14 10:51:01 INFO Jobs still running: ['5154903', '5154905', '5154906', '5154923', '5154928', '5154929', '5154930', '5154931', '5154932', '5154933', '5154934', '5154949', '5154970', '5154989', '5154994']. Waiting...
2026-07-14 10:51:16 INFO None 5154903: status RUNNING/PENDING
2026-07-14 10:51:16 INFO None 5154905: status RUNNING/PENDING
2026-07-14 10:51:16 INFO None 5154906: status RUNNING/PENDING
2026-07-14 10:51:16 INFO None 5154923: status RUNNING/PENDING
2026-07-14 10:51:16 INFO None 5154928: status RUNNING/PENDING
2026-07-14 10:51:16 INFO None 5154929: status RUNNING/PENDING
2026-07-14 10:51:16 INFO None 5154930: status RUNNING/PENDING
2026-07-14 10:51:16 INFO None 5154931: status FINISHED
2026-07-14 10:51:16 INFO None 5154932: status RUNNING/PENDING
2026-07-14 10:53:18 INFO None 5154933: status RUNNING/PENDING
2026-07-14 10:53:18 INFO None 5154934: status RUNNING/PENDING
2026-07-14 10:53:18 INFO None 5154949: status RUNNING/PENDING
2026-07-14 10:53:18 INFO None 5154970: status RUNNING/PENDING
2026-07-14 10:53:18 INFO None 5154989: status RUNNING/PENDING
2026-07-14 10:53:18 INFO None 5154994: status RUNNING/PENDING
2026-07-14 10:53:18 INFO Jobs still running: ['5154903', '5154905', '5154906', '5154923', '5154928', '5154929', '5154930', '5154932', '5154933', '5154934', '5154949', '5154970', '5154989', '5154994']. Waiting...
2026-07-14 10:53:34 INFO None 5154903: status RUNNING/PENDING
2026-07-14 10:53:34 INFO None 5154905: status RUNNING/PENDING
2026-07-14 10:53:34 INFO None 5154906: status RUNNING/PENDING
2026-07-14 10:53:34 INFO None 5154923: status RUNNING/PENDING
2026-07-14 10:53:34 INFO None 5154928: status RUNNING/PENDING
2026-07-14 10:53:34 INFO None 5154929: status RUNNING/PENDING
2026-07-14 10:53:34 INFO None 5154930: status RUNNING/PENDING
2026-07-14 10:53:34 INFO None 5154931: status FINISHED
2026-07-14 10:53:34 INFO None 5154932: status RUNNING/PENDING
2026-07-14 10:53:34 INFO None 5154933: status RUNNING/PENDING
2026-07-14 10:53:34 INFO None 5154934: status RUNNING/PENDING
2026-07-14 10:53:34 INFO None 5154949: status RUNNING/PENDING
2026-07-14 10:53:34 INFO None 5154970: status RUNNING/PENDING
2026-07-14 10:53:34 INFO None 5154989: status RUNNING/PENDING
2026-07-14 10:53:34 INFO None 5154994: status RUNNING/PENDING
2026-07-14 10:53:34 INFO Jobs still running: ['5154903', '5154905', '5154906', '5154923', '5154928', '5154929', '5154930', '5154932', '5154933', '5154934', '5154949', '5154970', '5154989', '5154994']. Waiting...
2026-07-14 10:53:49 INFO None 5154903: status RUNNING/PENDING
2026-07-14 10:53:49 INFO None 5154905: status RUNNING/PENDING
2026-07-14 10:53:49 INFO None 5154906: status RUNNING/PENDING
2026-07-14 10:53:49 INFO None 5154923: status RUNNING/PENDING
2026-07-14 10:53:49 INFO None 5154928: status RUNNING/PENDING
2026-07-14 10:53:49 INFO None 5154929: status RUNNING/PENDING
2026-07-14 10:53:49 INFO None 5154930: status RUNNING/PENDING
2026-07-14 10:53:49 INFO None 5154931: status FINISHED
2026-07-14 10:53:49 INFO None 5154932: status RUNNING/PENDING
2026-07-14 10:53:49 INFO None 5154933: status RUNNING/PENDING
2026-07-14 10:53:49 INFO None 5154934: status RUNNING/PENDING
2026-07-14 10:53:49 INFO None 5154949: status RUNNING/PENDING
2026-07-14 10:53:49 INFO None 5154970: status RUNNING/PENDING
2026-07-14 10:53:49 INFO None 5154989: status RUNNING/PENDING
2026-07-14 10:53:49 INFO None 5154994: status RUNNING/PENDING
2026-07-14 10:53:49 INFO Jobs still running: ['5154903', '5154905', '5154906', '5154923', '5154928', '5154929', '5154930', '5154932', '5154933', '5154934', '5154949', '5154970', '5154989', '5154994']. Waiting...
2026-07-14 10:54:04 INFO None 5154903: status RUNNING/PENDING
2026-07-14 10:54:04 INFO None 5154905: status RUNNING/PENDING
2026-07-14 10:54:04 INFO None 5154906: status RUNNING/PENDING
2026-07-14 10:54:04 INFO None 5154923: status RUNNING/PENDING
2026-07-14 10:54:04 INFO None 5154928: status RUNNING/PENDING
2026-07-14 10:54:04 INFO None 5154929: status RUNNING/PENDING
2026-07-14 10:54:04 INFO None 5154930: status RUNNING/PENDING
2026-07-14 10:54:04 INFO None 5154931: status FINISHED
2026-07-14 10:54:04 INFO None 5154932: status RUNNING/PENDING
2026-07-14 10:54:04 INFO None 5154933: status RUNNING/PENDING
2026-07-14 10:54:04 INFO None 5154934: status RUNNING/PENDING
2026-07-14 10:54:04 INFO None 5154949: status RUNNING/PENDING
2026-07-14 10:54:04 INFO None 5154970: status RUNNING/PENDING
2026-07-14 10:54:04 INFO None 5154989: status RUNNING/PENDING
2026-07-14 10:54:04 INFO None 5154994: status RUNNING/PENDING
2026-07-14 10:54:04 INFO Jobs still running: ['5154903', '5154905', '5154906', '5154923', '5154928', '5154929', '5154930', '5154932', '5154933', '5154934', '5154949', '5154970', '5154989', '5154994']. Waiting...
2026-07-14 10:54:19 INFO None 5154903: status RUNNING/PENDING
2026-07-14 10:54:19 INFO None 5154905: status RUNNING/PENDING
2026-07-14 10:54:19 INFO None 5154906: status RUNNING/PENDING
2026-07-14 10:54:19 INFO None 5154923: status RUNNING/PENDING
2026-07-14 10:54:19 INFO None 5154928: status RUNNING/PENDING
2026-07-14 10:54:19 INFO None 5154929: status RUNNING/PENDING
2026-07-14 10:54:19 INFO None 5154930: status RUNNING/PENDING
2026-07-14 10:54:20 INFO None 5154931: status FINISHED
2026-07-14 10:54:20 INFO None 5154932: status RUNNING/PENDING
2026-07-14 10:54:20 INFO None 5154933: status RUNNING/PENDING
2026-07-14 10:54:20 INFO None 5154934: status RUNNING/PENDING
2026-07-14 10:54:20 INFO None 5154949: status RUNNING/PENDING
2026-07-14 10:54:20 INFO None 5154970: status RUNNING/PENDING
2026-07-14 10:54:20 INFO None 5154989: status RUNNING/PENDING
2026-07-14 10:54:20 INFO None 5154994: status RUNNING/PENDING
2026-07-14 10:54:20 INFO Jobs still running: ['5154903', '5154905', '5154906', '5154923', '5154928', '5154929', '5154930', '5154932', '5154933', '5154934', '5154949', '5154970', '5154989', '5154994']. Waiting...
2026-07-14 10:54:35 INFO None 5154903: status RUNNING/PENDING
2026-07-14 10:54:35 INFO None 5154905: status RUNNING/PENDING
2026-07-14 10:54:35 INFO None 5154906: status RUNNING/PENDING
2026-07-14 10:54:35 INFO None 5154923: status RUNNING/PENDING
2026-07-14 10:54:35 INFO None 5154928: status RUNNING/PENDING
2026-07-14 10:54:35 INFO None 5154929: status RUNNING/PENDING
2026-07-14 10:54:35 INFO None 5154930: status RUNNING/PENDING
2026-07-14 10:54:35 INFO None 5154931: status FINISHED
2026-07-14 10:54:35 INFO None 5154932: status RUNNING/PENDING
2026-07-14 10:54:35 INFO None 5154933: status RUNNING/PENDING
2026-07-14 10:54:35 INFO None 5154934: status RUNNING/PENDING
2026-07-14 10:54:35 INFO None 5154949: status RUNNING/PENDING
2026-07-14 10:54:35 INFO None 5154970: status RUNNING/PENDING
2026-07-14 10:54:35 INFO None 5154989: status RUNNING/PENDING
2026-07-14 10:54:35 INFO None 5154994: status RUNNING/PENDING
2026-07-14 10:54:35 INFO Jobs still running: ['5154903', '5154905', '5154906', '5154923', '5154928', '5154929', '5154930', '5154932', '5154933', '5154934', '5154949', '5154970', '5154989', '5154994']. Waiting...
2026-07-14 10:54:50 INFO None 5154903: status RUNNING/PENDING
2026-07-14 10:54:50 INFO None 5154905: status RUNNING/PENDING
2026-07-14 10:54:50 INFO None 5154906: status RUNNING/PENDING
2026-07-14 10:54:50 INFO None 5154923: status RUNNING/PENDING
2026-07-14 10:54:50 INFO None 5154928: status RUNNING/PENDING
2026-07-14 10:54:50 INFO None 5154929: status RUNNING/PENDING
2026-07-14 10:54:50 INFO None 5154930: status RUNNING/PENDING
2026-07-14 10:54:50 INFO None 5154931: status FINISHED
2026-07-14 10:54:50 INFO None 5154932: status RUNNING/PENDING
2026-07-14 10:54:50 INFO None 5154933: status RUNNING/PENDING
2026-07-14 10:54:50 INFO None 5154934: status RUNNING/PENDING
2026-07-14 10:54:50 INFO None 5154949: status RUNNING/PENDING
2026-07-14 10:54:50 INFO None 5154970: status RUNNING/PENDING
2026-07-14 10:54:50 INFO None 5154989: status RUNNING/PENDING
2026-07-14 10:54:50 INFO None 5154994: status RUNNING/PENDING
2026-07-14 10:54:50 INFO Jobs still running: ['5154903', '5154905', '5154906', '5154923', '5154928', '5154929', '5154930', '5154932', '5154933', '5154934', '5154949', '5154970', '5154989', '5154994']. Waiting...
2026-07-14 10:55:05 INFO None 5154903: status RUNNING/PENDING
2026-07-14 10:55:05 INFO None 5154905: status RUNNING/PENDING
2026-07-14 10:55:05 INFO None 5154906: status RUNNING/PENDING
2026-07-14 10:55:05 INFO None 5154923: status RUNNING/PENDING
2026-07-14 10:55:05 INFO None 5154928: status RUNNING/PENDING
2026-07-14 10:55:05 INFO None 5154929: status RUNNING/PENDING
2026-07-14 10:55:05 INFO None 5154930: status RUNNING/PENDING
2026-07-14 10:55:05 INFO None 5154931: status FINISHED
2026-07-14 10:55:05 INFO None 5154932: status RUNNING/PENDING
2026-07-14 10:55:05 INFO None 5154933: status RUNNING/PENDING
2026-07-14 10:55:05 INFO None 5154934: status RUNNING/PENDING
2026-07-14 10:55:05 INFO None 5154949: status RUNNING/PENDING
2026-07-14 10:55:05 INFO None 5154970: status RUNNING/PENDING
2026-07-14 10:55:05 INFO None 5154989: status RUNNING/PENDING
2026-07-14 10:55:06 INFO None 5154994: status RUNNING/PENDING
2026-07-14 10:55:06 INFO Jobs still running: ['5154903', '5154905', '5154906', '5154923', '5154928', '5154929', '5154930', '5154932', '5154933', '5154934', '5154949', '5154970', '5154989', '5154994']. Waiting...
2026-07-14 10:55:21 INFO None 5154903: status RUNNING/PENDING
2026-07-14 10:55:21 INFO None 5154905: status RUNNING/PENDING
2026-07-14 10:55:21 INFO None 5154906: status RUNNING/PENDING
2026-07-14 10:55:21 INFO None 5154923: status RUNNING/PENDING
2026-07-14 10:55:21 INFO None 5154928: status RUNNING/PENDING
2026-07-14 10:55:21 INFO None 5154929: status RUNNING/PENDING
2026-07-14 10:55:21 INFO None 5154930: status RUNNING/PENDING
2026-07-14 10:55:21 INFO None 5154931: status FINISHED
2026-07-14 10:55:21 INFO None 5154932: status RUNNING/PENDING
2026-07-14 10:55:21 INFO None 5154933: status RUNNING/PENDING
2026-07-14 10:55:21 INFO None 5154934: status RUNNING/PENDING
2026-07-14 10:55:21 INFO None 5154949: status RUNNING/PENDING
2026-07-14 10:55:21 INFO None 5154970: status RUNNING/PENDING
2026-07-14 10:55:21 INFO None 5154989: status RUNNING/PENDING
2026-07-14 10:55:21 INFO None 5154994: status RUNNING/PENDING
2026-07-14 10:55:21 INFO Jobs still running: ['5154903', '5154905', '5154906', '5154923', '5154928', '5154929', '5154930', '5154932', '5154933', '5154934', '5154949', '5154970', '5154989', '5154994']. Waiting...
2026-07-14 10:55:36 INFO None 5154903: status RUNNING/PENDING
2026-07-14 10:55:36 INFO None 5154905: status RUNNING/PENDING
2026-07-14 10:55:36 INFO None 5154906: status RUNNING/PENDING
2026-07-14 10:55:36 INFO None 5154923: status RUNNING/PENDING
2026-07-14 10:55:36 INFO None 5154928: status RUNNING/PENDING
2026-07-14 10:55:36 INFO None 5154929: status RUNNING/PENDING
2026-07-14 10:55:36 INFO None 5154930: status RUNNING/PENDING
2026-07-14 10:55:36 INFO None 5154931: status FINISHED
2026-07-14 10:55:36 INFO None 5154932: status RUNNING/PENDING
2026-07-14 10:55:36 INFO None 5154933: status RUNNING/PENDING
2026-07-14 10:55:36 INFO None 5154934: status RUNNING/PENDING
2026-07-14 10:55:36 INFO None 5154949: status RUNNING/PENDING
2026-07-14 10:55:36 INFO None 5154970: status RUNNING/PENDING
2026-07-14 10:55:36 INFO None 5154989: status RUNNING/PENDING
2026-07-14 10:55:36 INFO None 5154994: status RUNNING/PENDING
2026-07-14 10:55:36 INFO Jobs still running: ['5154903', '5154905', '5154906', '5154923', '5154928', '5154929', '5154930', '5154932', '5154933', '5154934', '5154949', '5154970', '5154989', '5154994']. Waiting...
2026-07-14 10:55:51 INFO None 5154903: status RUNNING/PENDING
2026-07-14 10:55:51 INFO None 5154905: status RUNNING/PENDING
2026-07-14 10:55:51 INFO None 5154906: status RUNNING/PENDING
2026-07-14 10:55:51 INFO None 5154923: status RUNNING/PENDING
2026-07-14 10:55:51 INFO None 5154928: status RUNNING/PENDING
2026-07-14 10:55:51 INFO None 5154929: status RUNNING/PENDING
2026-07-14 10:55:51 INFO None 5154930: status RUNNING/PENDING
2026-07-14 10:55:51 INFO None 5154931: status FINISHED
2026-07-14 10:55:51 INFO None 5154932: status RUNNING/PENDING
2026-07-14 10:55:51 INFO None 5154933: status RUNNING/PENDING
2026-07-14 10:55:51 INFO None 5154934: status RUNNING/PENDING
2026-07-14 10:55:51 INFO None 5154949: status RUNNING/PENDING
2026-07-14 10:55:51 INFO None 5154970: status RUNNING/PENDING
2026-07-14 10:55:51 INFO None 5154989: status RUNNING/PENDING
2026-07-14 10:55:51 INFO None 5154994: status RUNNING/PENDING
2026-07-14 10:55:51 INFO Jobs still running: ['5154903', '5154905', '5154906', '5154923', '5154928', '5154929', '5154930', '5154932', '5154933', '5154934', '5154949', '5154970', '5154989', '5154994']. Waiting...
2026-07-14 10:56:06 INFO None 5154903: status RUNNING/PENDING
2026-07-14 10:56:06 INFO None 5154905: status RUNNING/PENDING
2026-07-14 10:56:06 INFO None 5154906: status RUNNING/PENDING
2026-07-14 10:56:06 INFO None 5154923: status RUNNING/PENDING
2026-07-14 10:56:06 INFO None 5154928: status RUNNING/PENDING
2026-07-14 10:56:06 INFO None 5154929: status RUNNING/PENDING
2026-07-14 10:56:07 INFO None 5154930: status RUNNING/PENDING
2026-07-14 10:56:07 INFO None 5154931: status FINISHED
2026-07-14 10:56:07 INFO None 5154932: status RUNNING/PENDING
2026-07-14 10:56:07 INFO None 5154933: status RUNNING/PENDING
2026-07-14 10:56:07 INFO None 5154934: status RUNNING/PENDING
2026-07-14 10:56:07 INFO None 5154949: status RUNNING/PENDING
2026-07-14 10:56:07 INFO None 5154970: status RUNNING/PENDING
2026-07-14 10:56:07 INFO None 5154989: status RUNNING/PENDING
2026-07-14 10:56:07 INFO None 5154994: status RUNNING/PENDING
2026-07-14 10:56:07 INFO Jobs still running: ['5154903', '5154905', '5154906', '5154923', '5154928', '5154929', '5154930', '5154932', '5154933', '5154934', '5154949', '5154970', '5154989', '5154994']. Waiting...
2026-07-14 10:56:22 INFO None 5154903: status RUNNING/PENDING
2026-07-14 10:56:22 INFO None 5154905: status RUNNING/PENDING
2026-07-14 10:56:22 INFO None 5154906: status RUNNING/PENDING
2026-07-14 10:56:22 INFO None 5154923: status RUNNING/PENDING
2026-07-14 10:56:22 INFO None 5154928: status RUNNING/PENDING
2026-07-14 10:56:22 INFO None 5154929: status RUNNING/PENDING
2026-07-14 10:56:22 INFO None 5154930: status RUNNING/PENDING
2026-07-14 10:56:22 INFO None 5154931: status FINISHED
2026-07-14 10:56:22 INFO None 5154932: status RUNNING/PENDING
2026-07-14 10:56:22 INFO None 5154933: status RUNNING/PENDING
2026-07-14 10:56:22 INFO None 5154934: status RUNNING/PENDING
2026-07-14 10:56:22 INFO None 5154949: status RUNNING/PENDING
2026-07-14 10:56:22 INFO None 5154970: status RUNNING/PENDING
2026-07-14 10:56:22 INFO None 5154989: status RUNNING/PENDING
2026-07-14 10:56:22 INFO None 5154994: status RUNNING/PENDING
2026-07-14 10:56:22 INFO Jobs still running: ['5154903', '5154905', '5154906', '5154923', '5154928', '5154929', '5154930', '5154932', '5154933', '5154934', '5154949', '5154970', '5154989', '5154994']. Waiting...
2026-07-14 10:56:37 INFO None 5154903: status RUNNING/PENDING
2026-07-14 10:58:35 INFO None 5154905: status RUNNING/PENDING
2026-07-14 10:58:35 INFO None 5154906: status RUNNING/PENDING
2026-07-14 10:58:35 INFO None 5154923: status RUNNING/PENDING
2026-07-14 10:58:35 INFO None 5154928: status RUNNING/PENDING
2026-07-14 10:58:35 INFO None 5154929: status RUNNING/PENDING
2026-07-14 10:58:35 INFO None 5154930: status RUNNING/PENDING
2026-07-14 10:58:35 INFO None 5154931: status FINISHED
2026-07-14 10:58:35 INFO None 5154932: status RUNNING/PENDING
2026-07-14 10:58:35 INFO None 5154933: status RUNNING/PENDING
2026-07-14 10:58:35 INFO None 5154934: status RUNNING/PENDING
2026-07-14 10:58:35 INFO None 5154949: status RUNNING/PENDING
2026-07-14 10:58:35 INFO None 5154970: status RUNNING/PENDING
2026-07-14 10:58:35 INFO None 5154989: status RUNNING/PENDING
2026-07-14 10:58:35 INFO None 5154994: status RUNNING/PENDING
2026-07-14 10:58:35 INFO Jobs still running: ['5154903', '5154905', '5154906', '5154923', '5154928', '5154929', '5154930', '5154932', '5154933', '5154934', '5154949', '5154970', '5154989', '5154994']. Waiting...
2026-07-14 10:58:50 INFO None 5154903: status RUNNING/PENDING
2026-07-14 10:58:50 INFO None 5154905: status RUNNING/PENDING
2026-07-14 10:58:50 INFO None 5154906: status RUNNING/PENDING
2026-07-14 10:58:50 INFO None 5154923: status RUNNING/PENDING
2026-07-14 10:58:50 INFO None 5154928: status RUNNING/PENDING
2026-07-14 10:58:50 INFO None 5154929: status RUNNING/PENDING
2026-07-14 10:58:50 INFO None 5154930: status RUNNING/PENDING
2026-07-14 10:58:50 INFO None 5154931: status FINISHED
2026-07-14 10:58:50 INFO None 5154932: status RUNNING/PENDING
2026-07-14 10:58:50 INFO None 5154933: status RUNNING/PENDING
2026-07-14 10:58:50 INFO None 5154934: status RUNNING/PENDING
2026-07-14 10:58:50 INFO None 5154949: status RUNNING/PENDING
2026-07-14 10:58:50 INFO None 5154970: status RUNNING/PENDING
2026-07-14 10:58:50 INFO None 5154989: status RUNNING/PENDING
2026-07-14 10:58:50 INFO None 5154994: status RUNNING/PENDING
2026-07-14 10:58:50 INFO Jobs still running: ['5154903', '5154905', '5154906', '5154923', '5154928', '5154929', '5154930', '5154932', '5154933', '5154934', '5154949', '5154970', '5154989', '5154994']. Waiting...
2026-07-14 10:59:05 INFO None 5154903: status RUNNING/PENDING
2026-07-14 10:59:05 INFO None 5154905: status RUNNING/PENDING
2026-07-14 10:59:05 INFO None 5154906: status RUNNING/PENDING
2026-07-14 10:59:05 INFO None 5154923: status RUNNING/PENDING
2026-07-14 10:59:05 INFO None 5154928: status FINISHED
2026-07-14 10:59:05 INFO None 5154929: status RUNNING/PENDING
2026-07-14 10:59:05 INFO None 5154930: status RUNNING/PENDING
2026-07-14 10:59:05 INFO None 5154931: status FINISHED
2026-07-14 10:59:05 INFO None 5154932: status RUNNING/PENDING
2026-07-14 10:59:06 INFO None 5154933: status RUNNING/PENDING
2026-07-14 10:59:06 INFO None 5154934: status RUNNING/PENDING
2026-07-14 10:59:06 INFO None 5154949: status RUNNING/PENDING
2026-07-14 10:59:06 INFO None 5154970: status RUNNING/PENDING
2026-07-14 10:59:06 INFO None 5154989: status RUNNING/PENDING
2026-07-14 10:59:06 INFO None 5154994: status RUNNING/PENDING
2026-07-14 10:59:06 INFO Jobs still running: ['5154903', '5154905', '5154906', '5154923', '5154929', '5154930', '5154932', '5154933', '5154934', '5154949', '5154970', '5154989', '5154994']. Waiting...
2026-07-14 10:59:21 INFO None 5154903: status RUNNING/PENDING
2026-07-14 10:59:21 INFO None 5154905: status RUNNING/PENDING
2026-07-14 10:59:21 INFO None 5154906: status RUNNING/PENDING
2026-07-14 10:59:21 INFO None 5154923: status RUNNING/PENDING
2026-07-14 10:59:21 INFO None 5154928: status FINISHED
2026-07-14 10:59:21 INFO None 5154929: status RUNNING/PENDING
2026-07-14 10:59:21 INFO None 5154930: status RUNNING/PENDING
2026-07-14 10:59:21 INFO None 5154931: status FINISHED
2026-07-14 10:59:21 INFO None 5154932: status RUNNING/PENDING
2026-07-14 10:59:21 INFO None 5154933: status RUNNING/PENDING
2026-07-14 10:59:21 INFO None 5154934: status RUNNING/PENDING
2026-07-14 10:59:21 INFO None 5154949: status RUNNING/PENDING
2026-07-14 10:59:21 INFO None 5154970: status RUNNING/PENDING
2026-07-14 10:59:21 INFO None 5154989: status RUNNING/PENDING
2026-07-14 10:59:21 INFO None 5154994: status RUNNING/PENDING
2026-07-14 10:59:21 INFO Jobs still running: ['5154903', '5154905', '5154906', '5154923', '5154929', '5154930', '5154932', '5154933', '5154934', '5154949', '5154970', '5154989', '5154994']. Waiting...
2026-07-14 10:59:36 INFO None 5154903: status FINISHED
2026-07-14 10:59:36 INFO None 5154905: status RUNNING/PENDING
2026-07-14 10:59:36 INFO None 5154906: status RUNNING/PENDING
2026-07-14 10:59:36 INFO None 5154923: status RUNNING/PENDING
2026-07-14 10:59:36 INFO None 5154928: status FINISHED
2026-07-14 10:59:36 INFO None 5154929: status RUNNING/PENDING
2026-07-14 10:59:36 INFO None 5154930: status RUNNING/PENDING
2026-07-14 10:59:36 INFO None 5154931: status FINISHED
2026-07-14 10:59:36 INFO None 5154932: status RUNNING/PENDING
2026-07-14 10:59:36 INFO None 5154933: status RUNNING/PENDING
2026-07-14 10:59:36 INFO None 5154934: status RUNNING/PENDING
2026-07-14 10:59:36 INFO None 5154949: status RUNNING/PENDING
2026-07-14 10:59:36 INFO None 5154970: status RUNNING/PENDING
2026-07-14 10:59:36 INFO None 5154989: status RUNNING/PENDING
2026-07-14 10:59:36 INFO None 5154994: status RUNNING/PENDING
2026-07-14 10:59:36 INFO Jobs still running: ['5154905', '5154906', '5154923', '5154929', '5154930', '5154932', '5154933', '5154934', '5154949', '5154970', '5154989', '5154994']. Waiting...
2026-07-14 10:59:51 INFO None 5154903: status FINISHED
2026-07-14 10:59:51 INFO None 5154905: status RUNNING/PENDING
2026-07-14 10:59:51 INFO None 5154906: status RUNNING/PENDING
2026-07-14 10:59:51 INFO None 5154923: status FINISHED
2026-07-14 10:59:51 INFO None 5154928: status FINISHED
2026-07-14 10:59:51 INFO None 5154929: status RUNNING/PENDING
2026-07-14 10:59:51 INFO None 5154930: status RUNNING/PENDING
2026-07-14 10:59:51 INFO None 5154931: status FINISHED
2026-07-14 10:59:51 INFO None 5154932: status FINISHED
2026-07-14 10:59:51 INFO None 5154933: status FINISHED
2026-07-14 10:59:51 INFO None 5154934: status FINISHED
2026-07-14 10:59:51 INFO None 5154949: status RUNNING/PENDING
2026-07-14 10:59:51 INFO None 5154970: status RUNNING/PENDING
2026-07-14 10:59:51 INFO None 5154989: status RUNNING/PENDING
2026-07-14 10:59:51 INFO None 5154994: status RUNNING/PENDING
2026-07-14 10:59:51 INFO Jobs still running: ['5154905', '5154906', '5154929', '5154930', '5154949', '5154970', '5154989', '5154994']. Waiting...
2026-07-14 11:00:06 INFO None 5154903: status FINISHED
2026-07-14 11:00:06 INFO None 5154905: status RUNNING/PENDING
2026-07-14 11:00:07 INFO None 5154906: status RUNNING/PENDING
2026-07-14 11:00:07 INFO None 5154923: status FINISHED
2026-07-14 11:00:07 INFO None 5154928: status FINISHED
2026-07-14 11:00:07 INFO None 5154929: status FINISHED
2026-07-14 11:00:07 INFO None 5154930: status FINISHED
2026-07-14 11:00:07 INFO None 5154931: status FINISHED
2026-07-14 11:00:07 INFO None 5154932: status FINISHED
2026-07-14 11:00:07 INFO None 5154933: status FINISHED
2026-07-14 11:00:07 INFO None 5154934: status FINISHED
2026-07-14 11:00:07 INFO None 5154949: status RUNNING/PENDING
2026-07-14 11:00:07 INFO None 5154970: status RUNNING/PENDING
2026-07-14 11:00:07 INFO None 5154989: status RUNNING/PENDING
2026-07-14 11:00:07 INFO None 5154994: status RUNNING/PENDING
2026-07-14 11:00:07 INFO Jobs still running: ['5154905', '5154906', '5154949', '5154970', '5154989', '5154994']. Waiting...
2026-07-14 11:00:22 INFO None 5154903: status FINISHED
2026-07-14 11:00:22 INFO None 5154905: status FINISHED
2026-07-14 11:00:22 INFO None 5154906: status FINISHED
2026-07-14 11:00:22 INFO None 5154923: status FINISHED
2026-07-14 11:00:22 INFO None 5154928: status FINISHED
2026-07-14 11:00:22 INFO None 5154929: status FINISHED
2026-07-14 11:00:22 INFO None 5154930: status FINISHED
2026-07-14 11:00:22 INFO None 5154931: status FINISHED
2026-07-14 11:00:22 INFO None 5154932: status FINISHED
2026-07-14 11:00:22 INFO None 5154933: status FINISHED
2026-07-14 11:00:22 INFO None 5154934: status FINISHED
2026-07-14 11:00:22 INFO None 5154949: status RUNNING/PENDING
2026-07-14 11:00:22 INFO None 5154970: status RUNNING/PENDING
2026-07-14 11:00:22 INFO None 5154989: status RUNNING/PENDING
2026-07-14 11:00:22 INFO None 5154994: status RUNNING/PENDING
2026-07-14 11:00:22 INFO Jobs still running: ['5154949', '5154970', '5154989', '5154994']. Waiting...
2026-07-14 11:00:37 INFO None 5154903: status FINISHED
2026-07-14 11:00:37 INFO None 5154905: status FINISHED
2026-07-14 11:00:37 INFO None 5154906: status FINISHED
2026-07-14 11:00:37 INFO None 5154923: status FINISHED
2026-07-14 11:00:37 INFO None 5154928: status FINISHED
2026-07-14 11:00:37 INFO None 5154929: status FINISHED
2026-07-14 11:00:37 INFO None 5154930: status FINISHED
2026-07-14 11:00:37 INFO None 5154931: status FINISHED
2026-07-14 11:00:37 INFO None 5154932: status FINISHED
2026-07-14 11:00:37 INFO None 5154933: status FINISHED
2026-07-14 11:00:37 INFO None 5154934: status FINISHED
2026-07-14 11:00:37 INFO None 5154949: status RUNNING/PENDING
2026-07-14 11:00:37 INFO None 5154970: status RUNNING/PENDING
2026-07-14 11:00:37 INFO None 5154989: status RUNNING/PENDING
2026-07-14 11:00:37 INFO None 5154994: status RUNNING/PENDING
2026-07-14 11:00:37 INFO Jobs still running: ['5154949', '5154970', '5154989', '5154994']. Waiting...
2026-07-14 11:00:52 INFO None 5154903: status FINISHED
2026-07-14 11:00:52 INFO None 5154905: status FINISHED
2026-07-14 11:00:52 INFO None 5154906: status FINISHED
2026-07-14 11:00:52 INFO None 5154923: status FINISHED
2026-07-14 11:00:52 INFO None 5154928: status FINISHED
2026-07-14 11:00:52 INFO None 5154929: status FINISHED
2026-07-14 11:00:52 INFO None 5154930: status FINISHED
2026-07-14 11:00:52 INFO None 5154931: status FINISHED
2026-07-14 11:00:52 INFO None 5154932: status FINISHED
2026-07-14 11:00:52 INFO None 5154933: status FINISHED
2026-07-14 11:00:53 INFO None 5154934: status FINISHED
2026-07-14 11:00:53 INFO None 5154949: status FINISHED
2026-07-14 11:00:53 INFO None 5154970: status RUNNING/PENDING
2026-07-14 11:00:53 INFO None 5154989: status RUNNING/PENDING
2026-07-14 11:00:53 INFO None 5154994: status RUNNING/PENDING
2026-07-14 11:00:53 INFO Jobs still running: ['5154970', '5154989', '5154994']. Waiting...
2026-07-14 11:01:08 INFO None 5154903: status FINISHED
2026-07-14 11:01:08 INFO None 5154905: status FINISHED
2026-07-14 11:01:08 INFO None 5154906: status FINISHED
2026-07-14 11:01:08 INFO None 5154923: status FINISHED
2026-07-14 11:01:08 INFO None 5154928: status FINISHED
2026-07-14 11:01:08 INFO None 5154929: status FINISHED
2026-07-14 11:01:08 INFO None 5154930: status FINISHED
2026-07-14 11:01:08 INFO None 5154931: status FINISHED
2026-07-14 11:01:08 INFO None 5154932: status FINISHED
2026-07-14 11:01:08 INFO None 5154933: status FINISHED
2026-07-14 11:01:08 INFO None 5154934: status FINISHED
2026-07-14 11:01:08 INFO None 5154949: status FINISHED
2026-07-14 11:01:08 INFO None 5154970: status RUNNING/PENDING
2026-07-14 11:01:08 INFO None 5154989: status RUNNING/PENDING
2026-07-14 11:01:08 INFO None 5154994: status RUNNING/PENDING
2026-07-14 11:01:08 INFO Jobs still running: ['5154970', '5154989', '5154994']. Waiting...
2026-07-14 11:01:23 INFO None 5154903: status FINISHED
2026-07-14 11:03:23 INFO None 5154905: status FINISHED
2026-07-14 11:03:23 INFO None 5154906: status FINISHED
2026-07-14 11:03:23 INFO None 5154923: status FINISHED
2026-07-14 11:03:23 INFO None 5154928: status FINISHED
2026-07-14 11:03:23 INFO None 5154929: status FINISHED
2026-07-14 11:03:23 INFO None 5154930: status FINISHED
2026-07-14 11:03:23 INFO None 5154931: status FINISHED (not in squeue)
2026-07-14 11:03:23 INFO None 5154932: status FINISHED
2026-07-14 11:03:23 INFO None 5154933: status FINISHED
2026-07-14 11:03:23 INFO None 5154934: status FINISHED
2026-07-14 11:03:23 INFO None 5154949: status FINISHED
2026-07-14 11:03:23 INFO None 5154970: status RUNNING/PENDING
2026-07-14 11:03:23 INFO None 5154989: status RUNNING/PENDING
2026-07-14 11:03:23 INFO None 5154994: status RUNNING/PENDING
2026-07-14 11:03:23 INFO Jobs still running: ['5154970', '5154989', '5154994']. Waiting...
2026-07-14 11:03:38 INFO None 5154903: status FINISHED
2026-07-14 11:03:38 INFO None 5154905: status FINISHED
2026-07-14 11:03:38 INFO None 5154906: status FINISHED
2026-07-14 11:03:38 INFO None 5154923: status FINISHED
2026-07-14 11:03:38 INFO None 5154928: status FINISHED
2026-07-14 11:03:38 INFO None 5154929: status FINISHED
2026-07-14 11:03:38 INFO None 5154930: status FINISHED
2026-07-14 11:03:38 INFO None 5154931: status FINISHED (not in squeue)
2026-07-14 11:03:38 INFO None 5154932: status FINISHED
2026-07-14 11:03:38 INFO None 5154933: status FINISHED
2026-07-14 11:03:38 INFO None 5154934: status FINISHED
2026-07-14 11:03:38 INFO None 5154949: status FINISHED
2026-07-14 11:03:38 INFO None 5154970: status RUNNING/PENDING
2026-07-14 11:03:38 INFO None 5154989: status RUNNING/PENDING
2026-07-14 11:03:38 INFO None 5154994: status RUNNING/PENDING
2026-07-14 11:03:38 INFO Jobs still running: ['5154970', '5154989', '5154994']. Waiting...
2026-07-14 11:03:53 INFO None 5154903: status FINISHED
2026-07-14 11:03:53 INFO None 5154905: status FINISHED
2026-07-14 11:03:53 INFO None 5154906: status FINISHED
2026-07-14 11:03:54 INFO None 5154923: status FINISHED
2026-07-14 11:03:54 INFO None 5154928: status FINISHED
2026-07-14 11:03:54 INFO None 5154929: status FINISHED
2026-07-14 11:03:54 INFO None 5154930: status FINISHED
2026-07-14 11:03:54 INFO None 5154931: status FINISHED (not in squeue)
2026-07-14 11:03:54 INFO None 5154932: status FINISHED
2026-07-14 11:03:54 INFO None 5154933: status FINISHED
2026-07-14 11:03:54 INFO None 5154934: status FINISHED
2026-07-14 11:03:54 INFO None 5154949: status FINISHED
2026-07-14 11:03:54 INFO None 5154970: status RUNNING/PENDING
2026-07-14 11:03:54 INFO None 5154989: status RUNNING/PENDING
2026-07-14 11:03:54 INFO None 5154994: status RUNNING/PENDING
2026-07-14 11:03:54 INFO Jobs still running: ['5154970', '5154989', '5154994']. Waiting...
2026-07-14 11:04:11 INFO None 5154903: status FINISHED
2026-07-14 11:04:11 INFO None 5154905: status FINISHED
2026-07-14 11:04:11 INFO None 5154906: status FINISHED
2026-07-14 11:04:11 INFO None 5154923: status FINISHED
2026-07-14 11:04:11 INFO None 5154928: status FINISHED
2026-07-14 11:04:11 INFO None 5154929: status FINISHED
2026-07-14 11:04:11 INFO None 5154930: status FINISHED
2026-07-14 11:04:11 INFO None 5154931: status FINISHED (not in squeue)
2026-07-14 11:04:11 INFO None 5154932: status FINISHED
2026-07-14 11:04:11 INFO None 5154933: status FINISHED
2026-07-14 11:04:11 INFO None 5154934: status FINISHED
2026-07-14 11:04:11 INFO None 5154949: status FINISHED
2026-07-14 11:04:11 INFO None 5154970: status RUNNING/PENDING
2026-07-14 11:04:11 INFO None 5154989: status RUNNING/PENDING
2026-07-14 11:04:11 INFO None 5154994: status RUNNING/PENDING
2026-07-14 11:04:11 INFO Jobs still running: ['5154970', '5154989', '5154994']. Waiting...
2026-07-14 11:04:26 INFO None 5154903: status FINISHED
2026-07-14 11:04:26 INFO None 5154905: status FINISHED
2026-07-14 11:04:26 INFO None 5154906: status FINISHED
2026-07-14 11:04:26 INFO None 5154923: status FINISHED
2026-07-14 11:04:26 INFO None 5154928: status FINISHED
2026-07-14 11:04:26 INFO None 5154929: status FINISHED
2026-07-14 11:04:26 INFO None 5154930: status FINISHED
2026-07-14 11:04:26 INFO None 5154931: status FINISHED (not in squeue)
2026-07-14 11:04:26 INFO None 5154932: status FINISHED
2026-07-14 11:04:26 INFO None 5154933: status FINISHED
2026-07-14 11:04:26 INFO None 5154934: status FINISHED
2026-07-14 11:04:26 INFO None 5154949: status FINISHED
2026-07-14 11:04:26 INFO None 5154970: status RUNNING/PENDING
2026-07-14 11:04:26 INFO None 5154989: status RUNNING/PENDING
2026-07-14 11:04:26 INFO None 5154994: status FINISHED
2026-07-14 11:04:26 INFO Jobs still running: ['5154970', '5154989']. Waiting...
2026-07-14 11:04:41 INFO None 5154903: status FINISHED
2026-07-14 11:04:41 INFO None 5154905: status FINISHED
2026-07-14 11:04:41 INFO None 5154906: status FINISHED
2026-07-14 11:04:41 INFO None 5154923: status FINISHED
2026-07-14 11:04:41 INFO None 5154928: status FINISHED
2026-07-14 11:04:41 INFO None 5154929: status FINISHED
2026-07-14 11:04:41 INFO None 5154930: status FINISHED
2026-07-14 11:04:41 INFO None 5154931: status FINISHED (not in squeue)
2026-07-14 11:04:41 INFO None 5154932: status FINISHED
2026-07-14 11:04:41 INFO None 5154933: status FINISHED
2026-07-14 11:04:41 INFO None 5154934: status FINISHED
2026-07-14 11:04:41 INFO None 5154949: status FINISHED
2026-07-14 11:04:41 INFO None 5154970: status RUNNING/PENDING
2026-07-14 11:04:41 INFO None 5154989: status FINISHED
2026-07-14 11:04:41 INFO None 5154994: status FINISHED
2026-07-14 11:04:41 INFO Jobs still running: ['5154970']. Waiting...
2026-07-14 11:04:58 INFO None 5154903: status FINISHED
2026-07-14 11:04:58 INFO None 5154905: status FINISHED
2026-07-14 11:04:58 INFO None 5154906: status FINISHED
2026-07-14 11:04:58 INFO None 5154923: status FINISHED
2026-07-14 11:04:58 INFO None 5154928: status FINISHED
2026-07-14 11:04:58 INFO None 5154929: status FINISHED
2026-07-14 11:04:58 INFO None 5154930: status FINISHED
2026-07-14 11:04:58 INFO None 5154931: status FINISHED (not in squeue)
2026-07-14 11:04:58 INFO None 5154932: status FINISHED
2026-07-14 11:04:58 INFO None 5154933: status FINISHED
2026-07-14 11:04:58 INFO None 5154934: status FINISHED
2026-07-14 11:04:58 INFO None 5154949: status FINISHED
2026-07-14 11:04:58 INFO None 5154970: status RUNNING/PENDING
2026-07-14 11:04:58 INFO None 5154989: status FINISHED
2026-07-14 11:04:58 INFO None 5154994: status FINISHED
2026-07-14 11:04:58 INFO Jobs still running: ['5154970']. Waiting...
2026-07-14 11:05:13 INFO None 5154903: status FINISHED
2026-07-14 11:05:13 INFO None 5154905: status FINISHED
2026-07-14 11:05:13 INFO None 5154906: status FINISHED
2026-07-14 11:05:14 INFO None 5154923: status FINISHED
2026-07-14 11:05:14 INFO None 5154928: status FINISHED
2026-07-14 11:05:14 INFO None 5154929: status FINISHED
2026-07-14 11:05:14 INFO None 5154930: status FINISHED
2026-07-14 11:05:14 INFO None 5154931: status FINISHED (not in squeue)
2026-07-14 11:05:14 INFO None 5154932: status FINISHED
2026-07-14 11:05:14 INFO None 5154933: status FINISHED
2026-07-14 11:05:14 INFO None 5154934: status FINISHED
2026-07-14 11:05:14 INFO None 5154949: status FINISHED
2026-07-14 11:05:14 INFO None 5154970: status RUNNING/PENDING
2026-07-14 11:05:14 INFO None 5154989: status FINISHED
2026-07-14 11:05:14 INFO None 5154994: status FINISHED
2026-07-14 11:05:14 INFO Jobs still running: ['5154970']. Waiting...
2026-07-14 11:05:29 INFO None 5154903: status FINISHED
2026-07-14 11:05:29 INFO None 5154905: status FINISHED
2026-07-14 11:05:29 INFO None 5154906: status FINISHED
2026-07-14 11:05:29 INFO None 5154923: status FINISHED
2026-07-14 11:05:29 INFO None 5154928: status FINISHED
2026-07-14 11:05:29 INFO None 5154929: status FINISHED
2026-07-14 11:05:29 INFO None 5154930: status FINISHED
2026-07-14 11:05:29 INFO None 5154931: status FINISHED (not in squeue)
2026-07-14 11:05:29 INFO None 5154932: status FINISHED
2026-07-14 11:05:29 INFO None 5154933: status FINISHED
2026-07-14 11:05:29 INFO None 5154934: status FINISHED
2026-07-14 11:05:29 INFO None 5154949: status FINISHED
2026-07-14 11:05:29 INFO None 5154970: status RUNNING/PENDING
2026-07-14 11:05:29 INFO None 5154989: status FINISHED
2026-07-14 11:05:29 INFO None 5154994: status FINISHED
2026-07-14 11:05:29 INFO Jobs still running: ['5154970']. Waiting...
2026-07-14 11:05:45 INFO None 5154903: status FINISHED
2026-07-14 11:05:45 INFO None 5154905: status FINISHED
2026-07-14 11:05:45 INFO None 5154906: status FINISHED
2026-07-14 11:05:45 INFO None 5154923: status FINISHED
2026-07-14 11:05:45 INFO None 5154928: status FINISHED
2026-07-14 11:05:45 INFO None 5154929: status FINISHED
2026-07-14 11:05:45 INFO None 5154930: status FINISHED
2026-07-14 11:05:45 INFO None 5154931: status FINISHED (not in squeue)
2026-07-14 11:05:45 INFO None 5154932: status FINISHED
2026-07-14 11:05:45 INFO None 5154933: status FINISHED
2026-07-14 11:05:45 INFO None 5154934: status FINISHED
2026-07-14 11:05:45 INFO None 5154949: status FINISHED
2026-07-14 11:05:45 INFO None 5154970: status RUNNING/PENDING
2026-07-14 11:05:45 INFO None 5154989: status FINISHED
2026-07-14 11:05:45 INFO None 5154994: status FINISHED
2026-07-14 11:05:45 INFO Jobs still running: ['5154970']. Waiting...
2026-07-14 11:06:00 INFO None 5154903: status FINISHED
2026-07-14 11:06:01 INFO None 5154905: status FINISHED
2026-07-14 11:06:01 INFO None 5154906: status FINISHED
2026-07-14 11:06:01 INFO None 5154923: status FINISHED
2026-07-14 11:06:01 INFO None 5154928: status FINISHED
2026-07-14 11:06:01 INFO None 5154929: status FINISHED
2026-07-14 11:06:01 INFO None 5154930: status FINISHED
2026-07-14 11:06:01 INFO None 5154931: status FINISHED (not in squeue)
2026-07-14 11:06:01 INFO None 5154932: status FINISHED
2026-07-14 11:06:01 INFO None 5154933: status FINISHED
2026-07-14 11:06:01 INFO None 5154934: status FINISHED
2026-07-14 11:06:01 INFO None 5154949: status FINISHED
2026-07-14 11:06:01 INFO None 5154970: status RUNNING/PENDING
2026-07-14 11:06:01 INFO None 5154989: status FINISHED
2026-07-14 11:06:01 INFO None 5154994: status FINISHED
2026-07-14 11:06:01 INFO Jobs still running: ['5154970']. Waiting...
2026-07-14 11:06:16 INFO None 5154903: status FINISHED
2026-07-14 11:08:57 INFO None 5154905: status FINISHED (not in squeue)
2026-07-14 11:08:57 INFO None 5154906: status FINISHED (not in squeue)
2026-07-14 11:08:57 INFO None 5154923: status FINISHED (not in squeue)
2026-07-14 11:08:57 INFO None 5154928: status FINISHED (not in squeue)
2026-07-14 11:08:57 INFO None 5154929: status FINISHED
2026-07-14 11:08:57 INFO None 5154930: status FINISHED
2026-07-14 11:08:57 INFO None 5154931: status FINISHED (not in squeue)
2026-07-14 11:08:57 INFO None 5154932: status FINISHED
2026-07-14 11:08:57 INFO None 5154933: status FINISHED
2026-07-14 11:08:57 INFO None 5154934: status FINISHED
2026-07-14 11:08:57 INFO None 5154949: status FINISHED
2026-07-14 11:08:57 INFO None 5154970: status FINISHED
2026-07-14 11:08:57 INFO None 5154989: status FINISHED
2026-07-14 11:08:57 INFO None 5154994: status FINISHED
2026-07-14 11:08:57 INFO Jobs ['5154903', '5154905', '5154906', '5154923', '5154928', '5154929', '5154930', '5154931', '5154932', '5154933', '5154934', '5154949', '5154970', '5154989', '5154994'] have finished
2026-07-14 11:08:57 INFO Checking restart files were created ...
2026-07-14 11:08:57 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020614_10_ENS1.nc(3673513755 bytes)
2026-07-14 11:08:57 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020614_10_ENS2.nc(3673513755 bytes)
2026-07-14 11:08:57 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020614_10_ENS3.nc(3673513755 bytes)
2026-07-14 11:08:57 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020614_10_ENS4.nc(3673513755 bytes)
2026-07-14 11:08:57 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020614_10_ENS5.nc(3673513755 bytes)
2026-07-14 11:08:57 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020614_10_ENS6.nc(3673513755 bytes)
2026-07-14 11:08:57 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020614_10_ENS7.nc(3673513755 bytes)
2026-07-14 11:08:57 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020614_10_ENS8.nc
2026-07-14 11:08:57 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020614_10_ENS9.nc(3673513755 bytes)
2026-07-14 11:08:57 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020614_10_ENS10.nc(3673513755 bytes)
2026-07-14 11:08:57 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020614_10_ENS11.nc(3673513755 bytes)
2026-07-14 11:08:57 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020614_10_ENS12.nc(3673513755 bytes)
2026-07-14 11:08:57 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020614_10_ENS13.nc(3673513755 bytes)
2026-07-14 11:08:57 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020614_10_ENS14.nc(3673513755 bytes)
2026-07-14 11:08:57 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020614_10_ENS15.nc(3673513755 bytes)
2026-07-14 11:08:57 INFO Check chimere log file at: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/ENS15_2020020614.out
2026-07-14 11:08:57 ERROR [PIPELINE] Error: The following chimere ENS run(s) failed (exit code 1): [8]
Traceback (most recent call last):
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/main.py", line 103, in <module>
    pipeline.run_pipeline()
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/pipelines/base_pipeline.py", line 138, in run_pipeline
    self.run_model()
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/pipelines/chimere2023/pipeline.py", line 334, in run_model
    raise ModelRunError(f"The following chimere ENS run(s) failed (exit code 1): {mems_to_rerun}")
pipeline_errors.ModelRunError: The following chimere ENS run(s) failed (exit code 1): [8]
+ exit 0
