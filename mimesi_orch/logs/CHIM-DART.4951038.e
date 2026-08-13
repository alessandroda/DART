+ SCRIPT_PID=2973331
+ /bin/bash -x /tmp/tmp.YcX5Gr6MZJ
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
2026-06-22 17:31:09 INFO 
███    ███ ██ ███    ███ ███████ ███████ ██
████  ████ ██ ████  ████ ██      ██      ██
██ ████ ██ ██ ██ ████ ██ █████   ███████ ██
██  ██  ██ ██ ██  ██  ██ ██           ██ ██
██      ██ ██ ██      ██ ███████ ███████ ██



2026-06-22 17:31:09 INFO [PIPELINE] =======================================
2026-06-22 17:31:09 INFO [PIPELINE] Starting chimere–DART orchestrator
2026-06-22 17:31:09 INFO [PIPELINE] Config file: config/config_irene_IM.yaml
2026-06-22 17:31:09 INFO [PIPELINE] Run dir: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart
2026-06-22 17:31:09 INFO [PIPELINE] Log file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/mimesi_orchestrator_logs/chimere_DART_20260622_173109.log
2026-06-22 17:31:09 INFO [PIPELINE] =======================================
2026-06-22 17:31:10 INFO Running assimilation with model_type=ModelType.CHIMERE
2026-06-22 17:31:10 INFO Using scheduler=Scheduler.SLURM, queue=rome
2026-06-22 17:31:10 INFO [STEP] ---- TIME LOOP START ----
2026-06-22 17:31:10 INFO [TIME] step_start current_time=2020-02-06 00:00:00 simulated_time=None dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-22 17:31:10 INFO [TIME] window start=2020-02-06 00:00:00 end=2020-02-06 01:00:00 run_hours=1 has_assimilation=False
2026-06-22 17:31:10 INFO Creating directories and links for ENS1 to run chimere's parallel part
2026-06-22 17:31:10 INFO Copying EMIS ...
2026-06-22 17:31:10 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens01.nc
2026-06-22 17:31:10 INFO Linking END ...
2026-06-22 17:31:10 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020500_24_ENS1.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-22 17:31:10 INFO >> Checking links...
2026-06-22 17:31:10 INFO >> All links are good for ENS1  ...
2026-06-22 17:31:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:31:18 INFO Hourly dataset computed and listing created
2026-06-22 17:31:28 INFO Hourly dataset computed
2026-06-22 17:31:28 INFO Creating directories and links for ENS2 to run chimere's parallel part
2026-06-22 17:31:28 INFO Copying EMIS ...
2026-06-22 17:31:28 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens02.nc
2026-06-22 17:31:28 INFO Linking END ...
2026-06-22 17:31:28 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020500_24_ENS2.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-22 17:31:28 INFO >> Checking links...
2026-06-22 17:31:28 INFO >> All links are good for ENS2  ...
2026-06-22 17:31:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:31:29 INFO Hourly dataset computed and listing created
2026-06-22 17:31:29 INFO Hourly dataset computed
2026-06-22 17:31:29 INFO Creating directories and links for ENS3 to run chimere's parallel part
2026-06-22 17:31:29 INFO Copying EMIS ...
2026-06-22 17:31:30 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens03.nc
2026-06-22 17:31:30 INFO Linking END ...
2026-06-22 17:31:30 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020500_24_ENS3.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-22 17:31:30 INFO >> Checking links...
2026-06-22 17:31:30 INFO >> All links are good for ENS3  ...
2026-06-22 17:31:30 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:31:31 INFO Hourly dataset computed and listing created
2026-06-22 17:31:31 INFO Hourly dataset computed
2026-06-22 17:31:31 INFO Creating directories and links for ENS4 to run chimere's parallel part
2026-06-22 17:31:31 INFO Copying EMIS ...
2026-06-22 17:31:32 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens04.nc
2026-06-22 17:31:32 INFO Linking END ...
2026-06-22 17:31:32 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020500_24_ENS4.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-22 17:31:32 INFO >> Checking links...
2026-06-22 17:31:32 INFO >> All links are good for ENS4  ...
2026-06-22 17:31:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:31:32 INFO Hourly dataset computed and listing created
2026-06-22 17:31:33 INFO Hourly dataset computed
2026-06-22 17:31:33 INFO Creating directories and links for ENS5 to run chimere's parallel part
2026-06-22 17:31:33 INFO Copying EMIS ...
2026-06-22 17:31:33 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens05.nc
2026-06-22 17:31:33 INFO Linking END ...
2026-06-22 17:31:33 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020500_24_ENS5.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-22 17:31:33 INFO >> Checking links...
2026-06-22 17:31:33 INFO >> All links are good for ENS5  ...
2026-06-22 17:31:33 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:31:34 INFO Hourly dataset computed and listing created
2026-06-22 17:31:35 INFO Hourly dataset computed
2026-06-22 17:31:35 INFO Creating directories and links for ENS6 to run chimere's parallel part
2026-06-22 17:31:35 INFO Copying EMIS ...
2026-06-22 17:31:35 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens06.nc
2026-06-22 17:31:35 INFO Linking END ...
2026-06-22 17:31:35 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020500_24_ENS6.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-22 17:31:35 INFO >> Checking links...
2026-06-22 17:31:35 INFO >> All links are good for ENS6  ...
2026-06-22 17:31:35 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:31:36 INFO Hourly dataset computed and listing created
2026-06-22 17:31:37 INFO Hourly dataset computed
2026-06-22 17:31:37 INFO Creating directories and links for ENS7 to run chimere's parallel part
2026-06-22 17:31:37 INFO Copying EMIS ...
2026-06-22 17:31:37 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens07.nc
2026-06-22 17:31:37 INFO Linking END ...
2026-06-22 17:31:37 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020500_24_ENS7.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-22 17:31:37 INFO >> Checking links...
2026-06-22 17:31:37 INFO >> All links are good for ENS7  ...
2026-06-22 17:31:37 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:31:38 INFO Hourly dataset computed and listing created
2026-06-22 17:31:39 INFO Hourly dataset computed
2026-06-22 17:31:39 INFO Creating directories and links for ENS8 to run chimere's parallel part
2026-06-22 17:31:39 INFO Copying EMIS ...
2026-06-22 17:31:39 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens08.nc
2026-06-22 17:31:39 INFO Linking END ...
2026-06-22 17:31:39 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020500_24_ENS8.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-22 17:31:39 INFO >> Checking links...
2026-06-22 17:31:39 INFO >> All links are good for ENS8  ...
2026-06-22 17:31:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:31:40 INFO Hourly dataset computed and listing created
2026-06-22 17:31:44 INFO Hourly dataset computed
2026-06-22 17:31:44 INFO Creating directories and links for ENS9 to run chimere's parallel part
2026-06-22 17:31:44 INFO Copying EMIS ...
2026-06-22 17:31:44 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens09.nc
2026-06-22 17:31:44 INFO Linking END ...
2026-06-22 17:31:44 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020500_24_ENS9.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-22 17:31:44 INFO >> Checking links...
2026-06-22 17:31:44 INFO >> All links are good for ENS9  ...
2026-06-22 17:31:44 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:31:45 INFO Hourly dataset computed and listing created
2026-06-22 17:31:46 INFO Hourly dataset computed
2026-06-22 17:31:46 INFO Creating directories and links for ENS10 to run chimere's parallel part
2026-06-22 17:31:46 INFO Copying EMIS ...
2026-06-22 17:31:46 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens10.nc
2026-06-22 17:31:46 INFO Linking END ...
2026-06-22 17:31:46 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020500_24_ENS10.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-22 17:31:46 INFO >> Checking links...
2026-06-22 17:31:46 INFO >> All links are good for ENS10  ...
2026-06-22 17:31:46 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:31:47 INFO Hourly dataset computed and listing created
2026-06-22 17:31:47 INFO Hourly dataset computed
2026-06-22 17:31:47 INFO Creating directories and links for ENS11 to run chimere's parallel part
2026-06-22 17:31:47 INFO Copying EMIS ...
2026-06-22 17:31:48 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens11.nc
2026-06-22 17:31:48 INFO Linking END ...
2026-06-22 17:31:48 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020500_24_ENS11.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-22 17:31:48 INFO >> Checking links...
2026-06-22 17:31:48 INFO >> All links are good for ENS11  ...
2026-06-22 17:31:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:31:49 INFO Hourly dataset computed and listing created
2026-06-22 17:31:49 INFO Hourly dataset computed
2026-06-22 17:31:49 INFO Creating directories and links for ENS12 to run chimere's parallel part
2026-06-22 17:31:49 INFO Copying EMIS ...
2026-06-22 17:31:50 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens12.nc
2026-06-22 17:31:50 INFO Linking END ...
2026-06-22 17:31:50 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020500_24_ENS12.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-22 17:31:50 INFO >> Checking links...
2026-06-22 17:31:50 INFO >> All links are good for ENS12  ...
2026-06-22 17:31:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:31:51 INFO Hourly dataset computed and listing created
2026-06-22 17:31:51 INFO Hourly dataset computed
2026-06-22 17:31:51 INFO Creating directories and links for ENS13 to run chimere's parallel part
2026-06-22 17:31:51 INFO Copying EMIS ...
2026-06-22 17:31:52 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens13.nc
2026-06-22 17:31:52 INFO Linking END ...
2026-06-22 17:31:52 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020500_24_ENS13.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-22 17:31:52 INFO >> Checking links...
2026-06-22 17:31:52 INFO >> All links are good for ENS13  ...
2026-06-22 17:31:52 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:31:53 INFO Hourly dataset computed and listing created
2026-06-22 17:31:53 INFO Hourly dataset computed
2026-06-22 17:31:53 INFO Creating directories and links for ENS14 to run chimere's parallel part
2026-06-22 17:31:53 INFO Copying EMIS ...
2026-06-22 17:31:54 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens14.nc
2026-06-22 17:31:54 INFO Linking END ...
2026-06-22 17:31:54 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020500_24_ENS14.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-22 17:31:54 INFO >> Checking links...
2026-06-22 17:31:54 INFO >> All links are good for ENS14  ...
2026-06-22 17:31:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:31:55 INFO Hourly dataset computed and listing created
2026-06-22 17:31:55 INFO Hourly dataset computed
2026-06-22 17:31:55 INFO Creating directories and links for ENS15 to run chimere's parallel part
2026-06-22 17:31:55 INFO Copying EMIS ...
2026-06-22 17:31:56 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens15.nc
2026-06-22 17:31:56 INFO Linking END ...
2026-06-22 17:31:56 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020500_24_ENS15.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-22 17:31:56 INFO >> Checking links...
2026-06-22 17:31:56 INFO >> All links are good for ENS15  ...
2026-06-22 17:31:56 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:31:56 INFO Hourly dataset computed and listing created
2026-06-22 17:31:57 INFO Hourly dataset computed
2026-06-22 17:31:57 INFO ---------->>> Running CHIMERE model from 2020-02-06 00:00:00 to 2020-02-06 01:00:00
2026-06-22 17:31:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:31:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1
2026-06-22 17:31:57 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020500_24_ENS1.nc
2026-06-22 17:31:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-06-22 17:31:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:31:57 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-06-22 17:31:57 INFO Queuing job for member 1...
2026-06-22 17:31:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:31:57 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-06-22 17:31:58 INFO Found: ['4951069']
2026-06-22 17:32:03 INFO [TGCC-IRENE] Submitted job with ID:['4951069']
2026-06-22 17:32:03 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:32:03 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2
2026-06-22 17:32:03 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020500_24_ENS2.nc
2026-06-22 17:32:03 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-06-22 17:32:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:32:03 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-06-22 17:32:03 INFO Queuing job for member 2...
2026-06-22 17:32:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:32:03 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-06-22 17:32:04 INFO Found: ['4951072']
2026-06-22 17:32:09 INFO [TGCC-IRENE] Submitted job with ID:['4951072']
2026-06-22 17:32:09 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:32:09 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3
2026-06-22 17:32:09 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020500_24_ENS3.nc
2026-06-22 17:32:09 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-06-22 17:32:09 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:32:09 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-06-22 17:32:09 INFO Queuing job for member 3...
2026-06-22 17:32:09 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:32:09 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-06-22 17:32:09 INFO Found: ['4951074']
2026-06-22 17:32:14 INFO [TGCC-IRENE] Submitted job with ID:['4951074']
2026-06-22 17:32:14 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:32:14 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4
2026-06-22 17:32:14 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020500_24_ENS4.nc
2026-06-22 17:32:14 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-06-22 17:32:14 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:32:14 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-06-22 17:32:14 INFO Queuing job for member 4...
2026-06-22 17:32:14 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:32:14 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-06-22 17:32:15 INFO Found: ['4951077']
2026-06-22 17:32:20 INFO [TGCC-IRENE] Submitted job with ID:['4951077']
2026-06-22 17:32:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:32:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5
2026-06-22 17:32:20 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020500_24_ENS5.nc
2026-06-22 17:32:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-06-22 17:32:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:32:20 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-06-22 17:32:20 INFO Queuing job for member 5...
2026-06-22 17:32:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:32:20 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-06-22 17:32:21 INFO Found: ['4951079']
2026-06-22 17:32:26 INFO [TGCC-IRENE] Submitted job with ID:['4951079']
2026-06-22 17:32:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:32:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6
2026-06-22 17:32:26 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020500_24_ENS6.nc
2026-06-22 17:32:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-06-22 17:32:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:32:26 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-06-22 17:32:26 INFO Queuing job for member 6...
2026-06-22 17:32:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:32:26 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-06-22 17:32:26 INFO Found: ['4951082']
2026-06-22 17:32:31 INFO [TGCC-IRENE] Submitted job with ID:['4951082']
2026-06-22 17:32:31 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:32:31 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7
2026-06-22 17:32:31 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020500_24_ENS7.nc
2026-06-22 17:32:31 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-06-22 17:32:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:32:32 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-06-22 17:32:32 INFO Queuing job for member 7...
2026-06-22 17:32:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:32:32 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-06-22 17:32:32 INFO Found: ['4951084']
2026-06-22 17:32:37 INFO [TGCC-IRENE] Submitted job with ID:['4951084']
2026-06-22 17:32:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:32:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8
2026-06-22 17:32:37 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020500_24_ENS8.nc
2026-06-22 17:32:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-06-22 17:32:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:32:37 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-06-22 17:32:37 INFO Queuing job for member 8...
2026-06-22 17:32:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:32:37 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-06-22 17:32:38 INFO Found: ['4951086']
2026-06-22 17:32:43 INFO [TGCC-IRENE] Submitted job with ID:['4951086']
2026-06-22 17:32:43 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:32:43 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9
2026-06-22 17:32:43 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020500_24_ENS9.nc
2026-06-22 17:32:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-06-22 17:32:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:32:43 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-06-22 17:32:43 INFO Queuing job for member 9...
2026-06-22 17:32:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:32:43 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-06-22 17:32:44 INFO Found: ['4951089']
2026-06-22 17:32:49 INFO [TGCC-IRENE] Submitted job with ID:['4951089']
2026-06-22 17:32:49 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:32:49 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10
2026-06-22 17:32:49 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020500_24_ENS10.nc
2026-06-22 17:32:49 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-06-22 17:32:49 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:32:49 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-06-22 17:32:49 INFO Queuing job for member 10...
2026-06-22 17:32:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:32:49 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-06-22 17:32:50 INFO Found: ['4951093']
2026-06-22 17:32:55 INFO [TGCC-IRENE] Submitted job with ID:['4951093']
2026-06-22 17:32:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:32:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11
2026-06-22 17:32:55 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020500_24_ENS11.nc
2026-06-22 17:32:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-06-22 17:32:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:32:55 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-06-22 17:32:55 INFO Queuing job for member 11...
2026-06-22 17:32:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:32:55 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-06-22 17:32:55 INFO Found: ['4951120']
2026-06-22 17:33:00 INFO [TGCC-IRENE] Submitted job with ID:['4951120']
2026-06-22 17:33:00 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:33:00 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12
2026-06-22 17:33:00 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020500_24_ENS12.nc
2026-06-22 17:33:00 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-06-22 17:33:00 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:33:00 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-06-22 17:33:00 INFO Queuing job for member 12...
2026-06-22 17:33:00 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:33:00 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-06-22 17:33:01 INFO Found: ['4951124']
2026-06-22 17:33:06 INFO [TGCC-IRENE] Submitted job with ID:['4951124']
2026-06-22 17:33:06 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:33:06 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13
2026-06-22 17:33:06 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020500_24_ENS13.nc
2026-06-22 17:33:06 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-06-22 17:33:06 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:33:06 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-06-22 17:33:06 INFO Queuing job for member 13...
2026-06-22 17:33:06 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:33:06 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-06-22 17:33:07 INFO Found: ['4951126']
2026-06-22 17:33:12 INFO [TGCC-IRENE] Submitted job with ID:['4951126']
2026-06-22 17:33:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:33:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14
2026-06-22 17:33:12 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020500_24_ENS14.nc
2026-06-22 17:33:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-06-22 17:33:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:33:12 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-06-22 17:33:12 INFO Queuing job for member 14...
2026-06-22 17:33:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:33:12 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-06-22 17:33:13 INFO Found: ['4951128']
2026-06-22 17:33:18 INFO [TGCC-IRENE] Submitted job with ID:['4951128']
2026-06-22 17:33:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:33:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15
2026-06-22 17:33:18 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020500_24_ENS15.nc
2026-06-22 17:33:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-06-22 17:33:18 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:33:18 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-06-22 17:33:18 INFO Queuing job for member 15...
2026-06-22 17:33:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:33:18 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-06-22 17:33:19 INFO Found: ['4951131']
2026-06-22 17:33:24 INFO [TGCC-IRENE] Submitted job with ID:['4951131']
2026-06-22 17:33:24 INFO Checking job status ...
2026-06-22 17:33:24 INFO None 4951069: status RUNNING/PENDING
2026-06-22 17:33:24 INFO None 4951072: status RUNNING/PENDING
2026-06-22 17:33:24 INFO None 4951074: status RUNNING/PENDING
2026-06-22 17:33:24 INFO None 4951077: status RUNNING/PENDING
2026-06-22 17:33:24 INFO None 4951079: status RUNNING/PENDING
2026-06-22 17:33:24 INFO None 4951082: status RUNNING/PENDING
2026-06-22 17:33:24 INFO None 4951084: status RUNNING/PENDING
2026-06-22 17:33:25 INFO None 4951086: status RUNNING/PENDING
2026-06-22 17:33:25 INFO None 4951089: status RUNNING/PENDING
2026-06-22 17:33:25 INFO None 4951093: status RUNNING/PENDING
2026-06-22 17:33:25 INFO None 4951120: status RUNNING/PENDING
2026-06-22 17:33:25 INFO None 4951124: status RUNNING/PENDING
2026-06-22 17:33:25 INFO None 4951126: status RUNNING/PENDING
2026-06-22 17:33:25 INFO None 4951128: status RUNNING/PENDING
2026-06-22 17:33:25 INFO None 4951131: status RUNNING/PENDING
2026-06-22 17:33:25 INFO Jobs still running: ['4951069', '4951072', '4951074', '4951077', '4951079', '4951082', '4951084', '4951086', '4951089', '4951093', '4951120', '4951124', '4951126', '4951128', '4951131']. Waiting...
2026-06-22 17:33:40 INFO None 4951069: status RUNNING/PENDING
2026-06-22 17:33:40 INFO None 4951072: status RUNNING/PENDING
2026-06-22 17:33:40 INFO None 4951074: status RUNNING/PENDING
2026-06-22 17:33:40 INFO None 4951077: status RUNNING/PENDING
2026-06-22 17:33:40 INFO None 4951079: status RUNNING/PENDING
2026-06-22 17:33:40 INFO None 4951082: status RUNNING/PENDING
2026-06-22 17:33:40 INFO None 4951084: status RUNNING/PENDING
2026-06-22 17:33:40 INFO None 4951086: status RUNNING/PENDING
2026-06-22 17:33:40 INFO None 4951089: status RUNNING/PENDING
2026-06-22 17:33:40 INFO None 4951093: status RUNNING/PENDING
2026-06-22 17:33:40 INFO None 4951120: status RUNNING/PENDING
2026-06-22 17:33:40 INFO None 4951124: status RUNNING/PENDING
2026-06-22 17:33:40 INFO None 4951126: status RUNNING/PENDING
2026-06-22 17:33:40 INFO None 4951128: status RUNNING/PENDING
2026-06-22 17:33:40 INFO None 4951131: status RUNNING/PENDING
2026-06-22 17:33:40 INFO Jobs still running: ['4951069', '4951072', '4951074', '4951077', '4951079', '4951082', '4951084', '4951086', '4951089', '4951093', '4951120', '4951124', '4951126', '4951128', '4951131']. Waiting...
2026-06-22 17:33:55 INFO None 4951069: status RUNNING/PENDING
2026-06-22 17:33:55 INFO None 4951072: status RUNNING/PENDING
2026-06-22 17:33:55 INFO None 4951074: status RUNNING/PENDING
2026-06-22 17:33:55 INFO None 4951077: status RUNNING/PENDING
2026-06-22 17:33:55 INFO None 4951079: status RUNNING/PENDING
2026-06-22 17:33:55 INFO None 4951082: status RUNNING/PENDING
2026-06-22 17:33:55 INFO None 4951084: status RUNNING/PENDING
2026-06-22 17:33:55 INFO None 4951086: status RUNNING/PENDING
2026-06-22 17:33:55 INFO None 4951089: status RUNNING/PENDING
2026-06-22 17:33:55 INFO None 4951093: status RUNNING/PENDING
2026-06-22 17:33:55 INFO None 4951120: status RUNNING/PENDING
2026-06-22 17:33:55 INFO None 4951124: status RUNNING/PENDING
2026-06-22 17:33:55 INFO None 4951126: status RUNNING/PENDING
2026-06-22 17:33:55 INFO None 4951128: status RUNNING/PENDING
2026-06-22 17:33:55 INFO None 4951131: status RUNNING/PENDING
2026-06-22 17:33:55 INFO Jobs still running: ['4951069', '4951072', '4951074', '4951077', '4951079', '4951082', '4951084', '4951086', '4951089', '4951093', '4951120', '4951124', '4951126', '4951128', '4951131']. Waiting...
2026-06-22 17:34:10 INFO None 4951069: status RUNNING/PENDING
2026-06-22 17:34:10 INFO None 4951072: status RUNNING/PENDING
2026-06-22 17:34:10 INFO None 4951074: status RUNNING/PENDING
2026-06-22 17:34:10 INFO None 4951077: status RUNNING/PENDING
2026-06-22 17:34:10 INFO None 4951079: status RUNNING/PENDING
2026-06-22 17:34:10 INFO None 4951082: status RUNNING/PENDING
2026-06-22 17:34:10 INFO None 4951084: status RUNNING/PENDING
2026-06-22 17:34:10 INFO None 4951086: status RUNNING/PENDING
2026-06-22 17:34:10 INFO None 4951089: status RUNNING/PENDING
2026-06-22 17:34:10 INFO None 4951093: status RUNNING/PENDING
2026-06-22 17:34:10 INFO None 4951120: status RUNNING/PENDING
2026-06-22 17:34:11 INFO None 4951124: status RUNNING/PENDING
2026-06-22 17:34:11 INFO None 4951126: status RUNNING/PENDING
2026-06-22 17:34:11 INFO None 4951128: status RUNNING/PENDING
2026-06-22 17:34:11 INFO None 4951131: status RUNNING/PENDING
2026-06-22 17:34:11 INFO Jobs still running: ['4951069', '4951072', '4951074', '4951077', '4951079', '4951082', '4951084', '4951086', '4951089', '4951093', '4951120', '4951124', '4951126', '4951128', '4951131']. Waiting...
2026-06-22 17:34:26 INFO None 4951069: status RUNNING/PENDING
2026-06-22 17:34:26 INFO None 4951072: status RUNNING/PENDING
2026-06-22 17:34:26 INFO None 4951074: status RUNNING/PENDING
2026-06-22 17:34:26 INFO None 4951077: status RUNNING/PENDING
2026-06-22 17:34:26 INFO None 4951079: status RUNNING/PENDING
2026-06-22 17:34:26 INFO None 4951082: status RUNNING/PENDING
2026-06-22 17:34:26 INFO None 4951084: status RUNNING/PENDING
2026-06-22 17:34:26 INFO None 4951086: status RUNNING/PENDING
2026-06-22 17:34:26 INFO None 4951089: status RUNNING/PENDING
2026-06-22 17:34:26 INFO None 4951093: status RUNNING/PENDING
2026-06-22 17:34:26 INFO None 4951120: status RUNNING/PENDING
2026-06-22 17:34:26 INFO None 4951124: status RUNNING/PENDING
2026-06-22 17:34:26 INFO None 4951126: status RUNNING/PENDING
2026-06-22 17:34:26 INFO None 4951128: status RUNNING/PENDING
2026-06-22 17:34:26 INFO None 4951131: status RUNNING/PENDING
2026-06-22 17:34:26 INFO Jobs still running: ['4951069', '4951072', '4951074', '4951077', '4951079', '4951082', '4951084', '4951086', '4951089', '4951093', '4951120', '4951124', '4951126', '4951128', '4951131']. Waiting...
2026-06-22 17:34:41 INFO None 4951069: status RUNNING/PENDING
2026-06-22 17:34:41 INFO None 4951072: status RUNNING/PENDING
2026-06-22 17:34:41 INFO None 4951074: status FINISHED
2026-06-22 17:34:41 INFO None 4951077: status FINISHED
2026-06-22 17:34:41 INFO None 4951079: status RUNNING/PENDING
2026-06-22 17:34:41 INFO None 4951082: status RUNNING/PENDING
2026-06-22 17:34:41 INFO None 4951084: status RUNNING/PENDING
2026-06-22 17:34:41 INFO None 4951086: status RUNNING/PENDING
2026-06-22 17:34:41 INFO None 4951089: status RUNNING/PENDING
2026-06-22 17:34:41 INFO None 4951093: status RUNNING/PENDING
2026-06-22 17:34:41 INFO None 4951120: status RUNNING/PENDING
2026-06-22 17:34:41 INFO None 4951124: status RUNNING/PENDING
2026-06-22 17:34:41 INFO None 4951126: status RUNNING/PENDING
2026-06-22 17:34:41 INFO None 4951128: status RUNNING/PENDING
2026-06-22 17:34:41 INFO None 4951131: status RUNNING/PENDING
2026-06-22 17:34:41 INFO Jobs still running: ['4951069', '4951072', '4951079', '4951082', '4951084', '4951086', '4951089', '4951093', '4951120', '4951124', '4951126', '4951128', '4951131']. Waiting...
2026-06-22 17:34:56 INFO None 4951069: status RUNNING/PENDING
2026-06-22 17:34:56 INFO None 4951072: status RUNNING/PENDING
2026-06-22 17:34:56 INFO None 4951074: status FINISHED
2026-06-22 17:34:56 INFO None 4951077: status FINISHED
2026-06-22 17:34:56 INFO None 4951079: status RUNNING/PENDING
2026-06-22 17:34:56 INFO None 4951082: status RUNNING/PENDING
2026-06-22 17:34:57 INFO None 4951084: status RUNNING/PENDING
2026-06-22 17:34:57 INFO None 4951086: status RUNNING/PENDING
2026-06-22 17:34:57 INFO None 4951089: status RUNNING/PENDING
2026-06-22 17:34:57 INFO None 4951093: status RUNNING/PENDING
2026-06-22 17:34:57 INFO None 4951120: status RUNNING/PENDING
2026-06-22 17:34:57 INFO None 4951124: status RUNNING/PENDING
2026-06-22 17:34:57 INFO None 4951126: status RUNNING/PENDING
2026-06-22 17:34:57 INFO None 4951128: status RUNNING/PENDING
2026-06-22 17:34:57 INFO None 4951131: status RUNNING/PENDING
2026-06-22 17:34:57 INFO Jobs still running: ['4951069', '4951072', '4951079', '4951082', '4951084', '4951086', '4951089', '4951093', '4951120', '4951124', '4951126', '4951128', '4951131']. Waiting...
2026-06-22 17:35:12 INFO None 4951069: status RUNNING/PENDING
2026-06-22 17:35:12 INFO None 4951072: status RUNNING/PENDING
2026-06-22 17:35:12 INFO None 4951074: status FINISHED
2026-06-22 17:35:12 INFO None 4951077: status FINISHED
2026-06-22 17:35:12 INFO None 4951079: status RUNNING/PENDING
2026-06-22 17:35:12 INFO None 4951082: status FINISHED
2026-06-22 17:35:12 INFO None 4951084: status RUNNING/PENDING
2026-06-22 17:35:12 INFO None 4951086: status RUNNING/PENDING
2026-06-22 17:35:12 INFO None 4951089: status RUNNING/PENDING
2026-06-22 17:35:12 INFO None 4951093: status RUNNING/PENDING
2026-06-22 17:35:12 INFO None 4951120: status RUNNING/PENDING
2026-06-22 17:35:12 INFO None 4951124: status RUNNING/PENDING
2026-06-22 17:35:12 INFO None 4951126: status RUNNING/PENDING
2026-06-22 17:35:12 INFO None 4951128: status RUNNING/PENDING
2026-06-22 17:35:12 INFO None 4951131: status RUNNING/PENDING
2026-06-22 17:35:12 INFO Jobs still running: ['4951069', '4951072', '4951079', '4951084', '4951086', '4951089', '4951093', '4951120', '4951124', '4951126', '4951128', '4951131']. Waiting...
2026-06-22 17:35:27 INFO None 4951069: status RUNNING/PENDING
2026-06-22 17:35:27 INFO None 4951072: status RUNNING/PENDING
2026-06-22 17:35:27 INFO None 4951074: status FINISHED
2026-06-22 17:35:27 INFO None 4951077: status FINISHED
2026-06-22 17:35:27 INFO None 4951079: status RUNNING/PENDING
2026-06-22 17:35:27 INFO None 4951082: status FINISHED
2026-06-22 17:35:27 INFO None 4951084: status RUNNING/PENDING
2026-06-22 17:35:27 INFO None 4951086: status RUNNING/PENDING
2026-06-22 17:35:27 INFO None 4951089: status RUNNING/PENDING
2026-06-22 17:35:27 INFO None 4951093: status RUNNING/PENDING
2026-06-22 17:35:27 INFO None 4951120: status RUNNING/PENDING
2026-06-22 17:35:27 INFO None 4951124: status RUNNING/PENDING
2026-06-22 17:35:27 INFO None 4951126: status FINISHED
2026-06-22 17:35:27 INFO None 4951128: status FINISHED
2026-06-22 17:35:27 INFO None 4951131: status RUNNING/PENDING
2026-06-22 17:35:27 INFO Jobs still running: ['4951069', '4951072', '4951079', '4951084', '4951086', '4951089', '4951093', '4951120', '4951124', '4951131']. Waiting...
2026-06-22 17:35:42 INFO None 4951069: status RUNNING/PENDING
2026-06-22 17:35:42 INFO None 4951072: status RUNNING/PENDING
2026-06-22 17:35:42 INFO None 4951074: status FINISHED
2026-06-22 17:35:42 INFO None 4951077: status FINISHED
2026-06-22 17:35:42 INFO None 4951079: status RUNNING/PENDING
2026-06-22 17:35:42 INFO None 4951082: status FINISHED
2026-06-22 17:35:43 INFO None 4951084: status RUNNING/PENDING
2026-06-22 17:35:43 INFO None 4951086: status RUNNING/PENDING
2026-06-22 17:35:43 INFO None 4951089: status RUNNING/PENDING
2026-06-22 17:35:43 INFO None 4951093: status RUNNING/PENDING
2026-06-22 17:35:43 INFO None 4951120: status RUNNING/PENDING
2026-06-22 17:35:43 INFO None 4951124: status RUNNING/PENDING
2026-06-22 17:35:43 INFO None 4951126: status FINISHED
2026-06-22 17:35:43 INFO None 4951128: status FINISHED
2026-06-22 17:35:43 INFO None 4951131: status RUNNING/PENDING
2026-06-22 17:35:43 INFO Jobs still running: ['4951069', '4951072', '4951079', '4951084', '4951086', '4951089', '4951093', '4951120', '4951124', '4951131']. Waiting...
2026-06-22 17:35:58 INFO None 4951069: status FINISHED
2026-06-22 17:35:58 INFO None 4951072: status RUNNING/PENDING
2026-06-22 17:35:58 INFO None 4951074: status FINISHED
2026-06-22 17:35:58 INFO None 4951077: status FINISHED
2026-06-22 17:35:58 INFO None 4951079: status FINISHED
2026-06-22 17:35:58 INFO None 4951082: status FINISHED
2026-06-22 17:35:58 INFO None 4951084: status RUNNING/PENDING
2026-06-22 17:35:58 INFO None 4951086: status RUNNING/PENDING
2026-06-22 17:35:58 INFO None 4951089: status RUNNING/PENDING
2026-06-22 17:35:58 INFO None 4951093: status RUNNING/PENDING
2026-06-22 17:35:58 INFO None 4951120: status RUNNING/PENDING
2026-06-22 17:35:58 INFO None 4951124: status FINISHED
2026-06-22 17:35:58 INFO None 4951126: status FINISHED
2026-06-22 17:35:58 INFO None 4951128: status FINISHED
2026-06-22 17:35:58 INFO None 4951131: status RUNNING/PENDING
2026-06-22 17:35:58 INFO Jobs still running: ['4951072', '4951084', '4951086', '4951089', '4951093', '4951120', '4951131']. Waiting...
2026-06-22 17:36:13 INFO None 4951069: status FINISHED
2026-06-22 17:36:13 INFO None 4951072: status RUNNING/PENDING
2026-06-22 17:36:13 INFO None 4951074: status FINISHED
2026-06-22 17:36:13 INFO None 4951077: status FINISHED
2026-06-22 17:36:13 INFO None 4951079: status FINISHED
2026-06-22 17:36:13 INFO None 4951082: status FINISHED
2026-06-22 17:36:13 INFO None 4951084: status RUNNING/PENDING
2026-06-22 17:36:13 INFO None 4951086: status RUNNING/PENDING
2026-06-22 17:36:13 INFO None 4951089: status RUNNING/PENDING
2026-06-22 17:36:13 INFO None 4951093: status RUNNING/PENDING
2026-06-22 17:36:13 INFO None 4951120: status RUNNING/PENDING
2026-06-22 17:36:13 INFO None 4951124: status FINISHED
2026-06-22 17:36:13 INFO None 4951126: status FINISHED
2026-06-22 17:36:13 INFO None 4951128: status FINISHED
2026-06-22 17:36:13 INFO None 4951131: status RUNNING/PENDING
2026-06-22 17:36:13 INFO Jobs still running: ['4951072', '4951084', '4951086', '4951089', '4951093', '4951120', '4951131']. Waiting...
2026-06-22 17:36:28 INFO None 4951069: status FINISHED
2026-06-22 17:36:28 INFO None 4951072: status RUNNING/PENDING
2026-06-22 17:36:28 INFO None 4951074: status FINISHED
2026-06-22 17:36:28 INFO None 4951077: status FINISHED
2026-06-22 17:36:28 INFO None 4951079: status FINISHED
2026-06-22 17:36:28 INFO None 4951082: status FINISHED
2026-06-22 17:36:28 INFO None 4951084: status RUNNING/PENDING
2026-06-22 17:36:28 INFO None 4951086: status RUNNING/PENDING
2026-06-22 17:36:29 INFO None 4951089: status RUNNING/PENDING
2026-06-22 17:36:29 INFO None 4951093: status RUNNING/PENDING
2026-06-22 17:36:29 INFO None 4951120: status RUNNING/PENDING
2026-06-22 17:36:29 INFO None 4951124: status FINISHED
2026-06-22 17:36:29 INFO None 4951126: status FINISHED
2026-06-22 17:36:29 INFO None 4951128: status FINISHED
2026-06-22 17:36:29 INFO None 4951131: status RUNNING/PENDING
2026-06-22 17:36:29 INFO Jobs still running: ['4951072', '4951084', '4951086', '4951089', '4951093', '4951120', '4951131']. Waiting...
2026-06-22 17:36:44 INFO None 4951069: status FINISHED
2026-06-22 17:36:44 INFO None 4951072: status FINISHED
2026-06-22 17:36:44 INFO None 4951074: status FINISHED
2026-06-22 17:36:44 INFO None 4951077: status FINISHED
2026-06-22 17:36:44 INFO None 4951079: status FINISHED
2026-06-22 17:36:44 INFO None 4951082: status FINISHED
2026-06-22 17:36:44 INFO None 4951084: status RUNNING/PENDING
2026-06-22 17:36:44 INFO None 4951086: status RUNNING/PENDING
2026-06-22 17:36:44 INFO None 4951089: status RUNNING/PENDING
2026-06-22 17:36:44 INFO None 4951093: status RUNNING/PENDING
2026-06-22 17:36:44 INFO None 4951120: status FINISHED
2026-06-22 17:36:44 INFO None 4951124: status FINISHED
2026-06-22 17:36:44 INFO None 4951126: status FINISHED
2026-06-22 17:36:44 INFO None 4951128: status FINISHED
2026-06-22 17:36:44 INFO None 4951131: status RUNNING/PENDING
2026-06-22 17:36:44 INFO Jobs still running: ['4951084', '4951086', '4951089', '4951093', '4951131']. Waiting...
2026-06-22 17:36:59 INFO None 4951069: status FINISHED
2026-06-22 17:36:59 INFO None 4951072: status FINISHED
2026-06-22 17:36:59 INFO None 4951074: status FINISHED
2026-06-22 17:36:59 INFO None 4951077: status FINISHED
2026-06-22 17:36:59 INFO None 4951079: status FINISHED
2026-06-22 17:36:59 INFO None 4951082: status FINISHED
2026-06-22 17:36:59 INFO None 4951084: status RUNNING/PENDING
2026-06-22 17:36:59 INFO None 4951086: status RUNNING/PENDING
2026-06-22 17:36:59 INFO None 4951089: status RUNNING/PENDING
2026-06-22 17:36:59 INFO None 4951093: status RUNNING/PENDING
2026-06-22 17:36:59 INFO None 4951120: status FINISHED
2026-06-22 17:36:59 INFO None 4951124: status FINISHED
2026-06-22 17:36:59 INFO None 4951126: status FINISHED
2026-06-22 17:36:59 INFO None 4951128: status FINISHED
2026-06-22 17:36:59 INFO None 4951131: status RUNNING/PENDING
2026-06-22 17:36:59 INFO Jobs still running: ['4951084', '4951086', '4951089', '4951093', '4951131']. Waiting...
2026-06-22 17:37:15 INFO None 4951069: status FINISHED
2026-06-22 17:37:15 INFO None 4951072: status FINISHED
2026-06-22 17:37:15 INFO None 4951074: status FINISHED
2026-06-22 17:37:15 INFO None 4951077: status FINISHED
2026-06-22 17:37:15 INFO None 4951079: status FINISHED
2026-06-22 17:37:15 INFO None 4951082: status FINISHED
2026-06-22 17:37:15 INFO None 4951084: status FINISHED
2026-06-22 17:37:15 INFO None 4951086: status FINISHED
2026-06-22 17:37:15 INFO None 4951089: status FINISHED
2026-06-22 17:37:15 INFO None 4951093: status FINISHED
2026-06-22 17:37:15 INFO None 4951120: status FINISHED
2026-06-22 17:37:15 INFO None 4951124: status FINISHED
2026-06-22 17:37:15 INFO None 4951126: status FINISHED
2026-06-22 17:37:15 INFO None 4951128: status FINISHED
2026-06-22 17:37:15 INFO None 4951131: status RUNNING/PENDING
2026-06-22 17:37:15 INFO Jobs still running: ['4951131']. Waiting...
2026-06-22 17:37:30 INFO None 4951069: status FINISHED
2026-06-22 17:37:30 INFO None 4951072: status FINISHED
2026-06-22 17:37:30 INFO None 4951074: status FINISHED
2026-06-22 17:37:30 INFO None 4951077: status FINISHED
2026-06-22 17:37:30 INFO None 4951079: status FINISHED
2026-06-22 17:37:30 INFO None 4951082: status FINISHED
2026-06-22 17:37:30 INFO None 4951084: status FINISHED
2026-06-22 17:37:30 INFO None 4951086: status FINISHED
2026-06-22 17:37:30 INFO None 4951089: status FINISHED
2026-06-22 17:37:30 INFO None 4951093: status FINISHED
2026-06-22 17:37:30 INFO None 4951120: status FINISHED
2026-06-22 17:37:30 INFO None 4951124: status FINISHED
2026-06-22 17:37:30 INFO None 4951126: status FINISHED
2026-06-22 17:37:30 INFO None 4951128: status FINISHED
2026-06-22 17:37:30 INFO None 4951131: status RUNNING/PENDING
2026-06-22 17:37:30 INFO Jobs still running: ['4951131']. Waiting...
2026-06-22 17:37:45 INFO None 4951069: status FINISHED
2026-06-22 17:37:45 INFO None 4951072: status FINISHED
2026-06-22 17:37:45 INFO None 4951074: status FINISHED
2026-06-22 17:37:45 INFO None 4951077: status FINISHED
2026-06-22 17:37:45 INFO None 4951079: status FINISHED
2026-06-22 17:37:45 INFO None 4951082: status FINISHED
2026-06-22 17:37:45 INFO None 4951084: status FINISHED
2026-06-22 17:37:45 INFO None 4951086: status FINISHED
2026-06-22 17:37:45 INFO None 4951089: status FINISHED
2026-06-22 17:37:45 INFO None 4951093: status FINISHED
2026-06-22 17:37:45 INFO None 4951120: status FINISHED
2026-06-22 17:37:45 INFO None 4951124: status FINISHED
2026-06-22 17:37:45 INFO None 4951126: status FINISHED
2026-06-22 17:37:45 INFO None 4951128: status FINISHED
2026-06-22 17:37:45 INFO None 4951131: status RUNNING/PENDING
2026-06-22 17:37:45 INFO Jobs still running: ['4951131']. Waiting...
2026-06-22 17:38:01 INFO None 4951069: status FINISHED
2026-06-22 17:38:01 INFO None 4951072: status FINISHED
2026-06-22 17:38:01 INFO None 4951074: status FINISHED
2026-06-22 17:38:01 INFO None 4951077: status FINISHED
2026-06-22 17:38:01 INFO None 4951079: status FINISHED
2026-06-22 17:38:01 INFO None 4951082: status FINISHED
2026-06-22 17:38:01 INFO None 4951084: status FINISHED
2026-06-22 17:38:01 INFO None 4951086: status FINISHED
2026-06-22 17:38:01 INFO None 4951089: status FINISHED
2026-06-22 17:38:01 INFO None 4951093: status FINISHED
2026-06-22 17:38:01 INFO None 4951120: status FINISHED
2026-06-22 17:38:01 INFO None 4951124: status FINISHED
2026-06-22 17:38:01 INFO None 4951126: status FINISHED
2026-06-22 17:38:01 INFO None 4951128: status FINISHED
2026-06-22 17:38:01 INFO None 4951131: status RUNNING/PENDING
2026-06-22 17:38:01 INFO Jobs still running: ['4951131']. Waiting...
2026-06-22 17:38:16 INFO None 4951069: status FINISHED
2026-06-22 17:38:16 INFO None 4951072: status FINISHED
2026-06-22 17:38:16 INFO None 4951074: status FINISHED
2026-06-22 17:38:16 INFO None 4951077: status FINISHED
2026-06-22 17:38:16 INFO None 4951079: status FINISHED
2026-06-22 17:38:16 INFO None 4951082: status FINISHED
2026-06-22 17:38:16 INFO None 4951084: status FINISHED
2026-06-22 17:38:16 INFO None 4951086: status FINISHED
2026-06-22 17:38:16 INFO None 4951089: status FINISHED
2026-06-22 17:38:16 INFO None 4951093: status FINISHED
2026-06-22 17:38:16 INFO None 4951120: status FINISHED
2026-06-22 17:38:16 INFO None 4951124: status FINISHED
2026-06-22 17:38:16 INFO None 4951126: status FINISHED
2026-06-22 17:38:16 INFO None 4951128: status FINISHED
2026-06-22 17:38:16 INFO None 4951131: status RUNNING/PENDING
2026-06-22 17:38:16 INFO Jobs still running: ['4951131']. Waiting...
2026-06-22 17:38:31 INFO None 4951069: status FINISHED
2026-06-22 17:38:31 INFO None 4951072: status FINISHED
2026-06-22 17:38:31 INFO None 4951074: status FINISHED
2026-06-22 17:38:31 INFO None 4951077: status FINISHED
2026-06-22 17:38:31 INFO None 4951079: status FINISHED
2026-06-22 17:38:31 INFO None 4951082: status FINISHED
2026-06-22 17:38:31 INFO None 4951084: status FINISHED
2026-06-22 17:38:31 INFO None 4951086: status FINISHED
2026-06-22 17:38:31 INFO None 4951089: status FINISHED
2026-06-22 17:38:31 INFO None 4951093: status FINISHED
2026-06-22 17:38:31 INFO None 4951120: status FINISHED
2026-06-22 17:38:31 INFO None 4951124: status FINISHED
2026-06-22 17:38:31 INFO None 4951126: status FINISHED
2026-06-22 17:38:31 INFO None 4951128: status FINISHED
2026-06-22 17:38:31 INFO None 4951131: status RUNNING/PENDING
2026-06-22 17:38:31 INFO Jobs still running: ['4951131']. Waiting...
2026-06-22 17:38:46 INFO None 4951069: status FINISHED
2026-06-22 17:38:46 INFO None 4951072: status FINISHED
2026-06-22 17:38:46 INFO None 4951074: status FINISHED
2026-06-22 17:38:46 INFO None 4951077: status FINISHED
2026-06-22 17:38:46 INFO None 4951079: status FINISHED
2026-06-22 17:38:47 INFO None 4951082: status FINISHED
2026-06-22 17:38:47 INFO None 4951084: status FINISHED
2026-06-22 17:38:47 INFO None 4951086: status FINISHED
2026-06-22 17:38:47 INFO None 4951089: status FINISHED
2026-06-22 17:38:47 INFO None 4951093: status FINISHED
2026-06-22 17:38:47 INFO None 4951120: status FINISHED
2026-06-22 17:38:47 INFO None 4951124: status FINISHED
2026-06-22 17:38:47 INFO None 4951126: status FINISHED
2026-06-22 17:38:47 INFO None 4951128: status FINISHED
2026-06-22 17:38:47 INFO None 4951131: status FINISHED
2026-06-22 17:38:47 INFO Jobs ['4951069', '4951072', '4951074', '4951077', '4951079', '4951082', '4951084', '4951086', '4951089', '4951093', '4951120', '4951124', '4951126', '4951128', '4951131'] have finished
2026-06-22 17:38:47 INFO Checking restart files were created ...
2026-06-22 17:38:47 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020600_1_ENS1.nc(668832435 bytes)
2026-06-22 17:38:47 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020600_1_ENS2.nc(668832435 bytes)
2026-06-22 17:38:47 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020600_1_ENS3.nc(668832435 bytes)
2026-06-22 17:38:47 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020600_1_ENS4.nc(668832435 bytes)
2026-06-22 17:38:47 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020600_1_ENS5.nc(668832435 bytes)
2026-06-22 17:38:47 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020600_1_ENS6.nc(668832435 bytes)
2026-06-22 17:38:47 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020600_1_ENS7.nc(668832435 bytes)
2026-06-22 17:38:47 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020600_1_ENS8.nc(668832435 bytes)
2026-06-22 17:38:47 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020600_1_ENS9.nc(668832435 bytes)
2026-06-22 17:38:47 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020600_1_ENS10.nc(668832435 bytes)
2026-06-22 17:38:47 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020600_1_ENS11.nc(668832435 bytes)
2026-06-22 17:38:47 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020600_1_ENS12.nc(668832435 bytes)
2026-06-22 17:38:47 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020600_1_ENS13.nc(668832435 bytes)
2026-06-22 17:38:47 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020600_1_ENS14.nc(668832435 bytes)
2026-06-22 17:38:47 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020600_1_ENS15.nc(668832435 bytes)
2026-06-22 17:38:47 INFO  Run_model() completed successfully.
2026-06-22 17:38:47 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 00:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-22 17:38:47 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 01:00:00 days=153072 seconds=3600
2026-06-22 17:38:47 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-06-22 17:38:47 INFO [TIME] increment current_time 2020-02-06 00:00:00 -> 2020-02-06 01:00:00
2026-06-22 17:38:47 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-22 17:38:47 INFO ---------->>> Running process_satellite_data()
2026-06-22 17:38:47 INFO [DART] No satellite data found, skipping assimilation
2026-06-22 17:38:47 INFO after_assimilation() skipped
2026-06-22 17:38:47 INFO Next run starts from 2020-02-06 01:00:00
2026-06-22 17:38:47 INFO Cycle is DONE; starting a new loop!
2026-06-22 17:38:47 INFO [TIME] step_end current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-22 17:38:47 INFO [TIME] step_start current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-22 17:38:47 INFO [TIME] window start=2020-02-06 01:00:00 end=2020-02-06 09:00:00 run_hours=8 has_assimilation=True
2026-06-22 17:38:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:38:48 INFO Hourly dataset computed and listing created
2026-06-22 17:39:04 INFO Hourly dataset computed
2026-06-22 17:39:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:39:05 INFO Hourly dataset computed and listing created
2026-06-22 17:39:07 INFO Hourly dataset computed
2026-06-22 17:39:07 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:39:08 INFO Hourly dataset computed and listing created
2026-06-22 17:39:10 INFO Hourly dataset computed
2026-06-22 17:39:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:39:11 INFO Hourly dataset computed and listing created
2026-06-22 17:39:13 INFO Hourly dataset computed
2026-06-22 17:39:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:39:14 INFO Hourly dataset computed and listing created
2026-06-22 17:39:16 INFO Hourly dataset computed
2026-06-22 17:39:16 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:39:18 INFO Hourly dataset computed and listing created
2026-06-22 17:39:20 INFO Hourly dataset computed
2026-06-22 17:39:20 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:39:21 INFO Hourly dataset computed and listing created
2026-06-22 17:39:23 INFO Hourly dataset computed
2026-06-22 17:39:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:39:24 INFO Hourly dataset computed and listing created
2026-06-22 17:39:26 INFO Hourly dataset computed
2026-06-22 17:39:26 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:39:27 INFO Hourly dataset computed and listing created
2026-06-22 17:39:29 INFO Hourly dataset computed
2026-06-22 17:39:29 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:39:30 INFO Hourly dataset computed and listing created
2026-06-22 17:39:32 INFO Hourly dataset computed
2026-06-22 17:39:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:39:33 INFO Hourly dataset computed and listing created
2026-06-22 17:39:35 INFO Hourly dataset computed
2026-06-22 17:39:35 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:39:36 INFO Hourly dataset computed and listing created
2026-06-22 17:39:39 INFO Hourly dataset computed
2026-06-22 17:39:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:39:40 INFO Hourly dataset computed and listing created
2026-06-22 17:39:42 INFO Hourly dataset computed
2026-06-22 17:39:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:39:43 INFO Hourly dataset computed and listing created
2026-06-22 17:39:45 INFO Hourly dataset computed
2026-06-22 17:39:45 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:39:46 INFO Hourly dataset computed and listing created
2026-06-22 17:39:48 INFO Hourly dataset computed
2026-06-22 17:39:48 INFO ---------->>> Running CHIMERE model from 2020-02-06 01:00:00 to 2020-02-06 09:00:00
2026-06-22 17:39:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:39:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1
2026-06-22 17:39:48 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020600_1_ENS1.nc
2026-06-22 17:39:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-06-22 17:39:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:39:48 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-06-22 17:39:48 INFO Queuing job for member 1...
2026-06-22 17:39:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:39:48 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-06-22 17:39:50 INFO Found: ['4951347']
2026-06-22 17:39:55 INFO [TGCC-IRENE] Submitted job with ID:['4951347']
2026-06-22 17:39:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:39:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2
2026-06-22 17:39:55 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020600_1_ENS2.nc
2026-06-22 17:39:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-06-22 17:39:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:39:55 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-06-22 17:39:56 INFO Queuing job for member 2...
2026-06-22 17:39:56 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:39:56 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-06-22 17:39:56 INFO Found: ['4951349']
2026-06-22 17:40:01 INFO [TGCC-IRENE] Submitted job with ID:['4951349']
2026-06-22 17:40:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:40:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3
2026-06-22 17:40:01 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020600_1_ENS3.nc
2026-06-22 17:40:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-06-22 17:40:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:40:01 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-06-22 17:40:01 INFO Queuing job for member 3...
2026-06-22 17:40:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:40:01 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-06-22 17:40:02 INFO Found: ['4951354']
2026-06-22 17:40:07 INFO [TGCC-IRENE] Submitted job with ID:['4951354']
2026-06-22 17:40:07 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:40:07 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4
2026-06-22 17:40:07 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020600_1_ENS4.nc
2026-06-22 17:40:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-06-22 17:40:07 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:40:07 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-06-22 17:40:07 INFO Queuing job for member 4...
2026-06-22 17:40:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:40:07 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-06-22 17:40:08 INFO Found: ['4951359']
2026-06-22 17:40:13 INFO [TGCC-IRENE] Submitted job with ID:['4951359']
2026-06-22 17:40:13 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:40:13 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5
2026-06-22 17:40:13 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020600_1_ENS5.nc
2026-06-22 17:40:13 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-06-22 17:40:13 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:40:13 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-06-22 17:40:13 INFO Queuing job for member 5...
2026-06-22 17:40:13 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:40:13 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-06-22 17:40:13 INFO Found: ['4951362']
2026-06-22 17:40:18 INFO [TGCC-IRENE] Submitted job with ID:['4951362']
2026-06-22 17:40:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:40:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6
2026-06-22 17:40:18 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020600_1_ENS6.nc
2026-06-22 17:40:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-06-22 17:40:19 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:40:19 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-06-22 17:40:19 INFO Queuing job for member 6...
2026-06-22 17:40:19 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:40:19 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-06-22 17:40:19 INFO Found: ['4951366']
2026-06-22 17:40:24 INFO [TGCC-IRENE] Submitted job with ID:['4951366']
2026-06-22 17:40:24 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:40:24 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7
2026-06-22 17:40:24 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020600_1_ENS7.nc
2026-06-22 17:40:24 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-06-22 17:40:24 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:40:24 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-06-22 17:40:24 INFO Queuing job for member 7...
2026-06-22 17:40:24 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:40:24 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-06-22 17:40:25 INFO Found: ['4951368']
2026-06-22 17:40:30 INFO [TGCC-IRENE] Submitted job with ID:['4951368']
2026-06-22 17:40:30 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:40:30 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8
2026-06-22 17:40:30 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020600_1_ENS8.nc
2026-06-22 17:40:30 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-06-22 17:40:30 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:40:30 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-06-22 17:40:30 INFO Queuing job for member 8...
2026-06-22 17:40:30 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:40:30 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-06-22 17:40:31 INFO Found: ['4951370']
2026-06-22 17:40:36 INFO [TGCC-IRENE] Submitted job with ID:['4951370']
2026-06-22 17:40:36 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:40:36 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9
2026-06-22 17:40:36 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020600_1_ENS9.nc
2026-06-22 17:40:36 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-06-22 17:40:36 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:40:36 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-06-22 17:40:36 INFO Queuing job for member 9...
2026-06-22 17:40:36 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:40:36 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-06-22 17:40:36 INFO Found: ['4951372']
2026-06-22 17:40:41 INFO [TGCC-IRENE] Submitted job with ID:['4951372']
2026-06-22 17:40:41 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:40:41 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10
2026-06-22 17:40:41 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020600_1_ENS10.nc
2026-06-22 17:40:41 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-06-22 17:40:41 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:40:41 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-06-22 17:40:42 INFO Queuing job for member 10...
2026-06-22 17:40:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:40:42 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-06-22 17:40:42 INFO Found: ['4951376']
2026-06-22 17:40:47 INFO [TGCC-IRENE] Submitted job with ID:['4951376']
2026-06-22 17:40:47 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:40:47 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11
2026-06-22 17:40:47 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020600_1_ENS11.nc
2026-06-22 17:40:47 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-06-22 17:40:47 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:40:47 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-06-22 17:40:47 INFO Queuing job for member 11...
2026-06-22 17:40:47 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:40:47 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-06-22 17:40:48 INFO Found: ['4951380']
2026-06-22 17:40:53 INFO [TGCC-IRENE] Submitted job with ID:['4951380']
2026-06-22 17:40:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:40:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12
2026-06-22 17:40:53 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020600_1_ENS12.nc
2026-06-22 17:40:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-06-22 17:40:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:40:53 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-06-22 17:40:53 INFO Queuing job for member 12...
2026-06-22 17:40:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:40:53 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-06-22 17:40:54 INFO Found: ['4951383']
2026-06-22 17:40:59 INFO [TGCC-IRENE] Submitted job with ID:['4951383']
2026-06-22 17:40:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:40:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13
2026-06-22 17:40:59 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020600_1_ENS13.nc
2026-06-22 17:40:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-06-22 17:40:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:40:59 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-06-22 17:40:59 INFO Queuing job for member 13...
2026-06-22 17:40:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:40:59 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-06-22 17:41:00 INFO Found: ['4951386']
2026-06-22 17:41:05 INFO [TGCC-IRENE] Submitted job with ID:['4951386']
2026-06-22 17:41:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:41:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14
2026-06-22 17:41:05 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020600_1_ENS14.nc
2026-06-22 17:41:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-06-22 17:41:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:41:05 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-06-22 17:41:05 INFO Queuing job for member 14...
2026-06-22 17:41:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:41:05 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-06-22 17:41:06 INFO Found: ['4951390']
2026-06-22 17:41:11 INFO [TGCC-IRENE] Submitted job with ID:['4951390']
2026-06-22 17:41:11 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:41:11 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15
2026-06-22 17:41:11 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020600_1_ENS15.nc
2026-06-22 17:41:11 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-06-22 17:41:11 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:41:11 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-06-22 17:41:11 INFO Queuing job for member 15...
2026-06-22 17:41:11 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:41:11 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-06-22 17:41:11 INFO Found: ['4951392']
2026-06-22 17:41:16 INFO [TGCC-IRENE] Submitted job with ID:['4951392']
2026-06-22 17:41:16 INFO Checking job status ...
2026-06-22 17:41:16 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:41:16 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:41:16 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:41:16 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:41:16 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:41:16 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:41:16 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:41:16 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:41:16 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:41:17 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:41:17 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:41:17 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:41:17 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:41:17 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:41:17 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:41:17 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:41:32 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:41:32 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:41:32 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:41:32 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:41:32 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:41:32 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:41:32 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:41:32 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:41:32 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:41:32 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:41:32 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:41:32 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:41:32 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:41:32 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:41:32 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:41:32 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:41:47 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:41:47 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:41:47 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:41:47 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:41:47 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:41:47 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:41:47 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:41:47 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:41:47 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:41:47 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:41:47 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:41:47 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:41:47 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:41:47 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:41:47 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:41:47 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:42:02 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:42:02 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:42:02 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:42:02 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:42:02 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:42:02 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:42:02 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:42:02 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:42:02 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:42:02 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:42:02 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:42:02 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:42:03 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:42:03 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:42:03 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:42:03 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:42:19 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:42:19 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:42:19 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:42:19 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:42:19 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:42:19 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:42:19 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:42:19 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:42:19 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:42:19 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:42:19 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:42:19 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:42:19 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:42:19 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:42:19 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:42:19 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:42:34 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:42:34 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:42:35 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:42:35 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:42:35 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:42:35 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:42:35 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:42:35 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:42:35 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:42:35 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:42:35 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:42:35 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:42:35 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:42:35 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:42:35 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:42:35 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:42:50 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:42:50 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:42:50 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:42:50 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:42:50 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:42:50 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:42:50 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:42:50 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:42:50 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:42:50 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:42:50 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:42:50 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:42:50 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:42:50 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:42:50 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:42:50 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:43:07 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:43:07 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:43:07 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:43:07 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:43:07 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:43:07 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:43:07 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:43:07 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:43:07 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:43:07 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:43:07 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:43:07 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:43:10 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:43:10 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:43:10 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:43:10 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:43:25 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:43:25 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:43:25 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:43:25 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:43:25 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:43:25 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:43:25 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:43:25 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:43:25 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:43:25 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:43:25 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:43:25 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:43:25 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:43:25 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:43:25 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:43:25 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:43:40 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:43:40 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:43:40 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:43:40 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:43:40 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:43:40 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:43:40 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:43:40 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:43:40 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:43:40 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:43:40 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:43:40 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:43:40 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:43:40 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:43:40 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:43:40 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:43:55 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:43:55 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:43:55 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:43:55 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:43:55 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:43:55 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:43:55 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:43:55 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:43:57 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:43:58 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:43:58 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:43:58 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:43:58 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:43:58 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:43:58 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:43:58 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:44:13 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:44:13 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:44:13 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:44:13 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:44:13 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:44:13 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:44:13 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:44:13 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:44:13 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:44:13 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:44:13 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:44:13 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:44:13 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:44:13 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:44:13 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:44:13 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:44:28 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:44:28 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:44:28 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:44:28 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:44:28 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:44:28 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:44:28 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:44:28 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:44:28 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:44:28 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:44:28 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:44:28 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:44:28 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:44:28 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:44:28 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:44:28 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:44:43 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:44:43 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:44:43 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:44:43 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:44:43 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:44:43 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:44:43 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:44:43 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:44:43 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:44:43 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:44:43 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:44:43 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:44:43 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:44:43 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:44:44 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:44:44 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:44:59 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:44:59 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:44:59 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:44:59 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:44:59 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:44:59 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:44:59 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:44:59 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:44:59 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:44:59 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:44:59 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:44:59 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:44:59 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:44:59 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:44:59 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:44:59 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:45:14 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:45:14 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:45:14 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:45:14 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:45:14 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:45:14 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:45:14 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:45:14 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:45:14 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:45:14 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:45:14 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:45:14 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:45:14 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:45:14 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:45:14 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:45:14 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:45:29 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:45:29 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:45:29 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:45:29 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:45:29 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:45:29 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:45:29 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:45:29 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:45:29 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:45:29 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:45:29 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:45:29 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:45:29 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:45:29 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:45:29 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:45:29 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:45:44 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:45:44 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:45:45 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:45:45 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:45:45 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:45:45 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:45:45 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:45:45 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:45:45 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:45:45 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:45:45 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:45:45 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:45:45 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:45:45 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:45:45 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:45:45 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:46:00 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:46:00 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:46:00 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:46:00 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:46:00 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:46:00 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:46:00 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:46:00 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:46:00 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:46:00 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:46:00 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:46:00 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:46:00 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:46:00 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:46:00 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:46:00 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:46:15 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:46:15 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:46:15 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:46:15 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:46:15 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:46:15 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:46:15 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:46:15 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:46:15 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:46:15 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:46:15 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:46:15 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:46:15 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:46:15 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:46:15 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:46:15 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:46:30 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:46:30 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:46:30 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:46:30 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:46:31 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:46:31 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:46:31 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:46:31 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:46:31 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:46:31 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:46:31 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:46:31 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:46:31 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:46:31 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:46:31 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:46:31 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:46:46 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:46:46 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:46:46 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:46:46 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:46:46 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:46:46 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:46:46 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:46:46 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:46:46 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:46:46 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:46:46 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:46:46 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:46:46 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:46:46 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:46:46 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:46:46 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:47:01 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:47:01 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:47:01 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:47:01 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:47:01 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:47:01 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:47:01 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:47:01 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:47:01 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:47:01 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:47:01 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:47:01 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:47:01 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:47:01 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:47:01 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:47:01 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:47:17 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:47:17 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:47:17 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:47:17 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:47:17 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:47:17 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:47:17 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:47:17 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:47:17 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:47:17 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:47:17 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:47:17 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:47:17 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:47:17 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:47:17 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:47:17 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:47:32 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:47:32 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:47:32 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:47:32 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:47:32 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:47:32 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:47:33 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:47:33 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:47:33 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:47:33 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:47:33 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:47:33 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:47:33 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:47:33 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:47:33 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:47:33 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:47:48 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:47:48 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:47:48 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:47:48 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:47:48 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:47:48 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:47:48 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:47:48 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:47:48 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:47:48 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:47:48 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:47:48 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:47:48 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:47:48 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:47:48 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:47:48 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:48:03 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:48:03 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:48:03 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:48:03 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:48:03 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:48:03 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:48:03 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:48:03 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:48:03 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:48:03 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:48:03 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:48:03 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:48:03 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:48:03 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:48:03 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:48:03 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:48:18 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:48:18 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:48:18 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:48:18 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:48:18 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:48:18 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:48:18 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:48:18 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:48:18 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:48:18 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:48:19 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:48:19 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:48:19 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:48:19 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:48:19 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:48:19 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:48:34 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:48:34 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:48:34 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:48:34 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:48:34 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:48:34 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:48:34 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:48:34 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:48:34 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:48:34 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:48:34 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:48:34 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:48:34 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:48:34 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:48:34 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:48:34 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:48:49 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:48:49 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:48:49 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:48:49 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:48:49 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:48:49 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:48:49 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:48:49 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:48:49 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:48:49 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:48:49 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:48:49 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:48:49 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:48:49 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:48:49 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:48:49 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:49:04 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:49:04 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:49:04 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:49:04 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:49:04 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:49:04 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:49:04 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:49:04 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:49:04 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:49:04 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:49:04 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:49:04 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:49:04 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:49:04 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:49:05 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:49:05 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:49:20 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:49:20 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:49:20 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:49:20 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:49:20 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:49:20 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:49:20 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:49:20 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:49:20 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:49:20 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:49:20 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:49:20 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:49:20 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:49:20 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:49:20 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:49:20 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:49:35 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:49:35 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:49:35 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:49:35 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:49:35 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:49:35 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:49:35 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:49:35 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:49:35 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:49:35 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:49:35 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:49:35 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:49:35 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:49:35 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:49:35 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:49:35 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:49:50 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:49:50 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:49:50 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:49:50 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:49:50 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:49:50 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:49:50 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:49:50 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:49:50 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:49:50 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:49:50 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:49:50 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:49:50 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:49:50 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:49:50 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:49:50 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:50:05 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:50:05 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:50:06 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:50:06 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:50:06 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:50:06 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:50:06 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:50:06 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:50:06 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:50:06 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:50:06 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:50:06 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:50:06 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:50:06 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:50:06 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:50:06 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:50:21 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:50:21 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:50:21 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:50:21 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:50:21 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:50:21 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:50:21 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:50:21 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:50:21 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:50:21 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:50:21 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:50:21 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:50:21 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:50:21 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:50:21 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:50:21 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:50:37 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:50:37 INFO None 4951349: status RUNNING/PENDING
2026-06-22 17:50:37 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:50:37 INFO None 4951359: status RUNNING/PENDING
2026-06-22 17:50:37 INFO None 4951362: status RUNNING/PENDING
2026-06-22 17:50:37 INFO None 4951366: status RUNNING/PENDING
2026-06-22 17:50:37 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:50:37 INFO None 4951370: status RUNNING/PENDING
2026-06-22 17:50:37 INFO None 4951372: status RUNNING/PENDING
2026-06-22 17:50:37 INFO None 4951376: status RUNNING/PENDING
2026-06-22 17:50:37 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:50:37 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:50:37 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:50:37 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:50:37 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:50:37 INFO Jobs still running: ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:50:52 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:50:53 INFO None 4951349: status FINISHED
2026-06-22 17:50:53 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:50:53 INFO None 4951359: status FINISHED
2026-06-22 17:50:53 INFO None 4951362: status FINISHED
2026-06-22 17:50:53 INFO None 4951366: status FINISHED
2026-06-22 17:50:53 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:50:53 INFO None 4951370: status FINISHED
2026-06-22 17:50:53 INFO None 4951372: status FINISHED
2026-06-22 17:50:53 INFO None 4951376: status FINISHED
2026-06-22 17:50:53 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:50:53 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:50:53 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:50:53 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:50:53 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:50:53 INFO Jobs still running: ['4951347', '4951354', '4951368', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:51:08 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:51:08 INFO None 4951349: status FINISHED
2026-06-22 17:51:08 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:51:08 INFO None 4951359: status FINISHED
2026-06-22 17:51:08 INFO None 4951362: status FINISHED
2026-06-22 17:51:08 INFO None 4951366: status FINISHED
2026-06-22 17:51:08 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:51:08 INFO None 4951370: status FINISHED
2026-06-22 17:51:08 INFO None 4951372: status FINISHED
2026-06-22 17:51:08 INFO None 4951376: status FINISHED
2026-06-22 17:51:08 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:51:08 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:51:08 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:51:08 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:51:08 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:51:08 INFO Jobs still running: ['4951347', '4951354', '4951368', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:51:23 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:51:23 INFO None 4951349: status FINISHED
2026-06-22 17:51:23 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:51:23 INFO None 4951359: status FINISHED
2026-06-22 17:51:23 INFO None 4951362: status FINISHED
2026-06-22 17:51:23 INFO None 4951366: status FINISHED
2026-06-22 17:51:23 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:51:23 INFO None 4951370: status FINISHED
2026-06-22 17:51:23 INFO None 4951372: status FINISHED
2026-06-22 17:51:23 INFO None 4951376: status FINISHED
2026-06-22 17:51:23 INFO None 4951380: status RUNNING/PENDING
2026-06-22 17:51:23 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:51:23 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:51:23 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:51:23 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:51:23 INFO Jobs still running: ['4951347', '4951354', '4951368', '4951380', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:51:38 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:51:38 INFO None 4951349: status FINISHED
2026-06-22 17:51:38 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:51:38 INFO None 4951359: status FINISHED
2026-06-22 17:51:38 INFO None 4951362: status FINISHED
2026-06-22 17:51:39 INFO None 4951366: status FINISHED
2026-06-22 17:51:39 INFO None 4951368: status RUNNING/PENDING
2026-06-22 17:51:39 INFO None 4951370: status FINISHED
2026-06-22 17:51:39 INFO None 4951372: status FINISHED
2026-06-22 17:51:39 INFO None 4951376: status FINISHED
2026-06-22 17:51:39 INFO None 4951380: status FINISHED
2026-06-22 17:51:39 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:51:39 INFO None 4951386: status RUNNING/PENDING
2026-06-22 17:51:39 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:51:39 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:51:39 INFO Jobs still running: ['4951347', '4951354', '4951368', '4951383', '4951386', '4951390', '4951392']. Waiting...
2026-06-22 17:51:54 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:51:54 INFO None 4951349: status FINISHED
2026-06-22 17:51:54 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:51:54 INFO None 4951359: status FINISHED
2026-06-22 17:51:54 INFO None 4951362: status FINISHED
2026-06-22 17:51:54 INFO None 4951366: status FINISHED
2026-06-22 17:51:54 INFO None 4951368: status FINISHED
2026-06-22 17:51:54 INFO None 4951370: status FINISHED
2026-06-22 17:51:54 INFO None 4951372: status FINISHED
2026-06-22 17:51:54 INFO None 4951376: status FINISHED
2026-06-22 17:51:54 INFO None 4951380: status FINISHED
2026-06-22 17:51:54 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:51:54 INFO None 4951386: status FINISHED
2026-06-22 17:51:54 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:51:54 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:51:54 INFO Jobs still running: ['4951347', '4951354', '4951383', '4951390', '4951392']. Waiting...
2026-06-22 17:52:09 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:52:09 INFO None 4951349: status FINISHED
2026-06-22 17:52:09 INFO None 4951354: status RUNNING/PENDING
2026-06-22 17:52:09 INFO None 4951359: status FINISHED
2026-06-22 17:52:09 INFO None 4951362: status FINISHED
2026-06-22 17:52:09 INFO None 4951366: status FINISHED
2026-06-22 17:52:09 INFO None 4951368: status FINISHED
2026-06-22 17:52:09 INFO None 4951370: status FINISHED
2026-06-22 17:52:09 INFO None 4951372: status FINISHED
2026-06-22 17:52:09 INFO None 4951376: status FINISHED
2026-06-22 17:52:09 INFO None 4951380: status FINISHED
2026-06-22 17:52:09 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:52:09 INFO None 4951386: status FINISHED
2026-06-22 17:52:09 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:52:09 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:52:09 INFO Jobs still running: ['4951347', '4951354', '4951383', '4951390', '4951392']. Waiting...
2026-06-22 17:52:24 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:52:24 INFO None 4951349: status FINISHED
2026-06-22 17:52:24 INFO None 4951354: status FINISHED
2026-06-22 17:52:24 INFO None 4951359: status FINISHED
2026-06-22 17:52:24 INFO None 4951362: status FINISHED
2026-06-22 17:52:25 INFO None 4951366: status FINISHED
2026-06-22 17:52:25 INFO None 4951368: status FINISHED
2026-06-22 17:52:25 INFO None 4951370: status FINISHED
2026-06-22 17:52:25 INFO None 4951372: status FINISHED
2026-06-22 17:52:25 INFO None 4951376: status FINISHED
2026-06-22 17:52:25 INFO None 4951380: status FINISHED
2026-06-22 17:52:25 INFO None 4951383: status RUNNING/PENDING
2026-06-22 17:52:25 INFO None 4951386: status FINISHED
2026-06-22 17:52:25 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:52:25 INFO None 4951392: status RUNNING/PENDING
2026-06-22 17:52:25 INFO Jobs still running: ['4951347', '4951383', '4951390', '4951392']. Waiting...
2026-06-22 17:52:40 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:52:40 INFO None 4951349: status FINISHED
2026-06-22 17:52:40 INFO None 4951354: status FINISHED
2026-06-22 17:52:40 INFO None 4951359: status FINISHED
2026-06-22 17:52:40 INFO None 4951362: status FINISHED
2026-06-22 17:52:40 INFO None 4951366: status FINISHED
2026-06-22 17:52:40 INFO None 4951368: status FINISHED
2026-06-22 17:52:40 INFO None 4951370: status FINISHED
2026-06-22 17:52:40 INFO None 4951372: status FINISHED
2026-06-22 17:52:40 INFO None 4951376: status FINISHED
2026-06-22 17:52:40 INFO None 4951380: status FINISHED
2026-06-22 17:52:40 INFO None 4951383: status FINISHED
2026-06-22 17:52:40 INFO None 4951386: status FINISHED
2026-06-22 17:52:40 INFO None 4951390: status RUNNING/PENDING
2026-06-22 17:52:40 INFO None 4951392: status FINISHED
2026-06-22 17:52:40 INFO Jobs still running: ['4951347', '4951390']. Waiting...
2026-06-22 17:52:55 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:52:55 INFO None 4951349: status FINISHED
2026-06-22 17:52:55 INFO None 4951354: status FINISHED
2026-06-22 17:52:55 INFO None 4951359: status FINISHED
2026-06-22 17:52:55 INFO None 4951362: status FINISHED
2026-06-22 17:52:55 INFO None 4951366: status FINISHED
2026-06-22 17:52:55 INFO None 4951368: status FINISHED
2026-06-22 17:52:55 INFO None 4951370: status FINISHED
2026-06-22 17:52:55 INFO None 4951372: status FINISHED
2026-06-22 17:52:55 INFO None 4951376: status FINISHED
2026-06-22 17:52:55 INFO None 4951380: status FINISHED
2026-06-22 17:52:55 INFO None 4951383: status FINISHED
2026-06-22 17:52:55 INFO None 4951386: status FINISHED
2026-06-22 17:52:55 INFO None 4951390: status FINISHED
2026-06-22 17:52:55 INFO None 4951392: status FINISHED
2026-06-22 17:52:55 INFO Jobs still running: ['4951347']. Waiting...
2026-06-22 17:53:10 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:53:10 INFO None 4951349: status FINISHED
2026-06-22 17:53:10 INFO None 4951354: status FINISHED
2026-06-22 17:53:10 INFO None 4951359: status FINISHED
2026-06-22 17:53:10 INFO None 4951362: status FINISHED
2026-06-22 17:53:10 INFO None 4951366: status FINISHED
2026-06-22 17:53:10 INFO None 4951368: status FINISHED
2026-06-22 17:53:10 INFO None 4951370: status FINISHED
2026-06-22 17:53:11 INFO None 4951372: status FINISHED
2026-06-22 17:53:11 INFO None 4951376: status FINISHED
2026-06-22 17:53:11 INFO None 4951380: status FINISHED
2026-06-22 17:53:11 INFO None 4951383: status FINISHED
2026-06-22 17:53:11 INFO None 4951386: status FINISHED
2026-06-22 17:53:11 INFO None 4951390: status FINISHED
2026-06-22 17:53:11 INFO None 4951392: status FINISHED
2026-06-22 17:53:11 INFO Jobs still running: ['4951347']. Waiting...
2026-06-22 17:53:26 INFO None 4951347: status RUNNING/PENDING
2026-06-22 17:53:26 INFO None 4951349: status FINISHED
2026-06-22 17:53:26 INFO None 4951354: status FINISHED
2026-06-22 17:53:26 INFO None 4951359: status FINISHED
2026-06-22 17:53:26 INFO None 4951362: status FINISHED
2026-06-22 17:53:26 INFO None 4951366: status FINISHED
2026-06-22 17:53:26 INFO None 4951368: status FINISHED
2026-06-22 17:53:26 INFO None 4951370: status FINISHED
2026-06-22 17:53:26 INFO None 4951372: status FINISHED
2026-06-22 17:53:26 INFO None 4951376: status FINISHED
2026-06-22 17:53:26 INFO None 4951380: status FINISHED
2026-06-22 17:53:26 INFO None 4951383: status FINISHED
2026-06-22 17:53:26 INFO None 4951386: status FINISHED
2026-06-22 17:53:26 INFO None 4951390: status FINISHED
2026-06-22 17:53:26 INFO None 4951392: status FINISHED
2026-06-22 17:53:26 INFO Jobs still running: ['4951347']. Waiting...
2026-06-22 17:53:41 INFO None 4951347: status FINISHED
2026-06-22 17:53:41 INFO None 4951349: status FINISHED
2026-06-22 17:53:41 INFO None 4951354: status FINISHED
2026-06-22 17:53:41 INFO None 4951359: status FINISHED
2026-06-22 17:53:41 INFO None 4951362: status FINISHED
2026-06-22 17:53:41 INFO None 4951366: status FINISHED
2026-06-22 17:53:41 INFO None 4951368: status FINISHED
2026-06-22 17:53:41 INFO None 4951370: status FINISHED
2026-06-22 17:53:41 INFO None 4951372: status FINISHED
2026-06-22 17:53:41 INFO None 4951376: status FINISHED
2026-06-22 17:53:41 INFO None 4951380: status FINISHED
2026-06-22 17:53:41 INFO None 4951383: status FINISHED
2026-06-22 17:53:41 INFO None 4951386: status FINISHED
2026-06-22 17:53:41 INFO None 4951390: status FINISHED
2026-06-22 17:53:41 INFO None 4951392: status FINISHED
2026-06-22 17:53:41 INFO Jobs ['4951347', '4951349', '4951354', '4951359', '4951362', '4951366', '4951368', '4951370', '4951372', '4951376', '4951380', '4951383', '4951386', '4951390', '4951392'] have finished
2026-06-22 17:53:41 INFO Checking restart files were created ...
2026-06-22 17:53:41 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020601_8_ENS1.nc(3005806795 bytes)
2026-06-22 17:53:41 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020601_8_ENS2.nc(3005806795 bytes)
2026-06-22 17:53:41 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020601_8_ENS3.nc(3005806795 bytes)
2026-06-22 17:53:41 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020601_8_ENS4.nc(3005806795 bytes)
2026-06-22 17:53:41 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020601_8_ENS5.nc(3005806795 bytes)
2026-06-22 17:53:41 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020601_8_ENS6.nc(3005806795 bytes)
2026-06-22 17:53:41 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020601_8_ENS7.nc(3005806795 bytes)
2026-06-22 17:53:41 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020601_8_ENS8.nc(3005806795 bytes)
2026-06-22 17:53:41 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020601_8_ENS9.nc(3005806795 bytes)
2026-06-22 17:53:41 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020601_8_ENS10.nc(3005806795 bytes)
2026-06-22 17:53:41 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020601_8_ENS11.nc(3005806795 bytes)
2026-06-22 17:53:41 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020601_8_ENS12.nc(3005806795 bytes)
2026-06-22 17:53:41 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020601_8_ENS13.nc(3005806795 bytes)
2026-06-22 17:53:41 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020601_8_ENS14.nc(3005806795 bytes)
2026-06-22 17:53:41 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020601_8_ENS15.nc(3005806795 bytes)
2026-06-22 17:53:41 INFO  Run_model() completed successfully.
2026-06-22 17:53:41 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-22 17:53:41 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 09:00:00 days=153072 seconds=32400
2026-06-22 17:53:41 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-06-22 17:53:41 INFO [TIME] increment current_time 2020-02-06 01:00:00 -> 2020-02-06 09:00:00
2026-06-22 17:53:41 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-22 17:53:41 INFO ---------->>> Running process_satellite_data()
2026-06-22 17:53:41 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12002.nc
2026-06-22 17:53:41 INFO ---------->>> Running run_obs_converter()
2026-06-22 17:53:41 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_33067_153072.out
2026-06-22 17:53:41 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_33067_153072.out
2026-06-22 17:53:41 INFO ---------->>> Running DART
2026-06-22 17:53:41 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-06-22 17:53:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/chim_ENS1_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/chim_ENS1_2020020609_1_out_toDART.nc
2026-06-22 17:53:42 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/chim_ENS2_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/chim_ENS2_2020020609_1_out_toDART.nc
2026-06-22 17:53:42 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/chim_ENS3_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/chim_ENS3_2020020609_1_out_toDART.nc
2026-06-22 17:53:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/chim_ENS4_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/chim_ENS4_2020020609_1_out_toDART.nc
2026-06-22 17:53:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/chim_ENS5_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/chim_ENS5_2020020609_1_out_toDART.nc
2026-06-22 17:53:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/chim_ENS6_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/chim_ENS6_2020020609_1_out_toDART.nc
2026-06-22 17:53:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/chim_ENS7_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/chim_ENS7_2020020609_1_out_toDART.nc
2026-06-22 17:53:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/chim_ENS8_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/chim_ENS8_2020020609_1_out_toDART.nc
2026-06-22 17:53:45 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/chim_ENS9_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/chim_ENS9_2020020609_1_out_toDART.nc
2026-06-22 17:53:45 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/chim_ENS10_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/chim_ENS10_2020020609_1_out_toDART.nc
2026-06-22 17:53:45 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/chim_ENS11_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/chim_ENS11_2020020609_1_out_toDART.nc
2026-06-22 17:53:46 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/chim_ENS12_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/chim_ENS12_2020020609_1_out_toDART.nc
2026-06-22 17:53:46 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/chim_ENS13_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/chim_ENS13_2020020609_1_out_toDART.nc
2026-06-22 17:53:47 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/chim_ENS14_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/chim_ENS14_2020020609_1_out_toDART.nc
2026-06-22 17:53:47 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/chim_ENS15_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/chim_ENS15_2020020609_1_out_toDART.nc
2026-06-22 17:53:48 INFO Replacement input_template.nml → input.nml completed successfully.
2026-06-22 17:53:48 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-06-22 17:53:48 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-06-22 17:53:48 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-06-22 17:53:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-06-22 17:53:48 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-06-22 17:53:55 INFO Found: []
2026-06-22 17:53:55 INFO No job id returned by command ./run_filter.bsh
2026-06-22 17:53:55 INFO No monitoring will be performed
2026-06-22 17:53:55 INFO Moving DART output files to analysis and preassim directories for date 2020020609 if present ...
2026-06-22 17:53:55 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-06-22 17:53:55 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-06-22 17:53:55 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-06-22 17:53:55 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-06-22 17:53:55 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-06-22 17:53:55 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-06-22 17:53:55 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-06-22 17:53:55 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-06-22 17:53:55 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-06-22 17:53:56 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-06-22 17:53:56 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-06-22 17:53:56 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-06-22 17:53:56 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-06-22 17:53:56 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-06-22 17:53:56 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-06-22 17:53:56 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-06-22 17:53:56 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-06-22 17:53:56 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-06-22 17:53:56 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-06-22 17:53:56 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-06-22 17:53:56 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-06-22 17:53:56 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-06-22 17:53:56 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-06-22 17:53:56 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-06-22 17:53:56 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-06-22 17:53:56 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-06-22 17:53:56 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-06-22 17:53:56 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-06-22 17:53:56 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-06-22 17:53:56 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-06-22 17:53:56 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-06-22 17:53:56 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020609'
2026-06-22 17:53:56 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-06-22 17:53:56 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020609'
2026-06-22 17:53:56 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-06-22 17:53:58 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-06-22 17:53:58 INFO run_dart() is DONE.
2026-06-22 17:53:58 INFO ---------->>> Running update_pollutant_in_end()
2026-06-22 17:53:59 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 17:53:59 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 17:54:00 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 17:54:00 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 17:54:00 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020601_8_ENS1.nc
2026-06-22 17:54:14 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 17:54:15 INFO [Posterior Diff ENS1] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS1_2020020609.nc
2026-06-22 17:54:15 INFO [Posterior Diff ENS1] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS1_2020020609.relative.nc
2026-06-22 17:54:15 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 17:54:15 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 17:54:16 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 17:54:16 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 17:54:17 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020601_8_ENS2.nc
2026-06-22 17:54:30 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 17:54:31 INFO [Posterior Diff ENS2] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS2_2020020609.nc
2026-06-22 17:54:31 INFO [Posterior Diff ENS2] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS2_2020020609.relative.nc
2026-06-22 17:54:31 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 17:54:31 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 17:54:32 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 17:54:32 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 17:54:33 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020601_8_ENS3.nc
2026-06-22 17:54:46 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 17:54:46 INFO [Posterior Diff ENS3] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS3_2020020609.nc
2026-06-22 17:54:46 INFO [Posterior Diff ENS3] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS3_2020020609.relative.nc
2026-06-22 17:54:47 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 17:54:47 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 17:54:47 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 17:54:47 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 17:54:48 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020601_8_ENS4.nc
2026-06-22 17:55:01 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 17:55:02 INFO [Posterior Diff ENS4] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS4_2020020609.nc
2026-06-22 17:55:02 INFO [Posterior Diff ENS4] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS4_2020020609.relative.nc
2026-06-22 17:55:02 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 17:55:02 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 17:55:03 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 17:55:03 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 17:55:04 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020601_8_ENS5.nc
2026-06-22 17:55:16 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 17:55:17 INFO [Posterior Diff ENS5] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS5_2020020609.nc
2026-06-22 17:55:17 INFO [Posterior Diff ENS5] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS5_2020020609.relative.nc
2026-06-22 17:55:17 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 17:55:17 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 17:55:18 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 17:55:18 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 17:55:19 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020601_8_ENS6.nc
2026-06-22 17:55:32 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 17:55:32 INFO [Posterior Diff ENS6] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS6_2020020609.nc
2026-06-22 17:55:32 INFO [Posterior Diff ENS6] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS6_2020020609.relative.nc
2026-06-22 17:55:33 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 17:55:33 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 17:55:34 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 17:55:34 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 17:55:34 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020601_8_ENS7.nc
2026-06-22 17:55:48 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 17:55:48 INFO [Posterior Diff ENS7] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS7_2020020609.nc
2026-06-22 17:55:48 INFO [Posterior Diff ENS7] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS7_2020020609.relative.nc
2026-06-22 17:55:48 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 17:55:48 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 17:55:49 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 17:55:49 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 17:55:50 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020601_8_ENS8.nc
2026-06-22 17:56:03 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 17:56:04 INFO [Posterior Diff ENS8] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS8_2020020609.nc
2026-06-22 17:56:04 INFO [Posterior Diff ENS8] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS8_2020020609.relative.nc
2026-06-22 17:56:04 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 17:56:04 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 17:56:05 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 17:56:05 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 17:56:06 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020601_8_ENS9.nc
2026-06-22 17:56:19 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 17:56:20 INFO [Posterior Diff ENS9] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS9_2020020609.nc
2026-06-22 17:56:20 INFO [Posterior Diff ENS9] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS9_2020020609.relative.nc
2026-06-22 17:56:20 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 17:56:20 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 17:56:21 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 17:56:21 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 17:56:22 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020601_8_ENS10.nc
2026-06-22 17:56:36 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 17:56:36 INFO [Posterior Diff ENS10] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS10_2020020609.nc
2026-06-22 17:56:36 INFO [Posterior Diff ENS10] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS10_2020020609.relative.nc
2026-06-22 17:56:37 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 17:56:37 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 17:56:38 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 17:56:38 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 17:56:39 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020601_8_ENS11.nc
2026-06-22 17:56:52 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 17:56:52 INFO [Posterior Diff ENS11] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS11_2020020609.nc
2026-06-22 17:56:52 INFO [Posterior Diff ENS11] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS11_2020020609.relative.nc
2026-06-22 17:56:53 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 17:56:53 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 17:56:54 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 17:56:54 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 17:56:55 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020601_8_ENS12.nc
2026-06-22 17:57:08 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 17:57:08 INFO [Posterior Diff ENS12] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS12_2020020609.nc
2026-06-22 17:57:08 INFO [Posterior Diff ENS12] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS12_2020020609.relative.nc
2026-06-22 17:57:08 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 17:57:08 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 17:57:10 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 17:57:10 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 17:57:10 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020601_8_ENS13.nc
2026-06-22 17:57:24 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 17:57:24 INFO [Posterior Diff ENS13] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS13_2020020609.nc
2026-06-22 17:57:24 INFO [Posterior Diff ENS13] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS13_2020020609.relative.nc
2026-06-22 17:57:25 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 17:57:25 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 17:57:25 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 17:57:25 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 17:57:26 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020601_8_ENS14.nc
2026-06-22 17:57:39 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 17:57:39 INFO [Posterior Diff ENS14] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS14_2020020609.nc
2026-06-22 17:57:39 INFO [Posterior Diff ENS14] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS14_2020020609.relative.nc
2026-06-22 17:57:40 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 17:57:40 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 17:57:41 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 17:57:41 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 17:57:41 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020601_8_ENS15.nc
2026-06-22 17:57:55 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 17:57:55 INFO [Posterior Diff ENS15] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS15_2020020609.nc
2026-06-22 17:57:55 INFO [Posterior Diff ENS15] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020609/diff_posterior_ENS15_2020020609.relative.nc
2026-06-22 17:57:55 INFO Next run starts from 2020-02-06 09:00:00
2026-06-22 17:57:55 INFO Cycle is DONE; starting a new loop!
2026-06-22 17:57:55 INFO [TIME] step_end current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-22 17:57:55 INFO [TIME] step_start current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-22 17:57:55 INFO [TIME] window start=2020-02-06 09:00:00 end=2020-02-06 11:00:00 run_hours=2 has_assimilation=True
2026-06-22 17:57:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:57:56 INFO Hourly dataset computed and listing created
2026-06-22 17:58:00 INFO Hourly dataset computed
2026-06-22 17:58:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:58:01 INFO Hourly dataset computed and listing created
2026-06-22 17:58:02 INFO Hourly dataset computed
2026-06-22 17:58:02 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:58:03 INFO Hourly dataset computed and listing created
2026-06-22 17:58:04 INFO Hourly dataset computed
2026-06-22 17:58:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:58:05 INFO Hourly dataset computed and listing created
2026-06-22 17:58:05 INFO Hourly dataset computed
2026-06-22 17:58:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:58:06 INFO Hourly dataset computed and listing created
2026-06-22 17:58:07 INFO Hourly dataset computed
2026-06-22 17:58:07 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:58:08 INFO Hourly dataset computed and listing created
2026-06-22 17:58:09 INFO Hourly dataset computed
2026-06-22 17:58:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:58:10 INFO Hourly dataset computed and listing created
2026-06-22 17:58:11 INFO Hourly dataset computed
2026-06-22 17:58:11 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:58:11 INFO Hourly dataset computed and listing created
2026-06-22 17:58:12 INFO Hourly dataset computed
2026-06-22 17:58:12 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:58:13 INFO Hourly dataset computed and listing created
2026-06-22 17:58:14 INFO Hourly dataset computed
2026-06-22 17:58:14 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:58:15 INFO Hourly dataset computed and listing created
2026-06-22 17:58:16 INFO Hourly dataset computed
2026-06-22 17:58:16 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:58:17 INFO Hourly dataset computed and listing created
2026-06-22 17:58:17 INFO Hourly dataset computed
2026-06-22 17:58:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:58:18 INFO Hourly dataset computed and listing created
2026-06-22 17:58:19 INFO Hourly dataset computed
2026-06-22 17:58:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:58:20 INFO Hourly dataset computed and listing created
2026-06-22 17:58:21 INFO Hourly dataset computed
2026-06-22 17:58:21 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:58:22 INFO Hourly dataset computed and listing created
2026-06-22 17:58:23 INFO Hourly dataset computed
2026-06-22 17:58:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 17:58:23 INFO Hourly dataset computed and listing created
2026-06-22 17:58:24 INFO Hourly dataset computed
2026-06-22 17:58:24 INFO ---------->>> Running CHIMERE model from 2020-02-06 09:00:00 to 2020-02-06 11:00:00
2026-06-22 17:58:24 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:58:24 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1
2026-06-22 17:58:24 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020601_8_ENS1.nc
2026-06-22 17:58:24 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-06-22 17:58:24 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:58:24 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-06-22 17:58:24 INFO Queuing job for member 1...
2026-06-22 17:58:24 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:58:24 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-06-22 17:58:25 INFO Found: ['4951831']
2026-06-22 17:58:30 INFO [TGCC-IRENE] Submitted job with ID:['4951831']
2026-06-22 17:58:30 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:58:30 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2
2026-06-22 17:58:30 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020601_8_ENS2.nc
2026-06-22 17:58:30 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-06-22 17:58:30 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:58:30 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-06-22 17:58:30 INFO Queuing job for member 2...
2026-06-22 17:58:30 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:58:30 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-06-22 17:58:31 INFO Found: ['4951837']
2026-06-22 17:58:36 INFO [TGCC-IRENE] Submitted job with ID:['4951837']
2026-06-22 17:58:36 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:58:36 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3
2026-06-22 17:58:36 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020601_8_ENS3.nc
2026-06-22 17:58:36 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-06-22 17:58:36 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:58:36 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-06-22 17:58:36 INFO Queuing job for member 3...
2026-06-22 17:58:36 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:58:36 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-06-22 17:58:37 INFO Found: ['4951840']
2026-06-22 17:58:42 INFO [TGCC-IRENE] Submitted job with ID:['4951840']
2026-06-22 17:58:42 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:58:42 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4
2026-06-22 17:58:42 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020601_8_ENS4.nc
2026-06-22 17:58:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-06-22 17:58:42 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:58:42 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-06-22 17:58:42 INFO Queuing job for member 4...
2026-06-22 17:58:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:58:42 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-06-22 17:58:42 INFO Found: ['4951841']
2026-06-22 17:58:47 INFO [TGCC-IRENE] Submitted job with ID:['4951841']
2026-06-22 17:58:47 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:58:47 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5
2026-06-22 17:58:47 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020601_8_ENS5.nc
2026-06-22 17:58:47 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-06-22 17:58:47 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:58:47 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-06-22 17:58:47 INFO Queuing job for member 5...
2026-06-22 17:58:47 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:58:47 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-06-22 17:58:48 INFO Found: ['4951843']
2026-06-22 17:58:53 INFO [TGCC-IRENE] Submitted job with ID:['4951843']
2026-06-22 17:58:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:58:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6
2026-06-22 17:58:53 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020601_8_ENS6.nc
2026-06-22 17:58:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-06-22 17:58:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:58:53 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-06-22 17:58:53 INFO Queuing job for member 6...
2026-06-22 17:58:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:58:53 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-06-22 17:58:54 INFO Found: ['4951847']
2026-06-22 17:58:59 INFO [TGCC-IRENE] Submitted job with ID:['4951847']
2026-06-22 17:58:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:58:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7
2026-06-22 17:58:59 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020601_8_ENS7.nc
2026-06-22 17:58:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-06-22 17:58:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:58:59 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-06-22 17:58:59 INFO Queuing job for member 7...
2026-06-22 17:58:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:58:59 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-06-22 17:59:00 INFO Found: ['4951850']
2026-06-22 17:59:05 INFO [TGCC-IRENE] Submitted job with ID:['4951850']
2026-06-22 17:59:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:59:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8
2026-06-22 17:59:05 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020601_8_ENS8.nc
2026-06-22 17:59:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-06-22 17:59:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:59:05 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-06-22 17:59:05 INFO Queuing job for member 8...
2026-06-22 17:59:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:59:05 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-06-22 17:59:05 INFO Found: ['4951854']
2026-06-22 17:59:10 INFO [TGCC-IRENE] Submitted job with ID:['4951854']
2026-06-22 17:59:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:59:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9
2026-06-22 17:59:10 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020601_8_ENS9.nc
2026-06-22 17:59:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-06-22 17:59:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:59:10 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-06-22 17:59:10 INFO Queuing job for member 9...
2026-06-22 17:59:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:59:10 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-06-22 17:59:11 INFO Found: ['4951856']
2026-06-22 17:59:16 INFO [TGCC-IRENE] Submitted job with ID:['4951856']
2026-06-22 17:59:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:59:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10
2026-06-22 17:59:16 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020601_8_ENS10.nc
2026-06-22 17:59:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-06-22 17:59:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:59:16 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-06-22 17:59:16 INFO Queuing job for member 10...
2026-06-22 17:59:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:59:16 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-06-22 17:59:17 INFO Found: ['4951858']
2026-06-22 17:59:22 INFO [TGCC-IRENE] Submitted job with ID:['4951858']
2026-06-22 17:59:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:59:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11
2026-06-22 17:59:22 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020601_8_ENS11.nc
2026-06-22 17:59:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-06-22 17:59:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:59:22 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-06-22 17:59:22 INFO Queuing job for member 11...
2026-06-22 17:59:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:59:22 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-06-22 17:59:22 INFO Found: ['4951859']
2026-06-22 17:59:27 INFO [TGCC-IRENE] Submitted job with ID:['4951859']
2026-06-22 17:59:27 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:59:27 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12
2026-06-22 17:59:27 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020601_8_ENS12.nc
2026-06-22 17:59:27 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-06-22 17:59:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:59:28 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-06-22 17:59:28 INFO Queuing job for member 12...
2026-06-22 17:59:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:59:28 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-06-22 17:59:28 INFO Found: ['4951863']
2026-06-22 17:59:33 INFO [TGCC-IRENE] Submitted job with ID:['4951863']
2026-06-22 17:59:33 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:59:33 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13
2026-06-22 17:59:33 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020601_8_ENS13.nc
2026-06-22 17:59:33 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-06-22 17:59:33 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:59:33 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-06-22 17:59:33 INFO Queuing job for member 13...
2026-06-22 17:59:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:59:33 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-06-22 17:59:34 INFO Found: ['4951866']
2026-06-22 17:59:39 INFO [TGCC-IRENE] Submitted job with ID:['4951866']
2026-06-22 17:59:39 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:59:39 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14
2026-06-22 17:59:39 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020601_8_ENS14.nc
2026-06-22 17:59:39 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-06-22 17:59:39 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:59:39 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-06-22 17:59:39 INFO Queuing job for member 14...
2026-06-22 17:59:39 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:59:39 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-06-22 17:59:40 INFO Found: ['4951869']
2026-06-22 17:59:45 INFO [TGCC-IRENE] Submitted job with ID:['4951869']
2026-06-22 17:59:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:59:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15
2026-06-22 17:59:45 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020601_8_ENS15.nc
2026-06-22 17:59:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-06-22 17:59:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:59:45 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-06-22 17:59:45 INFO Queuing job for member 15...
2026-06-22 17:59:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:59:45 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-06-22 17:59:45 INFO Found: ['4951870']
2026-06-22 17:59:50 INFO [TGCC-IRENE] Submitted job with ID:['4951870']
2026-06-22 17:59:50 INFO Checking job status ...
2026-06-22 17:59:50 INFO None 4951831: status RUNNING/PENDING
2026-06-22 17:59:50 INFO None 4951837: status RUNNING/PENDING
2026-06-22 17:59:50 INFO None 4951840: status RUNNING/PENDING
2026-06-22 17:59:50 INFO None 4951841: status RUNNING/PENDING
2026-06-22 17:59:51 INFO None 4951843: status RUNNING/PENDING
2026-06-22 17:59:51 INFO None 4951847: status RUNNING/PENDING
2026-06-22 17:59:51 INFO None 4951850: status RUNNING/PENDING
2026-06-22 17:59:51 INFO None 4951854: status RUNNING/PENDING
2026-06-22 17:59:51 INFO None 4951856: status RUNNING/PENDING
2026-06-22 17:59:51 INFO None 4951858: status RUNNING/PENDING
2026-06-22 17:59:51 INFO None 4951859: status RUNNING/PENDING
2026-06-22 17:59:51 INFO None 4951863: status RUNNING/PENDING
2026-06-22 17:59:51 INFO None 4951866: status RUNNING/PENDING
2026-06-22 17:59:51 INFO None 4951869: status RUNNING/PENDING
2026-06-22 17:59:51 INFO None 4951870: status RUNNING/PENDING
2026-06-22 17:59:51 INFO Jobs still running: ['4951831', '4951837', '4951840', '4951841', '4951843', '4951847', '4951850', '4951854', '4951856', '4951858', '4951859', '4951863', '4951866', '4951869', '4951870']. Waiting...
2026-06-22 18:00:06 INFO None 4951831: status RUNNING/PENDING
2026-06-22 18:00:06 INFO None 4951837: status RUNNING/PENDING
2026-06-22 18:00:06 INFO None 4951840: status RUNNING/PENDING
2026-06-22 18:00:06 INFO None 4951841: status RUNNING/PENDING
2026-06-22 18:00:06 INFO None 4951843: status RUNNING/PENDING
2026-06-22 18:00:06 INFO None 4951847: status RUNNING/PENDING
2026-06-22 18:00:06 INFO None 4951850: status RUNNING/PENDING
2026-06-22 18:00:06 INFO None 4951854: status RUNNING/PENDING
2026-06-22 18:00:06 INFO None 4951856: status RUNNING/PENDING
2026-06-22 18:00:06 INFO None 4951858: status RUNNING/PENDING
2026-06-22 18:00:06 INFO None 4951859: status RUNNING/PENDING
2026-06-22 18:00:06 INFO None 4951863: status RUNNING/PENDING
2026-06-22 18:00:06 INFO None 4951866: status RUNNING/PENDING
2026-06-22 18:00:06 INFO None 4951869: status RUNNING/PENDING
2026-06-22 18:00:06 INFO None 4951870: status RUNNING/PENDING
2026-06-22 18:00:06 INFO Jobs still running: ['4951831', '4951837', '4951840', '4951841', '4951843', '4951847', '4951850', '4951854', '4951856', '4951858', '4951859', '4951863', '4951866', '4951869', '4951870']. Waiting...
2026-06-22 18:00:21 INFO None 4951831: status RUNNING/PENDING
2026-06-22 18:00:21 INFO None 4951837: status RUNNING/PENDING
2026-06-22 18:00:21 INFO None 4951840: status RUNNING/PENDING
2026-06-22 18:00:21 INFO None 4951841: status RUNNING/PENDING
2026-06-22 18:00:21 INFO None 4951843: status RUNNING/PENDING
2026-06-22 18:00:21 INFO None 4951847: status RUNNING/PENDING
2026-06-22 18:00:21 INFO None 4951850: status RUNNING/PENDING
2026-06-22 18:00:21 INFO None 4951854: status RUNNING/PENDING
2026-06-22 18:00:21 INFO None 4951856: status RUNNING/PENDING
2026-06-22 18:00:21 INFO None 4951858: status RUNNING/PENDING
2026-06-22 18:00:21 INFO None 4951859: status RUNNING/PENDING
2026-06-22 18:00:21 INFO None 4951863: status RUNNING/PENDING
2026-06-22 18:00:21 INFO None 4951866: status RUNNING/PENDING
2026-06-22 18:00:21 INFO None 4951869: status RUNNING/PENDING
2026-06-22 18:00:21 INFO None 4951870: status RUNNING/PENDING
2026-06-22 18:00:21 INFO Jobs still running: ['4951831', '4951837', '4951840', '4951841', '4951843', '4951847', '4951850', '4951854', '4951856', '4951858', '4951859', '4951863', '4951866', '4951869', '4951870']. Waiting...
2026-06-22 18:00:36 INFO None 4951831: status RUNNING/PENDING
2026-06-22 18:00:36 INFO None 4951837: status RUNNING/PENDING
2026-06-22 18:00:36 INFO None 4951840: status RUNNING/PENDING
2026-06-22 18:00:36 INFO None 4951841: status RUNNING/PENDING
2026-06-22 18:00:36 INFO None 4951843: status RUNNING/PENDING
2026-06-22 18:00:36 INFO None 4951847: status RUNNING/PENDING
2026-06-22 18:00:36 INFO None 4951850: status RUNNING/PENDING
2026-06-22 18:00:36 INFO None 4951854: status RUNNING/PENDING
2026-06-22 18:00:36 INFO None 4951856: status RUNNING/PENDING
2026-06-22 18:00:37 INFO None 4951858: status RUNNING/PENDING
2026-06-22 18:00:37 INFO None 4951859: status RUNNING/PENDING
2026-06-22 18:00:37 INFO None 4951863: status RUNNING/PENDING
2026-06-22 18:00:37 INFO None 4951866: status RUNNING/PENDING
2026-06-22 18:00:37 INFO None 4951869: status RUNNING/PENDING
2026-06-22 18:00:37 INFO None 4951870: status RUNNING/PENDING
2026-06-22 18:00:37 INFO Jobs still running: ['4951831', '4951837', '4951840', '4951841', '4951843', '4951847', '4951850', '4951854', '4951856', '4951858', '4951859', '4951863', '4951866', '4951869', '4951870']. Waiting...
2026-06-22 18:00:52 INFO None 4951831: status RUNNING/PENDING
2026-06-22 18:00:52 INFO None 4951837: status RUNNING/PENDING
2026-06-22 18:00:52 INFO None 4951840: status RUNNING/PENDING
2026-06-22 18:00:52 INFO None 4951841: status RUNNING/PENDING
2026-06-22 18:00:52 INFO None 4951843: status RUNNING/PENDING
2026-06-22 18:00:52 INFO None 4951847: status RUNNING/PENDING
2026-06-22 18:00:52 INFO None 4951850: status RUNNING/PENDING
2026-06-22 18:00:52 INFO None 4951854: status RUNNING/PENDING
2026-06-22 18:00:52 INFO None 4951856: status RUNNING/PENDING
2026-06-22 18:00:52 INFO None 4951858: status RUNNING/PENDING
2026-06-22 18:00:52 INFO None 4951859: status RUNNING/PENDING
2026-06-22 18:00:52 INFO None 4951863: status RUNNING/PENDING
2026-06-22 18:00:52 INFO None 4951866: status RUNNING/PENDING
2026-06-22 18:00:52 INFO None 4951869: status RUNNING/PENDING
2026-06-22 18:00:52 INFO None 4951870: status RUNNING/PENDING
2026-06-22 18:00:52 INFO Jobs still running: ['4951831', '4951837', '4951840', '4951841', '4951843', '4951847', '4951850', '4951854', '4951856', '4951858', '4951859', '4951863', '4951866', '4951869', '4951870']. Waiting...
2026-06-22 18:01:07 INFO None 4951831: status RUNNING/PENDING
2026-06-22 18:01:07 INFO None 4951837: status RUNNING/PENDING
2026-06-22 18:01:07 INFO None 4951840: status RUNNING/PENDING
2026-06-22 18:01:07 INFO None 4951841: status RUNNING/PENDING
2026-06-22 18:01:07 INFO None 4951843: status RUNNING/PENDING
2026-06-22 18:01:07 INFO None 4951847: status RUNNING/PENDING
2026-06-22 18:01:07 INFO None 4951850: status RUNNING/PENDING
2026-06-22 18:01:07 INFO None 4951854: status RUNNING/PENDING
2026-06-22 18:01:07 INFO None 4951856: status RUNNING/PENDING
2026-06-22 18:01:07 INFO None 4951858: status RUNNING/PENDING
2026-06-22 18:01:07 INFO None 4951859: status RUNNING/PENDING
2026-06-22 18:01:07 INFO None 4951863: status RUNNING/PENDING
2026-06-22 18:01:07 INFO None 4951866: status RUNNING/PENDING
2026-06-22 18:01:07 INFO None 4951869: status RUNNING/PENDING
2026-06-22 18:01:07 INFO None 4951870: status RUNNING/PENDING
2026-06-22 18:01:07 INFO Jobs still running: ['4951831', '4951837', '4951840', '4951841', '4951843', '4951847', '4951850', '4951854', '4951856', '4951858', '4951859', '4951863', '4951866', '4951869', '4951870']. Waiting...
2026-06-22 18:01:22 INFO None 4951831: status RUNNING/PENDING
2026-06-22 18:01:22 INFO None 4951837: status RUNNING/PENDING
2026-06-22 18:01:22 INFO None 4951840: status RUNNING/PENDING
2026-06-22 18:01:22 INFO None 4951841: status RUNNING/PENDING
2026-06-22 18:01:22 INFO None 4951843: status RUNNING/PENDING
2026-06-22 18:01:22 INFO None 4951847: status RUNNING/PENDING
2026-06-22 18:01:22 INFO None 4951850: status RUNNING/PENDING
2026-06-22 18:01:22 INFO None 4951854: status RUNNING/PENDING
2026-06-22 18:01:22 INFO None 4951856: status RUNNING/PENDING
2026-06-22 18:01:22 INFO None 4951858: status RUNNING/PENDING
2026-06-22 18:01:22 INFO None 4951859: status RUNNING/PENDING
2026-06-22 18:01:22 INFO None 4951863: status RUNNING/PENDING
2026-06-22 18:01:22 INFO None 4951866: status RUNNING/PENDING
2026-06-22 18:01:23 INFO None 4951869: status RUNNING/PENDING
2026-06-22 18:01:23 INFO None 4951870: status RUNNING/PENDING
2026-06-22 18:01:23 INFO Jobs still running: ['4951831', '4951837', '4951840', '4951841', '4951843', '4951847', '4951850', '4951854', '4951856', '4951858', '4951859', '4951863', '4951866', '4951869', '4951870']. Waiting...
2026-06-22 18:01:38 INFO None 4951831: status RUNNING/PENDING
2026-06-22 18:01:38 INFO None 4951837: status RUNNING/PENDING
2026-06-22 18:01:38 INFO None 4951840: status RUNNING/PENDING
2026-06-22 18:01:38 INFO None 4951841: status RUNNING/PENDING
2026-06-22 18:01:38 INFO None 4951843: status RUNNING/PENDING
2026-06-22 18:01:38 INFO None 4951847: status RUNNING/PENDING
2026-06-22 18:01:38 INFO None 4951850: status RUNNING/PENDING
2026-06-22 18:01:38 INFO None 4951854: status RUNNING/PENDING
2026-06-22 18:01:38 INFO None 4951856: status RUNNING/PENDING
2026-06-22 18:01:38 INFO None 4951858: status RUNNING/PENDING
2026-06-22 18:01:38 INFO None 4951859: status RUNNING/PENDING
2026-06-22 18:01:38 INFO None 4951863: status RUNNING/PENDING
2026-06-22 18:01:38 INFO None 4951866: status RUNNING/PENDING
2026-06-22 18:01:38 INFO None 4951869: status RUNNING/PENDING
2026-06-22 18:01:38 INFO None 4951870: status RUNNING/PENDING
2026-06-22 18:01:38 INFO Jobs still running: ['4951831', '4951837', '4951840', '4951841', '4951843', '4951847', '4951850', '4951854', '4951856', '4951858', '4951859', '4951863', '4951866', '4951869', '4951870']. Waiting...
2026-06-22 18:01:53 INFO None 4951831: status RUNNING/PENDING
2026-06-22 18:01:53 INFO None 4951837: status RUNNING/PENDING
2026-06-22 18:01:53 INFO None 4951840: status RUNNING/PENDING
2026-06-22 18:01:53 INFO None 4951841: status RUNNING/PENDING
2026-06-22 18:01:53 INFO None 4951843: status RUNNING/PENDING
2026-06-22 18:01:53 INFO None 4951847: status RUNNING/PENDING
2026-06-22 18:01:53 INFO None 4951850: status RUNNING/PENDING
2026-06-22 18:01:53 INFO None 4951854: status RUNNING/PENDING
2026-06-22 18:01:53 INFO None 4951856: status RUNNING/PENDING
2026-06-22 18:01:53 INFO None 4951858: status RUNNING/PENDING
2026-06-22 18:01:53 INFO None 4951859: status RUNNING/PENDING
2026-06-22 18:01:53 INFO None 4951863: status RUNNING/PENDING
2026-06-22 18:01:53 INFO None 4951866: status RUNNING/PENDING
2026-06-22 18:01:53 INFO None 4951869: status RUNNING/PENDING
2026-06-22 18:01:53 INFO None 4951870: status RUNNING/PENDING
2026-06-22 18:01:53 INFO Jobs still running: ['4951831', '4951837', '4951840', '4951841', '4951843', '4951847', '4951850', '4951854', '4951856', '4951858', '4951859', '4951863', '4951866', '4951869', '4951870']. Waiting...
2026-06-22 18:02:08 INFO None 4951831: status RUNNING/PENDING
2026-06-22 18:02:08 INFO None 4951837: status RUNNING/PENDING
2026-06-22 18:02:08 INFO None 4951840: status RUNNING/PENDING
2026-06-22 18:02:08 INFO None 4951841: status RUNNING/PENDING
2026-06-22 18:02:08 INFO None 4951843: status RUNNING/PENDING
2026-06-22 18:02:08 INFO None 4951847: status RUNNING/PENDING
2026-06-22 18:02:08 INFO None 4951850: status RUNNING/PENDING
2026-06-22 18:02:08 INFO None 4951854: status RUNNING/PENDING
2026-06-22 18:02:08 INFO None 4951856: status RUNNING/PENDING
2026-06-22 18:02:08 INFO None 4951858: status RUNNING/PENDING
2026-06-22 18:02:08 INFO None 4951859: status RUNNING/PENDING
2026-06-22 18:02:08 INFO None 4951863: status RUNNING/PENDING
2026-06-22 18:02:08 INFO None 4951866: status RUNNING/PENDING
2026-06-22 18:02:08 INFO None 4951869: status RUNNING/PENDING
2026-06-22 18:02:08 INFO None 4951870: status RUNNING/PENDING
2026-06-22 18:02:08 INFO Jobs still running: ['4951831', '4951837', '4951840', '4951841', '4951843', '4951847', '4951850', '4951854', '4951856', '4951858', '4951859', '4951863', '4951866', '4951869', '4951870']. Waiting...
2026-06-22 18:02:23 INFO None 4951831: status RUNNING/PENDING
2026-06-22 18:02:23 INFO None 4951837: status RUNNING/PENDING
2026-06-22 18:02:24 INFO None 4951840: status RUNNING/PENDING
2026-06-22 18:02:24 INFO None 4951841: status RUNNING/PENDING
2026-06-22 18:02:24 INFO None 4951843: status RUNNING/PENDING
2026-06-22 18:02:24 INFO None 4951847: status FINISHED
2026-06-22 18:02:24 INFO None 4951850: status FINISHED
2026-06-22 18:02:24 INFO None 4951854: status RUNNING/PENDING
2026-06-22 18:02:24 INFO None 4951856: status RUNNING/PENDING
2026-06-22 18:02:24 INFO None 4951858: status RUNNING/PENDING
2026-06-22 18:02:24 INFO None 4951859: status RUNNING/PENDING
2026-06-22 18:02:24 INFO None 4951863: status RUNNING/PENDING
2026-06-22 18:02:24 INFO None 4951866: status RUNNING/PENDING
2026-06-22 18:02:24 INFO None 4951869: status RUNNING/PENDING
2026-06-22 18:02:24 INFO None 4951870: status RUNNING/PENDING
2026-06-22 18:02:24 INFO Jobs still running: ['4951831', '4951837', '4951840', '4951841', '4951843', '4951854', '4951856', '4951858', '4951859', '4951863', '4951866', '4951869', '4951870']. Waiting...
2026-06-22 18:02:39 INFO None 4951831: status RUNNING/PENDING
2026-06-22 18:02:39 INFO None 4951837: status RUNNING/PENDING
2026-06-22 18:02:39 INFO None 4951840: status RUNNING/PENDING
2026-06-22 18:02:39 INFO None 4951841: status RUNNING/PENDING
2026-06-22 18:02:39 INFO None 4951843: status RUNNING/PENDING
2026-06-22 18:02:39 INFO None 4951847: status FINISHED
2026-06-22 18:02:39 INFO None 4951850: status FINISHED
2026-06-22 18:02:39 INFO None 4951854: status RUNNING/PENDING
2026-06-22 18:02:39 INFO None 4951856: status RUNNING/PENDING
2026-06-22 18:02:39 INFO None 4951858: status RUNNING/PENDING
2026-06-22 18:02:39 INFO None 4951859: status RUNNING/PENDING
2026-06-22 18:02:39 INFO None 4951863: status RUNNING/PENDING
2026-06-22 18:02:39 INFO None 4951866: status RUNNING/PENDING
2026-06-22 18:02:39 INFO None 4951869: status RUNNING/PENDING
2026-06-22 18:02:39 INFO None 4951870: status RUNNING/PENDING
2026-06-22 18:02:39 INFO Jobs still running: ['4951831', '4951837', '4951840', '4951841', '4951843', '4951854', '4951856', '4951858', '4951859', '4951863', '4951866', '4951869', '4951870']. Waiting...
2026-06-22 18:02:54 INFO None 4951831: status RUNNING/PENDING
2026-06-22 18:02:54 INFO None 4951837: status RUNNING/PENDING
2026-06-22 18:02:54 INFO None 4951840: status RUNNING/PENDING
2026-06-22 18:02:54 INFO None 4951841: status RUNNING/PENDING
2026-06-22 18:02:54 INFO None 4951843: status RUNNING/PENDING
2026-06-22 18:02:54 INFO None 4951847: status FINISHED
2026-06-22 18:02:54 INFO None 4951850: status FINISHED
2026-06-22 18:02:54 INFO None 4951854: status FINISHED
2026-06-22 18:02:54 INFO None 4951856: status FINISHED
2026-06-22 18:02:54 INFO None 4951858: status RUNNING/PENDING
2026-06-22 18:02:54 INFO None 4951859: status FINISHED
2026-06-22 18:02:54 INFO None 4951863: status RUNNING/PENDING
2026-06-22 18:02:54 INFO None 4951866: status RUNNING/PENDING
2026-06-22 18:02:54 INFO None 4951869: status RUNNING/PENDING
2026-06-22 18:02:54 INFO None 4951870: status RUNNING/PENDING
2026-06-22 18:02:54 INFO Jobs still running: ['4951831', '4951837', '4951840', '4951841', '4951843', '4951858', '4951863', '4951866', '4951869', '4951870']. Waiting...
2026-06-22 18:03:09 INFO None 4951831: status RUNNING/PENDING
2026-06-22 18:03:09 INFO None 4951837: status RUNNING/PENDING
2026-06-22 18:03:09 INFO None 4951840: status FINISHED
2026-06-22 18:03:10 INFO None 4951841: status RUNNING/PENDING
2026-06-22 18:03:10 INFO None 4951843: status RUNNING/PENDING
2026-06-22 18:03:10 INFO None 4951847: status FINISHED
2026-06-22 18:03:10 INFO None 4951850: status FINISHED
2026-06-22 18:03:10 INFO None 4951854: status FINISHED
2026-06-22 18:03:10 INFO None 4951856: status FINISHED
2026-06-22 18:03:10 INFO None 4951858: status FINISHED
2026-06-22 18:03:10 INFO None 4951859: status FINISHED
2026-06-22 18:03:10 INFO None 4951863: status RUNNING/PENDING
2026-06-22 18:03:10 INFO None 4951866: status RUNNING/PENDING
2026-06-22 18:03:10 INFO None 4951869: status RUNNING/PENDING
2026-06-22 18:03:10 INFO None 4951870: status RUNNING/PENDING
2026-06-22 18:03:10 INFO Jobs still running: ['4951831', '4951837', '4951841', '4951843', '4951863', '4951866', '4951869', '4951870']. Waiting...
2026-06-22 18:03:25 INFO None 4951831: status RUNNING/PENDING
2026-06-22 18:03:25 INFO None 4951837: status RUNNING/PENDING
2026-06-22 18:03:25 INFO None 4951840: status FINISHED
2026-06-22 18:03:25 INFO None 4951841: status RUNNING/PENDING
2026-06-22 18:03:25 INFO None 4951843: status RUNNING/PENDING
2026-06-22 18:03:25 INFO None 4951847: status FINISHED
2026-06-22 18:03:25 INFO None 4951850: status FINISHED
2026-06-22 18:03:25 INFO None 4951854: status FINISHED
2026-06-22 18:03:25 INFO None 4951856: status FINISHED
2026-06-22 18:03:25 INFO None 4951858: status FINISHED
2026-06-22 18:03:25 INFO None 4951859: status FINISHED
2026-06-22 18:03:25 INFO None 4951863: status RUNNING/PENDING
2026-06-22 18:03:25 INFO None 4951866: status RUNNING/PENDING
2026-06-22 18:03:25 INFO None 4951869: status RUNNING/PENDING
2026-06-22 18:03:25 INFO None 4951870: status RUNNING/PENDING
2026-06-22 18:03:25 INFO Jobs still running: ['4951831', '4951837', '4951841', '4951843', '4951863', '4951866', '4951869', '4951870']. Waiting...
2026-06-22 18:03:41 INFO None 4951831: status RUNNING/PENDING
2026-06-22 18:03:41 INFO None 4951837: status RUNNING/PENDING
2026-06-22 18:03:41 INFO None 4951840: status FINISHED
2026-06-22 18:03:41 INFO None 4951841: status RUNNING/PENDING
2026-06-22 18:03:41 INFO None 4951843: status RUNNING/PENDING
2026-06-22 18:03:41 INFO None 4951847: status FINISHED
2026-06-22 18:03:41 INFO None 4951850: status FINISHED
2026-06-22 18:03:41 INFO None 4951854: status FINISHED
2026-06-22 18:03:41 INFO None 4951856: status FINISHED
2026-06-22 18:03:41 INFO None 4951858: status FINISHED
2026-06-22 18:03:41 INFO None 4951859: status FINISHED
2026-06-22 18:03:41 INFO None 4951863: status RUNNING/PENDING
2026-06-22 18:03:41 INFO None 4951866: status RUNNING/PENDING
2026-06-22 18:03:41 INFO None 4951869: status RUNNING/PENDING
2026-06-22 18:03:41 INFO None 4951870: status RUNNING/PENDING
2026-06-22 18:03:41 INFO Jobs still running: ['4951831', '4951837', '4951841', '4951843', '4951863', '4951866', '4951869', '4951870']. Waiting...
2026-06-22 18:03:56 INFO None 4951831: status RUNNING/PENDING
2026-06-22 18:03:56 INFO None 4951837: status RUNNING/PENDING
2026-06-22 18:03:56 INFO None 4951840: status FINISHED
2026-06-22 18:03:56 INFO None 4951841: status RUNNING/PENDING
2026-06-22 18:03:56 INFO None 4951843: status RUNNING/PENDING
2026-06-22 18:03:56 INFO None 4951847: status FINISHED
2026-06-22 18:03:56 INFO None 4951850: status FINISHED
2026-06-22 18:03:56 INFO None 4951854: status FINISHED
2026-06-22 18:03:56 INFO None 4951856: status FINISHED
2026-06-22 18:03:56 INFO None 4951858: status FINISHED
2026-06-22 18:03:56 INFO None 4951859: status FINISHED
2026-06-22 18:03:56 INFO None 4951863: status FINISHED
2026-06-22 18:03:56 INFO None 4951866: status FINISHED
2026-06-22 18:03:56 INFO None 4951869: status FINISHED
2026-06-22 18:03:56 INFO None 4951870: status FINISHED
2026-06-22 18:03:56 INFO Jobs still running: ['4951831', '4951837', '4951841', '4951843']. Waiting...
2026-06-22 18:04:11 INFO None 4951831: status FINISHED
2026-06-22 18:04:11 INFO None 4951837: status FINISHED
2026-06-22 18:04:11 INFO None 4951840: status FINISHED
2026-06-22 18:04:11 INFO None 4951841: status RUNNING/PENDING
2026-06-22 18:04:11 INFO None 4951843: status RUNNING/PENDING
2026-06-22 18:04:11 INFO None 4951847: status FINISHED
2026-06-22 18:04:11 INFO None 4951850: status FINISHED
2026-06-22 18:04:11 INFO None 4951854: status FINISHED
2026-06-22 18:04:11 INFO None 4951856: status FINISHED
2026-06-22 18:04:11 INFO None 4951858: status FINISHED
2026-06-22 18:04:11 INFO None 4951859: status FINISHED
2026-06-22 18:04:11 INFO None 4951863: status FINISHED
2026-06-22 18:04:11 INFO None 4951866: status FINISHED
2026-06-22 18:04:11 INFO None 4951869: status FINISHED
2026-06-22 18:04:11 INFO None 4951870: status FINISHED
2026-06-22 18:04:11 INFO Jobs still running: ['4951841', '4951843']. Waiting...
2026-06-22 18:04:26 INFO None 4951831: status FINISHED
2026-06-22 18:04:26 INFO None 4951837: status FINISHED
2026-06-22 18:04:26 INFO None 4951840: status FINISHED
2026-06-22 18:04:26 INFO None 4951841: status RUNNING/PENDING
2026-06-22 18:04:26 INFO None 4951843: status FINISHED
2026-06-22 18:04:26 INFO None 4951847: status FINISHED
2026-06-22 18:04:27 INFO None 4951850: status FINISHED
2026-06-22 18:04:27 INFO None 4951854: status FINISHED
2026-06-22 18:04:27 INFO None 4951856: status FINISHED
2026-06-22 18:04:27 INFO None 4951858: status FINISHED
2026-06-22 18:04:27 INFO None 4951859: status FINISHED
2026-06-22 18:04:27 INFO None 4951863: status FINISHED
2026-06-22 18:04:27 INFO None 4951866: status FINISHED
2026-06-22 18:04:27 INFO None 4951869: status FINISHED
2026-06-22 18:04:27 INFO None 4951870: status FINISHED
2026-06-22 18:04:27 INFO Jobs still running: ['4951841']. Waiting...
2026-06-22 18:04:42 INFO None 4951831: status FINISHED
2026-06-22 18:04:42 INFO None 4951837: status FINISHED
2026-06-22 18:04:42 INFO None 4951840: status FINISHED
2026-06-22 18:04:42 INFO None 4951841: status RUNNING/PENDING
2026-06-22 18:04:42 INFO None 4951843: status FINISHED
2026-06-22 18:04:42 INFO None 4951847: status FINISHED
2026-06-22 18:04:42 INFO None 4951850: status FINISHED
2026-06-22 18:04:42 INFO None 4951854: status FINISHED
2026-06-22 18:04:42 INFO None 4951856: status FINISHED
2026-06-22 18:04:42 INFO None 4951858: status FINISHED
2026-06-22 18:04:42 INFO None 4951859: status FINISHED
2026-06-22 18:04:42 INFO None 4951863: status FINISHED
2026-06-22 18:04:42 INFO None 4951866: status FINISHED
2026-06-22 18:04:42 INFO None 4951869: status FINISHED
2026-06-22 18:04:42 INFO None 4951870: status FINISHED
2026-06-22 18:04:42 INFO Jobs still running: ['4951841']. Waiting...
2026-06-22 18:04:57 INFO None 4951831: status FINISHED
2026-06-22 18:04:57 INFO None 4951837: status FINISHED
2026-06-22 18:04:57 INFO None 4951840: status FINISHED
2026-06-22 18:04:57 INFO None 4951841: status RUNNING/PENDING
2026-06-22 18:04:57 INFO None 4951843: status FINISHED
2026-06-22 18:04:57 INFO None 4951847: status FINISHED
2026-06-22 18:04:57 INFO None 4951850: status FINISHED
2026-06-22 18:04:57 INFO None 4951854: status FINISHED
2026-06-22 18:04:57 INFO None 4951856: status FINISHED
2026-06-22 18:04:57 INFO None 4951858: status FINISHED
2026-06-22 18:04:57 INFO None 4951859: status FINISHED
2026-06-22 18:04:57 INFO None 4951863: status FINISHED
2026-06-22 18:04:57 INFO None 4951866: status FINISHED
2026-06-22 18:04:57 INFO None 4951869: status FINISHED
2026-06-22 18:04:57 INFO None 4951870: status FINISHED
2026-06-22 18:04:57 INFO Jobs still running: ['4951841']. Waiting...
2026-06-22 18:05:12 INFO None 4951831: status FINISHED
2026-06-22 18:05:12 INFO None 4951837: status FINISHED
2026-06-22 18:05:12 INFO None 4951840: status FINISHED
2026-06-22 18:05:12 INFO None 4951841: status FINISHED
2026-06-22 18:05:12 INFO None 4951843: status FINISHED
2026-06-22 18:05:12 INFO None 4951847: status FINISHED
2026-06-22 18:05:12 INFO None 4951850: status FINISHED
2026-06-22 18:05:12 INFO None 4951854: status FINISHED
2026-06-22 18:05:12 INFO None 4951856: status FINISHED
2026-06-22 18:05:13 INFO None 4951858: status FINISHED
2026-06-22 18:05:13 INFO None 4951859: status FINISHED
2026-06-22 18:05:13 INFO None 4951863: status FINISHED
2026-06-22 18:05:13 INFO None 4951866: status FINISHED
2026-06-22 18:05:13 INFO None 4951869: status FINISHED
2026-06-22 18:05:13 INFO None 4951870: status FINISHED
2026-06-22 18:05:13 INFO Jobs ['4951831', '4951837', '4951840', '4951841', '4951843', '4951847', '4951850', '4951854', '4951856', '4951858', '4951859', '4951863', '4951866', '4951869', '4951870'] have finished
2026-06-22 18:05:13 INFO Checking restart files were created ...
2026-06-22 18:05:13 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc(1002685915 bytes)
2026-06-22 18:05:13 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc(1002685915 bytes)
2026-06-22 18:05:13 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc(1002685915 bytes)
2026-06-22 18:05:13 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc(1002685915 bytes)
2026-06-22 18:05:13 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc(1002685915 bytes)
2026-06-22 18:05:13 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc(1002685915 bytes)
2026-06-22 18:05:13 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc(1002685915 bytes)
2026-06-22 18:05:13 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc(1002685915 bytes)
2026-06-22 18:05:13 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc(1002685915 bytes)
2026-06-22 18:05:13 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc(1002685915 bytes)
2026-06-22 18:05:13 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc(1002685915 bytes)
2026-06-22 18:05:13 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc(1002685915 bytes)
2026-06-22 18:05:13 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc(1002685915 bytes)
2026-06-22 18:05:13 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc(1002685915 bytes)
2026-06-22 18:05:13 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc(1002685915 bytes)
2026-06-22 18:05:13 INFO  Run_model() completed successfully.
2026-06-22 18:05:13 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-22 18:05:13 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 11:00:00 days=153072 seconds=39600
2026-06-22 18:05:13 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-06-22 18:05:13 INFO [TIME] increment current_time 2020-02-06 09:00:00 -> 2020-02-06 11:00:00
2026-06-22 18:05:13 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-22 18:05:13 INFO ---------->>> Running process_satellite_data()
2026-06-22 18:05:13 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12003.nc
2026-06-22 18:05:13 INFO ---------->>> Running run_obs_converter()
2026-06-22 18:05:13 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_39045_153072.out
2026-06-22 18:05:13 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_39045_153072.out
2026-06-22 18:05:13 INFO ---------->>> Running DART
2026-06-22 18:05:13 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-06-22 18:05:13 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/chim_ENS1_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/chim_ENS1_2020020611_1_out_toDART.nc
2026-06-22 18:05:13 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/chim_ENS2_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/chim_ENS2_2020020611_1_out_toDART.nc
2026-06-22 18:05:13 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/chim_ENS3_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/chim_ENS3_2020020611_1_out_toDART.nc
2026-06-22 18:05:14 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/chim_ENS4_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/chim_ENS4_2020020611_1_out_toDART.nc
2026-06-22 18:05:14 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/chim_ENS5_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/chim_ENS5_2020020611_1_out_toDART.nc
2026-06-22 18:05:14 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/chim_ENS6_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/chim_ENS6_2020020611_1_out_toDART.nc
2026-06-22 18:05:15 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/chim_ENS7_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/chim_ENS7_2020020611_1_out_toDART.nc
2026-06-22 18:05:15 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/chim_ENS8_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/chim_ENS8_2020020611_1_out_toDART.nc
2026-06-22 18:05:15 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/chim_ENS9_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/chim_ENS9_2020020611_1_out_toDART.nc
2026-06-22 18:05:16 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/chim_ENS10_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/chim_ENS10_2020020611_1_out_toDART.nc
2026-06-22 18:05:16 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/chim_ENS11_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/chim_ENS11_2020020611_1_out_toDART.nc
2026-06-22 18:05:16 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/chim_ENS12_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/chim_ENS12_2020020611_1_out_toDART.nc
2026-06-22 18:05:17 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/chim_ENS13_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/chim_ENS13_2020020611_1_out_toDART.nc
2026-06-22 18:05:17 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/chim_ENS14_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/chim_ENS14_2020020611_1_out_toDART.nc
2026-06-22 18:05:17 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/chim_ENS15_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/chim_ENS15_2020020611_1_out_toDART.nc
2026-06-22 18:05:18 INFO Replacement input_template.nml → input.nml completed successfully.
2026-06-22 18:05:18 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-06-22 18:05:18 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-06-22 18:05:18 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-06-22 18:05:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-06-22 18:05:18 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-06-22 18:05:31 INFO Found: []
2026-06-22 18:05:31 INFO No job id returned by command ./run_filter.bsh
2026-06-22 18:05:31 INFO No monitoring will be performed
2026-06-22 18:05:31 INFO Moving DART output files to analysis and preassim directories for date 2020020611 if present ...
2026-06-22 18:05:31 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-06-22 18:05:31 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-06-22 18:05:31 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-06-22 18:05:31 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-06-22 18:05:31 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-06-22 18:05:31 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-06-22 18:05:31 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-06-22 18:05:31 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-06-22 18:05:31 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-06-22 18:05:31 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-06-22 18:05:31 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-06-22 18:05:31 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-06-22 18:05:31 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-06-22 18:05:31 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-06-22 18:05:31 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-06-22 18:05:31 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-06-22 18:05:31 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-06-22 18:05:31 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-06-22 18:05:31 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-06-22 18:05:31 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-06-22 18:05:31 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-06-22 18:05:31 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-06-22 18:05:31 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-06-22 18:05:31 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-06-22 18:05:31 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-06-22 18:05:31 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-06-22 18:05:31 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-06-22 18:05:31 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-06-22 18:05:31 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-06-22 18:05:31 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-06-22 18:05:31 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-06-22 18:05:31 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020611'
2026-06-22 18:05:31 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-06-22 18:05:31 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020611'
2026-06-22 18:05:31 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-06-22 18:05:31 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-06-22 18:05:31 INFO run_dart() is DONE.
2026-06-22 18:05:31 INFO ---------->>> Running update_pollutant_in_end()
2026-06-22 18:05:31 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:05:31 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:05:32 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:05:32 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:05:33 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc
2026-06-22 18:05:38 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:05:39 INFO [Posterior Diff ENS1] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS1_2020020611.nc
2026-06-22 18:05:39 INFO [Posterior Diff ENS1] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS1_2020020611.relative.nc
2026-06-22 18:05:39 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:05:39 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:05:40 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:05:40 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:05:41 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc
2026-06-22 18:05:46 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:05:46 INFO [Posterior Diff ENS2] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS2_2020020611.nc
2026-06-22 18:05:46 INFO [Posterior Diff ENS2] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS2_2020020611.relative.nc
2026-06-22 18:05:47 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:05:47 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:05:48 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:05:48 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:05:48 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc
2026-06-22 18:05:53 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:05:54 INFO [Posterior Diff ENS3] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS3_2020020611.nc
2026-06-22 18:05:54 INFO [Posterior Diff ENS3] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS3_2020020611.relative.nc
2026-06-22 18:05:54 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:05:54 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:05:55 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:05:55 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:05:56 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc
2026-06-22 18:06:01 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:06:02 INFO [Posterior Diff ENS4] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS4_2020020611.nc
2026-06-22 18:06:02 INFO [Posterior Diff ENS4] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS4_2020020611.relative.nc
2026-06-22 18:06:02 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:06:02 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:06:03 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:06:03 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:06:04 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc
2026-06-22 18:06:09 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:06:10 INFO [Posterior Diff ENS5] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS5_2020020611.nc
2026-06-22 18:06:10 INFO [Posterior Diff ENS5] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS5_2020020611.relative.nc
2026-06-22 18:06:10 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:06:10 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:06:11 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:06:11 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:06:12 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc
2026-06-22 18:06:17 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:06:17 INFO [Posterior Diff ENS6] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS6_2020020611.nc
2026-06-22 18:06:17 INFO [Posterior Diff ENS6] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS6_2020020611.relative.nc
2026-06-22 18:06:17 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:06:17 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:06:18 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:06:18 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:06:19 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc
2026-06-22 18:06:24 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:06:25 INFO [Posterior Diff ENS7] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS7_2020020611.nc
2026-06-22 18:06:25 INFO [Posterior Diff ENS7] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS7_2020020611.relative.nc
2026-06-22 18:06:25 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:06:25 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:06:26 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:06:26 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:06:27 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc
2026-06-22 18:06:32 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:06:33 INFO [Posterior Diff ENS8] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS8_2020020611.nc
2026-06-22 18:06:33 INFO [Posterior Diff ENS8] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS8_2020020611.relative.nc
2026-06-22 18:06:33 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:06:33 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:06:34 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:06:34 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:06:35 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc
2026-06-22 18:06:40 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:06:40 INFO [Posterior Diff ENS9] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS9_2020020611.nc
2026-06-22 18:06:40 INFO [Posterior Diff ENS9] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS9_2020020611.relative.nc
2026-06-22 18:06:41 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:06:41 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:06:42 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:06:42 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:06:42 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc
2026-06-22 18:06:48 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:06:48 INFO [Posterior Diff ENS10] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS10_2020020611.nc
2026-06-22 18:06:48 INFO [Posterior Diff ENS10] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS10_2020020611.relative.nc
2026-06-22 18:06:48 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:06:48 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:06:49 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:06:49 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:06:50 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc
2026-06-22 18:06:55 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:06:56 INFO [Posterior Diff ENS11] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS11_2020020611.nc
2026-06-22 18:06:56 INFO [Posterior Diff ENS11] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS11_2020020611.relative.nc
2026-06-22 18:06:56 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:06:56 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:06:57 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:06:57 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:06:58 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc
2026-06-22 18:07:03 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:07:03 INFO [Posterior Diff ENS12] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS12_2020020611.nc
2026-06-22 18:07:03 INFO [Posterior Diff ENS12] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS12_2020020611.relative.nc
2026-06-22 18:07:03 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:07:03 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:07:04 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:07:04 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:07:05 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc
2026-06-22 18:07:10 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:07:11 INFO [Posterior Diff ENS13] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS13_2020020611.nc
2026-06-22 18:07:11 INFO [Posterior Diff ENS13] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS13_2020020611.relative.nc
2026-06-22 18:07:11 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:07:11 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:07:12 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:07:12 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:07:13 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc
2026-06-22 18:07:18 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:07:18 INFO [Posterior Diff ENS14] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS14_2020020611.nc
2026-06-22 18:07:18 INFO [Posterior Diff ENS14] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS14_2020020611.relative.nc
2026-06-22 18:07:19 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:07:19 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:07:20 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:07:20 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:07:20 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc
2026-06-22 18:07:25 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:07:26 INFO [Posterior Diff ENS15] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS15_2020020611.nc
2026-06-22 18:07:26 INFO [Posterior Diff ENS15] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS15_2020020611.relative.nc
2026-06-22 18:07:26 INFO Next run starts from 2020-02-06 11:00:00
2026-06-22 18:07:26 INFO Cycle is DONE; starting a new loop!
2026-06-22 18:07:26 INFO [TIME] step_end current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-22 18:07:26 INFO [TIME] step_start current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-22 18:07:26 INFO [TIME] window start=2020-02-06 11:00:00 end=2020-02-06 13:00:00 run_hours=2 has_assimilation=True
2026-06-22 18:07:26 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:07:27 INFO Hourly dataset computed and listing created
2026-06-22 18:07:31 INFO Hourly dataset computed
2026-06-22 18:07:31 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:07:32 INFO Hourly dataset computed and listing created
2026-06-22 18:07:33 INFO Hourly dataset computed
2026-06-22 18:07:33 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:07:33 INFO Hourly dataset computed and listing created
2026-06-22 18:07:34 INFO Hourly dataset computed
2026-06-22 18:07:34 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:07:35 INFO Hourly dataset computed and listing created
2026-06-22 18:07:36 INFO Hourly dataset computed
2026-06-22 18:07:36 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:07:37 INFO Hourly dataset computed and listing created
2026-06-22 18:07:38 INFO Hourly dataset computed
2026-06-22 18:07:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:07:39 INFO Hourly dataset computed and listing created
2026-06-22 18:07:40 INFO Hourly dataset computed
2026-06-22 18:07:40 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:07:40 INFO Hourly dataset computed and listing created
2026-06-22 18:07:41 INFO Hourly dataset computed
2026-06-22 18:07:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:07:42 INFO Hourly dataset computed and listing created
2026-06-22 18:07:43 INFO Hourly dataset computed
2026-06-22 18:07:43 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:07:44 INFO Hourly dataset computed and listing created
2026-06-22 18:07:45 INFO Hourly dataset computed
2026-06-22 18:07:45 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:07:46 INFO Hourly dataset computed and listing created
2026-06-22 18:07:47 INFO Hourly dataset computed
2026-06-22 18:07:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:07:48 INFO Hourly dataset computed and listing created
2026-06-22 18:07:48 INFO Hourly dataset computed
2026-06-22 18:07:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:07:49 INFO Hourly dataset computed and listing created
2026-06-22 18:07:50 INFO Hourly dataset computed
2026-06-22 18:07:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:07:51 INFO Hourly dataset computed and listing created
2026-06-22 18:07:52 INFO Hourly dataset computed
2026-06-22 18:07:52 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:07:53 INFO Hourly dataset computed and listing created
2026-06-22 18:07:54 INFO Hourly dataset computed
2026-06-22 18:07:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:07:55 INFO Hourly dataset computed and listing created
2026-06-22 18:07:55 INFO Hourly dataset computed
2026-06-22 18:07:55 INFO ---------->>> Running CHIMERE model from 2020-02-06 11:00:00 to 2020-02-06 13:00:00
2026-06-22 18:07:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:07:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1
2026-06-22 18:07:55 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc
2026-06-22 18:07:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-06-22 18:07:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:07:55 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-06-22 18:07:55 INFO Queuing job for member 1...
2026-06-22 18:07:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:07:55 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-06-22 18:07:56 INFO Found: ['4952038']
2026-06-22 18:08:01 INFO [TGCC-IRENE] Submitted job with ID:['4952038']
2026-06-22 18:08:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:08:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2
2026-06-22 18:08:01 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc
2026-06-22 18:08:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-06-22 18:08:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:08:01 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-06-22 18:08:01 INFO Queuing job for member 2...
2026-06-22 18:08:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:08:01 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-06-22 18:08:02 INFO Found: ['4952041']
2026-06-22 18:08:07 INFO [TGCC-IRENE] Submitted job with ID:['4952041']
2026-06-22 18:08:07 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:08:07 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3
2026-06-22 18:08:07 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc
2026-06-22 18:08:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-06-22 18:08:07 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:08:07 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-06-22 18:08:07 INFO Queuing job for member 3...
2026-06-22 18:08:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:08:07 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-06-22 18:08:08 INFO Found: ['4952044']
2026-06-22 18:08:13 INFO [TGCC-IRENE] Submitted job with ID:['4952044']
2026-06-22 18:08:13 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:08:13 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4
2026-06-22 18:08:13 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc
2026-06-22 18:08:13 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-06-22 18:08:13 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:08:13 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-06-22 18:08:13 INFO Queuing job for member 4...
2026-06-22 18:08:13 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:08:13 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-06-22 18:08:13 INFO Found: ['4952049']
2026-06-22 18:08:18 INFO [TGCC-IRENE] Submitted job with ID:['4952049']
2026-06-22 18:08:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:08:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5
2026-06-22 18:08:18 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc
2026-06-22 18:08:19 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-06-22 18:08:19 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:08:19 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-06-22 18:08:19 INFO Queuing job for member 5...
2026-06-22 18:08:19 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:08:19 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-06-22 18:08:19 INFO Found: ['4952054']
2026-06-22 18:08:24 INFO [TGCC-IRENE] Submitted job with ID:['4952054']
2026-06-22 18:08:24 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:08:24 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6
2026-06-22 18:08:24 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc
2026-06-22 18:08:24 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-06-22 18:08:24 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:08:24 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-06-22 18:08:24 INFO Queuing job for member 6...
2026-06-22 18:08:24 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:08:24 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-06-22 18:08:25 INFO Found: ['4952058']
2026-06-22 18:08:30 INFO [TGCC-IRENE] Submitted job with ID:['4952058']
2026-06-22 18:08:30 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:08:30 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7
2026-06-22 18:08:30 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc
2026-06-22 18:08:30 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-06-22 18:08:30 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:08:30 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-06-22 18:08:30 INFO Queuing job for member 7...
2026-06-22 18:08:30 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:08:30 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-06-22 18:08:31 INFO Found: ['4952061']
2026-06-22 18:08:36 INFO [TGCC-IRENE] Submitted job with ID:['4952061']
2026-06-22 18:08:36 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:08:36 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8
2026-06-22 18:08:36 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc
2026-06-22 18:08:36 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-06-22 18:08:36 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:08:36 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-06-22 18:08:36 INFO Queuing job for member 8...
2026-06-22 18:08:36 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:08:36 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-06-22 18:08:36 INFO Found: ['4952064']
2026-06-22 18:08:41 INFO [TGCC-IRENE] Submitted job with ID:['4952064']
2026-06-22 18:08:41 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:08:41 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9
2026-06-22 18:08:41 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc
2026-06-22 18:08:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-06-22 18:08:42 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:08:42 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-06-22 18:08:42 INFO Queuing job for member 9...
2026-06-22 18:08:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:08:42 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-06-22 18:08:42 INFO Found: ['4952067']
2026-06-22 18:08:47 INFO [TGCC-IRENE] Submitted job with ID:['4952067']
2026-06-22 18:08:47 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:08:47 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10
2026-06-22 18:08:47 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc
2026-06-22 18:08:47 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-06-22 18:08:47 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:08:47 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-06-22 18:08:47 INFO Queuing job for member 10...
2026-06-22 18:08:47 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:08:47 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-06-22 18:08:48 INFO Found: ['4952072']
2026-06-22 18:08:53 INFO [TGCC-IRENE] Submitted job with ID:['4952072']
2026-06-22 18:08:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:08:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11
2026-06-22 18:08:53 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc
2026-06-22 18:08:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-06-22 18:08:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:08:53 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-06-22 18:08:53 INFO Queuing job for member 11...
2026-06-22 18:08:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:08:53 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-06-22 18:08:54 INFO Found: ['4952075']
2026-06-22 18:08:59 INFO [TGCC-IRENE] Submitted job with ID:['4952075']
2026-06-22 18:08:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:08:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12
2026-06-22 18:08:59 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc
2026-06-22 18:08:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-06-22 18:08:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:08:59 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-06-22 18:08:59 INFO Queuing job for member 12...
2026-06-22 18:08:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:08:59 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-06-22 18:09:00 INFO Found: ['4952078']
2026-06-22 18:09:05 INFO [TGCC-IRENE] Submitted job with ID:['4952078']
2026-06-22 18:09:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:09:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13
2026-06-22 18:09:05 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc
2026-06-22 18:09:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-06-22 18:09:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:09:05 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-06-22 18:09:05 INFO Queuing job for member 13...
2026-06-22 18:09:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:09:05 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-06-22 18:09:05 INFO Found: ['4952080']
2026-06-22 18:09:10 INFO [TGCC-IRENE] Submitted job with ID:['4952080']
2026-06-22 18:09:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:09:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14
2026-06-22 18:09:10 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc
2026-06-22 18:09:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-06-22 18:09:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:09:10 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-06-22 18:09:10 INFO Queuing job for member 14...
2026-06-22 18:09:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:09:10 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-06-22 18:09:11 INFO Found: ['4952084']
2026-06-22 18:09:16 INFO [TGCC-IRENE] Submitted job with ID:['4952084']
2026-06-22 18:09:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:09:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15
2026-06-22 18:09:16 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc
2026-06-22 18:09:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-06-22 18:09:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:09:16 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-06-22 18:09:16 INFO Queuing job for member 15...
2026-06-22 18:09:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:09:16 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-06-22 18:09:17 INFO Found: ['4952088']
2026-06-22 18:09:22 INFO [TGCC-IRENE] Submitted job with ID:['4952088']
2026-06-22 18:09:22 INFO Checking job status ...
2026-06-22 18:09:23 INFO None 4952038: status RUNNING/PENDING
2026-06-22 18:09:23 INFO None 4952041: status RUNNING/PENDING
2026-06-22 18:09:23 INFO None 4952044: status RUNNING/PENDING
2026-06-22 18:09:23 INFO None 4952049: status RUNNING/PENDING
2026-06-22 18:09:23 INFO None 4952054: status RUNNING/PENDING
2026-06-22 18:09:23 INFO None 4952058: status RUNNING/PENDING
2026-06-22 18:09:23 INFO None 4952061: status RUNNING/PENDING
2026-06-22 18:09:23 INFO None 4952064: status RUNNING/PENDING
2026-06-22 18:09:23 INFO None 4952067: status RUNNING/PENDING
2026-06-22 18:09:23 INFO None 4952072: status RUNNING/PENDING
2026-06-22 18:09:23 INFO None 4952075: status RUNNING/PENDING
2026-06-22 18:09:23 INFO None 4952078: status RUNNING/PENDING
2026-06-22 18:09:23 INFO None 4952080: status RUNNING/PENDING
2026-06-22 18:09:23 INFO None 4952084: status RUNNING/PENDING
2026-06-22 18:09:23 INFO None 4952088: status RUNNING/PENDING
2026-06-22 18:09:23 INFO Jobs still running: ['4952038', '4952041', '4952044', '4952049', '4952054', '4952058', '4952061', '4952064', '4952067', '4952072', '4952075', '4952078', '4952080', '4952084', '4952088']. Waiting...
2026-06-22 18:09:38 INFO None 4952038: status RUNNING/PENDING
2026-06-22 18:09:38 INFO None 4952041: status RUNNING/PENDING
2026-06-22 18:09:38 INFO None 4952044: status RUNNING/PENDING
2026-06-22 18:09:38 INFO None 4952049: status RUNNING/PENDING
2026-06-22 18:09:38 INFO None 4952054: status RUNNING/PENDING
2026-06-22 18:09:38 INFO None 4952058: status RUNNING/PENDING
2026-06-22 18:09:38 INFO None 4952061: status RUNNING/PENDING
2026-06-22 18:09:38 INFO None 4952064: status RUNNING/PENDING
2026-06-22 18:09:38 INFO None 4952067: status RUNNING/PENDING
2026-06-22 18:09:38 INFO None 4952072: status RUNNING/PENDING
2026-06-22 18:09:38 INFO None 4952075: status RUNNING/PENDING
2026-06-22 18:09:38 INFO None 4952078: status RUNNING/PENDING
2026-06-22 18:09:38 INFO None 4952080: status RUNNING/PENDING
2026-06-22 18:09:38 INFO None 4952084: status RUNNING/PENDING
2026-06-22 18:09:38 INFO None 4952088: status RUNNING/PENDING
2026-06-22 18:09:38 INFO Jobs still running: ['4952038', '4952041', '4952044', '4952049', '4952054', '4952058', '4952061', '4952064', '4952067', '4952072', '4952075', '4952078', '4952080', '4952084', '4952088']. Waiting...
2026-06-22 18:09:53 INFO None 4952038: status RUNNING/PENDING
2026-06-22 18:09:53 INFO None 4952041: status RUNNING/PENDING
2026-06-22 18:09:53 INFO None 4952044: status RUNNING/PENDING
2026-06-22 18:09:53 INFO None 4952049: status RUNNING/PENDING
2026-06-22 18:09:53 INFO None 4952054: status RUNNING/PENDING
2026-06-22 18:09:53 INFO None 4952058: status RUNNING/PENDING
2026-06-22 18:09:53 INFO None 4952061: status RUNNING/PENDING
2026-06-22 18:09:53 INFO None 4952064: status RUNNING/PENDING
2026-06-22 18:09:53 INFO None 4952067: status RUNNING/PENDING
2026-06-22 18:09:53 INFO None 4952072: status RUNNING/PENDING
2026-06-22 18:09:53 INFO None 4952075: status RUNNING/PENDING
2026-06-22 18:09:53 INFO None 4952078: status RUNNING/PENDING
2026-06-22 18:09:54 INFO None 4952080: status RUNNING/PENDING
2026-06-22 18:09:54 INFO None 4952084: status RUNNING/PENDING
2026-06-22 18:09:54 INFO None 4952088: status RUNNING/PENDING
2026-06-22 18:09:54 INFO Jobs still running: ['4952038', '4952041', '4952044', '4952049', '4952054', '4952058', '4952061', '4952064', '4952067', '4952072', '4952075', '4952078', '4952080', '4952084', '4952088']. Waiting...
2026-06-22 18:10:09 INFO None 4952038: status RUNNING/PENDING
2026-06-22 18:10:09 INFO None 4952041: status RUNNING/PENDING
2026-06-22 18:10:09 INFO None 4952044: status RUNNING/PENDING
2026-06-22 18:10:09 INFO None 4952049: status RUNNING/PENDING
2026-06-22 18:10:09 INFO None 4952054: status RUNNING/PENDING
2026-06-22 18:10:09 INFO None 4952058: status RUNNING/PENDING
2026-06-22 18:10:09 INFO None 4952061: status RUNNING/PENDING
2026-06-22 18:10:09 INFO None 4952064: status RUNNING/PENDING
2026-06-22 18:10:09 INFO None 4952067: status RUNNING/PENDING
2026-06-22 18:10:09 INFO None 4952072: status RUNNING/PENDING
2026-06-22 18:10:09 INFO None 4952075: status RUNNING/PENDING
2026-06-22 18:10:09 INFO None 4952078: status RUNNING/PENDING
2026-06-22 18:10:09 INFO None 4952080: status RUNNING/PENDING
2026-06-22 18:10:09 INFO None 4952084: status RUNNING/PENDING
2026-06-22 18:10:09 INFO None 4952088: status RUNNING/PENDING
2026-06-22 18:10:09 INFO Jobs still running: ['4952038', '4952041', '4952044', '4952049', '4952054', '4952058', '4952061', '4952064', '4952067', '4952072', '4952075', '4952078', '4952080', '4952084', '4952088']. Waiting...
2026-06-22 18:10:24 INFO None 4952038: status RUNNING/PENDING
2026-06-22 18:10:24 INFO None 4952041: status RUNNING/PENDING
2026-06-22 18:10:24 INFO None 4952044: status RUNNING/PENDING
2026-06-22 18:10:24 INFO None 4952049: status RUNNING/PENDING
2026-06-22 18:10:24 INFO None 4952054: status RUNNING/PENDING
2026-06-22 18:10:24 INFO None 4952058: status RUNNING/PENDING
2026-06-22 18:10:24 INFO None 4952061: status RUNNING/PENDING
2026-06-22 18:10:24 INFO None 4952064: status RUNNING/PENDING
2026-06-22 18:10:24 INFO None 4952067: status RUNNING/PENDING
2026-06-22 18:10:24 INFO None 4952072: status RUNNING/PENDING
2026-06-22 18:10:24 INFO None 4952075: status RUNNING/PENDING
2026-06-22 18:10:24 INFO None 4952078: status RUNNING/PENDING
2026-06-22 18:10:24 INFO None 4952080: status RUNNING/PENDING
2026-06-22 18:10:24 INFO None 4952084: status RUNNING/PENDING
2026-06-22 18:10:24 INFO None 4952088: status RUNNING/PENDING
2026-06-22 18:10:24 INFO Jobs still running: ['4952038', '4952041', '4952044', '4952049', '4952054', '4952058', '4952061', '4952064', '4952067', '4952072', '4952075', '4952078', '4952080', '4952084', '4952088']. Waiting...
2026-06-22 18:10:39 INFO None 4952038: status RUNNING/PENDING
2026-06-22 18:10:39 INFO None 4952041: status RUNNING/PENDING
2026-06-22 18:10:39 INFO None 4952044: status RUNNING/PENDING
2026-06-22 18:10:39 INFO None 4952049: status RUNNING/PENDING
2026-06-22 18:10:39 INFO None 4952054: status RUNNING/PENDING
2026-06-22 18:10:39 INFO None 4952058: status RUNNING/PENDING
2026-06-22 18:10:39 INFO None 4952061: status RUNNING/PENDING
2026-06-22 18:10:39 INFO None 4952064: status RUNNING/PENDING
2026-06-22 18:10:39 INFO None 4952067: status RUNNING/PENDING
2026-06-22 18:10:39 INFO None 4952072: status RUNNING/PENDING
2026-06-22 18:10:39 INFO None 4952075: status RUNNING/PENDING
2026-06-22 18:10:39 INFO None 4952078: status RUNNING/PENDING
2026-06-22 18:10:39 INFO None 4952080: status RUNNING/PENDING
2026-06-22 18:10:39 INFO None 4952084: status RUNNING/PENDING
2026-06-22 18:10:39 INFO None 4952088: status RUNNING/PENDING
2026-06-22 18:10:39 INFO Jobs still running: ['4952038', '4952041', '4952044', '4952049', '4952054', '4952058', '4952061', '4952064', '4952067', '4952072', '4952075', '4952078', '4952080', '4952084', '4952088']. Waiting...
2026-06-22 18:10:54 INFO None 4952038: status RUNNING/PENDING
2026-06-22 18:10:54 INFO None 4952041: status RUNNING/PENDING
2026-06-22 18:10:55 INFO None 4952044: status RUNNING/PENDING
2026-06-22 18:10:55 INFO None 4952049: status RUNNING/PENDING
2026-06-22 18:10:55 INFO None 4952054: status RUNNING/PENDING
2026-06-22 18:10:55 INFO None 4952058: status RUNNING/PENDING
2026-06-22 18:10:55 INFO None 4952061: status RUNNING/PENDING
2026-06-22 18:10:55 INFO None 4952064: status RUNNING/PENDING
2026-06-22 18:10:55 INFO None 4952067: status RUNNING/PENDING
2026-06-22 18:10:55 INFO None 4952072: status RUNNING/PENDING
2026-06-22 18:10:55 INFO None 4952075: status RUNNING/PENDING
2026-06-22 18:10:55 INFO None 4952078: status RUNNING/PENDING
2026-06-22 18:10:55 INFO None 4952080: status RUNNING/PENDING
2026-06-22 18:10:55 INFO None 4952084: status RUNNING/PENDING
2026-06-22 18:10:55 INFO None 4952088: status RUNNING/PENDING
2026-06-22 18:10:55 INFO Jobs still running: ['4952038', '4952041', '4952044', '4952049', '4952054', '4952058', '4952061', '4952064', '4952067', '4952072', '4952075', '4952078', '4952080', '4952084', '4952088']. Waiting...
2026-06-22 18:11:10 INFO None 4952038: status RUNNING/PENDING
2026-06-22 18:11:10 INFO None 4952041: status RUNNING/PENDING
2026-06-22 18:11:10 INFO None 4952044: status RUNNING/PENDING
2026-06-22 18:11:10 INFO None 4952049: status RUNNING/PENDING
2026-06-22 18:11:10 INFO None 4952054: status RUNNING/PENDING
2026-06-22 18:11:10 INFO None 4952058: status RUNNING/PENDING
2026-06-22 18:11:10 INFO None 4952061: status RUNNING/PENDING
2026-06-22 18:11:10 INFO None 4952064: status RUNNING/PENDING
2026-06-22 18:11:10 INFO None 4952067: status RUNNING/PENDING
2026-06-22 18:11:10 INFO None 4952072: status RUNNING/PENDING
2026-06-22 18:11:10 INFO None 4952075: status RUNNING/PENDING
2026-06-22 18:11:10 INFO None 4952078: status RUNNING/PENDING
2026-06-22 18:11:10 INFO None 4952080: status RUNNING/PENDING
2026-06-22 18:11:10 INFO None 4952084: status RUNNING/PENDING
2026-06-22 18:11:10 INFO None 4952088: status RUNNING/PENDING
2026-06-22 18:11:10 INFO Jobs still running: ['4952038', '4952041', '4952044', '4952049', '4952054', '4952058', '4952061', '4952064', '4952067', '4952072', '4952075', '4952078', '4952080', '4952084', '4952088']. Waiting...
2026-06-22 18:11:25 INFO None 4952038: status RUNNING/PENDING
2026-06-22 18:11:25 INFO None 4952041: status RUNNING/PENDING
2026-06-22 18:11:25 INFO None 4952044: status RUNNING/PENDING
2026-06-22 18:11:25 INFO None 4952049: status RUNNING/PENDING
2026-06-22 18:11:25 INFO None 4952054: status RUNNING/PENDING
2026-06-22 18:11:25 INFO None 4952058: status RUNNING/PENDING
2026-06-22 18:11:25 INFO None 4952061: status RUNNING/PENDING
2026-06-22 18:11:25 INFO None 4952064: status RUNNING/PENDING
2026-06-22 18:11:25 INFO None 4952067: status RUNNING/PENDING
2026-06-22 18:11:25 INFO None 4952072: status RUNNING/PENDING
2026-06-22 18:11:25 INFO None 4952075: status RUNNING/PENDING
2026-06-22 18:11:25 INFO None 4952078: status RUNNING/PENDING
2026-06-22 18:11:25 INFO None 4952080: status RUNNING/PENDING
2026-06-22 18:11:25 INFO None 4952084: status RUNNING/PENDING
2026-06-22 18:11:25 INFO None 4952088: status RUNNING/PENDING
2026-06-22 18:11:25 INFO Jobs still running: ['4952038', '4952041', '4952044', '4952049', '4952054', '4952058', '4952061', '4952064', '4952067', '4952072', '4952075', '4952078', '4952080', '4952084', '4952088']. Waiting...
2026-06-22 18:11:40 INFO None 4952038: status RUNNING/PENDING
2026-06-22 18:11:40 INFO None 4952041: status RUNNING/PENDING
2026-06-22 18:11:40 INFO None 4952044: status RUNNING/PENDING
2026-06-22 18:11:40 INFO None 4952049: status RUNNING/PENDING
2026-06-22 18:11:40 INFO None 4952054: status RUNNING/PENDING
2026-06-22 18:11:41 INFO None 4952058: status RUNNING/PENDING
2026-06-22 18:11:41 INFO None 4952061: status RUNNING/PENDING
2026-06-22 18:11:41 INFO None 4952064: status RUNNING/PENDING
2026-06-22 18:11:41 INFO None 4952067: status RUNNING/PENDING
2026-06-22 18:11:41 INFO None 4952072: status RUNNING/PENDING
2026-06-22 18:11:41 INFO None 4952075: status RUNNING/PENDING
2026-06-22 18:11:41 INFO None 4952078: status RUNNING/PENDING
2026-06-22 18:11:41 INFO None 4952080: status RUNNING/PENDING
2026-06-22 18:11:41 INFO None 4952084: status RUNNING/PENDING
2026-06-22 18:11:41 INFO None 4952088: status RUNNING/PENDING
2026-06-22 18:11:41 INFO Jobs still running: ['4952038', '4952041', '4952044', '4952049', '4952054', '4952058', '4952061', '4952064', '4952067', '4952072', '4952075', '4952078', '4952080', '4952084', '4952088']. Waiting...
2026-06-22 18:11:56 INFO None 4952038: status RUNNING/PENDING
2026-06-22 18:11:56 INFO None 4952041: status RUNNING/PENDING
2026-06-22 18:11:56 INFO None 4952044: status FINISHED
2026-06-22 18:11:56 INFO None 4952049: status FINISHED
2026-06-22 18:11:56 INFO None 4952054: status RUNNING/PENDING
2026-06-22 18:11:56 INFO None 4952058: status FINISHED
2026-06-22 18:11:56 INFO None 4952061: status FINISHED
2026-06-22 18:11:56 INFO None 4952064: status RUNNING/PENDING
2026-06-22 18:11:56 INFO None 4952067: status RUNNING/PENDING
2026-06-22 18:11:56 INFO None 4952072: status RUNNING/PENDING
2026-06-22 18:11:56 INFO None 4952075: status RUNNING/PENDING
2026-06-22 18:11:56 INFO None 4952078: status RUNNING/PENDING
2026-06-22 18:11:56 INFO None 4952080: status RUNNING/PENDING
2026-06-22 18:11:56 INFO None 4952084: status RUNNING/PENDING
2026-06-22 18:11:56 INFO None 4952088: status RUNNING/PENDING
2026-06-22 18:11:56 INFO Jobs still running: ['4952038', '4952041', '4952054', '4952064', '4952067', '4952072', '4952075', '4952078', '4952080', '4952084', '4952088']. Waiting...
2026-06-22 18:12:11 INFO None 4952038: status RUNNING/PENDING
2026-06-22 18:12:11 INFO None 4952041: status RUNNING/PENDING
2026-06-22 18:12:11 INFO None 4952044: status FINISHED
2026-06-22 18:12:11 INFO None 4952049: status FINISHED
2026-06-22 18:12:11 INFO None 4952054: status RUNNING/PENDING
2026-06-22 18:12:11 INFO None 4952058: status FINISHED
2026-06-22 18:12:11 INFO None 4952061: status FINISHED
2026-06-22 18:12:11 INFO None 4952064: status FINISHED
2026-06-22 18:12:11 INFO None 4952067: status FINISHED
2026-06-22 18:12:11 INFO None 4952072: status RUNNING/PENDING
2026-06-22 18:12:11 INFO None 4952075: status RUNNING/PENDING
2026-06-22 18:12:11 INFO None 4952078: status RUNNING/PENDING
2026-06-22 18:12:11 INFO None 4952080: status RUNNING/PENDING
2026-06-22 18:12:11 INFO None 4952084: status RUNNING/PENDING
2026-06-22 18:12:11 INFO None 4952088: status RUNNING/PENDING
2026-06-22 18:12:11 INFO Jobs still running: ['4952038', '4952041', '4952054', '4952072', '4952075', '4952078', '4952080', '4952084', '4952088']. Waiting...
2026-06-22 18:12:26 INFO None 4952038: status RUNNING/PENDING
2026-06-22 18:12:26 INFO None 4952041: status FINISHED
2026-06-22 18:12:26 INFO None 4952044: status FINISHED
2026-06-22 18:12:26 INFO None 4952049: status FINISHED
2026-06-22 18:12:26 INFO None 4952054: status RUNNING/PENDING
2026-06-22 18:12:26 INFO None 4952058: status FINISHED
2026-06-22 18:12:26 INFO None 4952061: status FINISHED
2026-06-22 18:12:26 INFO None 4952064: status FINISHED
2026-06-22 18:12:26 INFO None 4952067: status FINISHED
2026-06-22 18:12:26 INFO None 4952072: status RUNNING/PENDING
2026-06-22 18:12:27 INFO None 4952075: status RUNNING/PENDING
2026-06-22 18:12:27 INFO None 4952078: status RUNNING/PENDING
2026-06-22 18:12:27 INFO None 4952080: status RUNNING/PENDING
2026-06-22 18:12:27 INFO None 4952084: status RUNNING/PENDING
2026-06-22 18:12:27 INFO None 4952088: status RUNNING/PENDING
2026-06-22 18:12:27 INFO Jobs still running: ['4952038', '4952054', '4952072', '4952075', '4952078', '4952080', '4952084', '4952088']. Waiting...
2026-06-22 18:12:42 INFO None 4952038: status RUNNING/PENDING
2026-06-22 18:12:42 INFO None 4952041: status FINISHED
2026-06-22 18:12:42 INFO None 4952044: status FINISHED
2026-06-22 18:12:42 INFO None 4952049: status FINISHED
2026-06-22 18:12:42 INFO None 4952054: status RUNNING/PENDING
2026-06-22 18:12:42 INFO None 4952058: status FINISHED
2026-06-22 18:12:42 INFO None 4952061: status FINISHED
2026-06-22 18:12:42 INFO None 4952064: status FINISHED
2026-06-22 18:12:42 INFO None 4952067: status FINISHED
2026-06-22 18:12:42 INFO None 4952072: status FINISHED
2026-06-22 18:12:42 INFO None 4952075: status FINISHED
2026-06-22 18:12:42 INFO None 4952078: status RUNNING/PENDING
2026-06-22 18:12:42 INFO None 4952080: status FINISHED
2026-06-22 18:12:42 INFO None 4952084: status RUNNING/PENDING
2026-06-22 18:12:42 INFO None 4952088: status RUNNING/PENDING
2026-06-22 18:12:42 INFO Jobs still running: ['4952038', '4952054', '4952078', '4952084', '4952088']. Waiting...
2026-06-22 18:12:57 INFO None 4952038: status RUNNING/PENDING
2026-06-22 18:12:57 INFO None 4952041: status FINISHED
2026-06-22 18:12:57 INFO None 4952044: status FINISHED
2026-06-22 18:12:57 INFO None 4952049: status FINISHED
2026-06-22 18:12:57 INFO None 4952054: status RUNNING/PENDING
2026-06-22 18:12:57 INFO None 4952058: status FINISHED
2026-06-22 18:12:57 INFO None 4952061: status FINISHED
2026-06-22 18:12:57 INFO None 4952064: status FINISHED
2026-06-22 18:12:57 INFO None 4952067: status FINISHED
2026-06-22 18:12:57 INFO None 4952072: status FINISHED
2026-06-22 18:12:57 INFO None 4952075: status FINISHED
2026-06-22 18:12:57 INFO None 4952078: status FINISHED
2026-06-22 18:12:57 INFO None 4952080: status FINISHED
2026-06-22 18:12:57 INFO None 4952084: status RUNNING/PENDING
2026-06-22 18:12:57 INFO None 4952088: status RUNNING/PENDING
2026-06-22 18:12:57 INFO Jobs still running: ['4952038', '4952054', '4952084', '4952088']. Waiting...
2026-06-22 18:13:12 INFO None 4952038: status RUNNING/PENDING
2026-06-22 18:13:12 INFO None 4952041: status FINISHED
2026-06-22 18:13:12 INFO None 4952044: status FINISHED
2026-06-22 18:13:12 INFO None 4952049: status FINISHED
2026-06-22 18:13:12 INFO None 4952054: status RUNNING/PENDING
2026-06-22 18:13:12 INFO None 4952058: status FINISHED
2026-06-22 18:13:12 INFO None 4952061: status FINISHED
2026-06-22 18:13:12 INFO None 4952064: status FINISHED
2026-06-22 18:13:12 INFO None 4952067: status FINISHED
2026-06-22 18:13:12 INFO None 4952072: status FINISHED
2026-06-22 18:13:12 INFO None 4952075: status FINISHED
2026-06-22 18:13:13 INFO None 4952078: status FINISHED
2026-06-22 18:13:13 INFO None 4952080: status FINISHED
2026-06-22 18:13:13 INFO None 4952084: status RUNNING/PENDING
2026-06-22 18:13:13 INFO None 4952088: status RUNNING/PENDING
2026-06-22 18:13:13 INFO Jobs still running: ['4952038', '4952054', '4952084', '4952088']. Waiting...
2026-06-22 18:13:28 INFO None 4952038: status RUNNING/PENDING
2026-06-22 18:13:28 INFO None 4952041: status FINISHED
2026-06-22 18:13:28 INFO None 4952044: status FINISHED
2026-06-22 18:13:28 INFO None 4952049: status FINISHED
2026-06-22 18:13:28 INFO None 4952054: status FINISHED
2026-06-22 18:13:28 INFO None 4952058: status FINISHED
2026-06-22 18:13:28 INFO None 4952061: status FINISHED
2026-06-22 18:13:28 INFO None 4952064: status FINISHED
2026-06-22 18:13:28 INFO None 4952067: status FINISHED
2026-06-22 18:13:28 INFO None 4952072: status FINISHED
2026-06-22 18:13:28 INFO None 4952075: status FINISHED
2026-06-22 18:13:28 INFO None 4952078: status FINISHED
2026-06-22 18:13:28 INFO None 4952080: status FINISHED
2026-06-22 18:13:28 INFO None 4952084: status FINISHED
2026-06-22 18:13:29 INFO None 4952088: status FINISHED
2026-06-22 18:13:29 INFO Jobs still running: ['4952038']. Waiting...
2026-06-22 18:13:44 INFO None 4952038: status RUNNING/PENDING
2026-06-22 18:13:44 INFO None 4952041: status FINISHED
2026-06-22 18:13:44 INFO None 4952044: status FINISHED
2026-06-22 18:13:44 INFO None 4952049: status FINISHED
2026-06-22 18:13:44 INFO None 4952054: status FINISHED
2026-06-22 18:13:44 INFO None 4952058: status FINISHED
2026-06-22 18:13:44 INFO None 4952061: status FINISHED
2026-06-22 18:13:44 INFO None 4952064: status FINISHED
2026-06-22 18:13:44 INFO None 4952067: status FINISHED
2026-06-22 18:13:44 INFO None 4952072: status FINISHED
2026-06-22 18:13:44 INFO None 4952075: status FINISHED
2026-06-22 18:13:44 INFO None 4952078: status FINISHED
2026-06-22 18:13:44 INFO None 4952080: status FINISHED
2026-06-22 18:13:44 INFO None 4952084: status FINISHED
2026-06-22 18:13:44 INFO None 4952088: status FINISHED
2026-06-22 18:13:44 INFO Jobs still running: ['4952038']. Waiting...
2026-06-22 18:13:59 INFO None 4952038: status RUNNING/PENDING
2026-06-22 18:13:59 INFO None 4952041: status FINISHED
2026-06-22 18:13:59 INFO None 4952044: status FINISHED
2026-06-22 18:13:59 INFO None 4952049: status FINISHED
2026-06-22 18:13:59 INFO None 4952054: status FINISHED
2026-06-22 18:13:59 INFO None 4952058: status FINISHED
2026-06-22 18:13:59 INFO None 4952061: status FINISHED
2026-06-22 18:13:59 INFO None 4952064: status FINISHED
2026-06-22 18:13:59 INFO None 4952067: status FINISHED
2026-06-22 18:13:59 INFO None 4952072: status FINISHED
2026-06-22 18:13:59 INFO None 4952075: status FINISHED
2026-06-22 18:13:59 INFO None 4952078: status FINISHED
2026-06-22 18:13:59 INFO None 4952080: status FINISHED
2026-06-22 18:13:59 INFO None 4952084: status FINISHED
2026-06-22 18:13:59 INFO None 4952088: status FINISHED
2026-06-22 18:13:59 INFO Jobs still running: ['4952038']. Waiting...
2026-06-22 18:14:15 INFO None 4952038: status RUNNING/PENDING
2026-06-22 18:14:15 INFO None 4952041: status FINISHED
2026-06-22 18:14:15 INFO None 4952044: status FINISHED
2026-06-22 18:14:15 INFO None 4952049: status FINISHED
2026-06-22 18:14:15 INFO None 4952054: status FINISHED
2026-06-22 18:14:15 INFO None 4952058: status FINISHED
2026-06-22 18:14:15 INFO None 4952061: status FINISHED
2026-06-22 18:14:15 INFO None 4952064: status FINISHED
2026-06-22 18:14:15 INFO None 4952067: status FINISHED
2026-06-22 18:14:15 INFO None 4952072: status FINISHED
2026-06-22 18:14:15 INFO None 4952075: status FINISHED
2026-06-22 18:14:15 INFO None 4952078: status FINISHED
2026-06-22 18:14:16 INFO None 4952080: status FINISHED
2026-06-22 18:14:16 INFO None 4952084: status FINISHED
2026-06-22 18:14:16 INFO None 4952088: status FINISHED
2026-06-22 18:14:16 INFO Jobs still running: ['4952038']. Waiting...
2026-06-22 18:14:31 INFO None 4952038: status RUNNING/PENDING
2026-06-22 18:14:31 INFO None 4952041: status FINISHED
2026-06-22 18:14:31 INFO None 4952044: status FINISHED
2026-06-22 18:14:31 INFO None 4952049: status FINISHED
2026-06-22 18:14:31 INFO None 4952054: status FINISHED
2026-06-22 18:14:31 INFO None 4952058: status FINISHED
2026-06-22 18:14:31 INFO None 4952061: status FINISHED
2026-06-22 18:14:31 INFO None 4952064: status FINISHED
2026-06-22 18:14:31 INFO None 4952067: status FINISHED
2026-06-22 18:14:31 INFO None 4952072: status FINISHED
2026-06-22 18:14:31 INFO None 4952075: status FINISHED
2026-06-22 18:14:31 INFO None 4952078: status FINISHED
2026-06-22 18:14:31 INFO None 4952080: status FINISHED
2026-06-22 18:14:31 INFO None 4952084: status FINISHED
2026-06-22 18:14:31 INFO None 4952088: status FINISHED
2026-06-22 18:14:31 INFO Jobs still running: ['4952038']. Waiting...
2026-06-22 18:14:46 INFO None 4952038: status FINISHED
2026-06-22 18:14:46 INFO None 4952041: status FINISHED
2026-06-22 18:14:46 INFO None 4952044: status FINISHED
2026-06-22 18:14:46 INFO None 4952049: status FINISHED
2026-06-22 18:14:46 INFO None 4952054: status FINISHED
2026-06-22 18:14:46 INFO None 4952058: status FINISHED
2026-06-22 18:14:46 INFO None 4952061: status FINISHED
2026-06-22 18:14:46 INFO None 4952064: status FINISHED
2026-06-22 18:14:46 INFO None 4952067: status FINISHED
2026-06-22 18:14:46 INFO None 4952072: status FINISHED
2026-06-22 18:14:46 INFO None 4952075: status FINISHED
2026-06-22 18:14:46 INFO None 4952078: status FINISHED
2026-06-22 18:14:46 INFO None 4952080: status FINISHED
2026-06-22 18:14:46 INFO None 4952084: status FINISHED
2026-06-22 18:14:46 INFO None 4952088: status FINISHED
2026-06-22 18:14:46 INFO Jobs ['4952038', '4952041', '4952044', '4952049', '4952054', '4952058', '4952061', '4952064', '4952067', '4952072', '4952075', '4952078', '4952080', '4952084', '4952088'] have finished
2026-06-22 18:14:46 INFO Checking restart files were created ...
2026-06-22 18:14:46 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020611_2_ENS1.nc(1002685915 bytes)
2026-06-22 18:14:46 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020611_2_ENS2.nc(1002685915 bytes)
2026-06-22 18:14:46 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020611_2_ENS3.nc(1002685915 bytes)
2026-06-22 18:14:46 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020611_2_ENS4.nc(1002685915 bytes)
2026-06-22 18:14:46 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020611_2_ENS5.nc(1002685915 bytes)
2026-06-22 18:14:46 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020611_2_ENS6.nc(1002685915 bytes)
2026-06-22 18:14:46 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020611_2_ENS7.nc(1002685915 bytes)
2026-06-22 18:14:46 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020611_2_ENS8.nc(1002685915 bytes)
2026-06-22 18:14:46 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020611_2_ENS9.nc(1002685915 bytes)
2026-06-22 18:14:46 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020611_2_ENS10.nc(1002685915 bytes)
2026-06-22 18:14:46 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020611_2_ENS11.nc(1002685915 bytes)
2026-06-22 18:14:46 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020611_2_ENS12.nc(1002685915 bytes)
2026-06-22 18:14:46 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020611_2_ENS13.nc(1002685915 bytes)
2026-06-22 18:14:46 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020611_2_ENS14.nc(1002685915 bytes)
2026-06-22 18:14:46 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020611_2_ENS15.nc(1002685915 bytes)
2026-06-22 18:14:46 INFO  Run_model() completed successfully.
2026-06-22 18:14:46 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-22 18:14:46 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 13:00:00 days=153072 seconds=46800
2026-06-22 18:14:46 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-06-22 18:14:46 INFO [TIME] increment current_time 2020-02-06 11:00:00 -> 2020-02-06 13:00:00
2026-06-22 18:14:46 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-22 18:14:46 INFO ---------->>> Running process_satellite_data()
2026-06-22 18:14:46 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12004.nc
2026-06-22 18:14:46 INFO ---------->>> Running run_obs_converter()
2026-06-22 18:14:46 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_45135_153072.out
2026-06-22 18:14:46 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_45135_153072.out
2026-06-22 18:14:46 INFO ---------->>> Running DART
2026-06-22 18:14:46 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-06-22 18:14:46 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/chim_ENS1_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/chim_ENS1_2020020613_1_out_toDART.nc
2026-06-22 18:14:47 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/chim_ENS2_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/chim_ENS2_2020020613_1_out_toDART.nc
2026-06-22 18:14:47 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/chim_ENS3_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/chim_ENS3_2020020613_1_out_toDART.nc
2026-06-22 18:14:47 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/chim_ENS4_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/chim_ENS4_2020020613_1_out_toDART.nc
2026-06-22 18:14:48 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/chim_ENS5_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/chim_ENS5_2020020613_1_out_toDART.nc
2026-06-22 18:14:48 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/chim_ENS6_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/chim_ENS6_2020020613_1_out_toDART.nc
2026-06-22 18:14:48 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/chim_ENS7_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/chim_ENS7_2020020613_1_out_toDART.nc
2026-06-22 18:14:48 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/chim_ENS8_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/chim_ENS8_2020020613_1_out_toDART.nc
2026-06-22 18:14:49 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/chim_ENS9_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/chim_ENS9_2020020613_1_out_toDART.nc
2026-06-22 18:14:49 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/chim_ENS10_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/chim_ENS10_2020020613_1_out_toDART.nc
2026-06-22 18:14:49 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/chim_ENS11_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/chim_ENS11_2020020613_1_out_toDART.nc
2026-06-22 18:14:50 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/chim_ENS12_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/chim_ENS12_2020020613_1_out_toDART.nc
2026-06-22 18:14:50 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/chim_ENS13_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/chim_ENS13_2020020613_1_out_toDART.nc
2026-06-22 18:14:50 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/chim_ENS14_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/chim_ENS14_2020020613_1_out_toDART.nc
2026-06-22 18:14:51 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/chim_ENS15_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/chim_ENS15_2020020613_1_out_toDART.nc
2026-06-22 18:14:51 INFO Replacement input_template.nml → input.nml completed successfully.
2026-06-22 18:14:51 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-06-22 18:14:51 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-06-22 18:14:51 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-06-22 18:14:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-06-22 18:14:51 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-06-22 18:15:05 INFO Found: []
2026-06-22 18:15:05 INFO No job id returned by command ./run_filter.bsh
2026-06-22 18:15:05 INFO No monitoring will be performed
2026-06-22 18:15:05 INFO Moving DART output files to analysis and preassim directories for date 2020020613 if present ...
2026-06-22 18:15:05 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-06-22 18:15:05 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-06-22 18:15:05 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-06-22 18:15:05 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-06-22 18:15:05 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-06-22 18:15:05 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-06-22 18:15:05 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-06-22 18:15:05 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-06-22 18:15:05 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-06-22 18:15:05 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-06-22 18:15:05 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-06-22 18:15:05 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-06-22 18:15:05 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-06-22 18:15:05 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-06-22 18:15:05 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-06-22 18:15:05 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-06-22 18:15:05 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-06-22 18:15:05 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-06-22 18:15:05 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-06-22 18:15:05 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-06-22 18:15:05 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-06-22 18:15:05 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-06-22 18:15:05 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-06-22 18:15:05 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-06-22 18:15:05 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-06-22 18:15:05 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-06-22 18:15:05 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-06-22 18:15:05 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-06-22 18:15:05 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-06-22 18:15:05 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-06-22 18:15:05 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-06-22 18:15:05 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020613'
2026-06-22 18:15:05 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-06-22 18:15:05 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020613'
2026-06-22 18:15:05 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-06-22 18:15:05 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-06-22 18:15:05 INFO run_dart() is DONE.
2026-06-22 18:15:05 INFO ---------->>> Running update_pollutant_in_end()
2026-06-22 18:15:06 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:15:06 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:15:07 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:15:07 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:15:07 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020611_2_ENS1.nc
2026-06-22 18:15:12 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:15:13 INFO [Posterior Diff ENS1] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS1_2020020613.nc
2026-06-22 18:15:13 INFO [Posterior Diff ENS1] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS1_2020020613.relative.nc
2026-06-22 18:15:13 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:15:13 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:15:14 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:15:14 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:15:15 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020611_2_ENS2.nc
2026-06-22 18:15:20 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:15:20 INFO [Posterior Diff ENS2] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS2_2020020613.nc
2026-06-22 18:15:20 INFO [Posterior Diff ENS2] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS2_2020020613.relative.nc
2026-06-22 18:15:21 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:15:21 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:15:21 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:15:21 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:15:22 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020611_2_ENS3.nc
2026-06-22 18:15:27 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:15:28 INFO [Posterior Diff ENS3] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS3_2020020613.nc
2026-06-22 18:15:28 INFO [Posterior Diff ENS3] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS3_2020020613.relative.nc
2026-06-22 18:15:28 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:15:28 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:15:29 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:15:29 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:15:30 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020611_2_ENS4.nc
2026-06-22 18:15:35 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:15:35 INFO [Posterior Diff ENS4] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS4_2020020613.nc
2026-06-22 18:15:35 INFO [Posterior Diff ENS4] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS4_2020020613.relative.nc
2026-06-22 18:15:35 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:15:35 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:15:36 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:15:36 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:15:37 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020611_2_ENS5.nc
2026-06-22 18:15:42 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:15:43 INFO [Posterior Diff ENS5] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS5_2020020613.nc
2026-06-22 18:15:43 INFO [Posterior Diff ENS5] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS5_2020020613.relative.nc
2026-06-22 18:15:43 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:15:43 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:15:44 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:15:44 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:15:45 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020611_2_ENS6.nc
2026-06-22 18:15:50 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:15:50 INFO [Posterior Diff ENS6] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS6_2020020613.nc
2026-06-22 18:15:50 INFO [Posterior Diff ENS6] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS6_2020020613.relative.nc
2026-06-22 18:15:51 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:15:51 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:15:52 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:15:52 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:15:53 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020611_2_ENS7.nc
2026-06-22 18:15:58 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:15:58 INFO [Posterior Diff ENS7] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS7_2020020613.nc
2026-06-22 18:15:58 INFO [Posterior Diff ENS7] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS7_2020020613.relative.nc
2026-06-22 18:15:59 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:15:59 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:16:00 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:16:00 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:16:01 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020611_2_ENS8.nc
2026-06-22 18:16:06 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:16:06 INFO [Posterior Diff ENS8] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS8_2020020613.nc
2026-06-22 18:16:06 INFO [Posterior Diff ENS8] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS8_2020020613.relative.nc
2026-06-22 18:16:06 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:16:06 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:16:07 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:16:07 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:16:08 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020611_2_ENS9.nc
2026-06-22 18:16:14 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:16:14 INFO [Posterior Diff ENS9] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS9_2020020613.nc
2026-06-22 18:16:14 INFO [Posterior Diff ENS9] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS9_2020020613.relative.nc
2026-06-22 18:16:14 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:16:14 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:16:15 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:16:15 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:16:16 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020611_2_ENS10.nc
2026-06-22 18:16:21 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:16:22 INFO [Posterior Diff ENS10] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS10_2020020613.nc
2026-06-22 18:16:22 INFO [Posterior Diff ENS10] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS10_2020020613.relative.nc
2026-06-22 18:16:22 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:16:22 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:16:23 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:16:23 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:16:24 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020611_2_ENS11.nc
2026-06-22 18:16:29 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:16:29 INFO [Posterior Diff ENS11] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS11_2020020613.nc
2026-06-22 18:16:29 INFO [Posterior Diff ENS11] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS11_2020020613.relative.nc
2026-06-22 18:16:30 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:16:30 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:16:30 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:16:30 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:16:31 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020611_2_ENS12.nc
2026-06-22 18:16:36 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:16:36 INFO [Posterior Diff ENS12] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS12_2020020613.nc
2026-06-22 18:16:36 INFO [Posterior Diff ENS12] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS12_2020020613.relative.nc
2026-06-22 18:16:37 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:16:37 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:16:38 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:16:38 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:16:39 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020611_2_ENS13.nc
2026-06-22 18:16:44 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:16:44 INFO [Posterior Diff ENS13] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS13_2020020613.nc
2026-06-22 18:16:44 INFO [Posterior Diff ENS13] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS13_2020020613.relative.nc
2026-06-22 18:16:45 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:16:45 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:16:45 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:16:45 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:16:46 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020611_2_ENS14.nc
2026-06-22 18:16:51 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:16:52 INFO [Posterior Diff ENS14] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS14_2020020613.nc
2026-06-22 18:16:52 INFO [Posterior Diff ENS14] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS14_2020020613.relative.nc
2026-06-22 18:16:52 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:16:52 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:16:53 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:16:53 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:16:54 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020611_2_ENS15.nc
2026-06-22 18:16:59 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:16:59 INFO [Posterior Diff ENS15] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS15_2020020613.nc
2026-06-22 18:16:59 INFO [Posterior Diff ENS15] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020613/diff_posterior_ENS15_2020020613.relative.nc
2026-06-22 18:17:00 INFO Next run starts from 2020-02-06 13:00:00
2026-06-22 18:17:00 INFO Cycle is DONE; starting a new loop!
2026-06-22 18:17:00 INFO [TIME] step_end current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-22 18:17:00 INFO [TIME] step_start current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-22 18:17:00 INFO [TIME] window start=2020-02-06 13:00:00 end=2020-02-06 14:00:00 run_hours=1 has_assimilation=True
2026-06-22 18:17:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:17:01 INFO Hourly dataset computed and listing created
2026-06-22 18:17:03 INFO Hourly dataset computed
2026-06-22 18:17:03 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:17:05 INFO Hourly dataset computed and listing created
2026-06-22 18:17:05 INFO Hourly dataset computed
2026-06-22 18:17:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:17:06 INFO Hourly dataset computed and listing created
2026-06-22 18:17:06 INFO Hourly dataset computed
2026-06-22 18:17:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:17:07 INFO Hourly dataset computed and listing created
2026-06-22 18:17:08 INFO Hourly dataset computed
2026-06-22 18:17:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:17:09 INFO Hourly dataset computed and listing created
2026-06-22 18:17:09 INFO Hourly dataset computed
2026-06-22 18:17:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:17:10 INFO Hourly dataset computed and listing created
2026-06-22 18:17:11 INFO Hourly dataset computed
2026-06-22 18:17:11 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:17:12 INFO Hourly dataset computed and listing created
2026-06-22 18:17:12 INFO Hourly dataset computed
2026-06-22 18:17:12 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:17:13 INFO Hourly dataset computed and listing created
2026-06-22 18:17:14 INFO Hourly dataset computed
2026-06-22 18:17:14 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:17:15 INFO Hourly dataset computed and listing created
2026-06-22 18:17:15 INFO Hourly dataset computed
2026-06-22 18:17:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:17:16 INFO Hourly dataset computed and listing created
2026-06-22 18:17:17 INFO Hourly dataset computed
2026-06-22 18:17:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:17:18 INFO Hourly dataset computed and listing created
2026-06-22 18:17:18 INFO Hourly dataset computed
2026-06-22 18:17:18 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:17:19 INFO Hourly dataset computed and listing created
2026-06-22 18:17:20 INFO Hourly dataset computed
2026-06-22 18:17:20 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:17:20 INFO Hourly dataset computed and listing created
2026-06-22 18:17:21 INFO Hourly dataset computed
2026-06-22 18:17:21 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:17:22 INFO Hourly dataset computed and listing created
2026-06-22 18:17:22 INFO Hourly dataset computed
2026-06-22 18:17:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:17:23 INFO Hourly dataset computed and listing created
2026-06-22 18:17:24 INFO Hourly dataset computed
2026-06-22 18:17:24 INFO ---------->>> Running CHIMERE model from 2020-02-06 13:00:00 to 2020-02-06 14:00:00
2026-06-22 18:17:24 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:17:24 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1
2026-06-22 18:17:24 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020611_2_ENS1.nc
2026-06-22 18:17:24 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-06-22 18:17:24 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:17:24 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-06-22 18:17:24 INFO Queuing job for member 1...
2026-06-22 18:17:24 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:17:24 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-06-22 18:17:25 INFO Found: ['4952234']
2026-06-22 18:17:30 INFO [TGCC-IRENE] Submitted job with ID:['4952234']
2026-06-22 18:17:30 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:17:30 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2
2026-06-22 18:17:30 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020611_2_ENS2.nc
2026-06-22 18:17:30 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-06-22 18:17:30 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:17:30 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-06-22 18:17:30 INFO Queuing job for member 2...
2026-06-22 18:17:30 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:17:30 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-06-22 18:17:31 INFO Found: ['4952237']
2026-06-22 18:17:36 INFO [TGCC-IRENE] Submitted job with ID:['4952237']
2026-06-22 18:17:36 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:17:36 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3
2026-06-22 18:17:36 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020611_2_ENS3.nc
2026-06-22 18:17:36 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-06-22 18:17:36 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:17:36 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-06-22 18:17:36 INFO Queuing job for member 3...
2026-06-22 18:17:36 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:17:36 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-06-22 18:17:37 INFO Found: ['4952240']
2026-06-22 18:17:42 INFO [TGCC-IRENE] Submitted job with ID:['4952240']
2026-06-22 18:17:42 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:17:42 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4
2026-06-22 18:17:42 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020611_2_ENS4.nc
2026-06-22 18:17:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-06-22 18:17:42 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:17:42 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-06-22 18:17:42 INFO Queuing job for member 4...
2026-06-22 18:17:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:17:42 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-06-22 18:17:42 INFO Found: ['4952241']
2026-06-22 18:17:47 INFO [TGCC-IRENE] Submitted job with ID:['4952241']
2026-06-22 18:17:47 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:17:47 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5
2026-06-22 18:17:47 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020611_2_ENS5.nc
2026-06-22 18:17:47 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-06-22 18:17:47 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:17:47 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-06-22 18:17:47 INFO Queuing job for member 5...
2026-06-22 18:17:47 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:17:47 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-06-22 18:17:48 INFO Found: ['4952244']
2026-06-22 18:17:53 INFO [TGCC-IRENE] Submitted job with ID:['4952244']
2026-06-22 18:17:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:17:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6
2026-06-22 18:17:53 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020611_2_ENS6.nc
2026-06-22 18:17:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-06-22 18:17:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:17:53 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-06-22 18:17:53 INFO Queuing job for member 6...
2026-06-22 18:17:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:17:53 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-06-22 18:17:54 INFO Found: ['4952248']
2026-06-22 18:17:59 INFO [TGCC-IRENE] Submitted job with ID:['4952248']
2026-06-22 18:17:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:17:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7
2026-06-22 18:17:59 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020611_2_ENS7.nc
2026-06-22 18:17:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-06-22 18:17:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:17:59 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-06-22 18:17:59 INFO Queuing job for member 7...
2026-06-22 18:17:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:17:59 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-06-22 18:18:00 INFO Found: ['4952250']
2026-06-22 18:18:05 INFO [TGCC-IRENE] Submitted job with ID:['4952250']
2026-06-22 18:18:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:18:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8
2026-06-22 18:18:05 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020611_2_ENS8.nc
2026-06-22 18:18:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-06-22 18:18:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:18:05 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-06-22 18:18:05 INFO Queuing job for member 8...
2026-06-22 18:18:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:18:05 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-06-22 18:18:05 INFO Found: ['4952252']
2026-06-22 18:18:10 INFO [TGCC-IRENE] Submitted job with ID:['4952252']
2026-06-22 18:18:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:18:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9
2026-06-22 18:18:10 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020611_2_ENS9.nc
2026-06-22 18:18:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-06-22 18:18:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:18:10 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-06-22 18:18:10 INFO Queuing job for member 9...
2026-06-22 18:18:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:18:10 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-06-22 18:18:11 INFO Found: ['4952253']
2026-06-22 18:18:16 INFO [TGCC-IRENE] Submitted job with ID:['4952253']
2026-06-22 18:18:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:18:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10
2026-06-22 18:18:16 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020611_2_ENS10.nc
2026-06-22 18:18:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-06-22 18:18:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:18:16 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-06-22 18:18:16 INFO Queuing job for member 10...
2026-06-22 18:18:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:18:16 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-06-22 18:18:17 INFO Found: ['4952255']
2026-06-22 18:18:22 INFO [TGCC-IRENE] Submitted job with ID:['4952255']
2026-06-22 18:18:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:18:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11
2026-06-22 18:18:22 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020611_2_ENS11.nc
2026-06-22 18:18:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-06-22 18:18:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:18:22 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-06-22 18:18:22 INFO Queuing job for member 11...
2026-06-22 18:18:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:18:22 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-06-22 18:18:23 INFO Found: ['4952259']
2026-06-22 18:18:28 INFO [TGCC-IRENE] Submitted job with ID:['4952259']
2026-06-22 18:18:28 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:18:28 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12
2026-06-22 18:18:28 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020611_2_ENS12.nc
2026-06-22 18:18:28 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-06-22 18:18:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:18:28 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-06-22 18:18:28 INFO Queuing job for member 12...
2026-06-22 18:18:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:18:28 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-06-22 18:18:28 INFO Found: ['4952261']
2026-06-22 18:18:33 INFO [TGCC-IRENE] Submitted job with ID:['4952261']
2026-06-22 18:18:33 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:18:33 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13
2026-06-22 18:18:33 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020611_2_ENS13.nc
2026-06-22 18:18:33 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-06-22 18:18:33 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:18:33 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-06-22 18:18:33 INFO Queuing job for member 13...
2026-06-22 18:18:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:18:33 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-06-22 18:18:34 INFO Found: ['4952266']
2026-06-22 18:18:39 INFO [TGCC-IRENE] Submitted job with ID:['4952266']
2026-06-22 18:18:39 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:18:39 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14
2026-06-22 18:18:39 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020611_2_ENS14.nc
2026-06-22 18:18:39 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-06-22 18:18:39 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:18:39 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-06-22 18:18:39 INFO Queuing job for member 14...
2026-06-22 18:18:39 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:18:39 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-06-22 18:18:40 INFO Found: ['4952270']
2026-06-22 18:18:45 INFO [TGCC-IRENE] Submitted job with ID:['4952270']
2026-06-22 18:18:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:18:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15
2026-06-22 18:18:45 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020611_2_ENS15.nc
2026-06-22 18:18:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-06-22 18:18:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:18:45 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-06-22 18:18:45 INFO Queuing job for member 15...
2026-06-22 18:18:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:18:45 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-06-22 18:18:46 INFO Found: ['4952273']
2026-06-22 18:18:51 INFO [TGCC-IRENE] Submitted job with ID:['4952273']
2026-06-22 18:18:51 INFO Checking job status ...
2026-06-22 18:18:51 INFO None 4952234: status RUNNING/PENDING
2026-06-22 18:18:51 INFO None 4952237: status RUNNING/PENDING
2026-06-22 18:18:51 INFO None 4952240: status RUNNING/PENDING
2026-06-22 18:18:51 INFO None 4952241: status RUNNING/PENDING
2026-06-22 18:18:51 INFO None 4952244: status RUNNING/PENDING
2026-06-22 18:18:51 INFO None 4952248: status RUNNING/PENDING
2026-06-22 18:18:51 INFO None 4952250: status RUNNING/PENDING
2026-06-22 18:18:51 INFO None 4952252: status RUNNING/PENDING
2026-06-22 18:18:51 INFO None 4952253: status RUNNING/PENDING
2026-06-22 18:18:51 INFO None 4952255: status RUNNING/PENDING
2026-06-22 18:18:51 INFO None 4952259: status RUNNING/PENDING
2026-06-22 18:18:51 INFO None 4952261: status RUNNING/PENDING
2026-06-22 18:18:51 INFO None 4952266: status RUNNING/PENDING
2026-06-22 18:18:51 INFO None 4952270: status RUNNING/PENDING
2026-06-22 18:18:51 INFO None 4952273: status RUNNING/PENDING
2026-06-22 18:18:51 INFO Jobs still running: ['4952234', '4952237', '4952240', '4952241', '4952244', '4952248', '4952250', '4952252', '4952253', '4952255', '4952259', '4952261', '4952266', '4952270', '4952273']. Waiting...
2026-06-22 18:19:07 INFO None 4952234: status RUNNING/PENDING
2026-06-22 18:19:07 INFO None 4952237: status RUNNING/PENDING
2026-06-22 18:19:07 INFO None 4952240: status RUNNING/PENDING
2026-06-22 18:19:07 INFO None 4952241: status RUNNING/PENDING
2026-06-22 18:19:07 INFO None 4952244: status RUNNING/PENDING
2026-06-22 18:19:07 INFO None 4952248: status RUNNING/PENDING
2026-06-22 18:19:07 INFO None 4952250: status RUNNING/PENDING
2026-06-22 18:19:07 INFO None 4952252: status RUNNING/PENDING
2026-06-22 18:19:07 INFO None 4952253: status RUNNING/PENDING
2026-06-22 18:19:07 INFO None 4952255: status RUNNING/PENDING
2026-06-22 18:19:07 INFO None 4952259: status RUNNING/PENDING
2026-06-22 18:19:07 INFO None 4952261: status RUNNING/PENDING
2026-06-22 18:19:07 INFO None 4952266: status RUNNING/PENDING
2026-06-22 18:19:07 INFO None 4952270: status RUNNING/PENDING
2026-06-22 18:19:07 INFO None 4952273: status RUNNING/PENDING
2026-06-22 18:19:07 INFO Jobs still running: ['4952234', '4952237', '4952240', '4952241', '4952244', '4952248', '4952250', '4952252', '4952253', '4952255', '4952259', '4952261', '4952266', '4952270', '4952273']. Waiting...
2026-06-22 18:19:23 INFO None 4952234: status RUNNING/PENDING
2026-06-22 18:19:23 INFO None 4952237: status RUNNING/PENDING
2026-06-22 18:19:23 INFO None 4952240: status RUNNING/PENDING
2026-06-22 18:19:23 INFO None 4952241: status RUNNING/PENDING
2026-06-22 18:19:23 INFO None 4952244: status RUNNING/PENDING
2026-06-22 18:19:23 INFO None 4952248: status RUNNING/PENDING
2026-06-22 18:19:23 INFO None 4952250: status RUNNING/PENDING
2026-06-22 18:19:23 INFO None 4952252: status RUNNING/PENDING
2026-06-22 18:19:23 INFO None 4952253: status RUNNING/PENDING
2026-06-22 18:19:23 INFO None 4952255: status RUNNING/PENDING
2026-06-22 18:19:23 INFO None 4952259: status RUNNING/PENDING
2026-06-22 18:19:23 INFO None 4952261: status RUNNING/PENDING
2026-06-22 18:19:23 INFO None 4952266: status RUNNING/PENDING
2026-06-22 18:19:23 INFO None 4952270: status RUNNING/PENDING
2026-06-22 18:19:23 INFO None 4952273: status RUNNING/PENDING
2026-06-22 18:19:23 INFO Jobs still running: ['4952234', '4952237', '4952240', '4952241', '4952244', '4952248', '4952250', '4952252', '4952253', '4952255', '4952259', '4952261', '4952266', '4952270', '4952273']. Waiting...
2026-06-22 18:19:38 INFO None 4952234: status RUNNING/PENDING
2026-06-22 18:19:38 INFO None 4952237: status RUNNING/PENDING
2026-06-22 18:19:38 INFO None 4952240: status RUNNING/PENDING
2026-06-22 18:19:38 INFO None 4952241: status RUNNING/PENDING
2026-06-22 18:19:38 INFO None 4952244: status RUNNING/PENDING
2026-06-22 18:19:38 INFO None 4952248: status RUNNING/PENDING
2026-06-22 18:19:38 INFO None 4952250: status RUNNING/PENDING
2026-06-22 18:19:38 INFO None 4952252: status RUNNING/PENDING
2026-06-22 18:19:38 INFO None 4952253: status RUNNING/PENDING
2026-06-22 18:19:38 INFO None 4952255: status RUNNING/PENDING
2026-06-22 18:19:38 INFO None 4952259: status RUNNING/PENDING
2026-06-22 18:19:38 INFO None 4952261: status RUNNING/PENDING
2026-06-22 18:19:38 INFO None 4952266: status RUNNING/PENDING
2026-06-22 18:19:38 INFO None 4952270: status RUNNING/PENDING
2026-06-22 18:19:38 INFO None 4952273: status RUNNING/PENDING
2026-06-22 18:19:38 INFO Jobs still running: ['4952234', '4952237', '4952240', '4952241', '4952244', '4952248', '4952250', '4952252', '4952253', '4952255', '4952259', '4952261', '4952266', '4952270', '4952273']. Waiting...
2026-06-22 18:19:53 INFO None 4952234: status RUNNING/PENDING
2026-06-22 18:19:53 INFO None 4952237: status RUNNING/PENDING
2026-06-22 18:19:53 INFO None 4952240: status RUNNING/PENDING
2026-06-22 18:19:53 INFO None 4952241: status RUNNING/PENDING
2026-06-22 18:19:53 INFO None 4952244: status RUNNING/PENDING
2026-06-22 18:19:53 INFO None 4952248: status RUNNING/PENDING
2026-06-22 18:19:53 INFO None 4952250: status RUNNING/PENDING
2026-06-22 18:19:53 INFO None 4952252: status RUNNING/PENDING
2026-06-22 18:19:55 INFO None 4952253: status RUNNING/PENDING
2026-06-22 18:19:55 INFO None 4952255: status FINISHED
2026-06-22 18:19:55 INFO None 4952259: status RUNNING/PENDING
2026-06-22 18:19:55 INFO None 4952261: status RUNNING/PENDING
2026-06-22 18:19:55 INFO None 4952266: status RUNNING/PENDING
2026-06-22 18:19:55 INFO None 4952270: status RUNNING/PENDING
2026-06-22 18:19:55 INFO None 4952273: status RUNNING/PENDING
2026-06-22 18:19:55 INFO Jobs still running: ['4952234', '4952237', '4952240', '4952241', '4952244', '4952248', '4952250', '4952252', '4952253', '4952259', '4952261', '4952266', '4952270', '4952273']. Waiting...
2026-06-22 18:20:11 INFO None 4952234: status RUNNING/PENDING
2026-06-22 18:20:11 INFO None 4952237: status RUNNING/PENDING
2026-06-22 18:20:11 INFO None 4952240: status RUNNING/PENDING
2026-06-22 18:20:11 INFO None 4952241: status RUNNING/PENDING
2026-06-22 18:20:11 INFO None 4952244: status RUNNING/PENDING
2026-06-22 18:20:11 INFO None 4952248: status RUNNING/PENDING
2026-06-22 18:20:11 INFO None 4952250: status FINISHED
2026-06-22 18:20:11 INFO None 4952252: status RUNNING/PENDING
2026-06-22 18:20:11 INFO None 4952253: status RUNNING/PENDING
2026-06-22 18:20:11 INFO None 4952255: status FINISHED
2026-06-22 18:20:11 INFO None 4952259: status RUNNING/PENDING
2026-06-22 18:20:11 INFO None 4952261: status RUNNING/PENDING
2026-06-22 18:20:11 INFO None 4952266: status RUNNING/PENDING
2026-06-22 18:20:11 INFO None 4952270: status RUNNING/PENDING
2026-06-22 18:20:11 INFO None 4952273: status RUNNING/PENDING
2026-06-22 18:20:11 INFO Jobs still running: ['4952234', '4952237', '4952240', '4952241', '4952244', '4952248', '4952252', '4952253', '4952259', '4952261', '4952266', '4952270', '4952273']. Waiting...
2026-06-22 18:20:26 INFO None 4952234: status RUNNING/PENDING
2026-06-22 18:20:26 INFO None 4952237: status RUNNING/PENDING
2026-06-22 18:20:26 INFO None 4952240: status RUNNING/PENDING
2026-06-22 18:20:26 INFO None 4952241: status FINISHED
2026-06-22 18:20:26 INFO None 4952244: status RUNNING/PENDING
2026-06-22 18:20:26 INFO None 4952248: status RUNNING/PENDING
2026-06-22 18:20:26 INFO None 4952250: status FINISHED
2026-06-22 18:20:26 INFO None 4952252: status FINISHED
2026-06-22 18:20:26 INFO None 4952253: status RUNNING/PENDING
2026-06-22 18:20:26 INFO None 4952255: status FINISHED
2026-06-22 18:20:26 INFO None 4952259: status RUNNING/PENDING
2026-06-22 18:20:26 INFO None 4952261: status RUNNING/PENDING
2026-06-22 18:20:26 INFO None 4952266: status RUNNING/PENDING
2026-06-22 18:20:26 INFO None 4952270: status RUNNING/PENDING
2026-06-22 18:20:26 INFO None 4952273: status RUNNING/PENDING
2026-06-22 18:20:26 INFO Jobs still running: ['4952234', '4952237', '4952240', '4952244', '4952248', '4952253', '4952259', '4952261', '4952266', '4952270', '4952273']. Waiting...
2026-06-22 18:20:41 INFO None 4952234: status RUNNING/PENDING
2026-06-22 18:20:41 INFO None 4952237: status RUNNING/PENDING
2026-06-22 18:20:41 INFO None 4952240: status RUNNING/PENDING
2026-06-22 18:20:41 INFO None 4952241: status FINISHED
2026-06-22 18:20:41 INFO None 4952244: status RUNNING/PENDING
2026-06-22 18:20:41 INFO None 4952248: status RUNNING/PENDING
2026-06-22 18:20:41 INFO None 4952250: status FINISHED
2026-06-22 18:20:41 INFO None 4952252: status FINISHED
2026-06-22 18:20:41 INFO None 4952253: status FINISHED
2026-06-22 18:20:41 INFO None 4952255: status FINISHED
2026-06-22 18:20:41 INFO None 4952259: status RUNNING/PENDING
2026-06-22 18:20:41 INFO None 4952261: status RUNNING/PENDING
2026-06-22 18:20:41 INFO None 4952266: status RUNNING/PENDING
2026-06-22 18:20:41 INFO None 4952270: status RUNNING/PENDING
2026-06-22 18:20:41 INFO None 4952273: status RUNNING/PENDING
2026-06-22 18:20:41 INFO Jobs still running: ['4952234', '4952237', '4952240', '4952244', '4952248', '4952259', '4952261', '4952266', '4952270', '4952273']. Waiting...
2026-06-22 18:20:56 INFO None 4952234: status RUNNING/PENDING
2026-06-22 18:20:56 INFO None 4952237: status RUNNING/PENDING
2026-06-22 18:20:56 INFO None 4952240: status FINISHED
2026-06-22 18:20:56 INFO None 4952241: status FINISHED
2026-06-22 18:20:56 INFO None 4952244: status FINISHED
2026-06-22 18:20:57 INFO None 4952248: status RUNNING/PENDING
2026-06-22 18:20:57 INFO None 4952250: status FINISHED
2026-06-22 18:20:57 INFO None 4952252: status FINISHED
2026-06-22 18:20:57 INFO None 4952253: status FINISHED
2026-06-22 18:20:57 INFO None 4952255: status FINISHED
2026-06-22 18:20:57 INFO None 4952259: status RUNNING/PENDING
2026-06-22 18:20:57 INFO None 4952261: status FINISHED
2026-06-22 18:20:57 INFO None 4952266: status RUNNING/PENDING
2026-06-22 18:20:57 INFO None 4952270: status RUNNING/PENDING
2026-06-22 18:20:57 INFO None 4952273: status RUNNING/PENDING
2026-06-22 18:20:57 INFO Jobs still running: ['4952234', '4952237', '4952248', '4952259', '4952266', '4952270', '4952273']. Waiting...
2026-06-22 18:21:12 INFO None 4952234: status RUNNING/PENDING
2026-06-22 18:21:12 INFO None 4952237: status FINISHED
2026-06-22 18:21:12 INFO None 4952240: status FINISHED
2026-06-22 18:21:12 INFO None 4952241: status FINISHED
2026-06-22 18:21:12 INFO None 4952244: status FINISHED
2026-06-22 18:21:12 INFO None 4952248: status FINISHED
2026-06-22 18:21:12 INFO None 4952250: status FINISHED
2026-06-22 18:21:12 INFO None 4952252: status FINISHED
2026-06-22 18:21:12 INFO None 4952253: status FINISHED
2026-06-22 18:21:12 INFO None 4952255: status FINISHED
2026-06-22 18:21:12 INFO None 4952259: status RUNNING/PENDING
2026-06-22 18:21:12 INFO None 4952261: status FINISHED
2026-06-22 18:21:12 INFO None 4952266: status RUNNING/PENDING
2026-06-22 18:21:12 INFO None 4952270: status RUNNING/PENDING
2026-06-22 18:21:12 INFO None 4952273: status RUNNING/PENDING
2026-06-22 18:21:12 INFO Jobs still running: ['4952234', '4952259', '4952266', '4952270', '4952273']. Waiting...
2026-06-22 18:21:27 INFO None 4952234: status RUNNING/PENDING
2026-06-22 18:21:27 INFO None 4952237: status FINISHED
2026-06-22 18:21:27 INFO None 4952240: status FINISHED
2026-06-22 18:21:27 INFO None 4952241: status FINISHED
2026-06-22 18:21:27 INFO None 4952244: status FINISHED
2026-06-22 18:21:27 INFO None 4952248: status FINISHED
2026-06-22 18:21:27 INFO None 4952250: status FINISHED
2026-06-22 18:21:27 INFO None 4952252: status FINISHED
2026-06-22 18:21:27 INFO None 4952253: status FINISHED
2026-06-22 18:21:27 INFO None 4952255: status FINISHED
2026-06-22 18:21:27 INFO None 4952259: status FINISHED
2026-06-22 18:21:27 INFO None 4952261: status FINISHED
2026-06-22 18:21:27 INFO None 4952266: status FINISHED
2026-06-22 18:21:27 INFO None 4952270: status FINISHED
2026-06-22 18:21:27 INFO None 4952273: status RUNNING/PENDING
2026-06-22 18:21:27 INFO Jobs still running: ['4952234', '4952273']. Waiting...
2026-06-22 18:21:42 INFO None 4952234: status RUNNING/PENDING
2026-06-22 18:21:42 INFO None 4952237: status FINISHED
2026-06-22 18:21:42 INFO None 4952240: status FINISHED
2026-06-22 18:21:42 INFO None 4952241: status FINISHED
2026-06-22 18:21:42 INFO None 4952244: status FINISHED
2026-06-22 18:21:42 INFO None 4952248: status FINISHED
2026-06-22 18:21:42 INFO None 4952250: status FINISHED
2026-06-22 18:21:42 INFO None 4952252: status FINISHED
2026-06-22 18:21:42 INFO None 4952253: status FINISHED
2026-06-22 18:21:42 INFO None 4952255: status FINISHED
2026-06-22 18:21:43 INFO None 4952259: status FINISHED
2026-06-22 18:21:43 INFO None 4952261: status FINISHED
2026-06-22 18:21:43 INFO None 4952266: status FINISHED
2026-06-22 18:21:43 INFO None 4952270: status FINISHED
2026-06-22 18:21:43 INFO None 4952273: status FINISHED
2026-06-22 18:21:43 INFO Jobs still running: ['4952234']. Waiting...
2026-06-22 18:21:58 INFO None 4952234: status RUNNING/PENDING
2026-06-22 18:21:58 INFO None 4952237: status FINISHED
2026-06-22 18:21:58 INFO None 4952240: status FINISHED
2026-06-22 18:21:58 INFO None 4952241: status FINISHED
2026-06-22 18:21:58 INFO None 4952244: status FINISHED
2026-06-22 18:21:58 INFO None 4952248: status FINISHED
2026-06-22 18:21:58 INFO None 4952250: status FINISHED
2026-06-22 18:21:58 INFO None 4952252: status FINISHED
2026-06-22 18:21:58 INFO None 4952253: status FINISHED
2026-06-22 18:21:58 INFO None 4952255: status FINISHED
2026-06-22 18:21:58 INFO None 4952259: status FINISHED
2026-06-22 18:21:58 INFO None 4952261: status FINISHED
2026-06-22 18:21:58 INFO None 4952266: status FINISHED
2026-06-22 18:21:58 INFO None 4952270: status FINISHED
2026-06-22 18:21:58 INFO None 4952273: status FINISHED
2026-06-22 18:21:58 INFO Jobs still running: ['4952234']. Waiting...
2026-06-22 18:22:13 INFO None 4952234: status RUNNING/PENDING
2026-06-22 18:22:13 INFO None 4952237: status FINISHED
2026-06-22 18:22:13 INFO None 4952240: status FINISHED
2026-06-22 18:22:13 INFO None 4952241: status FINISHED
2026-06-22 18:22:13 INFO None 4952244: status FINISHED
2026-06-22 18:22:13 INFO None 4952248: status FINISHED
2026-06-22 18:22:13 INFO None 4952250: status FINISHED
2026-06-22 18:22:13 INFO None 4952252: status FINISHED
2026-06-22 18:22:13 INFO None 4952253: status FINISHED
2026-06-22 18:22:13 INFO None 4952255: status FINISHED
2026-06-22 18:22:13 INFO None 4952259: status FINISHED
2026-06-22 18:22:13 INFO None 4952261: status FINISHED
2026-06-22 18:22:13 INFO None 4952266: status FINISHED
2026-06-22 18:22:13 INFO None 4952270: status FINISHED
2026-06-22 18:22:13 INFO None 4952273: status FINISHED
2026-06-22 18:22:13 INFO Jobs still running: ['4952234']. Waiting...
2026-06-22 18:22:28 INFO None 4952234: status FINISHED
2026-06-22 18:22:28 INFO None 4952237: status FINISHED
2026-06-22 18:22:28 INFO None 4952240: status FINISHED
2026-06-22 18:22:28 INFO None 4952241: status FINISHED
2026-06-22 18:22:28 INFO None 4952244: status FINISHED
2026-06-22 18:22:28 INFO None 4952248: status FINISHED
2026-06-22 18:22:28 INFO None 4952250: status FINISHED
2026-06-22 18:22:28 INFO None 4952252: status FINISHED
2026-06-22 18:22:28 INFO None 4952253: status FINISHED
2026-06-22 18:22:28 INFO None 4952255: status FINISHED
2026-06-22 18:22:28 INFO None 4952259: status FINISHED
2026-06-22 18:22:28 INFO None 4952261: status FINISHED
2026-06-22 18:22:28 INFO None 4952266: status FINISHED
2026-06-22 18:22:28 INFO None 4952270: status FINISHED
2026-06-22 18:22:28 INFO None 4952273: status FINISHED
2026-06-22 18:22:28 INFO Jobs ['4952234', '4952237', '4952240', '4952241', '4952244', '4952248', '4952250', '4952252', '4952253', '4952255', '4952259', '4952261', '4952266', '4952270', '4952273'] have finished
2026-06-22 18:22:28 INFO Checking restart files were created ...
2026-06-22 18:22:28 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020613_1_ENS1.nc(668832435 bytes)
2026-06-22 18:22:28 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020613_1_ENS2.nc(668832435 bytes)
2026-06-22 18:22:28 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020613_1_ENS3.nc(668832435 bytes)
2026-06-22 18:22:28 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020613_1_ENS4.nc(668832435 bytes)
2026-06-22 18:22:28 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020613_1_ENS5.nc(668832435 bytes)
2026-06-22 18:22:28 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020613_1_ENS6.nc(668832435 bytes)
2026-06-22 18:22:28 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020613_1_ENS7.nc(668832435 bytes)
2026-06-22 18:22:28 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020613_1_ENS8.nc(668832435 bytes)
2026-06-22 18:22:28 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020613_1_ENS9.nc(668832435 bytes)
2026-06-22 18:22:28 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020613_1_ENS10.nc(668832435 bytes)
2026-06-22 18:22:28 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020613_1_ENS11.nc(668832435 bytes)
2026-06-22 18:22:28 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020613_1_ENS12.nc(668832435 bytes)
2026-06-22 18:22:28 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020613_1_ENS13.nc(668832435 bytes)
2026-06-22 18:22:28 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020613_1_ENS14.nc(668832435 bytes)
2026-06-22 18:22:28 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020613_1_ENS15.nc(668832435 bytes)
2026-06-22 18:22:28 INFO  Run_model() completed successfully.
2026-06-22 18:22:28 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-22 18:22:28 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 14:00:00 days=153072 seconds=50400
2026-06-22 18:22:28 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-06-22 18:22:28 INFO [TIME] increment current_time 2020-02-06 13:00:00 -> 2020-02-06 14:00:00
2026-06-22 18:22:28 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-22 18:22:28 INFO ---------->>> Running process_satellite_data()
2026-06-22 18:22:29 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12005.nc
2026-06-22 18:22:29 INFO ---------->>> Running run_obs_converter()
2026-06-22 18:22:29 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_51225_153072.out
2026-06-22 18:22:29 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_51225_153072.out
2026-06-22 18:22:29 INFO ---------->>> Running DART
2026-06-22 18:22:29 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-06-22 18:22:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/chim_ENS1_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/chim_ENS1_2020020614_1_out_toDART.nc
2026-06-22 18:22:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/chim_ENS2_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/chim_ENS2_2020020614_1_out_toDART.nc
2026-06-22 18:22:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/chim_ENS3_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/chim_ENS3_2020020614_1_out_toDART.nc
2026-06-22 18:22:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/chim_ENS4_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/chim_ENS4_2020020614_1_out_toDART.nc
2026-06-22 18:22:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/chim_ENS5_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/chim_ENS5_2020020614_1_out_toDART.nc
2026-06-22 18:22:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/chim_ENS6_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/chim_ENS6_2020020614_1_out_toDART.nc
2026-06-22 18:22:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/chim_ENS7_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/chim_ENS7_2020020614_1_out_toDART.nc
2026-06-22 18:22:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/chim_ENS8_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/chim_ENS8_2020020614_1_out_toDART.nc
2026-06-22 18:22:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/chim_ENS9_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/chim_ENS9_2020020614_1_out_toDART.nc
2026-06-22 18:22:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/chim_ENS10_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/chim_ENS10_2020020614_1_out_toDART.nc
2026-06-22 18:22:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/chim_ENS11_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/chim_ENS11_2020020614_1_out_toDART.nc
2026-06-22 18:22:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/chim_ENS12_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/chim_ENS12_2020020614_1_out_toDART.nc
2026-06-22 18:22:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/chim_ENS13_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/chim_ENS13_2020020614_1_out_toDART.nc
2026-06-22 18:22:33 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/chim_ENS14_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/chim_ENS14_2020020614_1_out_toDART.nc
2026-06-22 18:22:33 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/chim_ENS15_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/chim_ENS15_2020020614_1_out_toDART.nc
2026-06-22 18:22:33 INFO Replacement input_template.nml → input.nml completed successfully.
2026-06-22 18:22:33 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-06-22 18:22:33 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-06-22 18:22:33 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-06-22 18:22:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-06-22 18:22:33 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-06-22 18:22:41 INFO Found: []
2026-06-22 18:22:41 INFO No job id returned by command ./run_filter.bsh
2026-06-22 18:22:41 INFO No monitoring will be performed
2026-06-22 18:22:41 INFO Moving DART output files to analysis and preassim directories for date 2020020614 if present ...
2026-06-22 18:22:41 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-06-22 18:22:41 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-06-22 18:22:41 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-06-22 18:22:42 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-06-22 18:22:42 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-06-22 18:22:42 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-06-22 18:22:42 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-06-22 18:22:42 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-06-22 18:22:42 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-06-22 18:22:42 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-06-22 18:22:42 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-06-22 18:22:42 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-06-22 18:22:42 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-06-22 18:22:42 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-06-22 18:22:42 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-06-22 18:22:42 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-06-22 18:22:42 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-06-22 18:22:42 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-06-22 18:22:42 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-06-22 18:22:42 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-06-22 18:22:42 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-06-22 18:22:42 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-06-22 18:22:42 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-06-22 18:22:42 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-06-22 18:22:42 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-06-22 18:22:42 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-06-22 18:22:42 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-06-22 18:22:42 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-06-22 18:22:42 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-06-22 18:22:42 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-06-22 18:22:42 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-06-22 18:22:42 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/analysis/2020020614'
2026-06-22 18:22:42 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-06-22 18:22:42 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/preassim/2020020614'
2026-06-22 18:22:42 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-06-22 18:22:42 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-06-22 18:22:42 INFO run_dart() is DONE.
2026-06-22 18:22:42 INFO ---------->>> Running update_pollutant_in_end()
2026-06-22 18:22:42 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:22:42 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:22:43 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:22:43 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:22:44 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020613_1_ENS1.nc
2026-06-22 18:22:47 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:22:47 INFO [Posterior Diff ENS1] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS1_2020020614.nc
2026-06-22 18:22:47 INFO [Posterior Diff ENS1] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS1_2020020614.relative.nc
2026-06-22 18:22:48 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:22:48 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:22:49 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:22:49 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:22:49 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020613_1_ENS2.nc
2026-06-22 18:22:53 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:22:53 INFO [Posterior Diff ENS2] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS2_2020020614.nc
2026-06-22 18:22:53 INFO [Posterior Diff ENS2] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS2_2020020614.relative.nc
2026-06-22 18:22:53 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:22:53 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:22:54 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:22:54 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:22:55 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020613_1_ENS3.nc
2026-06-22 18:22:59 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:22:59 INFO [Posterior Diff ENS3] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS3_2020020614.nc
2026-06-22 18:22:59 INFO [Posterior Diff ENS3] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS3_2020020614.relative.nc
2026-06-22 18:22:59 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:22:59 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:23:00 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:23:00 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:23:01 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020613_1_ENS4.nc
2026-06-22 18:23:05 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:23:05 INFO [Posterior Diff ENS4] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS4_2020020614.nc
2026-06-22 18:23:05 INFO [Posterior Diff ENS4] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS4_2020020614.relative.nc
2026-06-22 18:23:05 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:23:05 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:23:06 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:23:06 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:23:07 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020613_1_ENS5.nc
2026-06-22 18:23:10 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:23:11 INFO [Posterior Diff ENS5] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS5_2020020614.nc
2026-06-22 18:23:11 INFO [Posterior Diff ENS5] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS5_2020020614.relative.nc
2026-06-22 18:23:11 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:23:11 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:23:12 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:23:12 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:23:13 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020613_1_ENS6.nc
2026-06-22 18:23:16 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:23:16 INFO [Posterior Diff ENS6] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS6_2020020614.nc
2026-06-22 18:23:16 INFO [Posterior Diff ENS6] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS6_2020020614.relative.nc
2026-06-22 18:23:17 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:23:17 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:23:18 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:23:18 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:23:18 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020613_1_ENS7.nc
2026-06-22 18:23:22 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:23:22 INFO [Posterior Diff ENS7] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS7_2020020614.nc
2026-06-22 18:23:22 INFO [Posterior Diff ENS7] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS7_2020020614.relative.nc
2026-06-22 18:23:22 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:23:22 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:23:23 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:23:23 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:23:24 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020613_1_ENS8.nc
2026-06-22 18:23:27 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:23:28 INFO [Posterior Diff ENS8] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS8_2020020614.nc
2026-06-22 18:23:28 INFO [Posterior Diff ENS8] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS8_2020020614.relative.nc
2026-06-22 18:23:28 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:23:28 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:23:29 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:23:29 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:23:30 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020613_1_ENS9.nc
2026-06-22 18:23:33 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:23:34 INFO [Posterior Diff ENS9] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS9_2020020614.nc
2026-06-22 18:23:34 INFO [Posterior Diff ENS9] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS9_2020020614.relative.nc
2026-06-22 18:23:34 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:23:34 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:23:35 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:23:35 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:23:36 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020613_1_ENS10.nc
2026-06-22 18:23:39 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:23:39 INFO [Posterior Diff ENS10] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS10_2020020614.nc
2026-06-22 18:23:39 INFO [Posterior Diff ENS10] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS10_2020020614.relative.nc
2026-06-22 18:23:40 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:23:40 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:23:41 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:23:41 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:23:41 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020613_1_ENS11.nc
2026-06-22 18:23:45 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:23:45 INFO [Posterior Diff ENS11] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS11_2020020614.nc
2026-06-22 18:23:45 INFO [Posterior Diff ENS11] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS11_2020020614.relative.nc
2026-06-22 18:23:46 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:23:46 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:23:47 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:23:47 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:23:47 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020613_1_ENS12.nc
2026-06-22 18:23:51 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:23:51 INFO [Posterior Diff ENS12] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS12_2020020614.nc
2026-06-22 18:23:51 INFO [Posterior Diff ENS12] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS12_2020020614.relative.nc
2026-06-22 18:23:51 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:23:51 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:23:52 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:23:52 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:23:53 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020613_1_ENS13.nc
2026-06-22 18:23:57 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:23:57 INFO [Posterior Diff ENS13] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS13_2020020614.nc
2026-06-22 18:23:57 INFO [Posterior Diff ENS13] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS13_2020020614.relative.nc
2026-06-22 18:23:57 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:23:57 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:23:58 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:23:58 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:23:59 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020613_1_ENS14.nc
2026-06-22 18:24:02 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:24:03 INFO [Posterior Diff ENS14] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS14_2020020614.nc
2026-06-22 18:24:03 INFO [Posterior Diff ENS14] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS14_2020020614.relative.nc
2026-06-22 18:24:03 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:24:03 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:24:04 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 18:24:04 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 18:24:05 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020613_1_ENS15.nc
2026-06-22 18:24:08 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 18:24:09 INFO [Posterior Diff ENS15] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS15_2020020614.nc
2026-06-22 18:24:09 INFO [Posterior Diff ENS15] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmpemis_0615_15m_low_v2/posteriors/2020020614/diff_posterior_ENS15_2020020614.relative.nc
2026-06-22 18:24:09 INFO Next run starts from 2020-02-06 14:00:00
2026-06-22 18:24:09 INFO Cycle is DONE; starting a new loop!
2026-06-22 18:24:09 INFO [TIME] step_end current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-22 18:24:09 INFO [TIME] step_start current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-22 18:24:09 INFO [TIME] window start=2020-02-06 14:00:00 end=2020-02-07 00:00:00 run_hours=10 has_assimilation=False
2026-06-22 18:24:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:24:10 INFO Hourly dataset computed and listing created
2026-06-22 18:24:23 INFO Hourly dataset computed
2026-06-22 18:24:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:24:25 INFO Hourly dataset computed and listing created
2026-06-22 18:24:28 INFO Hourly dataset computed
2026-06-22 18:24:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:24:29 INFO Hourly dataset computed and listing created
2026-06-22 18:24:32 INFO Hourly dataset computed
2026-06-22 18:24:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:24:33 INFO Hourly dataset computed and listing created
2026-06-22 18:24:35 INFO Hourly dataset computed
2026-06-22 18:24:35 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:24:36 INFO Hourly dataset computed and listing created
2026-06-22 18:24:39 INFO Hourly dataset computed
2026-06-22 18:24:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:24:40 INFO Hourly dataset computed and listing created
2026-06-22 18:24:43 INFO Hourly dataset computed
2026-06-22 18:24:43 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:24:44 INFO Hourly dataset computed and listing created
2026-06-22 18:24:46 INFO Hourly dataset computed
2026-06-22 18:24:46 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:24:47 INFO Hourly dataset computed and listing created
2026-06-22 18:24:49 INFO Hourly dataset computed
2026-06-22 18:24:49 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:24:51 INFO Hourly dataset computed and listing created
2026-06-22 18:24:53 INFO Hourly dataset computed
2026-06-22 18:24:53 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:24:54 INFO Hourly dataset computed and listing created
2026-06-22 18:24:57 INFO Hourly dataset computed
2026-06-22 18:24:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:24:58 INFO Hourly dataset computed and listing created
2026-06-22 18:25:00 INFO Hourly dataset computed
2026-06-22 18:25:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:25:02 INFO Hourly dataset computed and listing created
2026-06-22 18:25:04 INFO Hourly dataset computed
2026-06-22 18:25:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:25:05 INFO Hourly dataset computed and listing created
2026-06-22 18:25:08 INFO Hourly dataset computed
2026-06-22 18:25:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:25:09 INFO Hourly dataset computed and listing created
2026-06-22 18:25:11 INFO Hourly dataset computed
2026-06-22 18:25:12 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-22 18:25:13 INFO Hourly dataset computed and listing created
2026-06-22 18:25:15 INFO Hourly dataset computed
2026-06-22 18:25:15 INFO ---------->>> Running CHIMERE model from 2020-02-06 14:00:00 to 2020-02-07 00:00:00
2026-06-22 18:25:15 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:25:15 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1
2026-06-22 18:25:15 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020613_1_ENS1.nc
2026-06-22 18:25:15 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-06-22 18:25:15 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:25:15 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-06-22 18:25:15 INFO Queuing job for member 1...
2026-06-22 18:25:15 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:25:15 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-06-22 18:25:16 INFO Found: ['4952423']
2026-06-22 18:25:21 INFO [TGCC-IRENE] Submitted job with ID:['4952423']
2026-06-22 18:25:21 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:25:21 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2
2026-06-22 18:25:21 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020613_1_ENS2.nc
2026-06-22 18:25:21 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-06-22 18:25:21 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:25:21 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-06-22 18:25:21 INFO Queuing job for member 2...
2026-06-22 18:25:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:25:21 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-06-22 18:25:22 INFO Found: ['4952425']
2026-06-22 18:25:27 INFO [TGCC-IRENE] Submitted job with ID:['4952425']
2026-06-22 18:25:27 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:25:27 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3
2026-06-22 18:25:27 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020613_1_ENS3.nc
2026-06-22 18:25:27 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-06-22 18:25:27 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:25:27 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-06-22 18:25:27 INFO Queuing job for member 3...
2026-06-22 18:25:27 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:25:27 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-06-22 18:25:27 INFO Found: ['4952433']
2026-06-22 18:25:32 INFO [TGCC-IRENE] Submitted job with ID:['4952433']
2026-06-22 18:25:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:25:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4
2026-06-22 18:25:32 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020613_1_ENS4.nc
2026-06-22 18:25:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-06-22 18:25:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:25:32 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-06-22 18:25:32 INFO Queuing job for member 4...
2026-06-22 18:25:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:25:32 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-06-22 18:25:33 INFO Found: ['4952436']
2026-06-22 18:25:38 INFO [TGCC-IRENE] Submitted job with ID:['4952436']
2026-06-22 18:25:38 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:25:38 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5
2026-06-22 18:25:38 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020613_1_ENS5.nc
2026-06-22 18:25:38 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-06-22 18:25:38 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:25:38 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-06-22 18:25:38 INFO Queuing job for member 5...
2026-06-22 18:25:38 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:25:38 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-06-22 18:25:39 INFO Found: ['4952442']
2026-06-22 18:25:44 INFO [TGCC-IRENE] Submitted job with ID:['4952442']
2026-06-22 18:25:44 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:25:44 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6
2026-06-22 18:25:44 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020613_1_ENS6.nc
2026-06-22 18:25:44 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-06-22 18:25:44 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:25:44 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-06-22 18:25:44 INFO Queuing job for member 6...
2026-06-22 18:25:44 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:25:44 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-06-22 18:25:45 INFO Found: ['4952449']
2026-06-22 18:25:50 INFO [TGCC-IRENE] Submitted job with ID:['4952449']
2026-06-22 18:25:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:25:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7
2026-06-22 18:25:50 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020613_1_ENS7.nc
2026-06-22 18:25:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-06-22 18:25:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:25:50 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-06-22 18:25:50 INFO Queuing job for member 7...
2026-06-22 18:25:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:25:50 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-06-22 18:25:50 INFO Found: ['4952458']
2026-06-22 18:25:55 INFO [TGCC-IRENE] Submitted job with ID:['4952458']
2026-06-22 18:25:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:25:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8
2026-06-22 18:25:55 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020613_1_ENS8.nc
2026-06-22 18:25:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-06-22 18:25:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:25:55 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-06-22 18:25:55 INFO Queuing job for member 8...
2026-06-22 18:25:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:25:55 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-06-22 18:25:56 INFO Found: ['4952466']
2026-06-22 18:26:01 INFO [TGCC-IRENE] Submitted job with ID:['4952466']
2026-06-22 18:26:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:26:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9
2026-06-22 18:26:01 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020613_1_ENS9.nc
2026-06-22 18:26:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-06-22 18:26:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:26:01 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-06-22 18:26:01 INFO Queuing job for member 9...
2026-06-22 18:26:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:26:01 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-06-22 18:26:02 INFO Found: ['4952472']
2026-06-22 18:26:07 INFO [TGCC-IRENE] Submitted job with ID:['4952472']
2026-06-22 18:26:07 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:26:07 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10
2026-06-22 18:26:07 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020613_1_ENS10.nc
2026-06-22 18:26:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-06-22 18:26:07 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:26:07 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-06-22 18:26:07 INFO Queuing job for member 10...
2026-06-22 18:26:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:26:07 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-06-22 18:26:08 INFO Found: ['4952481']
2026-06-22 18:26:13 INFO [TGCC-IRENE] Submitted job with ID:['4952481']
2026-06-22 18:26:13 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:26:13 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11
2026-06-22 18:26:13 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020613_1_ENS11.nc
2026-06-22 18:26:13 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-06-22 18:26:13 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:26:13 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-06-22 18:26:13 INFO Queuing job for member 11...
2026-06-22 18:26:13 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:26:13 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-06-22 18:26:13 INFO Found: ['4952487']
2026-06-22 18:26:18 INFO [TGCC-IRENE] Submitted job with ID:['4952487']
2026-06-22 18:26:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:26:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12
2026-06-22 18:26:18 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020613_1_ENS12.nc
2026-06-22 18:26:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-06-22 18:26:18 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:26:18 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-06-22 18:26:18 INFO Queuing job for member 12...
2026-06-22 18:26:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:26:18 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-06-22 18:26:20 INFO Found: ['4952493']
2026-06-22 18:26:25 INFO [TGCC-IRENE] Submitted job with ID:['4952493']
2026-06-22 18:26:25 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:26:25 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13
2026-06-22 18:26:25 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020613_1_ENS13.nc
2026-06-22 18:26:25 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-06-22 18:26:25 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:26:25 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-06-22 18:26:25 INFO Queuing job for member 13...
2026-06-22 18:26:25 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:26:25 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-06-22 18:26:26 INFO Found: ['4952501']
2026-06-22 18:26:31 INFO [TGCC-IRENE] Submitted job with ID:['4952501']
2026-06-22 18:26:31 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:26:31 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14
2026-06-22 18:26:31 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020613_1_ENS14.nc
2026-06-22 18:26:31 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-06-22 18:26:31 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:26:31 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-06-22 18:26:31 INFO Queuing job for member 14...
2026-06-22 18:26:31 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:26:31 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-06-22 18:26:32 INFO Found: ['4952504']
2026-06-22 18:26:37 INFO [TGCC-IRENE] Submitted job with ID:['4952504']
2026-06-22 18:26:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 18:26:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15
2026-06-22 18:26:37 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020613_1_ENS15.nc
2026-06-22 18:26:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-06-22 18:26:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 18:26:37 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-06-22 18:26:37 INFO Queuing job for member 15...
2026-06-22 18:26:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 18:26:37 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-06-22 18:26:37 INFO Found: ['4952511']
2026-06-22 18:26:42 INFO [TGCC-IRENE] Submitted job with ID:['4952511']
2026-06-22 18:26:42 INFO Checking job status ...
2026-06-22 18:26:42 INFO None 4952423: status FINISHED
2026-06-22 18:26:42 INFO None 4952425: status RUNNING/PENDING
2026-06-22 18:26:42 INFO None 4952433: status RUNNING/PENDING
2026-06-22 18:26:42 INFO None 4952436: status RUNNING/PENDING
2026-06-22 18:26:42 INFO None 4952442: status RUNNING/PENDING
2026-06-22 18:26:42 INFO None 4952449: status RUNNING/PENDING
2026-06-22 18:26:42 INFO None 4952458: status RUNNING/PENDING
2026-06-22 18:26:43 INFO None 4952466: status RUNNING/PENDING
2026-06-22 18:26:43 INFO None 4952472: status RUNNING/PENDING
2026-06-22 18:26:43 INFO None 4952481: status RUNNING/PENDING
2026-06-22 18:26:43 INFO None 4952487: status RUNNING/PENDING
2026-06-22 18:26:43 INFO None 4952493: status RUNNING/PENDING
2026-06-22 18:26:43 INFO None 4952501: status RUNNING/PENDING
2026-06-22 18:26:43 INFO None 4952504: status RUNNING/PENDING
2026-06-22 18:26:43 INFO None 4952511: status RUNNING/PENDING
2026-06-22 18:26:43 INFO Jobs still running: ['4952425', '4952433', '4952436', '4952442', '4952449', '4952458', '4952466', '4952472', '4952481', '4952487', '4952493', '4952501', '4952504', '4952511']. Waiting...
2026-06-22 18:26:58 INFO None 4952423: status FINISHED
2026-06-22 18:26:58 INFO None 4952425: status FINISHED
2026-06-22 18:26:58 INFO None 4952433: status RUNNING/PENDING
2026-06-22 18:26:58 INFO None 4952436: status RUNNING/PENDING
2026-06-22 18:26:58 INFO None 4952442: status RUNNING/PENDING
2026-06-22 18:26:58 INFO None 4952449: status RUNNING/PENDING
2026-06-22 18:26:58 INFO None 4952458: status FINISHED
2026-06-22 18:26:58 INFO None 4952466: status FINISHED
2026-06-22 18:26:58 INFO None 4952472: status FINISHED
2026-06-22 18:26:58 INFO None 4952481: status RUNNING/PENDING
2026-06-22 18:26:58 INFO None 4952487: status RUNNING/PENDING
2026-06-22 18:26:58 INFO None 4952493: status RUNNING/PENDING
2026-06-22 18:26:58 INFO None 4952501: status RUNNING/PENDING
2026-06-22 18:26:58 INFO None 4952504: status RUNNING/PENDING
2026-06-22 18:26:58 INFO None 4952511: status RUNNING/PENDING
2026-06-22 18:26:58 INFO Jobs still running: ['4952433', '4952436', '4952442', '4952449', '4952481', '4952487', '4952493', '4952501', '4952504', '4952511']. Waiting...
2026-06-22 18:27:13 INFO None 4952423: status FINISHED
2026-06-22 18:27:13 INFO None 4952425: status FINISHED
2026-06-22 18:27:13 INFO None 4952433: status RUNNING/PENDING
2026-06-22 18:27:13 INFO None 4952436: status FINISHED
2026-06-22 18:27:13 INFO None 4952442: status RUNNING/PENDING
2026-06-22 18:27:13 INFO None 4952449: status RUNNING/PENDING
2026-06-22 18:27:13 INFO None 4952458: status FINISHED
2026-06-22 18:27:13 INFO None 4952466: status FINISHED
2026-06-22 18:27:13 INFO None 4952472: status FINISHED
2026-06-22 18:27:13 INFO None 4952481: status RUNNING/PENDING
2026-06-22 18:27:13 INFO None 4952487: status RUNNING/PENDING
2026-06-22 18:27:13 INFO None 4952493: status RUNNING/PENDING
2026-06-22 18:27:13 INFO None 4952501: status RUNNING/PENDING
2026-06-22 18:27:13 INFO None 4952504: status RUNNING/PENDING
2026-06-22 18:27:13 INFO None 4952511: status RUNNING/PENDING
2026-06-22 18:27:13 INFO Jobs still running: ['4952433', '4952442', '4952449', '4952481', '4952487', '4952493', '4952501', '4952504', '4952511']. Waiting...
2026-06-22 18:27:28 INFO None 4952423: status FINISHED
2026-06-22 18:27:28 INFO None 4952425: status FINISHED
2026-06-22 18:27:28 INFO None 4952433: status RUNNING/PENDING
2026-06-22 18:27:28 INFO None 4952436: status FINISHED
2026-06-22 18:27:28 INFO None 4952442: status FINISHED
2026-06-22 18:27:28 INFO None 4952449: status FINISHED
2026-06-22 18:27:28 INFO None 4952458: status FINISHED
2026-06-22 18:27:28 INFO None 4952466: status FINISHED
2026-06-22 18:27:28 INFO None 4952472: status FINISHED
2026-06-22 18:27:29 INFO None 4952481: status RUNNING/PENDING
2026-06-22 18:27:29 INFO None 4952487: status RUNNING/PENDING
2026-06-22 18:27:29 INFO None 4952493: status RUNNING/PENDING
2026-06-22 18:27:29 INFO None 4952501: status RUNNING/PENDING
2026-06-22 18:27:29 INFO None 4952504: status RUNNING/PENDING
2026-06-22 18:27:29 INFO None 4952511: status RUNNING/PENDING
2026-06-22 18:27:29 INFO Jobs still running: ['4952433', '4952481', '4952487', '4952493', '4952501', '4952504', '4952511']. Waiting...
2026-06-22 18:27:44 INFO None 4952423: status FINISHED
2026-06-22 18:27:44 INFO None 4952425: status FINISHED
2026-06-22 18:27:44 INFO None 4952433: status RUNNING/PENDING
2026-06-22 18:27:44 INFO None 4952436: status FINISHED
2026-06-22 18:27:44 INFO None 4952442: status FINISHED
2026-06-22 18:27:44 INFO None 4952449: status FINISHED
2026-06-22 18:27:44 INFO None 4952458: status FINISHED
2026-06-22 18:27:44 INFO None 4952466: status FINISHED
2026-06-22 18:27:44 INFO None 4952472: status FINISHED
2026-06-22 18:27:44 INFO None 4952481: status RUNNING/PENDING
2026-06-22 18:27:44 INFO None 4952487: status RUNNING/PENDING
2026-06-22 18:27:44 INFO None 4952493: status RUNNING/PENDING
2026-06-22 18:27:44 INFO None 4952501: status RUNNING/PENDING
2026-06-22 18:27:44 INFO None 4952504: status FINISHED
2026-06-22 18:27:44 INFO None 4952511: status RUNNING/PENDING
2026-06-22 18:27:44 INFO Jobs still running: ['4952433', '4952481', '4952487', '4952493', '4952501', '4952511']. Waiting...
2026-06-22 18:27:59 INFO None 4952423: status FINISHED
2026-06-22 18:27:59 INFO None 4952425: status FINISHED
2026-06-22 18:27:59 INFO None 4952433: status FINISHED
2026-06-22 18:27:59 INFO None 4952436: status FINISHED
2026-06-22 18:27:59 INFO None 4952442: status FINISHED
2026-06-22 18:27:59 INFO None 4952449: status FINISHED
2026-06-22 18:27:59 INFO None 4952458: status FINISHED
2026-06-22 18:27:59 INFO None 4952466: status FINISHED
2026-06-22 18:27:59 INFO None 4952472: status FINISHED
2026-06-22 18:27:59 INFO None 4952481: status RUNNING/PENDING
2026-06-22 18:27:59 INFO None 4952487: status FINISHED
2026-06-22 18:27:59 INFO None 4952493: status RUNNING/PENDING
2026-06-22 18:27:59 INFO None 4952501: status RUNNING/PENDING
2026-06-22 18:27:59 INFO None 4952504: status FINISHED
2026-06-22 18:27:59 INFO None 4952511: status RUNNING/PENDING
2026-06-22 18:27:59 INFO Jobs still running: ['4952481', '4952493', '4952501', '4952511']. Waiting...
2026-06-22 18:28:14 INFO None 4952423: status FINISHED
2026-06-22 18:28:14 INFO None 4952425: status FINISHED
2026-06-22 18:28:14 INFO None 4952433: status FINISHED
2026-06-22 18:28:14 INFO None 4952436: status FINISHED
2026-06-22 18:28:14 INFO None 4952442: status FINISHED
2026-06-22 18:28:14 INFO None 4952449: status FINISHED
2026-06-22 18:28:14 INFO None 4952458: status FINISHED
2026-06-22 18:28:14 INFO None 4952466: status FINISHED
2026-06-22 18:28:14 INFO None 4952472: status FINISHED
2026-06-22 18:28:14 INFO None 4952481: status RUNNING/PENDING
2026-06-22 18:28:15 INFO None 4952487: status FINISHED
2026-06-22 18:28:15 INFO None 4952493: status FINISHED
2026-06-22 18:28:15 INFO None 4952501: status RUNNING/PENDING
2026-06-22 18:28:15 INFO None 4952504: status FINISHED
2026-06-22 18:28:15 INFO None 4952511: status RUNNING/PENDING
2026-06-22 18:28:15 INFO Jobs still running: ['4952481', '4952501', '4952511']. Waiting...
2026-06-22 18:28:30 INFO None 4952423: status FINISHED
2026-06-22 18:28:30 INFO None 4952425: status FINISHED
2026-06-22 18:28:30 INFO None 4952433: status FINISHED
2026-06-22 18:28:30 INFO None 4952436: status FINISHED
2026-06-22 18:28:30 INFO None 4952442: status FINISHED
2026-06-22 18:28:30 INFO None 4952449: status FINISHED
2026-06-22 18:28:30 INFO None 4952458: status FINISHED
2026-06-22 18:28:30 INFO None 4952466: status FINISHED
2026-06-22 18:28:30 INFO None 4952472: status FINISHED
2026-06-22 18:28:30 INFO None 4952481: status FINISHED
2026-06-22 18:28:30 INFO None 4952487: status FINISHED
2026-06-22 18:28:30 INFO None 4952493: status FINISHED
2026-06-22 18:28:30 INFO None 4952501: status RUNNING/PENDING
2026-06-22 18:28:30 INFO None 4952504: status FINISHED
2026-06-22 18:28:30 INFO None 4952511: status RUNNING/PENDING
2026-06-22 18:28:30 INFO Jobs still running: ['4952501', '4952511']. Waiting...
2026-06-22 18:28:45 INFO None 4952423: status FINISHED
2026-06-22 18:28:45 INFO None 4952425: status FINISHED
2026-06-22 18:28:45 INFO None 4952433: status FINISHED
2026-06-22 18:28:45 INFO None 4952436: status FINISHED
2026-06-22 18:28:45 INFO None 4952442: status FINISHED
2026-06-22 18:28:45 INFO None 4952449: status FINISHED
2026-06-22 18:28:45 INFO None 4952458: status FINISHED
2026-06-22 18:28:45 INFO None 4952466: status FINISHED
2026-06-22 18:28:45 INFO None 4952472: status FINISHED
2026-06-22 18:28:45 INFO None 4952481: status FINISHED
2026-06-22 18:28:45 INFO None 4952487: status FINISHED
2026-06-22 18:28:45 INFO None 4952493: status FINISHED
2026-06-22 18:28:45 INFO None 4952501: status FINISHED
2026-06-22 18:28:45 INFO None 4952504: status FINISHED
2026-06-22 18:28:45 INFO None 4952511: status RUNNING/PENDING
2026-06-22 18:28:45 INFO Jobs still running: ['4952511']. Waiting...
2026-06-22 18:29:00 INFO None 4952423: status FINISHED
2026-06-22 18:29:00 INFO None 4952425: status FINISHED
2026-06-22 18:29:01 INFO None 4952433: status FINISHED
2026-06-22 18:29:01 INFO None 4952436: status FINISHED
2026-06-22 18:29:01 INFO None 4952442: status FINISHED
2026-06-22 18:29:01 INFO None 4952449: status FINISHED
2026-06-22 18:29:01 INFO None 4952458: status FINISHED
2026-06-22 18:29:01 INFO None 4952466: status FINISHED
2026-06-22 18:29:01 INFO None 4952472: status FINISHED
2026-06-22 18:29:01 INFO None 4952481: status FINISHED
2026-06-22 18:29:01 INFO None 4952487: status FINISHED
2026-06-22 18:29:01 INFO None 4952493: status FINISHED
2026-06-22 18:29:01 INFO None 4952501: status FINISHED
2026-06-22 18:29:01 INFO None 4952504: status FINISHED
2026-06-22 18:29:01 INFO None 4952511: status RUNNING/PENDING
2026-06-22 18:29:01 INFO Jobs still running: ['4952511']. Waiting...
2026-06-22 18:29:16 INFO None 4952423: status FINISHED
2026-06-22 18:29:16 INFO None 4952425: status FINISHED
2026-06-22 18:29:16 INFO None 4952433: status FINISHED
2026-06-22 18:29:16 INFO None 4952436: status FINISHED
2026-06-22 18:29:16 INFO None 4952442: status FINISHED
2026-06-22 18:29:16 INFO None 4952449: status FINISHED
2026-06-22 18:29:16 INFO None 4952458: status FINISHED
2026-06-22 18:29:16 INFO None 4952466: status FINISHED
2026-06-22 18:29:16 INFO None 4952472: status FINISHED
2026-06-22 18:29:16 INFO None 4952481: status FINISHED
2026-06-22 18:29:16 INFO None 4952487: status FINISHED
2026-06-22 18:29:16 INFO None 4952493: status FINISHED
2026-06-22 18:29:16 INFO None 4952501: status FINISHED
2026-06-22 18:29:16 INFO None 4952504: status FINISHED
2026-06-22 18:29:16 INFO None 4952511: status FINISHED
2026-06-22 18:29:16 INFO Jobs ['4952423', '4952425', '4952433', '4952436', '4952442', '4952449', '4952458', '4952466', '4952472', '4952481', '4952487', '4952493', '4952501', '4952504', '4952511'] have finished
2026-06-22 18:29:16 INFO Checking restart files were created ...
2026-06-22 18:29:16 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020614_10_ENS1.nc
2026-06-22 18:29:16 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020614_10_ENS2.nc
2026-06-22 18:29:16 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020614_10_ENS3.nc
2026-06-22 18:29:16 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020614_10_ENS4.nc
2026-06-22 18:29:16 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020614_10_ENS5.nc
2026-06-22 18:29:16 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020614_10_ENS6.nc
2026-06-22 18:29:16 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020614_10_ENS7.nc
2026-06-22 18:29:16 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020614_10_ENS8.nc
2026-06-22 18:29:16 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020614_10_ENS9.nc
2026-06-22 18:29:16 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020614_10_ENS10.nc
2026-06-22 18:29:16 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020614_10_ENS11.nc
2026-06-22 18:29:16 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020614_10_ENS12.nc
2026-06-22 18:29:16 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020614_10_ENS13.nc
2026-06-22 18:29:16 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020614_10_ENS14.nc
2026-06-22 18:29:16 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020614_10_ENS15.nc
2026-06-22 18:29:16 INFO Check chimere log file at: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/ENS15_2020020614.out
2026-06-22 18:29:16 ERROR [PIPELINE] Error: The following chimere ENS run(s) failed (exit code 1): [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
Traceback (most recent call last):
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/main.py", line 103, in <module>
    pipeline.run_pipeline()
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/pipelines/base_pipeline.py", line 138, in run_pipeline
    self.run_model()
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/pipelines/chimere2023/pipeline.py", line 327, in run_model
    raise ModelRunError(f"The following chimere ENS run(s) failed (exit code 1): {mems_to_rerun}")
pipeline_errors.ModelRunError: The following chimere ENS run(s) failed (exit code 1): [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
+ exit 0
