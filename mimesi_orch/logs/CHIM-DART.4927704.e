+ /bin/bash -x /tmp/tmp.dx50BjlJJu
+ SCRIPT_PID=569744
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
2026-06-20 16:26:58 INFO 
███    ███ ██ ███    ███ ███████ ███████ ██
████  ████ ██ ████  ████ ██      ██      ██
██ ████ ██ ██ ██ ████ ██ █████   ███████ ██
██  ██  ██ ██ ██  ██  ██ ██           ██ ██
██      ██ ██ ██      ██ ███████ ███████ ██



2026-06-20 16:26:58 INFO [PIPELINE] =======================================
2026-06-20 16:26:58 INFO [PIPELINE] Starting chimere–DART orchestrator
2026-06-20 16:26:58 INFO [PIPELINE] Config file: config/config_irene_IM.yaml
2026-06-20 16:26:58 INFO [PIPELINE] Run dir: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart
2026-06-20 16:26:58 INFO [PIPELINE] Log file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/mimesi_orchestrator_logs/chimere_DART_20260620_162658.log
2026-06-20 16:26:58 INFO [PIPELINE] =======================================
2026-06-20 16:26:59 INFO Running assimilation with model_type=ModelType.CHIMERE
2026-06-20 16:26:59 INFO Using scheduler=Scheduler.SLURM, queue=rome
2026-06-20 16:26:59 INFO [STEP] ---- TIME LOOP START ----
2026-06-20 16:26:59 INFO [TIME] step_start current_time=2020-02-06 00:00:00 simulated_time=None dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-20 16:26:59 INFO [TIME] window start=2020-02-06 00:00:00 end=2020-02-06 01:00:00 run_hours=1 has_assimilation=False
2026-06-20 16:26:59 INFO Creating directories and links for ENS1 to run chimere's parallel part
2026-06-20 16:26:59 INFO Linking EMIS ...
2026-06-20 16:26:59 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens01.nc
2026-06-20 16:26:59 INFO Linking END ...
2026-06-20 16:26:59 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS1/end.2020020500_24_ENS1.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-20 16:26:59 INFO >> Checking links...
2026-06-20 16:26:59 INFO >> All links are good for ENS1  ...
2026-06-20 16:26:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-20 16:27:07 INFO Hourly dataset computed and listing created
2026-06-20 16:27:11 INFO Hourly dataset computed
2026-06-20 16:27:11 INFO Creating directories and links for ENS2 to run chimere's parallel part
2026-06-20 16:27:11 INFO Linking EMIS ...
2026-06-20 16:27:11 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens02.nc
2026-06-20 16:27:11 INFO Linking END ...
2026-06-20 16:27:11 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS2/end.2020020500_24_ENS2.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-20 16:27:11 INFO >> Checking links...
2026-06-20 16:27:11 INFO >> All links are good for ENS2  ...
2026-06-20 16:27:11 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-20 16:27:12 INFO Hourly dataset computed and listing created
2026-06-20 16:27:13 INFO Hourly dataset computed
2026-06-20 16:27:13 INFO Creating directories and links for ENS3 to run chimere's parallel part
2026-06-20 16:27:13 INFO Linking EMIS ...
2026-06-20 16:27:13 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens03.nc
2026-06-20 16:27:13 INFO Linking END ...
2026-06-20 16:27:13 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS3/end.2020020500_24_ENS3.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-20 16:27:13 INFO >> Checking links...
2026-06-20 16:27:13 INFO >> All links are good for ENS3  ...
2026-06-20 16:27:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-20 16:27:14 INFO Hourly dataset computed and listing created
2026-06-20 16:27:15 INFO Hourly dataset computed
2026-06-20 16:27:15 INFO Creating directories and links for ENS4 to run chimere's parallel part
2026-06-20 16:27:15 INFO Linking EMIS ...
2026-06-20 16:27:15 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens04.nc
2026-06-20 16:27:15 INFO Linking END ...
2026-06-20 16:27:15 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS4/end.2020020500_24_ENS4.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-20 16:27:15 INFO >> Checking links...
2026-06-20 16:27:15 INFO >> All links are good for ENS4  ...
2026-06-20 16:27:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-20 16:27:16 INFO Hourly dataset computed and listing created
2026-06-20 16:27:17 INFO Hourly dataset computed
2026-06-20 16:27:17 INFO Creating directories and links for ENS5 to run chimere's parallel part
2026-06-20 16:27:17 INFO Linking EMIS ...
2026-06-20 16:27:17 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens05.nc
2026-06-20 16:27:18 INFO Linking END ...
2026-06-20 16:27:18 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS5/end.2020020500_24_ENS5.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-20 16:27:18 INFO >> Checking links...
2026-06-20 16:27:18 INFO >> All links are good for ENS5  ...
2026-06-20 16:27:18 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-20 16:27:18 INFO Hourly dataset computed and listing created
2026-06-20 16:27:19 INFO Hourly dataset computed
2026-06-20 16:27:19 INFO Creating directories and links for ENS6 to run chimere's parallel part
2026-06-20 16:27:19 INFO Linking EMIS ...
2026-06-20 16:27:20 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens06.nc
2026-06-20 16:27:20 INFO Linking END ...
2026-06-20 16:27:20 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS6/end.2020020500_24_ENS6.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-20 16:27:20 INFO >> Checking links...
2026-06-20 16:27:20 INFO >> All links are good for ENS6  ...
2026-06-20 16:27:20 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-20 16:27:21 INFO Hourly dataset computed and listing created
2026-06-20 16:27:21 INFO Hourly dataset computed
2026-06-20 16:27:21 INFO Creating directories and links for ENS7 to run chimere's parallel part
2026-06-20 16:27:21 INFO Linking EMIS ...
2026-06-20 16:27:22 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens07.nc
2026-06-20 16:27:22 INFO Linking END ...
2026-06-20 16:27:22 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS7/end.2020020500_24_ENS7.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-20 16:27:22 INFO >> Checking links...
2026-06-20 16:27:22 INFO >> All links are good for ENS7  ...
2026-06-20 16:27:22 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-20 16:27:23 INFO Hourly dataset computed and listing created
2026-06-20 16:27:23 INFO Hourly dataset computed
2026-06-20 16:27:24 INFO Creating directories and links for ENS8 to run chimere's parallel part
2026-06-20 16:27:24 INFO Linking EMIS ...
2026-06-20 16:27:24 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens08.nc
2026-06-20 16:27:24 INFO Linking END ...
2026-06-20 16:27:24 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS8/end.2020020500_24_ENS8.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-20 16:27:24 INFO >> Checking links...
2026-06-20 16:27:24 INFO >> All links are good for ENS8  ...
2026-06-20 16:27:24 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-20 16:27:25 INFO Hourly dataset computed and listing created
2026-06-20 16:27:28 INFO Hourly dataset computed
2026-06-20 16:27:28 INFO Creating directories and links for ENS9 to run chimere's parallel part
2026-06-20 16:27:28 INFO Linking EMIS ...
2026-06-20 16:27:29 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens09.nc
2026-06-20 16:27:29 INFO Linking END ...
2026-06-20 16:27:29 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS9/end.2020020500_24_ENS9.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-20 16:27:29 INFO >> Checking links...
2026-06-20 16:27:29 INFO >> All links are good for ENS9  ...
2026-06-20 16:27:29 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-20 16:27:30 INFO Hourly dataset computed and listing created
2026-06-20 16:27:33 INFO Hourly dataset computed
2026-06-20 16:27:33 INFO Creating directories and links for ENS10 to run chimere's parallel part
2026-06-20 16:27:33 INFO Linking EMIS ...
2026-06-20 16:27:34 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens10.nc
2026-06-20 16:27:34 INFO Linking END ...
2026-06-20 16:27:34 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS10/end.2020020500_24_ENS10.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-20 16:27:34 INFO >> Checking links...
2026-06-20 16:27:34 INFO >> All links are good for ENS10  ...
2026-06-20 16:27:34 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-20 16:27:36 INFO Hourly dataset computed and listing created
2026-06-20 16:27:38 INFO Hourly dataset computed
2026-06-20 16:27:38 INFO Creating directories and links for ENS11 to run chimere's parallel part
2026-06-20 16:27:38 INFO Linking EMIS ...
2026-06-20 16:27:39 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens11.nc
2026-06-20 16:27:39 INFO Linking END ...
2026-06-20 16:27:39 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS11/end.2020020500_24_ENS11.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-20 16:27:39 INFO >> Checking links...
2026-06-20 16:27:39 INFO >> All links are good for ENS11  ...
2026-06-20 16:27:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-20 16:27:40 INFO Hourly dataset computed and listing created
2026-06-20 16:27:40 INFO Hourly dataset computed
2026-06-20 16:27:40 INFO Creating directories and links for ENS12 to run chimere's parallel part
2026-06-20 16:27:40 INFO Linking EMIS ...
2026-06-20 16:27:41 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens12.nc
2026-06-20 16:27:41 INFO Linking END ...
2026-06-20 16:27:41 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS12/end.2020020500_24_ENS12.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-20 16:27:41 INFO >> Checking links...
2026-06-20 16:27:41 INFO >> All links are good for ENS12  ...
2026-06-20 16:27:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-20 16:27:42 INFO Hourly dataset computed and listing created
2026-06-20 16:27:42 INFO Hourly dataset computed
2026-06-20 16:27:43 INFO Creating directories and links for ENS13 to run chimere's parallel part
2026-06-20 16:27:43 INFO Linking EMIS ...
2026-06-20 16:27:43 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens13.nc
2026-06-20 16:27:43 INFO Linking END ...
2026-06-20 16:27:43 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS13/end.2020020500_24_ENS13.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-20 16:27:43 INFO >> Checking links...
2026-06-20 16:27:43 INFO >> All links are good for ENS13  ...
2026-06-20 16:27:43 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-20 16:27:44 INFO Hourly dataset computed and listing created
2026-06-20 16:27:45 INFO Hourly dataset computed
2026-06-20 16:27:45 INFO Creating directories and links for ENS14 to run chimere's parallel part
2026-06-20 16:27:45 INFO Linking EMIS ...
2026-06-20 16:27:45 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens14.nc
2026-06-20 16:27:45 INFO Linking END ...
2026-06-20 16:27:45 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS14/end.2020020500_24_ENS14.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-20 16:27:45 INFO >> Checking links...
2026-06-20 16:27:45 INFO >> All links are good for ENS14  ...
2026-06-20 16:27:45 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-20 16:27:46 INFO Hourly dataset computed and listing created
2026-06-20 16:27:47 INFO Hourly dataset computed
2026-06-20 16:27:47 INFO Creating directories and links for ENS15 to run chimere's parallel part
2026-06-20 16:27:47 INFO Linking EMIS ...
2026-06-20 16:27:47 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens15.nc
2026-06-20 16:27:47 INFO Linking END ...
2026-06-20 16:27:47 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS15/end.2020020500_24_ENS15.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-20 16:27:47 INFO >> Checking links...
2026-06-20 16:27:47 INFO >> All links are good for ENS15  ...
2026-06-20 16:27:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-20 16:27:48 INFO Hourly dataset computed and listing created
2026-06-20 16:27:49 INFO Hourly dataset computed
2026-06-20 16:27:49 INFO ---------->>> Running CHIMERE model from 2020-02-06 00:00:00 to 2020-02-06 01:00:00
2026-06-20 16:27:49 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-20 16:27:49 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS1
2026-06-20 16:27:49 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS1/end.2020020500_24_ENS1.nc
2026-06-20 16:27:49 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-06-20 16:27:49 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-20 16:27:49 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-06-20 16:27:49 INFO Queuing job for member 1...
2026-06-20 16:27:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-20 16:27:49 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-06-20 16:27:50 INFO Found: ['4927708']
2026-06-20 16:27:55 INFO [TGCC-IRENE] Submitted job with ID:['4927708']
2026-06-20 16:27:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-20 16:27:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS2
2026-06-20 16:27:55 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS2/end.2020020500_24_ENS2.nc
2026-06-20 16:27:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-06-20 16:27:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-20 16:27:55 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-06-20 16:27:55 INFO Queuing job for member 2...
2026-06-20 16:27:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-20 16:27:55 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-06-20 16:27:56 INFO Found: ['4927709']
2026-06-20 16:28:01 INFO [TGCC-IRENE] Submitted job with ID:['4927709']
2026-06-20 16:28:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-20 16:28:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS3
2026-06-20 16:28:01 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS3/end.2020020500_24_ENS3.nc
2026-06-20 16:28:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-06-20 16:28:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-20 16:28:01 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-06-20 16:28:01 INFO Queuing job for member 3...
2026-06-20 16:28:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-20 16:28:01 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-06-20 16:28:01 INFO Found: ['4927710']
2026-06-20 16:28:06 INFO [TGCC-IRENE] Submitted job with ID:['4927710']
2026-06-20 16:28:06 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-20 16:28:06 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS4
2026-06-20 16:28:06 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS4/end.2020020500_24_ENS4.nc
2026-06-20 16:28:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-06-20 16:28:07 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-20 16:28:07 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-06-20 16:28:07 INFO Queuing job for member 4...
2026-06-20 16:28:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-20 16:28:07 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-06-20 16:28:07 INFO Found: ['4927712']
2026-06-20 16:28:12 INFO [TGCC-IRENE] Submitted job with ID:['4927712']
2026-06-20 16:28:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-20 16:28:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS5
2026-06-20 16:28:12 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS5/end.2020020500_24_ENS5.nc
2026-06-20 16:28:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-06-20 16:28:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-20 16:28:12 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-06-20 16:28:12 INFO Queuing job for member 5...
2026-06-20 16:28:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-20 16:28:12 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-06-20 16:28:13 INFO Found: ['4927713']
2026-06-20 16:28:18 INFO [TGCC-IRENE] Submitted job with ID:['4927713']
2026-06-20 16:28:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-20 16:28:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS6
2026-06-20 16:28:18 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS6/end.2020020500_24_ENS6.nc
2026-06-20 16:28:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-06-20 16:28:18 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-20 16:28:18 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-06-20 16:28:18 INFO Queuing job for member 6...
2026-06-20 16:28:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-20 16:28:18 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-06-20 16:28:19 INFO Found: ['4927715']
2026-06-20 16:28:24 INFO [TGCC-IRENE] Submitted job with ID:['4927715']
2026-06-20 16:28:24 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-20 16:28:24 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS7
2026-06-20 16:28:24 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS7/end.2020020500_24_ENS7.nc
2026-06-20 16:28:24 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-06-20 16:28:24 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-20 16:28:24 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-06-20 16:28:24 INFO Queuing job for member 7...
2026-06-20 16:28:24 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-20 16:28:24 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-06-20 16:28:25 INFO Found: ['4927716']
2026-06-20 16:28:30 INFO [TGCC-IRENE] Submitted job with ID:['4927716']
2026-06-20 16:28:30 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-20 16:28:30 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS8
2026-06-20 16:28:30 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS8/end.2020020500_24_ENS8.nc
2026-06-20 16:28:30 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-06-20 16:28:30 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-20 16:28:30 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-06-20 16:28:30 INFO Queuing job for member 8...
2026-06-20 16:28:30 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-20 16:28:30 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-06-20 16:28:31 INFO Found: ['4927719']
2026-06-20 16:28:36 INFO [TGCC-IRENE] Submitted job with ID:['4927719']
2026-06-20 16:28:36 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-20 16:28:36 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS9
2026-06-20 16:28:36 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS9/end.2020020500_24_ENS9.nc
2026-06-20 16:28:36 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-06-20 16:28:36 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-20 16:28:36 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-06-20 16:28:36 INFO Queuing job for member 9...
2026-06-20 16:28:36 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-20 16:28:36 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-06-20 16:28:37 INFO Found: ['4927720']
2026-06-20 16:28:42 INFO [TGCC-IRENE] Submitted job with ID:['4927720']
2026-06-20 16:28:42 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-20 16:28:42 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS10
2026-06-20 16:28:42 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS10/end.2020020500_24_ENS10.nc
2026-06-20 16:28:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-06-20 16:28:42 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-20 16:28:42 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-06-20 16:28:42 INFO Queuing job for member 10...
2026-06-20 16:28:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-20 16:28:42 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-06-20 16:28:42 INFO Found: ['4927721']
2026-06-20 16:28:47 INFO [TGCC-IRENE] Submitted job with ID:['4927721']
2026-06-20 16:28:47 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-20 16:28:47 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS11
2026-06-20 16:28:47 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS11/end.2020020500_24_ENS11.nc
2026-06-20 16:28:47 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-06-20 16:28:47 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-20 16:28:47 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-06-20 16:28:47 INFO Queuing job for member 11...
2026-06-20 16:28:47 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-20 16:28:47 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-06-20 16:28:48 INFO Found: ['4927722']
2026-06-20 16:28:53 INFO [TGCC-IRENE] Submitted job with ID:['4927722']
2026-06-20 16:28:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-20 16:28:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS12
2026-06-20 16:28:53 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS12/end.2020020500_24_ENS12.nc
2026-06-20 16:28:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-06-20 16:28:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-20 16:28:53 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-06-20 16:28:53 INFO Queuing job for member 12...
2026-06-20 16:28:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-20 16:28:53 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-06-20 16:28:54 INFO Found: ['4927723']
2026-06-20 16:28:59 INFO [TGCC-IRENE] Submitted job with ID:['4927723']
2026-06-20 16:28:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-20 16:28:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS13
2026-06-20 16:28:59 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS13/end.2020020500_24_ENS13.nc
2026-06-20 16:28:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-06-20 16:28:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-20 16:28:59 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-06-20 16:28:59 INFO Queuing job for member 13...
2026-06-20 16:28:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-20 16:28:59 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-06-20 16:29:00 INFO Found: ['4927724']
2026-06-20 16:29:05 INFO [TGCC-IRENE] Submitted job with ID:['4927724']
2026-06-20 16:29:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-20 16:29:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS14
2026-06-20 16:29:05 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS14/end.2020020500_24_ENS14.nc
2026-06-20 16:29:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-06-20 16:29:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-20 16:29:05 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-06-20 16:29:05 INFO Queuing job for member 14...
2026-06-20 16:29:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-20 16:29:05 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-06-20 16:29:06 INFO Found: ['4927726']
2026-06-20 16:29:11 INFO [TGCC-IRENE] Submitted job with ID:['4927726']
2026-06-20 16:29:11 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-20 16:29:11 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS15
2026-06-20 16:29:11 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS15/end.2020020500_24_ENS15.nc
2026-06-20 16:29:11 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-06-20 16:29:11 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-20 16:29:11 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-06-20 16:29:11 INFO Queuing job for member 15...
2026-06-20 16:29:11 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-20 16:29:11 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-06-20 16:29:12 INFO Found: ['4927727']
2026-06-20 16:29:17 INFO [TGCC-IRENE] Submitted job with ID:['4927727']
2026-06-20 16:29:17 INFO Checking job status ...
2026-06-20 16:29:17 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:29:17 INFO None 4927709: status RUNNING/PENDING
2026-06-20 16:29:17 INFO None 4927710: status RUNNING/PENDING
2026-06-20 16:29:17 INFO None 4927712: status RUNNING/PENDING
2026-06-20 16:29:17 INFO None 4927713: status RUNNING/PENDING
2026-06-20 16:29:17 INFO None 4927715: status RUNNING/PENDING
2026-06-20 16:29:17 INFO None 4927716: status RUNNING/PENDING
2026-06-20 16:29:17 INFO None 4927719: status RUNNING/PENDING
2026-06-20 16:29:17 INFO None 4927720: status RUNNING/PENDING
2026-06-20 16:29:17 INFO None 4927721: status RUNNING/PENDING
2026-06-20 16:29:17 INFO None 4927722: status RUNNING/PENDING
2026-06-20 16:29:17 INFO None 4927723: status RUNNING/PENDING
2026-06-20 16:29:17 INFO None 4927724: status RUNNING/PENDING
2026-06-20 16:29:17 INFO None 4927726: status RUNNING/PENDING
2026-06-20 16:29:17 INFO None 4927727: status RUNNING/PENDING
2026-06-20 16:29:17 INFO Jobs still running: ['4927708', '4927709', '4927710', '4927712', '4927713', '4927715', '4927716', '4927719', '4927720', '4927721', '4927722', '4927723', '4927724', '4927726', '4927727']. Waiting...
2026-06-20 16:29:32 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:29:32 INFO None 4927709: status RUNNING/PENDING
2026-06-20 16:29:32 INFO None 4927710: status RUNNING/PENDING
2026-06-20 16:29:32 INFO None 4927712: status RUNNING/PENDING
2026-06-20 16:29:32 INFO None 4927713: status RUNNING/PENDING
2026-06-20 16:29:32 INFO None 4927715: status RUNNING/PENDING
2026-06-20 16:29:32 INFO None 4927716: status RUNNING/PENDING
2026-06-20 16:29:32 INFO None 4927719: status RUNNING/PENDING
2026-06-20 16:29:32 INFO None 4927720: status RUNNING/PENDING
2026-06-20 16:29:32 INFO None 4927721: status RUNNING/PENDING
2026-06-20 16:29:32 INFO None 4927722: status RUNNING/PENDING
2026-06-20 16:29:32 INFO None 4927723: status RUNNING/PENDING
2026-06-20 16:29:32 INFO None 4927724: status RUNNING/PENDING
2026-06-20 16:29:32 INFO None 4927726: status RUNNING/PENDING
2026-06-20 16:29:32 INFO None 4927727: status RUNNING/PENDING
2026-06-20 16:29:32 INFO Jobs still running: ['4927708', '4927709', '4927710', '4927712', '4927713', '4927715', '4927716', '4927719', '4927720', '4927721', '4927722', '4927723', '4927724', '4927726', '4927727']. Waiting...
2026-06-20 16:29:47 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:29:47 INFO None 4927709: status RUNNING/PENDING
2026-06-20 16:29:47 INFO None 4927710: status RUNNING/PENDING
2026-06-20 16:29:47 INFO None 4927712: status RUNNING/PENDING
2026-06-20 16:29:47 INFO None 4927713: status RUNNING/PENDING
2026-06-20 16:29:47 INFO None 4927715: status RUNNING/PENDING
2026-06-20 16:29:47 INFO None 4927716: status RUNNING/PENDING
2026-06-20 16:29:47 INFO None 4927719: status RUNNING/PENDING
2026-06-20 16:29:47 INFO None 4927720: status RUNNING/PENDING
2026-06-20 16:29:47 INFO None 4927721: status RUNNING/PENDING
2026-06-20 16:29:48 INFO None 4927722: status RUNNING/PENDING
2026-06-20 16:29:48 INFO None 4927723: status RUNNING/PENDING
2026-06-20 16:29:48 INFO None 4927724: status RUNNING/PENDING
2026-06-20 16:29:48 INFO None 4927726: status RUNNING/PENDING
2026-06-20 16:29:48 INFO None 4927727: status RUNNING/PENDING
2026-06-20 16:29:48 INFO Jobs still running: ['4927708', '4927709', '4927710', '4927712', '4927713', '4927715', '4927716', '4927719', '4927720', '4927721', '4927722', '4927723', '4927724', '4927726', '4927727']. Waiting...
2026-06-20 16:30:03 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:30:03 INFO None 4927709: status RUNNING/PENDING
2026-06-20 16:30:03 INFO None 4927710: status RUNNING/PENDING
2026-06-20 16:30:03 INFO None 4927712: status RUNNING/PENDING
2026-06-20 16:30:03 INFO None 4927713: status RUNNING/PENDING
2026-06-20 16:30:03 INFO None 4927715: status RUNNING/PENDING
2026-06-20 16:30:03 INFO None 4927716: status RUNNING/PENDING
2026-06-20 16:30:03 INFO None 4927719: status RUNNING/PENDING
2026-06-20 16:30:03 INFO None 4927720: status RUNNING/PENDING
2026-06-20 16:30:03 INFO None 4927721: status RUNNING/PENDING
2026-06-20 16:30:03 INFO None 4927722: status RUNNING/PENDING
2026-06-20 16:30:03 INFO None 4927723: status RUNNING/PENDING
2026-06-20 16:30:03 INFO None 4927724: status RUNNING/PENDING
2026-06-20 16:30:03 INFO None 4927726: status RUNNING/PENDING
2026-06-20 16:30:03 INFO None 4927727: status RUNNING/PENDING
2026-06-20 16:30:03 INFO Jobs still running: ['4927708', '4927709', '4927710', '4927712', '4927713', '4927715', '4927716', '4927719', '4927720', '4927721', '4927722', '4927723', '4927724', '4927726', '4927727']. Waiting...
2026-06-20 16:30:18 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:30:18 INFO None 4927709: status RUNNING/PENDING
2026-06-20 16:30:18 INFO None 4927710: status RUNNING/PENDING
2026-06-20 16:30:18 INFO None 4927712: status RUNNING/PENDING
2026-06-20 16:30:18 INFO None 4927713: status RUNNING/PENDING
2026-06-20 16:30:18 INFO None 4927715: status RUNNING/PENDING
2026-06-20 16:30:18 INFO None 4927716: status RUNNING/PENDING
2026-06-20 16:30:18 INFO None 4927719: status RUNNING/PENDING
2026-06-20 16:30:18 INFO None 4927720: status RUNNING/PENDING
2026-06-20 16:30:18 INFO None 4927721: status RUNNING/PENDING
2026-06-20 16:30:18 INFO None 4927722: status RUNNING/PENDING
2026-06-20 16:30:18 INFO None 4927723: status RUNNING/PENDING
2026-06-20 16:30:18 INFO None 4927724: status RUNNING/PENDING
2026-06-20 16:30:18 INFO None 4927726: status RUNNING/PENDING
2026-06-20 16:30:18 INFO None 4927727: status RUNNING/PENDING
2026-06-20 16:30:18 INFO Jobs still running: ['4927708', '4927709', '4927710', '4927712', '4927713', '4927715', '4927716', '4927719', '4927720', '4927721', '4927722', '4927723', '4927724', '4927726', '4927727']. Waiting...
2026-06-20 16:30:33 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:30:33 INFO None 4927709: status RUNNING/PENDING
2026-06-20 16:30:33 INFO None 4927710: status RUNNING/PENDING
2026-06-20 16:30:33 INFO None 4927712: status FINISHED
2026-06-20 16:30:33 INFO None 4927713: status RUNNING/PENDING
2026-06-20 16:30:33 INFO None 4927715: status RUNNING/PENDING
2026-06-20 16:30:33 INFO None 4927716: status RUNNING/PENDING
2026-06-20 16:30:33 INFO None 4927719: status RUNNING/PENDING
2026-06-20 16:30:33 INFO None 4927720: status RUNNING/PENDING
2026-06-20 16:30:33 INFO None 4927721: status RUNNING/PENDING
2026-06-20 16:30:34 INFO None 4927722: status RUNNING/PENDING
2026-06-20 16:30:34 INFO None 4927723: status RUNNING/PENDING
2026-06-20 16:30:34 INFO None 4927724: status RUNNING/PENDING
2026-06-20 16:30:34 INFO None 4927726: status RUNNING/PENDING
2026-06-20 16:30:34 INFO None 4927727: status RUNNING/PENDING
2026-06-20 16:30:34 INFO Jobs still running: ['4927708', '4927709', '4927710', '4927713', '4927715', '4927716', '4927719', '4927720', '4927721', '4927722', '4927723', '4927724', '4927726', '4927727']. Waiting...
2026-06-20 16:30:49 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:30:49 INFO None 4927709: status RUNNING/PENDING
2026-06-20 16:30:49 INFO None 4927710: status RUNNING/PENDING
2026-06-20 16:30:49 INFO None 4927712: status FINISHED
2026-06-20 16:30:49 INFO None 4927713: status FINISHED
2026-06-20 16:30:49 INFO None 4927715: status RUNNING/PENDING
2026-06-20 16:30:49 INFO None 4927716: status RUNNING/PENDING
2026-06-20 16:30:49 INFO None 4927719: status RUNNING/PENDING
2026-06-20 16:30:49 INFO None 4927720: status RUNNING/PENDING
2026-06-20 16:30:49 INFO None 4927721: status RUNNING/PENDING
2026-06-20 16:30:49 INFO None 4927722: status RUNNING/PENDING
2026-06-20 16:30:49 INFO None 4927723: status RUNNING/PENDING
2026-06-20 16:30:49 INFO None 4927724: status RUNNING/PENDING
2026-06-20 16:30:49 INFO None 4927726: status RUNNING/PENDING
2026-06-20 16:30:49 INFO None 4927727: status RUNNING/PENDING
2026-06-20 16:30:49 INFO Jobs still running: ['4927708', '4927709', '4927710', '4927715', '4927716', '4927719', '4927720', '4927721', '4927722', '4927723', '4927724', '4927726', '4927727']. Waiting...
2026-06-20 16:31:04 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:31:04 INFO None 4927709: status RUNNING/PENDING
2026-06-20 16:31:04 INFO None 4927710: status RUNNING/PENDING
2026-06-20 16:31:04 INFO None 4927712: status FINISHED
2026-06-20 16:31:04 INFO None 4927713: status FINISHED
2026-06-20 16:31:04 INFO None 4927715: status RUNNING/PENDING
2026-06-20 16:31:04 INFO None 4927716: status RUNNING/PENDING
2026-06-20 16:31:04 INFO None 4927719: status RUNNING/PENDING
2026-06-20 16:31:04 INFO None 4927720: status RUNNING/PENDING
2026-06-20 16:31:04 INFO None 4927721: status RUNNING/PENDING
2026-06-20 16:31:04 INFO None 4927722: status RUNNING/PENDING
2026-06-20 16:31:04 INFO None 4927723: status RUNNING/PENDING
2026-06-20 16:31:04 INFO None 4927724: status RUNNING/PENDING
2026-06-20 16:31:04 INFO None 4927726: status FINISHED
2026-06-20 16:31:04 INFO None 4927727: status RUNNING/PENDING
2026-06-20 16:31:04 INFO Jobs still running: ['4927708', '4927709', '4927710', '4927715', '4927716', '4927719', '4927720', '4927721', '4927722', '4927723', '4927724', '4927727']. Waiting...
2026-06-20 16:31:19 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:31:19 INFO None 4927709: status RUNNING/PENDING
2026-06-20 16:31:19 INFO None 4927710: status RUNNING/PENDING
2026-06-20 16:31:19 INFO None 4927712: status FINISHED
2026-06-20 16:31:19 INFO None 4927713: status FINISHED
2026-06-20 16:31:19 INFO None 4927715: status RUNNING/PENDING
2026-06-20 16:31:19 INFO None 4927716: status RUNNING/PENDING
2026-06-20 16:31:19 INFO None 4927719: status RUNNING/PENDING
2026-06-20 16:31:19 INFO None 4927720: status RUNNING/PENDING
2026-06-20 16:31:19 INFO None 4927721: status FINISHED
2026-06-20 16:31:19 INFO None 4927722: status FINISHED
2026-06-20 16:31:20 INFO None 4927723: status FINISHED
2026-06-20 16:31:20 INFO None 4927724: status FINISHED
2026-06-20 16:31:20 INFO None 4927726: status FINISHED
2026-06-20 16:31:20 INFO None 4927727: status RUNNING/PENDING
2026-06-20 16:31:20 INFO Jobs still running: ['4927708', '4927709', '4927710', '4927715', '4927716', '4927719', '4927720', '4927727']. Waiting...
2026-06-20 16:31:35 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:31:35 INFO None 4927709: status RUNNING/PENDING
2026-06-20 16:31:35 INFO None 4927710: status FINISHED
2026-06-20 16:31:35 INFO None 4927712: status FINISHED
2026-06-20 16:31:35 INFO None 4927713: status FINISHED
2026-06-20 16:31:35 INFO None 4927715: status RUNNING/PENDING
2026-06-20 16:31:35 INFO None 4927716: status FINISHED
2026-06-20 16:31:35 INFO None 4927719: status FINISHED
2026-06-20 16:31:35 INFO None 4927720: status FINISHED
2026-06-20 16:31:35 INFO None 4927721: status FINISHED
2026-06-20 16:31:35 INFO None 4927722: status FINISHED
2026-06-20 16:31:35 INFO None 4927723: status FINISHED
2026-06-20 16:31:35 INFO None 4927724: status FINISHED
2026-06-20 16:31:35 INFO None 4927726: status FINISHED
2026-06-20 16:31:35 INFO None 4927727: status RUNNING/PENDING
2026-06-20 16:31:35 INFO Jobs still running: ['4927708', '4927709', '4927715', '4927727']. Waiting...
2026-06-20 16:31:50 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:31:50 INFO None 4927709: status FINISHED
2026-06-20 16:31:50 INFO None 4927710: status FINISHED
2026-06-20 16:31:50 INFO None 4927712: status FINISHED
2026-06-20 16:31:50 INFO None 4927713: status FINISHED
2026-06-20 16:31:50 INFO None 4927715: status FINISHED
2026-06-20 16:31:50 INFO None 4927716: status FINISHED
2026-06-20 16:31:50 INFO None 4927719: status FINISHED
2026-06-20 16:31:50 INFO None 4927720: status FINISHED
2026-06-20 16:31:50 INFO None 4927721: status FINISHED
2026-06-20 16:31:50 INFO None 4927722: status FINISHED
2026-06-20 16:31:50 INFO None 4927723: status FINISHED
2026-06-20 16:31:50 INFO None 4927724: status FINISHED
2026-06-20 16:31:50 INFO None 4927726: status FINISHED
2026-06-20 16:31:50 INFO None 4927727: status RUNNING/PENDING
2026-06-20 16:31:50 INFO Jobs still running: ['4927708', '4927727']. Waiting...
2026-06-20 16:32:05 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:32:05 INFO None 4927709: status FINISHED
2026-06-20 16:32:05 INFO None 4927710: status FINISHED
2026-06-20 16:32:05 INFO None 4927712: status FINISHED
2026-06-20 16:32:05 INFO None 4927713: status FINISHED
2026-06-20 16:32:05 INFO None 4927715: status FINISHED
2026-06-20 16:32:05 INFO None 4927716: status FINISHED
2026-06-20 16:32:05 INFO None 4927719: status FINISHED
2026-06-20 16:32:05 INFO None 4927720: status FINISHED
2026-06-20 16:32:05 INFO None 4927721: status FINISHED
2026-06-20 16:32:05 INFO None 4927722: status FINISHED
2026-06-20 16:32:05 INFO None 4927723: status FINISHED
2026-06-20 16:32:06 INFO None 4927724: status FINISHED
2026-06-20 16:32:06 INFO None 4927726: status FINISHED
2026-06-20 16:32:06 INFO None 4927727: status RUNNING/PENDING
2026-06-20 16:32:06 INFO Jobs still running: ['4927708', '4927727']. Waiting...
2026-06-20 16:32:21 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:32:21 INFO None 4927709: status FINISHED
2026-06-20 16:32:21 INFO None 4927710: status FINISHED
2026-06-20 16:32:21 INFO None 4927712: status FINISHED
2026-06-20 16:32:21 INFO None 4927713: status FINISHED
2026-06-20 16:32:21 INFO None 4927715: status FINISHED
2026-06-20 16:32:21 INFO None 4927716: status FINISHED
2026-06-20 16:32:21 INFO None 4927719: status FINISHED
2026-06-20 16:32:21 INFO None 4927720: status FINISHED
2026-06-20 16:32:21 INFO None 4927721: status FINISHED
2026-06-20 16:32:21 INFO None 4927722: status FINISHED
2026-06-20 16:32:21 INFO None 4927723: status FINISHED
2026-06-20 16:32:21 INFO None 4927724: status FINISHED
2026-06-20 16:32:21 INFO None 4927726: status FINISHED
2026-06-20 16:32:21 INFO None 4927727: status FINISHED
2026-06-20 16:32:21 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:32:36 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:32:36 INFO None 4927709: status FINISHED
2026-06-20 16:32:36 INFO None 4927710: status FINISHED
2026-06-20 16:32:36 INFO None 4927712: status FINISHED
2026-06-20 16:32:36 INFO None 4927713: status FINISHED
2026-06-20 16:32:36 INFO None 4927715: status FINISHED
2026-06-20 16:32:36 INFO None 4927716: status FINISHED
2026-06-20 16:32:36 INFO None 4927719: status FINISHED
2026-06-20 16:32:36 INFO None 4927720: status FINISHED
2026-06-20 16:32:36 INFO None 4927721: status FINISHED
2026-06-20 16:32:36 INFO None 4927722: status FINISHED
2026-06-20 16:32:36 INFO None 4927723: status FINISHED
2026-06-20 16:32:36 INFO None 4927724: status FINISHED
2026-06-20 16:32:36 INFO None 4927726: status FINISHED
2026-06-20 16:32:36 INFO None 4927727: status FINISHED
2026-06-20 16:32:36 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:32:51 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:32:51 INFO None 4927709: status FINISHED
2026-06-20 16:32:51 INFO None 4927710: status FINISHED
2026-06-20 16:32:51 INFO None 4927712: status FINISHED
2026-06-20 16:32:51 INFO None 4927713: status FINISHED
2026-06-20 16:32:51 INFO None 4927715: status FINISHED
2026-06-20 16:32:51 INFO None 4927716: status FINISHED
2026-06-20 16:32:51 INFO None 4927719: status FINISHED
2026-06-20 16:32:51 INFO None 4927720: status FINISHED
2026-06-20 16:32:51 INFO None 4927721: status FINISHED
2026-06-20 16:32:51 INFO None 4927722: status FINISHED
2026-06-20 16:32:51 INFO None 4927723: status FINISHED
2026-06-20 16:32:51 INFO None 4927724: status FINISHED
2026-06-20 16:32:52 INFO None 4927726: status FINISHED
2026-06-20 16:32:52 INFO None 4927727: status FINISHED
2026-06-20 16:32:52 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:33:07 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:33:07 INFO None 4927709: status FINISHED
2026-06-20 16:33:07 INFO None 4927710: status FINISHED
2026-06-20 16:33:07 INFO None 4927712: status FINISHED
2026-06-20 16:33:07 INFO None 4927713: status FINISHED
2026-06-20 16:33:07 INFO None 4927715: status FINISHED
2026-06-20 16:33:07 INFO None 4927716: status FINISHED
2026-06-20 16:33:07 INFO None 4927719: status FINISHED
2026-06-20 16:33:07 INFO None 4927720: status FINISHED
2026-06-20 16:33:07 INFO None 4927721: status FINISHED
2026-06-20 16:33:07 INFO None 4927722: status FINISHED
2026-06-20 16:33:07 INFO None 4927723: status FINISHED
2026-06-20 16:33:07 INFO None 4927724: status FINISHED
2026-06-20 16:33:07 INFO None 4927726: status FINISHED
2026-06-20 16:33:07 INFO None 4927727: status FINISHED
2026-06-20 16:33:07 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:33:22 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:33:22 INFO None 4927709: status FINISHED
2026-06-20 16:33:22 INFO None 4927710: status FINISHED
2026-06-20 16:33:22 INFO None 4927712: status FINISHED
2026-06-20 16:33:22 INFO None 4927713: status FINISHED
2026-06-20 16:33:22 INFO None 4927715: status FINISHED
2026-06-20 16:33:22 INFO None 4927716: status FINISHED
2026-06-20 16:33:22 INFO None 4927719: status FINISHED
2026-06-20 16:33:22 INFO None 4927720: status FINISHED
2026-06-20 16:33:22 INFO None 4927721: status FINISHED
2026-06-20 16:33:22 INFO None 4927722: status FINISHED
2026-06-20 16:33:22 INFO None 4927723: status FINISHED
2026-06-20 16:33:22 INFO None 4927724: status FINISHED
2026-06-20 16:33:22 INFO None 4927726: status FINISHED
2026-06-20 16:33:22 INFO None 4927727: status FINISHED
2026-06-20 16:33:22 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:33:37 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:33:37 INFO None 4927709: status FINISHED
2026-06-20 16:33:37 INFO None 4927710: status FINISHED
2026-06-20 16:33:37 INFO None 4927712: status FINISHED
2026-06-20 16:33:37 INFO None 4927713: status FINISHED
2026-06-20 16:33:37 INFO None 4927715: status FINISHED
2026-06-20 16:33:37 INFO None 4927716: status FINISHED
2026-06-20 16:33:37 INFO None 4927719: status FINISHED
2026-06-20 16:33:37 INFO None 4927720: status FINISHED
2026-06-20 16:33:37 INFO None 4927721: status FINISHED
2026-06-20 16:33:37 INFO None 4927722: status FINISHED
2026-06-20 16:33:37 INFO None 4927723: status FINISHED
2026-06-20 16:33:37 INFO None 4927724: status FINISHED
2026-06-20 16:33:37 INFO None 4927726: status FINISHED
2026-06-20 16:33:37 INFO None 4927727: status FINISHED
2026-06-20 16:33:37 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:33:53 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:33:53 INFO None 4927709: status FINISHED
2026-06-20 16:33:53 INFO None 4927710: status FINISHED
2026-06-20 16:33:53 INFO None 4927712: status FINISHED
2026-06-20 16:33:53 INFO None 4927713: status FINISHED
2026-06-20 16:33:53 INFO None 4927715: status FINISHED
2026-06-20 16:33:53 INFO None 4927716: status FINISHED
2026-06-20 16:33:53 INFO None 4927719: status FINISHED
2026-06-20 16:33:53 INFO None 4927720: status FINISHED
2026-06-20 16:33:53 INFO None 4927721: status FINISHED
2026-06-20 16:33:53 INFO None 4927722: status FINISHED
2026-06-20 16:33:53 INFO None 4927723: status FINISHED
2026-06-20 16:33:53 INFO None 4927724: status FINISHED
2026-06-20 16:33:53 INFO None 4927726: status FINISHED
2026-06-20 16:33:53 INFO None 4927727: status FINISHED
2026-06-20 16:33:53 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:34:08 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:34:08 INFO None 4927709: status FINISHED
2026-06-20 16:34:08 INFO None 4927710: status FINISHED
2026-06-20 16:34:08 INFO None 4927712: status FINISHED
2026-06-20 16:34:08 INFO None 4927713: status FINISHED
2026-06-20 16:34:08 INFO None 4927715: status FINISHED
2026-06-20 16:34:08 INFO None 4927716: status FINISHED
2026-06-20 16:34:08 INFO None 4927719: status FINISHED
2026-06-20 16:34:08 INFO None 4927720: status FINISHED
2026-06-20 16:34:08 INFO None 4927721: status FINISHED
2026-06-20 16:34:08 INFO None 4927722: status FINISHED
2026-06-20 16:34:08 INFO None 4927723: status FINISHED
2026-06-20 16:34:08 INFO None 4927724: status FINISHED
2026-06-20 16:34:08 INFO None 4927726: status FINISHED
2026-06-20 16:34:08 INFO None 4927727: status FINISHED
2026-06-20 16:34:08 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:34:23 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:34:23 INFO None 4927709: status FINISHED
2026-06-20 16:34:23 INFO None 4927710: status FINISHED
2026-06-20 16:34:23 INFO None 4927712: status FINISHED
2026-06-20 16:34:23 INFO None 4927713: status FINISHED
2026-06-20 16:34:23 INFO None 4927715: status FINISHED
2026-06-20 16:34:23 INFO None 4927716: status FINISHED
2026-06-20 16:34:23 INFO None 4927719: status FINISHED
2026-06-20 16:34:23 INFO None 4927720: status FINISHED
2026-06-20 16:34:23 INFO None 4927721: status FINISHED
2026-06-20 16:34:23 INFO None 4927722: status FINISHED
2026-06-20 16:34:23 INFO None 4927723: status FINISHED
2026-06-20 16:34:23 INFO None 4927724: status FINISHED
2026-06-20 16:34:23 INFO None 4927726: status FINISHED
2026-06-20 16:34:23 INFO None 4927727: status FINISHED
2026-06-20 16:34:23 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:34:38 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:34:38 INFO None 4927709: status FINISHED
2026-06-20 16:34:39 INFO None 4927710: status FINISHED
2026-06-20 16:34:39 INFO None 4927712: status FINISHED
2026-06-20 16:34:39 INFO None 4927713: status FINISHED
2026-06-20 16:34:39 INFO None 4927715: status FINISHED
2026-06-20 16:34:39 INFO None 4927716: status FINISHED
2026-06-20 16:34:39 INFO None 4927719: status FINISHED
2026-06-20 16:34:39 INFO None 4927720: status FINISHED
2026-06-20 16:34:39 INFO None 4927721: status FINISHED
2026-06-20 16:34:39 INFO None 4927722: status FINISHED
2026-06-20 16:34:39 INFO None 4927723: status FINISHED
2026-06-20 16:34:39 INFO None 4927724: status FINISHED
2026-06-20 16:34:39 INFO None 4927726: status FINISHED
2026-06-20 16:34:39 INFO None 4927727: status FINISHED
2026-06-20 16:34:39 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:34:54 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:34:54 INFO None 4927709: status FINISHED
2026-06-20 16:34:54 INFO None 4927710: status FINISHED
2026-06-20 16:34:54 INFO None 4927712: status FINISHED
2026-06-20 16:34:54 INFO None 4927713: status FINISHED
2026-06-20 16:34:54 INFO None 4927715: status FINISHED
2026-06-20 16:34:54 INFO None 4927716: status FINISHED
2026-06-20 16:34:54 INFO None 4927719: status FINISHED
2026-06-20 16:34:54 INFO None 4927720: status FINISHED
2026-06-20 16:34:54 INFO None 4927721: status FINISHED
2026-06-20 16:34:54 INFO None 4927722: status FINISHED
2026-06-20 16:34:54 INFO None 4927723: status FINISHED
2026-06-20 16:34:54 INFO None 4927724: status FINISHED
2026-06-20 16:34:54 INFO None 4927726: status FINISHED
2026-06-20 16:34:54 INFO None 4927727: status FINISHED
2026-06-20 16:34:54 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:35:09 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:35:09 INFO None 4927709: status FINISHED
2026-06-20 16:35:09 INFO None 4927710: status FINISHED
2026-06-20 16:35:09 INFO None 4927712: status FINISHED
2026-06-20 16:35:09 INFO None 4927713: status FINISHED
2026-06-20 16:35:09 INFO None 4927715: status FINISHED
2026-06-20 16:35:09 INFO None 4927716: status FINISHED
2026-06-20 16:35:09 INFO None 4927719: status FINISHED
2026-06-20 16:35:09 INFO None 4927720: status FINISHED
2026-06-20 16:35:09 INFO None 4927721: status FINISHED
2026-06-20 16:35:09 INFO None 4927722: status FINISHED
2026-06-20 16:35:09 INFO None 4927723: status FINISHED
2026-06-20 16:35:09 INFO None 4927724: status FINISHED
2026-06-20 16:35:09 INFO None 4927726: status FINISHED
2026-06-20 16:35:09 INFO None 4927727: status FINISHED
2026-06-20 16:35:09 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:35:24 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:35:24 INFO None 4927709: status FINISHED
2026-06-20 16:35:24 INFO None 4927710: status FINISHED
2026-06-20 16:35:24 INFO None 4927712: status FINISHED
2026-06-20 16:35:24 INFO None 4927713: status FINISHED
2026-06-20 16:35:24 INFO None 4927715: status FINISHED
2026-06-20 16:35:25 INFO None 4927716: status FINISHED
2026-06-20 16:35:25 INFO None 4927719: status FINISHED
2026-06-20 16:35:25 INFO None 4927720: status FINISHED
2026-06-20 16:35:25 INFO None 4927721: status FINISHED
2026-06-20 16:35:25 INFO None 4927722: status FINISHED
2026-06-20 16:35:25 INFO None 4927723: status FINISHED
2026-06-20 16:35:25 INFO None 4927724: status FINISHED
2026-06-20 16:35:25 INFO None 4927726: status FINISHED
2026-06-20 16:35:25 INFO None 4927727: status FINISHED
2026-06-20 16:35:25 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:35:40 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:35:40 INFO None 4927709: status FINISHED
2026-06-20 16:35:40 INFO None 4927710: status FINISHED
2026-06-20 16:35:40 INFO None 4927712: status FINISHED
2026-06-20 16:35:40 INFO None 4927713: status FINISHED
2026-06-20 16:35:40 INFO None 4927715: status FINISHED
2026-06-20 16:35:40 INFO None 4927716: status FINISHED
2026-06-20 16:35:40 INFO None 4927719: status FINISHED
2026-06-20 16:35:40 INFO None 4927720: status FINISHED
2026-06-20 16:35:40 INFO None 4927721: status FINISHED
2026-06-20 16:35:40 INFO None 4927722: status FINISHED
2026-06-20 16:35:40 INFO None 4927723: status FINISHED
2026-06-20 16:35:40 INFO None 4927724: status FINISHED
2026-06-20 16:35:40 INFO None 4927726: status FINISHED
2026-06-20 16:35:40 INFO None 4927727: status FINISHED
2026-06-20 16:35:40 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:35:55 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:35:55 INFO None 4927709: status FINISHED
2026-06-20 16:35:55 INFO None 4927710: status FINISHED
2026-06-20 16:35:55 INFO None 4927712: status FINISHED
2026-06-20 16:35:55 INFO None 4927713: status FINISHED
2026-06-20 16:35:55 INFO None 4927715: status FINISHED
2026-06-20 16:35:55 INFO None 4927716: status FINISHED
2026-06-20 16:35:55 INFO None 4927719: status FINISHED
2026-06-20 16:35:55 INFO None 4927720: status FINISHED
2026-06-20 16:35:55 INFO None 4927721: status FINISHED
2026-06-20 16:35:55 INFO None 4927722: status FINISHED
2026-06-20 16:35:55 INFO None 4927723: status FINISHED
2026-06-20 16:35:55 INFO None 4927724: status FINISHED
2026-06-20 16:35:55 INFO None 4927726: status FINISHED
2026-06-20 16:35:55 INFO None 4927727: status FINISHED
2026-06-20 16:35:55 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:36:10 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:36:10 INFO None 4927709: status FINISHED
2026-06-20 16:36:10 INFO None 4927710: status FINISHED
2026-06-20 16:36:10 INFO None 4927712: status FINISHED
2026-06-20 16:36:10 INFO None 4927713: status FINISHED
2026-06-20 16:36:10 INFO None 4927715: status FINISHED
2026-06-20 16:36:10 INFO None 4927716: status FINISHED
2026-06-20 16:36:11 INFO None 4927719: status FINISHED
2026-06-20 16:36:11 INFO None 4927720: status FINISHED
2026-06-20 16:36:11 INFO None 4927721: status FINISHED
2026-06-20 16:36:11 INFO None 4927722: status FINISHED
2026-06-20 16:36:11 INFO None 4927723: status FINISHED
2026-06-20 16:36:11 INFO None 4927724: status FINISHED
2026-06-20 16:36:11 INFO None 4927726: status FINISHED
2026-06-20 16:36:11 INFO None 4927727: status FINISHED
2026-06-20 16:36:11 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:36:26 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:36:26 INFO None 4927709: status FINISHED
2026-06-20 16:36:26 INFO None 4927710: status FINISHED
2026-06-20 16:36:26 INFO None 4927712: status FINISHED
2026-06-20 16:36:26 INFO None 4927713: status FINISHED
2026-06-20 16:36:26 INFO None 4927715: status FINISHED
2026-06-20 16:36:26 INFO None 4927716: status FINISHED
2026-06-20 16:36:26 INFO None 4927719: status FINISHED
2026-06-20 16:36:26 INFO None 4927720: status FINISHED
2026-06-20 16:36:26 INFO None 4927721: status FINISHED
2026-06-20 16:36:26 INFO None 4927722: status FINISHED
2026-06-20 16:36:26 INFO None 4927723: status FINISHED
2026-06-20 16:36:26 INFO None 4927724: status FINISHED
2026-06-20 16:36:26 INFO None 4927726: status FINISHED
2026-06-20 16:36:26 INFO None 4927727: status FINISHED
2026-06-20 16:36:26 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:36:41 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:36:41 INFO None 4927709: status FINISHED
2026-06-20 16:36:41 INFO None 4927710: status FINISHED
2026-06-20 16:36:41 INFO None 4927712: status FINISHED
2026-06-20 16:36:41 INFO None 4927713: status FINISHED
2026-06-20 16:36:41 INFO None 4927715: status FINISHED
2026-06-20 16:36:41 INFO None 4927716: status FINISHED
2026-06-20 16:36:41 INFO None 4927719: status FINISHED
2026-06-20 16:36:41 INFO None 4927720: status FINISHED
2026-06-20 16:36:41 INFO None 4927721: status FINISHED
2026-06-20 16:36:41 INFO None 4927722: status FINISHED
2026-06-20 16:36:41 INFO None 4927723: status FINISHED
2026-06-20 16:36:41 INFO None 4927724: status FINISHED
2026-06-20 16:36:41 INFO None 4927726: status FINISHED
2026-06-20 16:36:41 INFO None 4927727: status FINISHED
2026-06-20 16:36:41 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:36:56 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:36:56 INFO None 4927709: status FINISHED
2026-06-20 16:36:56 INFO None 4927710: status FINISHED
2026-06-20 16:36:56 INFO None 4927712: status FINISHED
2026-06-20 16:36:56 INFO None 4927713: status FINISHED
2026-06-20 16:36:56 INFO None 4927715: status FINISHED
2026-06-20 16:36:57 INFO None 4927716: status FINISHED
2026-06-20 16:36:57 INFO None 4927719: status FINISHED
2026-06-20 16:36:57 INFO None 4927720: status FINISHED
2026-06-20 16:36:57 INFO None 4927721: status FINISHED
2026-06-20 16:36:57 INFO None 4927722: status FINISHED
2026-06-20 16:36:57 INFO None 4927723: status FINISHED
2026-06-20 16:36:57 INFO None 4927724: status FINISHED
2026-06-20 16:36:57 INFO None 4927726: status FINISHED
2026-06-20 16:36:57 INFO None 4927727: status FINISHED
2026-06-20 16:36:57 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:37:12 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:37:12 INFO None 4927709: status FINISHED
2026-06-20 16:37:12 INFO None 4927710: status FINISHED
2026-06-20 16:37:12 INFO None 4927712: status FINISHED
2026-06-20 16:37:12 INFO None 4927713: status FINISHED
2026-06-20 16:37:12 INFO None 4927715: status FINISHED
2026-06-20 16:37:12 INFO None 4927716: status FINISHED
2026-06-20 16:37:12 INFO None 4927719: status FINISHED
2026-06-20 16:37:12 INFO None 4927720: status FINISHED
2026-06-20 16:37:12 INFO None 4927721: status FINISHED
2026-06-20 16:37:12 INFO None 4927722: status FINISHED
2026-06-20 16:37:12 INFO None 4927723: status FINISHED
2026-06-20 16:37:12 INFO None 4927724: status FINISHED
2026-06-20 16:37:12 INFO None 4927726: status FINISHED
2026-06-20 16:37:12 INFO None 4927727: status FINISHED
2026-06-20 16:37:12 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:37:27 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:37:27 INFO None 4927709: status FINISHED
2026-06-20 16:37:27 INFO None 4927710: status FINISHED
2026-06-20 16:37:27 INFO None 4927712: status FINISHED
2026-06-20 16:37:27 INFO None 4927713: status FINISHED
2026-06-20 16:37:27 INFO None 4927715: status FINISHED
2026-06-20 16:37:27 INFO None 4927716: status FINISHED
2026-06-20 16:37:27 INFO None 4927719: status FINISHED
2026-06-20 16:37:27 INFO None 4927720: status FINISHED
2026-06-20 16:37:27 INFO None 4927721: status FINISHED
2026-06-20 16:37:27 INFO None 4927722: status FINISHED
2026-06-20 16:37:27 INFO None 4927723: status FINISHED
2026-06-20 16:37:27 INFO None 4927724: status FINISHED
2026-06-20 16:37:27 INFO None 4927726: status FINISHED
2026-06-20 16:37:27 INFO None 4927727: status FINISHED
2026-06-20 16:37:27 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:37:42 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:37:42 INFO None 4927709: status FINISHED
2026-06-20 16:37:42 INFO None 4927710: status FINISHED
2026-06-20 16:37:42 INFO None 4927712: status FINISHED
2026-06-20 16:37:42 INFO None 4927713: status FINISHED
2026-06-20 16:37:42 INFO None 4927715: status FINISHED
2026-06-20 16:37:42 INFO None 4927716: status FINISHED
2026-06-20 16:37:42 INFO None 4927719: status FINISHED
2026-06-20 16:37:43 INFO None 4927720: status FINISHED
2026-06-20 16:37:43 INFO None 4927721: status FINISHED
2026-06-20 16:37:43 INFO None 4927722: status FINISHED
2026-06-20 16:37:43 INFO None 4927723: status FINISHED
2026-06-20 16:37:43 INFO None 4927724: status FINISHED
2026-06-20 16:37:43 INFO None 4927726: status FINISHED
2026-06-20 16:37:43 INFO None 4927727: status FINISHED
2026-06-20 16:37:43 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:37:58 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:37:58 INFO None 4927709: status FINISHED
2026-06-20 16:37:58 INFO None 4927710: status FINISHED
2026-06-20 16:37:58 INFO None 4927712: status FINISHED
2026-06-20 16:37:58 INFO None 4927713: status FINISHED
2026-06-20 16:37:58 INFO None 4927715: status FINISHED
2026-06-20 16:37:58 INFO None 4927716: status FINISHED
2026-06-20 16:37:58 INFO None 4927719: status FINISHED
2026-06-20 16:37:58 INFO None 4927720: status FINISHED
2026-06-20 16:37:58 INFO None 4927721: status FINISHED
2026-06-20 16:37:58 INFO None 4927722: status FINISHED
2026-06-20 16:37:58 INFO None 4927723: status FINISHED
2026-06-20 16:37:58 INFO None 4927724: status FINISHED
2026-06-20 16:37:58 INFO None 4927726: status FINISHED
2026-06-20 16:37:58 INFO None 4927727: status FINISHED
2026-06-20 16:37:58 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:38:13 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:38:13 INFO None 4927709: status FINISHED
2026-06-20 16:38:13 INFO None 4927710: status FINISHED
2026-06-20 16:38:13 INFO None 4927712: status FINISHED
2026-06-20 16:38:13 INFO None 4927713: status FINISHED
2026-06-20 16:38:13 INFO None 4927715: status FINISHED
2026-06-20 16:38:13 INFO None 4927716: status FINISHED
2026-06-20 16:38:13 INFO None 4927719: status FINISHED
2026-06-20 16:38:13 INFO None 4927720: status FINISHED
2026-06-20 16:38:13 INFO None 4927721: status FINISHED
2026-06-20 16:38:13 INFO None 4927722: status FINISHED
2026-06-20 16:38:13 INFO None 4927723: status FINISHED
2026-06-20 16:38:13 INFO None 4927724: status FINISHED
2026-06-20 16:38:13 INFO None 4927726: status FINISHED
2026-06-20 16:38:13 INFO None 4927727: status FINISHED
2026-06-20 16:38:13 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:38:28 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:38:28 INFO None 4927709: status FINISHED
2026-06-20 16:38:28 INFO None 4927710: status FINISHED
2026-06-20 16:38:28 INFO None 4927712: status FINISHED
2026-06-20 16:38:28 INFO None 4927713: status FINISHED
2026-06-20 16:38:28 INFO None 4927715: status FINISHED
2026-06-20 16:38:28 INFO None 4927716: status FINISHED
2026-06-20 16:38:28 INFO None 4927719: status FINISHED
2026-06-20 16:38:28 INFO None 4927720: status FINISHED
2026-06-20 16:38:28 INFO None 4927721: status FINISHED
2026-06-20 16:38:28 INFO None 4927722: status FINISHED
2026-06-20 16:38:28 INFO None 4927723: status FINISHED
2026-06-20 16:38:28 INFO None 4927724: status FINISHED
2026-06-20 16:38:29 INFO None 4927726: status FINISHED
2026-06-20 16:38:29 INFO None 4927727: status FINISHED
2026-06-20 16:38:29 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:38:44 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:38:44 INFO None 4927709: status FINISHED
2026-06-20 16:38:44 INFO None 4927710: status FINISHED
2026-06-20 16:38:44 INFO None 4927712: status FINISHED
2026-06-20 16:38:44 INFO None 4927713: status FINISHED
2026-06-20 16:38:44 INFO None 4927715: status FINISHED
2026-06-20 16:38:44 INFO None 4927716: status FINISHED
2026-06-20 16:38:44 INFO None 4927719: status FINISHED
2026-06-20 16:38:44 INFO None 4927720: status FINISHED
2026-06-20 16:38:44 INFO None 4927721: status FINISHED
2026-06-20 16:38:44 INFO None 4927722: status FINISHED
2026-06-20 16:38:44 INFO None 4927723: status FINISHED
2026-06-20 16:38:44 INFO None 4927724: status FINISHED
2026-06-20 16:38:44 INFO None 4927726: status FINISHED
2026-06-20 16:38:44 INFO None 4927727: status FINISHED
2026-06-20 16:38:44 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:38:59 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:38:59 INFO None 4927709: status FINISHED
2026-06-20 16:38:59 INFO None 4927710: status FINISHED
2026-06-20 16:38:59 INFO None 4927712: status FINISHED
2026-06-20 16:38:59 INFO None 4927713: status FINISHED
2026-06-20 16:38:59 INFO None 4927715: status FINISHED
2026-06-20 16:38:59 INFO None 4927716: status FINISHED
2026-06-20 16:38:59 INFO None 4927719: status FINISHED
2026-06-20 16:38:59 INFO None 4927720: status FINISHED
2026-06-20 16:38:59 INFO None 4927721: status FINISHED
2026-06-20 16:38:59 INFO None 4927722: status FINISHED
2026-06-20 16:38:59 INFO None 4927723: status FINISHED
2026-06-20 16:38:59 INFO None 4927724: status FINISHED
2026-06-20 16:38:59 INFO None 4927726: status FINISHED
2026-06-20 16:38:59 INFO None 4927727: status FINISHED
2026-06-20 16:38:59 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:39:14 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:39:14 INFO None 4927709: status FINISHED
2026-06-20 16:39:14 INFO None 4927710: status FINISHED
2026-06-20 16:39:14 INFO None 4927712: status FINISHED
2026-06-20 16:39:14 INFO None 4927713: status FINISHED
2026-06-20 16:39:14 INFO None 4927715: status FINISHED
2026-06-20 16:39:14 INFO None 4927716: status FINISHED
2026-06-20 16:39:14 INFO None 4927719: status FINISHED
2026-06-20 16:39:14 INFO None 4927720: status FINISHED
2026-06-20 16:39:14 INFO None 4927721: status FINISHED
2026-06-20 16:39:14 INFO None 4927722: status FINISHED
2026-06-20 16:39:14 INFO None 4927723: status FINISHED
2026-06-20 16:39:14 INFO None 4927724: status FINISHED
2026-06-20 16:39:14 INFO None 4927726: status FINISHED
2026-06-20 16:39:14 INFO None 4927727: status FINISHED
2026-06-20 16:39:14 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:39:29 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:39:30 INFO None 4927709: status FINISHED
2026-06-20 16:39:30 INFO None 4927710: status FINISHED
2026-06-20 16:39:30 INFO None 4927712: status FINISHED
2026-06-20 16:39:30 INFO None 4927713: status FINISHED
2026-06-20 16:39:30 INFO None 4927715: status FINISHED
2026-06-20 16:39:30 INFO None 4927716: status FINISHED
2026-06-20 16:39:30 INFO None 4927719: status FINISHED
2026-06-20 16:39:30 INFO None 4927720: status FINISHED
2026-06-20 16:39:30 INFO None 4927721: status FINISHED
2026-06-20 16:39:30 INFO None 4927722: status FINISHED
2026-06-20 16:39:30 INFO None 4927723: status FINISHED
2026-06-20 16:39:30 INFO None 4927724: status FINISHED
2026-06-20 16:39:30 INFO None 4927726: status FINISHED
2026-06-20 16:39:30 INFO None 4927727: status FINISHED
2026-06-20 16:39:30 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:39:45 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:39:45 INFO None 4927709: status FINISHED
2026-06-20 16:39:45 INFO None 4927710: status FINISHED
2026-06-20 16:39:45 INFO None 4927712: status FINISHED
2026-06-20 16:39:45 INFO None 4927713: status FINISHED
2026-06-20 16:39:45 INFO None 4927715: status FINISHED
2026-06-20 16:39:45 INFO None 4927716: status FINISHED
2026-06-20 16:39:45 INFO None 4927719: status FINISHED
2026-06-20 16:39:45 INFO None 4927720: status FINISHED
2026-06-20 16:39:45 INFO None 4927721: status FINISHED
2026-06-20 16:39:45 INFO None 4927722: status FINISHED
2026-06-20 16:39:45 INFO None 4927723: status FINISHED
2026-06-20 16:39:45 INFO None 4927724: status FINISHED
2026-06-20 16:39:45 INFO None 4927726: status FINISHED
2026-06-20 16:39:45 INFO None 4927727: status FINISHED
2026-06-20 16:39:45 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:40:00 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:40:00 INFO None 4927709: status FINISHED
2026-06-20 16:40:00 INFO None 4927710: status FINISHED
2026-06-20 16:40:00 INFO None 4927712: status FINISHED
2026-06-20 16:40:00 INFO None 4927713: status FINISHED
2026-06-20 16:40:00 INFO None 4927715: status FINISHED
2026-06-20 16:40:00 INFO None 4927716: status FINISHED
2026-06-20 16:40:00 INFO None 4927719: status FINISHED
2026-06-20 16:40:00 INFO None 4927720: status FINISHED
2026-06-20 16:40:00 INFO None 4927721: status FINISHED
2026-06-20 16:40:00 INFO None 4927722: status FINISHED
2026-06-20 16:40:00 INFO None 4927723: status FINISHED
2026-06-20 16:40:00 INFO None 4927724: status FINISHED
2026-06-20 16:40:00 INFO None 4927726: status FINISHED
2026-06-20 16:40:00 INFO None 4927727: status FINISHED
2026-06-20 16:40:00 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:40:15 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:40:15 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:40:15 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:40:15 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:40:15 INFO None 4927713: status FINISHED
2026-06-20 16:40:16 INFO None 4927715: status FINISHED
2026-06-20 16:40:16 INFO None 4927716: status FINISHED
2026-06-20 16:40:16 INFO None 4927719: status FINISHED
2026-06-20 16:40:16 INFO None 4927720: status FINISHED
2026-06-20 16:40:16 INFO None 4927721: status FINISHED
2026-06-20 16:40:16 INFO None 4927722: status FINISHED
2026-06-20 16:40:16 INFO None 4927723: status FINISHED
2026-06-20 16:40:16 INFO None 4927724: status FINISHED
2026-06-20 16:40:16 INFO None 4927726: status FINISHED
2026-06-20 16:40:16 INFO None 4927727: status FINISHED
2026-06-20 16:40:16 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:40:31 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:40:31 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:40:31 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:40:31 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:40:31 INFO None 4927713: status FINISHED
2026-06-20 16:40:31 INFO None 4927715: status FINISHED
2026-06-20 16:40:31 INFO None 4927716: status FINISHED
2026-06-20 16:40:31 INFO None 4927719: status FINISHED
2026-06-20 16:40:31 INFO None 4927720: status FINISHED
2026-06-20 16:40:31 INFO None 4927721: status FINISHED
2026-06-20 16:40:31 INFO None 4927722: status FINISHED
2026-06-20 16:40:31 INFO None 4927723: status FINISHED
2026-06-20 16:40:31 INFO None 4927724: status FINISHED
2026-06-20 16:40:31 INFO None 4927726: status FINISHED
2026-06-20 16:40:31 INFO None 4927727: status FINISHED
2026-06-20 16:40:31 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:40:46 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:40:46 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:40:46 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:40:46 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:40:46 INFO None 4927713: status FINISHED
2026-06-20 16:40:46 INFO None 4927715: status FINISHED
2026-06-20 16:40:46 INFO None 4927716: status FINISHED
2026-06-20 16:40:46 INFO None 4927719: status FINISHED
2026-06-20 16:40:46 INFO None 4927720: status FINISHED
2026-06-20 16:40:46 INFO None 4927721: status FINISHED
2026-06-20 16:40:46 INFO None 4927722: status FINISHED
2026-06-20 16:40:46 INFO None 4927723: status FINISHED
2026-06-20 16:40:46 INFO None 4927724: status FINISHED
2026-06-20 16:40:47 INFO None 4927726: status FINISHED
2026-06-20 16:40:47 INFO None 4927727: status FINISHED
2026-06-20 16:40:47 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:41:02 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:41:02 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:41:02 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:41:02 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:41:02 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:41:02 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:41:02 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:41:02 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:41:02 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:41:02 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:41:02 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:41:02 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:41:02 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:41:02 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:41:02 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:41:02 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:41:17 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:41:17 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:41:17 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:41:17 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:41:17 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:41:17 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:41:17 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:41:17 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:41:17 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:41:17 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:41:17 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:41:17 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:41:17 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:41:17 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:41:17 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:41:17 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:41:32 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:41:32 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:41:32 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:41:32 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:41:32 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:41:32 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:41:32 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:41:32 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:41:32 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:41:32 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:41:32 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:41:32 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:41:32 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:41:32 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:41:32 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:41:32 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:41:47 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:41:47 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:41:47 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:41:47 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:41:48 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:41:48 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:41:48 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:41:48 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:41:48 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:41:48 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:41:48 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:41:48 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:41:48 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:41:48 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:41:48 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:41:48 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:42:03 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:42:03 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:42:03 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:42:03 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:42:03 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:42:03 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:42:03 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:42:03 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:42:03 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:42:03 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:42:03 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:42:03 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:42:03 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:42:03 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:42:03 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:42:03 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:42:18 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:42:18 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:42:18 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:42:18 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:42:18 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:42:18 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:42:18 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:42:18 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:42:18 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:42:18 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:42:18 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:42:18 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:42:18 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:42:18 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:42:18 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:42:18 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:42:33 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:42:33 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:42:33 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:42:33 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:42:33 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:42:33 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:42:33 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:42:33 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:42:34 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:42:34 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:42:34 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:42:34 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:42:34 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:42:34 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:42:34 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:42:34 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:42:49 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:42:49 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:42:49 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:42:49 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:42:49 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:42:49 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:42:49 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:42:49 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:42:49 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:42:49 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:42:49 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:42:49 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:42:49 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:42:49 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:42:49 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:42:49 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:43:04 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:43:04 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:43:04 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:43:04 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:43:04 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:43:04 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:43:04 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:43:04 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:43:04 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:43:04 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:43:04 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:43:04 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:43:04 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:43:04 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:43:04 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:43:04 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:43:19 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:43:19 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:43:19 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:43:19 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:43:19 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:43:19 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:43:19 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:43:19 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:43:19 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:43:19 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:43:19 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:43:19 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:43:19 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:43:20 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:43:20 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:43:20 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:43:35 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:43:35 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:43:35 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:43:35 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:43:35 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:43:35 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:43:35 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:43:35 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:43:35 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:43:35 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:43:35 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:43:35 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:43:35 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:43:35 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:43:35 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:43:35 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:43:50 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:43:50 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:43:50 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:43:50 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:43:50 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:43:50 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:43:50 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:43:50 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:43:50 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:43:50 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:43:50 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:43:50 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:43:50 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:43:50 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:43:50 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:43:50 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:44:05 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:44:05 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:44:05 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:44:05 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:44:05 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:44:05 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:44:05 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:44:05 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:44:05 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:44:05 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:44:05 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:44:05 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:44:05 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:44:05 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:44:05 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:44:05 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:44:20 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:44:20 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:44:20 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:44:21 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:44:21 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:44:21 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:44:21 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:44:21 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:44:21 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:44:21 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:44:21 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:44:21 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:44:21 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:44:21 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:44:21 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:44:21 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:44:36 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:44:36 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:44:36 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:44:36 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:44:36 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:44:36 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:44:36 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:44:36 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:44:36 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:44:36 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:44:36 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:44:36 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:44:36 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:44:36 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:44:36 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:44:36 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:44:51 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:44:51 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:44:51 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:44:51 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:44:51 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:44:51 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:44:51 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:44:51 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:44:51 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:44:51 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:44:51 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:44:51 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:44:51 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:44:51 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:44:51 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:44:51 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:45:06 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:45:06 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:45:06 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:45:06 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:45:06 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:45:06 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:45:06 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:45:06 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:45:06 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:45:07 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:45:07 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:45:07 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:45:07 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:45:07 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:45:07 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:45:07 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:45:22 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:45:22 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:45:22 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:45:22 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:45:22 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:45:22 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:45:22 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:45:22 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:45:22 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:45:22 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:45:22 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:45:22 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:45:22 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:45:22 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:45:22 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:45:22 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:45:37 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:45:37 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:45:37 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:45:37 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:45:37 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:45:37 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:45:37 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:45:37 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:45:37 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:45:37 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:45:37 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:45:37 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:45:37 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:45:37 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:45:37 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:45:37 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:45:52 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:45:52 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:45:52 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:45:52 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:45:52 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:45:52 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:45:52 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:45:52 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:45:52 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:45:52 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:45:52 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:45:52 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:45:52 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:45:52 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:45:53 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:45:53 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:46:08 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:46:08 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:46:08 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:46:08 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:46:08 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:46:08 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:46:08 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:46:08 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:46:08 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:46:08 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:46:08 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:46:08 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:46:08 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:46:08 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:46:08 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:46:08 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:46:23 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:46:23 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:46:23 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:46:23 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:46:23 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:46:23 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:46:23 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:46:23 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:46:23 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:46:23 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:46:23 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:46:23 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:46:23 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:46:23 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:46:23 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:46:23 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:46:38 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:46:38 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:46:38 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:46:38 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:46:38 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:46:38 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:46:38 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:46:38 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:46:38 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:46:38 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:46:38 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:46:38 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:46:38 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:46:38 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:46:38 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:46:38 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:46:53 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:46:53 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:46:54 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:46:54 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:46:54 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:46:54 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:46:54 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:46:54 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:46:54 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:46:54 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:46:54 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:46:54 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:46:54 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:46:54 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:46:54 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:46:54 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:47:09 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:47:09 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:47:09 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:47:09 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:47:09 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:47:09 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:47:09 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:47:09 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:47:09 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:47:09 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:47:09 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:47:09 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:47:09 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:47:09 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:47:09 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:47:09 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:47:24 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:47:24 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:47:24 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:47:24 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:47:24 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:47:24 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:47:24 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:47:24 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:47:24 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:47:24 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:47:24 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:47:24 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:47:24 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:47:24 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:47:24 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:47:24 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:47:39 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:47:39 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:47:39 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:47:39 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:47:39 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:47:39 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:47:39 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:47:40 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:47:40 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:47:40 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:47:40 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:47:40 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:47:40 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:47:40 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:47:40 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:47:40 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:47:55 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:47:55 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:47:55 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:47:55 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:47:55 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:47:55 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:47:55 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:47:55 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:47:55 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:47:55 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:47:55 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:47:55 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:47:55 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:47:55 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:47:55 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:47:55 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:48:10 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:48:10 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:48:10 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:48:10 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:48:10 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:48:10 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:48:10 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:48:10 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:48:10 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:48:10 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:48:10 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:48:10 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:48:10 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:48:10 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:48:10 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:48:10 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:48:25 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:48:25 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:48:25 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:48:25 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:48:25 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:48:25 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:48:25 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:48:25 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:48:25 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:48:25 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:48:25 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:48:25 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:48:25 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:48:25 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:48:26 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:48:26 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:48:41 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:48:41 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:48:41 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:48:41 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:48:41 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:48:41 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:48:41 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:48:41 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:48:41 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:48:41 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:48:41 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:48:41 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:48:41 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:48:41 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:48:41 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:48:41 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:48:56 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:48:56 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:48:56 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:48:56 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:48:56 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:48:56 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:48:56 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:48:56 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:48:56 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:48:56 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:48:56 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:48:56 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:48:56 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:48:56 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:48:56 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:48:56 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:49:11 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:49:11 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:49:11 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:49:11 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:49:11 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:49:11 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:49:11 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:49:11 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:49:11 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:49:11 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:49:11 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:49:11 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:49:11 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:49:11 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:49:11 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:49:11 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:49:26 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:49:26 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:49:26 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:49:26 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:49:26 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:49:26 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:49:27 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:49:27 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:49:27 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:49:27 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:49:27 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:49:27 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:49:27 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:49:27 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:49:27 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:49:27 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:49:42 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:49:42 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:49:42 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:49:42 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:49:42 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:49:42 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:49:42 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:49:42 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:49:42 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:49:42 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:49:42 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:49:42 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:49:42 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:49:42 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:49:42 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:49:42 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:49:57 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:49:57 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:49:57 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:49:57 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:49:57 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:49:57 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:49:57 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:49:57 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:49:57 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:49:57 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:49:57 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:49:57 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:49:57 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:49:57 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:49:57 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:49:57 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:50:12 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:50:12 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:50:12 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:50:12 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:50:12 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:50:12 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:50:12 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:50:12 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:50:12 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:50:12 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:50:12 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:50:12 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:50:12 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:50:13 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:50:13 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:50:13 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:50:28 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:50:28 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:50:28 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:50:28 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:50:28 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:50:28 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:50:28 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:50:28 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:50:28 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:50:28 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:50:28 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:50:28 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:50:28 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:50:28 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:50:28 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:50:28 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:50:43 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:50:43 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:50:43 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:50:43 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:50:43 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:50:43 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:50:43 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:50:43 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:50:43 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:50:43 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:50:43 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:50:43 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:50:43 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:50:43 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:50:43 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:50:43 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:50:58 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:50:58 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:50:58 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:50:58 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:50:58 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:50:58 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:50:58 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:50:58 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:50:58 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:50:58 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:50:58 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:50:58 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:50:58 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:50:58 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:50:58 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:50:58 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:51:13 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:51:13 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:51:13 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:51:13 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:51:13 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:51:14 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:51:14 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:51:14 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:51:14 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:51:14 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:51:14 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:51:14 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:51:14 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:51:14 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:51:14 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:51:14 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:51:29 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:51:29 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:51:29 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:51:29 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:51:29 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:51:29 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:51:29 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:51:29 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:51:29 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:51:29 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:51:29 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:51:29 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:51:29 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:51:29 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:51:29 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:51:29 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:51:44 INFO None 4927708: status RUNNING/PENDING
2026-06-20 16:51:44 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:51:44 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:51:44 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:51:44 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:51:44 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:51:44 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:51:44 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:51:44 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:51:44 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:51:44 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:51:44 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:51:44 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:51:44 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:51:44 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:51:44 INFO Jobs still running: ['4927708']. Waiting...
2026-06-20 16:51:59 INFO None 4927708: status FINISHED
2026-06-20 16:51:59 INFO None 4927709: status FINISHED (not in squeue)
2026-06-20 16:51:59 INFO None 4927710: status FINISHED (not in squeue)
2026-06-20 16:51:59 INFO None 4927712: status FINISHED (not in squeue)
2026-06-20 16:51:59 INFO None 4927713: status FINISHED (not in squeue)
2026-06-20 16:51:59 INFO None 4927715: status FINISHED (not in squeue)
2026-06-20 16:51:59 INFO None 4927716: status FINISHED (not in squeue)
2026-06-20 16:51:59 INFO None 4927719: status FINISHED (not in squeue)
2026-06-20 16:51:59 INFO None 4927720: status FINISHED (not in squeue)
2026-06-20 16:51:59 INFO None 4927721: status FINISHED (not in squeue)
2026-06-20 16:51:59 INFO None 4927722: status FINISHED (not in squeue)
2026-06-20 16:51:59 INFO None 4927723: status FINISHED (not in squeue)
2026-06-20 16:52:00 INFO None 4927724: status FINISHED (not in squeue)
2026-06-20 16:52:00 INFO None 4927726: status FINISHED (not in squeue)
2026-06-20 16:52:00 INFO None 4927727: status FINISHED (not in squeue)
2026-06-20 16:52:00 INFO Jobs ['4927708', '4927709', '4927710', '4927712', '4927713', '4927715', '4927716', '4927719', '4927720', '4927721', '4927722', '4927723', '4927724', '4927726', '4927727'] have finished
2026-06-20 16:52:00 INFO Checking restart files were created ...
2026-06-20 16:52:00 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS1/end.2020020600_1_ENS1.nc(334978955 bytes)
2026-06-20 16:52:00 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS2/end.2020020600_1_ENS2.nc(668832435 bytes)
2026-06-20 16:52:00 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS3/end.2020020600_1_ENS3.nc(668832435 bytes)
2026-06-20 16:52:00 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS4/end.2020020600_1_ENS4.nc(668832435 bytes)
2026-06-20 16:52:00 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS5/end.2020020600_1_ENS5.nc(668832435 bytes)
2026-06-20 16:52:00 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS6/end.2020020600_1_ENS6.nc(668832435 bytes)
2026-06-20 16:52:00 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS7/end.2020020600_1_ENS7.nc(668832435 bytes)
2026-06-20 16:52:00 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS8/end.2020020600_1_ENS8.nc(668832435 bytes)
2026-06-20 16:52:00 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS9/end.2020020600_1_ENS9.nc(668832435 bytes)
2026-06-20 16:52:00 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS10/end.2020020600_1_ENS10.nc(668832435 bytes)
2026-06-20 16:52:00 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS11/end.2020020600_1_ENS11.nc(668832435 bytes)
2026-06-20 16:52:00 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS12/end.2020020600_1_ENS12.nc(668832435 bytes)
2026-06-20 16:52:00 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS13/end.2020020600_1_ENS13.nc(668832435 bytes)
2026-06-20 16:52:00 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS14/end.2020020600_1_ENS14.nc(668832435 bytes)
2026-06-20 16:52:00 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS15/end.2020020600_1_ENS15.nc(668832435 bytes)
2026-06-20 16:52:00 INFO  Run_model() completed successfully.
2026-06-20 16:52:00 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 00:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-20 16:52:00 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 01:00:00 days=153072 seconds=3600
2026-06-20 16:52:00 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-06-20 16:52:00 INFO [TIME] increment current_time 2020-02-06 00:00:00 -> 2020-02-06 01:00:00
2026-06-20 16:52:00 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-20 16:52:00 INFO ---------->>> Running process_satellite_data()
2026-06-20 16:52:00 INFO [DART] No satellite data found, skipping assimilation
2026-06-20 16:52:00 INFO after_assimilation() skipped
2026-06-20 16:52:00 INFO Next run starts from 2020-02-06 01:00:00
2026-06-20 16:52:00 INFO Cycle is DONE; starting a new loop!
2026-06-20 16:52:00 INFO [TIME] step_end current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-20 16:52:00 INFO [TIME] step_start current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-20 16:52:00 INFO [TIME] window start=2020-02-06 01:00:00 end=2020-02-06 09:00:00 run_hours=8 has_assimilation=True
2026-06-20 16:52:00 INFO Creating directories and links for ENS1 to run chimere's parallel part
2026-06-20 16:52:00 INFO Linking EMIS ...
2026-06-20 16:52:00 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens01.nc
2026-06-20 16:52:00 INFO >> Checking links...
2026-06-20 16:52:00 INFO >> All links are good for ENS1  ...
2026-06-20 16:52:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-20 16:52:02 INFO Hourly dataset computed and listing created
2026-06-20 16:52:14 INFO Hourly dataset computed
2026-06-20 16:52:14 INFO Creating directories and links for ENS2 to run chimere's parallel part
2026-06-20 16:52:14 INFO Linking EMIS ...
2026-06-20 16:52:15 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens02.nc
2026-06-20 16:52:15 INFO >> Checking links...
2026-06-20 16:52:15 INFO >> All links are good for ENS2  ...
2026-06-20 16:52:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-20 16:52:16 INFO Hourly dataset computed and listing created
2026-06-20 16:52:18 INFO Hourly dataset computed
2026-06-20 16:52:18 INFO Creating directories and links for ENS3 to run chimere's parallel part
2026-06-20 16:52:18 INFO Linking EMIS ...
2026-06-20 16:52:19 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens03.nc
2026-06-20 16:52:19 INFO >> Checking links...
2026-06-20 16:52:19 INFO >> All links are good for ENS3  ...
2026-06-20 16:52:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-20 16:52:20 INFO Hourly dataset computed and listing created
2026-06-20 16:52:22 INFO Hourly dataset computed
2026-06-20 16:52:23 INFO Creating directories and links for ENS4 to run chimere's parallel part
2026-06-20 16:52:23 INFO Linking EMIS ...
2026-06-20 16:52:23 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens04.nc
2026-06-20 16:52:23 INFO >> Checking links...
2026-06-20 16:52:23 INFO >> All links are good for ENS4  ...
2026-06-20 16:52:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-20 16:52:24 INFO Hourly dataset computed and listing created
2026-06-20 16:52:26 INFO Hourly dataset computed
2026-06-20 16:52:27 INFO Creating directories and links for ENS5 to run chimere's parallel part
2026-06-20 16:52:27 INFO Linking EMIS ...
2026-06-20 16:52:27 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens05.nc
2026-06-20 16:52:27 INFO >> Checking links...
2026-06-20 16:52:27 INFO >> All links are good for ENS5  ...
2026-06-20 16:52:27 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-20 16:52:28 INFO Hourly dataset computed and listing created
2026-06-20 16:52:31 INFO Hourly dataset computed
2026-06-20 16:52:31 INFO Creating directories and links for ENS6 to run chimere's parallel part
2026-06-20 16:52:31 INFO Linking EMIS ...
2026-06-20 16:52:31 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens06.nc
2026-06-20 16:52:31 INFO >> Checking links...
2026-06-20 16:52:31 INFO >> All links are good for ENS6  ...
2026-06-20 16:52:31 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-20 16:52:33 INFO Hourly dataset computed and listing created
2026-06-20 16:52:35 INFO Hourly dataset computed
2026-06-20 16:52:35 INFO Creating directories and links for ENS7 to run chimere's parallel part
2026-06-20 16:52:35 INFO Linking EMIS ...
2026-06-20 16:52:35 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens07.nc
2026-06-20 16:52:36 INFO >> Checking links...
2026-06-20 16:52:36 INFO >> All links are good for ENS7  ...
2026-06-20 16:52:36 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-20 16:52:37 INFO Hourly dataset computed and listing created
2026-06-20 16:52:39 INFO Hourly dataset computed
2026-06-20 16:52:39 INFO Creating directories and links for ENS8 to run chimere's parallel part
2026-06-20 16:52:39 INFO Linking EMIS ...
2026-06-20 16:52:40 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdamp_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens08.nc
2026-06-20 16:52:40 INFO >> Checking links...
2026-06-20 16:52:40 INFO >> All links are good for ENS8  ...
2026-06-20 16:52:40 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-20 16:52:42 INFO Hourly dataset computed and listing created
[2026-06-20T16:52:48.921] error: *** JOB 4927704 ON irene4580 CANCELLED AT 2026-06-20T16:52:48 DUE to SIGNAL Terminated ***
