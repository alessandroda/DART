+ /bin/bash -x /tmp/tmp.1In6lNA0O1
+ SCRIPT_PID=2958223
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
2026-06-22 17:24:39 INFO 
███    ███ ██ ███    ███ ███████ ███████ ██
████  ████ ██ ████  ████ ██      ██      ██
██ ████ ██ ██ ██ ████ ██ █████   ███████ ██
██  ██  ██ ██ ██  ██  ██ ██           ██ ██
██      ██ ██ ██      ██ ███████ ███████ ██



2026-06-22 17:24:39 INFO [PIPELINE] =======================================
2026-06-22 17:24:39 INFO [PIPELINE] Starting chimere–DART orchestrator
2026-06-22 17:24:39 INFO [PIPELINE] Config file: config/config_irene_IM.yaml
2026-06-22 17:24:39 INFO [PIPELINE] Run dir: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart
2026-06-22 17:24:39 INFO [PIPELINE] Log file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/mimesi_orchestrator_logs/chimere_DART_20260622_172439.log
2026-06-22 17:24:39 INFO [PIPELINE] =======================================
2026-06-22 17:24:39 INFO Running assimilation with model_type=ModelType.CHIMERE
2026-06-22 17:24:39 INFO Using scheduler=Scheduler.SLURM, queue=rome
2026-06-22 17:24:39 INFO [STEP] ---- TIME LOOP START ----
2026-06-22 17:24:39 INFO [TIME] step_start current_time=2020-02-06 00:00:00 simulated_time=None dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-22 17:24:39 INFO [TIME] window start=2020-02-06 00:00:00 end=2020-02-06 01:00:00 run_hours=1 has_assimilation=False
2026-06-22 17:24:39 INFO Creating directories and links for ENS1 to run chimere's parallel part
2026-06-22 17:24:39 INFO Copying EMIS ...
2026-06-22 17:24:40 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens01.nc
2026-06-22 17:24:40 INFO Linking END ...
2026-06-22 17:24:40 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020500_24_ENS1.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-22 17:24:40 INFO >> Checking links...
2026-06-22 17:24:40 INFO >> All links are good for ENS1  ...
2026-06-22 17:24:40 INFO Creating directories and links for ENS2 to run chimere's parallel part
2026-06-22 17:24:40 INFO Copying EMIS ...
2026-06-22 17:24:40 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens02.nc
2026-06-22 17:24:40 INFO Linking END ...
2026-06-22 17:24:40 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020500_24_ENS2.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-22 17:24:40 INFO >> Checking links...
2026-06-22 17:24:40 INFO >> All links are good for ENS2  ...
2026-06-22 17:24:40 INFO Creating directories and links for ENS3 to run chimere's parallel part
2026-06-22 17:24:40 INFO Copying EMIS ...
2026-06-22 17:24:40 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens03.nc
2026-06-22 17:24:40 INFO Linking END ...
2026-06-22 17:24:40 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020500_24_ENS3.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-22 17:24:40 INFO >> Checking links...
2026-06-22 17:24:40 INFO >> All links are good for ENS3  ...
2026-06-22 17:24:40 INFO Creating directories and links for ENS4 to run chimere's parallel part
2026-06-22 17:24:40 INFO Copying EMIS ...
2026-06-22 17:24:41 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens04.nc
2026-06-22 17:24:41 INFO Linking END ...
2026-06-22 17:24:41 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020500_24_ENS4.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-22 17:24:41 INFO >> Checking links...
2026-06-22 17:24:41 INFO >> All links are good for ENS4  ...
2026-06-22 17:24:41 INFO Creating directories and links for ENS5 to run chimere's parallel part
2026-06-22 17:24:41 INFO Copying EMIS ...
2026-06-22 17:24:41 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens05.nc
2026-06-22 17:24:41 INFO Linking END ...
2026-06-22 17:24:41 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020500_24_ENS5.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-22 17:24:41 INFO >> Checking links...
2026-06-22 17:24:41 INFO >> All links are good for ENS5  ...
2026-06-22 17:24:41 INFO Creating directories and links for ENS6 to run chimere's parallel part
2026-06-22 17:24:41 INFO Copying EMIS ...
2026-06-22 17:24:42 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens06.nc
2026-06-22 17:24:42 INFO Linking END ...
2026-06-22 17:24:42 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020500_24_ENS6.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-22 17:24:42 INFO >> Checking links...
2026-06-22 17:24:42 INFO >> All links are good for ENS6  ...
2026-06-22 17:24:42 INFO Creating directories and links for ENS7 to run chimere's parallel part
2026-06-22 17:24:42 INFO Copying EMIS ...
2026-06-22 17:24:42 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens07.nc
2026-06-22 17:24:42 INFO Linking END ...
2026-06-22 17:24:42 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020500_24_ENS7.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-22 17:24:42 INFO >> Checking links...
2026-06-22 17:24:42 INFO >> All links are good for ENS7  ...
2026-06-22 17:24:42 INFO Creating directories and links for ENS8 to run chimere's parallel part
2026-06-22 17:24:42 INFO Copying EMIS ...
2026-06-22 17:24:42 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens08.nc
2026-06-22 17:24:42 INFO Linking END ...
2026-06-22 17:24:42 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020500_24_ENS8.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-22 17:24:42 INFO >> Checking links...
2026-06-22 17:24:42 INFO >> All links are good for ENS8  ...
2026-06-22 17:24:42 INFO Creating directories and links for ENS9 to run chimere's parallel part
2026-06-22 17:24:42 INFO Copying EMIS ...
2026-06-22 17:24:43 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens09.nc
2026-06-22 17:24:43 INFO Linking END ...
2026-06-22 17:24:43 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020500_24_ENS9.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-22 17:24:43 INFO >> Checking links...
2026-06-22 17:24:43 INFO >> All links are good for ENS9  ...
2026-06-22 17:24:43 INFO Creating directories and links for ENS10 to run chimere's parallel part
2026-06-22 17:24:43 INFO Copying EMIS ...
2026-06-22 17:24:43 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens10.nc
2026-06-22 17:24:43 INFO Linking END ...
2026-06-22 17:24:43 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020500_24_ENS10.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-22 17:24:43 INFO >> Checking links...
2026-06-22 17:24:43 INFO >> All links are good for ENS10  ...
2026-06-22 17:24:43 INFO Creating directories and links for ENS11 to run chimere's parallel part
2026-06-22 17:24:43 INFO Copying EMIS ...
2026-06-22 17:24:44 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens11.nc
2026-06-22 17:24:44 INFO Linking END ...
2026-06-22 17:24:44 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020500_24_ENS11.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-22 17:24:44 INFO >> Checking links...
2026-06-22 17:24:44 INFO >> All links are good for ENS11  ...
2026-06-22 17:24:44 INFO Creating directories and links for ENS12 to run chimere's parallel part
2026-06-22 17:24:44 INFO Copying EMIS ...
2026-06-22 17:24:44 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens12.nc
2026-06-22 17:24:44 INFO Linking END ...
2026-06-22 17:24:44 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020500_24_ENS12.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-22 17:24:44 INFO >> Checking links...
2026-06-22 17:24:44 INFO >> All links are good for ENS12  ...
2026-06-22 17:24:44 INFO Creating directories and links for ENS13 to run chimere's parallel part
2026-06-22 17:24:44 INFO Copying EMIS ...
2026-06-22 17:24:45 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens13.nc
2026-06-22 17:24:45 INFO Linking END ...
2026-06-22 17:24:45 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020500_24_ENS13.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-22 17:24:45 INFO >> Checking links...
2026-06-22 17:24:45 INFO >> All links are good for ENS13  ...
2026-06-22 17:24:45 INFO Creating directories and links for ENS14 to run chimere's parallel part
2026-06-22 17:24:45 INFO Copying EMIS ...
2026-06-22 17:24:45 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens14.nc
2026-06-22 17:24:45 INFO Linking END ...
2026-06-22 17:24:45 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020500_24_ENS14.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-22 17:24:45 INFO >> Checking links...
2026-06-22 17:24:45 INFO >> All links are good for ENS14  ...
2026-06-22 17:24:45 INFO Creating directories and links for ENS15 to run chimere's parallel part
2026-06-22 17:24:45 INFO Copying EMIS ...
2026-06-22 17:24:45 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens15.nc
2026-06-22 17:24:45 INFO Linking END ...
2026-06-22 17:24:45 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020500_24_ENS15.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-22 17:24:45 INFO >> Checking links...
2026-06-22 17:24:45 INFO >> All links are good for ENS15  ...
2026-06-22 17:24:45 INFO ---------->>> Running CHIMERE model from 2020-02-06 00:00:00 to 2020-02-06 01:00:00
2026-06-22 17:24:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:24:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1
2026-06-22 17:24:45 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS1/end.2020020500_24_ENS1.nc
2026-06-22 17:24:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-06-22 17:24:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:24:45 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-06-22 17:24:45 INFO Queuing job for member 1...
2026-06-22 17:24:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:24:45 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-06-22 17:24:46 INFO Found: ['4950624']
2026-06-22 17:24:51 INFO [TGCC-IRENE] Submitted job with ID:['4950624']
2026-06-22 17:24:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:24:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2
2026-06-22 17:24:51 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS2/end.2020020500_24_ENS2.nc
2026-06-22 17:24:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-06-22 17:24:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:24:51 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-06-22 17:24:51 INFO Queuing job for member 2...
2026-06-22 17:24:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:24:51 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-06-22 17:24:52 INFO Found: ['4950627']
2026-06-22 17:24:57 INFO [TGCC-IRENE] Submitted job with ID:['4950627']
2026-06-22 17:24:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:24:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3
2026-06-22 17:24:57 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS3/end.2020020500_24_ENS3.nc
2026-06-22 17:24:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-06-22 17:24:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:24:57 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-06-22 17:24:57 INFO Queuing job for member 3...
2026-06-22 17:24:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:24:57 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-06-22 17:24:58 INFO Found: ['4950630']
2026-06-22 17:25:03 INFO [TGCC-IRENE] Submitted job with ID:['4950630']
2026-06-22 17:25:03 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:25:03 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4
2026-06-22 17:25:03 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS4/end.2020020500_24_ENS4.nc
2026-06-22 17:25:03 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-06-22 17:25:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:25:03 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-06-22 17:25:03 INFO Queuing job for member 4...
2026-06-22 17:25:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:25:03 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-06-22 17:25:03 INFO Found: ['4950634']
2026-06-22 17:25:08 INFO [TGCC-IRENE] Submitted job with ID:['4950634']
2026-06-22 17:25:08 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:25:08 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5
2026-06-22 17:25:08 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS5/end.2020020500_24_ENS5.nc
2026-06-22 17:25:08 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-06-22 17:25:08 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:25:08 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-06-22 17:25:08 INFO Queuing job for member 5...
2026-06-22 17:25:08 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:25:08 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-06-22 17:25:11 INFO Found: ['4950637']
2026-06-22 17:25:16 INFO [TGCC-IRENE] Submitted job with ID:['4950637']
2026-06-22 17:25:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:25:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6
2026-06-22 17:25:16 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS6/end.2020020500_24_ENS6.nc
2026-06-22 17:25:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-06-22 17:25:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:25:16 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-06-22 17:25:16 INFO Queuing job for member 6...
2026-06-22 17:25:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:25:16 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-06-22 17:25:17 INFO Found: ['4950641']
2026-06-22 17:25:22 INFO [TGCC-IRENE] Submitted job with ID:['4950641']
2026-06-22 17:25:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:25:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7
2026-06-22 17:25:22 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS7/end.2020020500_24_ENS7.nc
2026-06-22 17:25:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-06-22 17:25:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:25:22 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-06-22 17:25:22 INFO Queuing job for member 7...
2026-06-22 17:25:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:25:22 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-06-22 17:25:22 INFO Found: ['4950644']
2026-06-22 17:25:27 INFO [TGCC-IRENE] Submitted job with ID:['4950644']
2026-06-22 17:25:27 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:25:27 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8
2026-06-22 17:25:27 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS8/end.2020020500_24_ENS8.nc
2026-06-22 17:25:27 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-06-22 17:25:27 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:25:27 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-06-22 17:25:27 INFO Queuing job for member 8...
2026-06-22 17:25:27 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:25:27 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-06-22 17:25:28 INFO Found: ['4950649']
2026-06-22 17:25:33 INFO [TGCC-IRENE] Submitted job with ID:['4950649']
2026-06-22 17:25:33 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:25:33 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9
2026-06-22 17:25:33 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS9/end.2020020500_24_ENS9.nc
2026-06-22 17:25:33 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-06-22 17:25:33 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:25:33 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-06-22 17:25:33 INFO Queuing job for member 9...
2026-06-22 17:25:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:25:33 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-06-22 17:25:34 INFO Found: ['4950654']
2026-06-22 17:25:39 INFO [TGCC-IRENE] Submitted job with ID:['4950654']
2026-06-22 17:25:39 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:25:39 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10
2026-06-22 17:25:39 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS10/end.2020020500_24_ENS10.nc
2026-06-22 17:25:39 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-06-22 17:25:39 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:25:39 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-06-22 17:25:39 INFO Queuing job for member 10...
2026-06-22 17:25:39 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:25:39 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-06-22 17:25:39 INFO Found: ['4950660']
2026-06-22 17:25:44 INFO [TGCC-IRENE] Submitted job with ID:['4950660']
2026-06-22 17:25:44 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:25:44 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11
2026-06-22 17:25:44 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS11/end.2020020500_24_ENS11.nc
2026-06-22 17:25:44 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-06-22 17:25:44 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:25:44 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-06-22 17:25:44 INFO Queuing job for member 11...
2026-06-22 17:25:44 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:25:44 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-06-22 17:25:45 INFO Found: ['4950665']
2026-06-22 17:25:50 INFO [TGCC-IRENE] Submitted job with ID:['4950665']
2026-06-22 17:25:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:25:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12
2026-06-22 17:25:50 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS12/end.2020020500_24_ENS12.nc
2026-06-22 17:25:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-06-22 17:25:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:25:50 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-06-22 17:25:50 INFO Queuing job for member 12...
2026-06-22 17:25:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:25:50 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-06-22 17:25:51 INFO Found: ['4950669']
2026-06-22 17:25:56 INFO [TGCC-IRENE] Submitted job with ID:['4950669']
2026-06-22 17:25:56 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:25:56 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13
2026-06-22 17:25:56 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS13/end.2020020500_24_ENS13.nc
2026-06-22 17:25:56 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-06-22 17:25:56 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:25:56 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-06-22 17:25:56 INFO Queuing job for member 13...
2026-06-22 17:25:56 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:25:56 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-06-22 17:25:57 INFO Found: ['4950671']
2026-06-22 17:26:02 INFO [TGCC-IRENE] Submitted job with ID:['4950671']
2026-06-22 17:26:02 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:26:02 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14
2026-06-22 17:26:02 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS14/end.2020020500_24_ENS14.nc
2026-06-22 17:26:02 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-06-22 17:26:02 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:26:02 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-06-22 17:26:02 INFO Queuing job for member 14...
2026-06-22 17:26:02 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:26:02 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-06-22 17:26:04 INFO Found: ['4950678']
2026-06-22 17:26:09 INFO [TGCC-IRENE] Submitted job with ID:['4950678']
2026-06-22 17:26:09 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:26:09 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15
2026-06-22 17:26:09 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmpemis_0615_15m_low_v2/ENS15/end.2020020500_24_ENS15.nc
2026-06-22 17:26:09 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-06-22 17:26:09 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:26:09 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-06-22 17:26:09 INFO Queuing job for member 15...
2026-06-22 17:26:09 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:26:09 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-06-22 17:26:10 INFO Found: ['4950682']
2026-06-22 17:26:15 INFO [TGCC-IRENE] Submitted job with ID:['4950682']
2026-06-22 17:26:15 INFO Checking job status ...
2026-06-22 17:26:15 INFO None 4950624: status RUNNING/PENDING
2026-06-22 17:26:15 INFO None 4950627: status RUNNING/PENDING
2026-06-22 17:26:15 INFO None 4950630: status RUNNING/PENDING
2026-06-22 17:26:15 INFO None 4950634: status RUNNING/PENDING
2026-06-22 17:26:15 INFO None 4950637: status RUNNING/PENDING
2026-06-22 17:26:15 INFO None 4950641: status RUNNING/PENDING
2026-06-22 17:26:15 INFO None 4950644: status RUNNING/PENDING
2026-06-22 17:26:15 INFO None 4950649: status RUNNING/PENDING
2026-06-22 17:26:15 INFO None 4950654: status RUNNING/PENDING
2026-06-22 17:26:15 INFO None 4950660: status RUNNING/PENDING
2026-06-22 17:26:15 INFO None 4950665: status RUNNING/PENDING
2026-06-22 17:26:15 INFO None 4950669: status RUNNING/PENDING
2026-06-22 17:26:15 INFO None 4950671: status RUNNING/PENDING
2026-06-22 17:26:15 INFO None 4950678: status RUNNING/PENDING
2026-06-22 17:26:15 INFO None 4950682: status RUNNING/PENDING
2026-06-22 17:26:15 INFO Jobs still running: ['4950624', '4950627', '4950630', '4950634', '4950637', '4950641', '4950644', '4950649', '4950654', '4950660', '4950665', '4950669', '4950671', '4950678', '4950682']. Waiting...
2026-06-22 17:26:30 INFO None 4950624: status RUNNING/PENDING
2026-06-22 17:26:30 INFO None 4950627: status RUNNING/PENDING
2026-06-22 17:26:30 INFO None 4950630: status RUNNING/PENDING
2026-06-22 17:26:30 INFO None 4950634: status FINISHED
2026-06-22 17:26:30 INFO None 4950637: status RUNNING/PENDING
2026-06-22 17:26:31 INFO None 4950641: status RUNNING/PENDING
2026-06-22 17:26:31 INFO None 4950644: status RUNNING/PENDING
2026-06-22 17:26:31 INFO None 4950649: status RUNNING/PENDING
2026-06-22 17:26:31 INFO None 4950654: status RUNNING/PENDING
2026-06-22 17:26:31 INFO None 4950660: status RUNNING/PENDING
2026-06-22 17:26:31 INFO None 4950665: status RUNNING/PENDING
2026-06-22 17:26:31 INFO None 4950669: status RUNNING/PENDING
2026-06-22 17:26:31 INFO None 4950671: status RUNNING/PENDING
2026-06-22 17:26:31 INFO None 4950678: status RUNNING/PENDING
2026-06-22 17:26:31 INFO None 4950682: status RUNNING/PENDING
2026-06-22 17:26:31 INFO Jobs still running: ['4950624', '4950627', '4950630', '4950637', '4950641', '4950644', '4950649', '4950654', '4950660', '4950665', '4950669', '4950671', '4950678', '4950682']. Waiting...
2026-06-22 17:26:46 INFO None 4950624: status RUNNING/PENDING
2026-06-22 17:26:46 INFO None 4950627: status RUNNING/PENDING
2026-06-22 17:26:46 INFO None 4950630: status RUNNING/PENDING
2026-06-22 17:26:46 INFO None 4950634: status FINISHED
2026-06-22 17:26:46 INFO None 4950637: status RUNNING/PENDING
2026-06-22 17:26:46 INFO None 4950641: status RUNNING/PENDING
2026-06-22 17:26:46 INFO None 4950644: status RUNNING/PENDING
2026-06-22 17:26:46 INFO None 4950649: status RUNNING/PENDING
2026-06-22 17:26:46 INFO None 4950654: status RUNNING/PENDING
2026-06-22 17:26:46 INFO None 4950660: status RUNNING/PENDING
2026-06-22 17:26:46 INFO None 4950665: status RUNNING/PENDING
2026-06-22 17:26:46 INFO None 4950669: status RUNNING/PENDING
2026-06-22 17:26:46 INFO None 4950671: status RUNNING/PENDING
2026-06-22 17:26:46 INFO None 4950678: status RUNNING/PENDING
2026-06-22 17:26:46 INFO None 4950682: status RUNNING/PENDING
2026-06-22 17:26:46 INFO Jobs still running: ['4950624', '4950627', '4950630', '4950637', '4950641', '4950644', '4950649', '4950654', '4950660', '4950665', '4950669', '4950671', '4950678', '4950682']. Waiting...
2026-06-22 17:27:01 INFO None 4950624: status RUNNING/PENDING
2026-06-22 17:27:01 INFO None 4950627: status RUNNING/PENDING
2026-06-22 17:27:01 INFO None 4950630: status RUNNING/PENDING
2026-06-22 17:27:01 INFO None 4950634: status FINISHED
2026-06-22 17:27:01 INFO None 4950637: status FINISHED
2026-06-22 17:27:01 INFO None 4950641: status FINISHED
2026-06-22 17:27:01 INFO None 4950644: status RUNNING/PENDING
2026-06-22 17:27:01 INFO None 4950649: status RUNNING/PENDING
2026-06-22 17:27:01 INFO None 4950654: status RUNNING/PENDING
2026-06-22 17:27:01 INFO None 4950660: status RUNNING/PENDING
2026-06-22 17:27:01 INFO None 4950665: status RUNNING/PENDING
2026-06-22 17:27:01 INFO None 4950669: status RUNNING/PENDING
2026-06-22 17:27:01 INFO None 4950671: status FINISHED
2026-06-22 17:27:01 INFO None 4950678: status RUNNING/PENDING
2026-06-22 17:27:01 INFO None 4950682: status RUNNING/PENDING
2026-06-22 17:27:01 INFO Jobs still running: ['4950624', '4950627', '4950630', '4950644', '4950649', '4950654', '4950660', '4950665', '4950669', '4950678', '4950682']. Waiting...
2026-06-22 17:27:16 INFO None 4950624: status FINISHED
2026-06-22 17:27:16 INFO None 4950627: status RUNNING/PENDING
2026-06-22 17:27:16 INFO None 4950630: status RUNNING/PENDING
2026-06-22 17:27:16 INFO None 4950634: status FINISHED
2026-06-22 17:27:17 INFO None 4950637: status FINISHED
2026-06-22 17:27:17 INFO None 4950641: status FINISHED
2026-06-22 17:27:17 INFO None 4950644: status RUNNING/PENDING
2026-06-22 17:27:17 INFO None 4950649: status RUNNING/PENDING
2026-06-22 17:27:17 INFO None 4950654: status RUNNING/PENDING
2026-06-22 17:27:17 INFO None 4950660: status RUNNING/PENDING
2026-06-22 17:27:17 INFO None 4950665: status RUNNING/PENDING
2026-06-22 17:27:17 INFO None 4950669: status RUNNING/PENDING
2026-06-22 17:27:17 INFO None 4950671: status FINISHED
2026-06-22 17:27:17 INFO None 4950678: status RUNNING/PENDING
2026-06-22 17:27:17 INFO None 4950682: status RUNNING/PENDING
2026-06-22 17:27:17 INFO Jobs still running: ['4950627', '4950630', '4950644', '4950649', '4950654', '4950660', '4950665', '4950669', '4950678', '4950682']. Waiting...
2026-06-22 17:27:32 INFO None 4950624: status FINISHED
2026-06-22 17:27:32 INFO None 4950627: status RUNNING/PENDING
2026-06-22 17:27:32 INFO None 4950630: status RUNNING/PENDING
2026-06-22 17:27:32 INFO None 4950634: status FINISHED
2026-06-22 17:27:32 INFO None 4950637: status FINISHED
2026-06-22 17:27:32 INFO None 4950641: status FINISHED
2026-06-22 17:27:32 INFO None 4950644: status FINISHED
2026-06-22 17:27:32 INFO None 4950649: status FINISHED
2026-06-22 17:27:32 INFO None 4950654: status RUNNING/PENDING
2026-06-22 17:27:32 INFO None 4950660: status RUNNING/PENDING
2026-06-22 17:27:32 INFO None 4950665: status RUNNING/PENDING
2026-06-22 17:27:32 INFO None 4950669: status RUNNING/PENDING
2026-06-22 17:27:32 INFO None 4950671: status FINISHED
2026-06-22 17:27:32 INFO None 4950678: status RUNNING/PENDING
2026-06-22 17:27:32 INFO None 4950682: status RUNNING/PENDING
2026-06-22 17:27:32 INFO Jobs still running: ['4950627', '4950630', '4950654', '4950660', '4950665', '4950669', '4950678', '4950682']. Waiting...
2026-06-22 17:27:47 INFO None 4950624: status FINISHED
2026-06-22 17:27:47 INFO None 4950627: status RUNNING/PENDING
2026-06-22 17:27:47 INFO None 4950630: status RUNNING/PENDING
2026-06-22 17:27:47 INFO None 4950634: status FINISHED
2026-06-22 17:27:47 INFO None 4950637: status FINISHED
2026-06-22 17:27:47 INFO None 4950641: status FINISHED
2026-06-22 17:27:47 INFO None 4950644: status FINISHED
2026-06-22 17:27:47 INFO None 4950649: status FINISHED
2026-06-22 17:27:47 INFO None 4950654: status FINISHED
2026-06-22 17:27:47 INFO None 4950660: status FINISHED
2026-06-22 17:27:47 INFO None 4950665: status FINISHED
2026-06-22 17:27:47 INFO None 4950669: status RUNNING/PENDING
2026-06-22 17:27:47 INFO None 4950671: status FINISHED
2026-06-22 17:27:47 INFO None 4950678: status RUNNING/PENDING
2026-06-22 17:27:47 INFO None 4950682: status RUNNING/PENDING
2026-06-22 17:27:47 INFO Jobs still running: ['4950627', '4950630', '4950669', '4950678', '4950682']. Waiting...
2026-06-22 17:28:02 INFO None 4950624: status FINISHED
2026-06-22 17:28:02 INFO None 4950627: status RUNNING/PENDING
2026-06-22 17:28:02 INFO None 4950630: status RUNNING/PENDING
2026-06-22 17:28:02 INFO None 4950634: status FINISHED
2026-06-22 17:28:02 INFO None 4950637: status FINISHED
2026-06-22 17:28:02 INFO None 4950641: status FINISHED
2026-06-22 17:28:02 INFO None 4950644: status FINISHED
2026-06-22 17:28:02 INFO None 4950649: status FINISHED
2026-06-22 17:28:03 INFO None 4950654: status FINISHED
2026-06-22 17:28:03 INFO None 4950660: status FINISHED
2026-06-22 17:28:03 INFO None 4950665: status FINISHED
2026-06-22 17:28:03 INFO None 4950669: status RUNNING/PENDING
2026-06-22 17:28:03 INFO None 4950671: status FINISHED
2026-06-22 17:28:03 INFO None 4950678: status RUNNING/PENDING
2026-06-22 17:28:03 INFO None 4950682: status RUNNING/PENDING
2026-06-22 17:28:03 INFO Jobs still running: ['4950627', '4950630', '4950669', '4950678', '4950682']. Waiting...
2026-06-22 17:28:18 INFO None 4950624: status FINISHED
2026-06-22 17:28:18 INFO None 4950627: status RUNNING/PENDING
2026-06-22 17:28:18 INFO None 4950630: status RUNNING/PENDING
2026-06-22 17:28:18 INFO None 4950634: status FINISHED
2026-06-22 17:28:18 INFO None 4950637: status FINISHED
2026-06-22 17:28:18 INFO None 4950641: status FINISHED
2026-06-22 17:28:18 INFO None 4950644: status FINISHED
2026-06-22 17:28:18 INFO None 4950649: status FINISHED
2026-06-22 17:28:18 INFO None 4950654: status FINISHED
2026-06-22 17:28:18 INFO None 4950660: status FINISHED
2026-06-22 17:28:18 INFO None 4950665: status FINISHED
2026-06-22 17:28:18 INFO None 4950669: status RUNNING/PENDING
2026-06-22 17:28:18 INFO None 4950671: status FINISHED
2026-06-22 17:28:18 INFO None 4950678: status RUNNING/PENDING
2026-06-22 17:28:18 INFO None 4950682: status RUNNING/PENDING
2026-06-22 17:28:18 INFO Jobs still running: ['4950627', '4950630', '4950669', '4950678', '4950682']. Waiting...
2026-06-22 17:28:33 INFO None 4950624: status FINISHED
2026-06-22 17:28:33 INFO None 4950627: status RUNNING/PENDING
2026-06-22 17:28:33 INFO None 4950630: status RUNNING/PENDING
2026-06-22 17:28:33 INFO None 4950634: status FINISHED
2026-06-22 17:28:33 INFO None 4950637: status FINISHED
2026-06-22 17:28:33 INFO None 4950641: status FINISHED
2026-06-22 17:28:33 INFO None 4950644: status FINISHED
2026-06-22 17:28:33 INFO None 4950649: status FINISHED
2026-06-22 17:28:33 INFO None 4950654: status FINISHED
2026-06-22 17:28:33 INFO None 4950660: status FINISHED
2026-06-22 17:28:33 INFO None 4950665: status FINISHED
2026-06-22 17:28:33 INFO None 4950669: status FINISHED
2026-06-22 17:28:33 INFO None 4950671: status FINISHED
2026-06-22 17:28:33 INFO None 4950678: status RUNNING/PENDING
2026-06-22 17:28:33 INFO None 4950682: status RUNNING/PENDING
2026-06-22 17:28:33 INFO Jobs still running: ['4950627', '4950630', '4950678', '4950682']. Waiting...
2026-06-22 17:28:48 INFO None 4950624: status FINISHED
2026-06-22 17:28:48 INFO None 4950627: status FINISHED
2026-06-22 17:28:48 INFO None 4950630: status FINISHED
2026-06-22 17:28:48 INFO None 4950634: status FINISHED
2026-06-22 17:28:48 INFO None 4950637: status FINISHED
2026-06-22 17:28:48 INFO None 4950641: status FINISHED
2026-06-22 17:28:48 INFO None 4950644: status FINISHED
2026-06-22 17:28:48 INFO None 4950649: status FINISHED
2026-06-22 17:28:48 INFO None 4950654: status FINISHED
2026-06-22 17:28:48 INFO None 4950660: status FINISHED
2026-06-22 17:28:48 INFO None 4950665: status FINISHED
2026-06-22 17:28:48 INFO None 4950669: status FINISHED
2026-06-22 17:28:48 INFO None 4950671: status FINISHED
2026-06-22 17:28:48 INFO None 4950678: status RUNNING/PENDING
2026-06-22 17:28:48 INFO None 4950682: status RUNNING/PENDING
2026-06-22 17:28:48 INFO Jobs still running: ['4950678', '4950682']. Waiting...
2026-06-22 17:29:04 INFO None 4950624: status FINISHED
2026-06-22 17:29:04 INFO None 4950627: status FINISHED
2026-06-22 17:29:04 INFO None 4950630: status FINISHED
2026-06-22 17:29:04 INFO None 4950634: status FINISHED
2026-06-22 17:29:04 INFO None 4950637: status FINISHED
2026-06-22 17:29:04 INFO None 4950641: status FINISHED
2026-06-22 17:29:04 INFO None 4950644: status FINISHED
2026-06-22 17:29:04 INFO None 4950649: status FINISHED
2026-06-22 17:29:04 INFO None 4950654: status FINISHED
2026-06-22 17:29:04 INFO None 4950660: status FINISHED
2026-06-22 17:29:04 INFO None 4950665: status FINISHED
2026-06-22 17:29:04 INFO None 4950669: status FINISHED
2026-06-22 17:29:04 INFO None 4950671: status FINISHED
2026-06-22 17:29:04 INFO None 4950678: status FINISHED
2026-06-22 17:29:04 INFO None 4950682: status RUNNING/PENDING
2026-06-22 17:29:04 INFO Jobs still running: ['4950682']. Waiting...
2026-06-22 17:29:19 INFO None 4950624: status FINISHED
2026-06-22 17:29:19 INFO None 4950627: status FINISHED
2026-06-22 17:29:19 INFO None 4950630: status FINISHED
2026-06-22 17:29:19 INFO None 4950634: status FINISHED
2026-06-22 17:29:19 INFO None 4950637: status FINISHED
2026-06-22 17:29:19 INFO None 4950641: status FINISHED
2026-06-22 17:29:19 INFO None 4950644: status FINISHED
2026-06-22 17:29:19 INFO None 4950649: status FINISHED
2026-06-22 17:29:19 INFO None 4950654: status FINISHED
2026-06-22 17:29:19 INFO None 4950660: status FINISHED
2026-06-22 17:29:19 INFO None 4950665: status FINISHED
2026-06-22 17:29:19 INFO None 4950669: status FINISHED
2026-06-22 17:29:19 INFO None 4950671: status FINISHED
2026-06-22 17:29:19 INFO None 4950678: status FINISHED
2026-06-22 17:29:19 INFO None 4950682: status RUNNING/PENDING
2026-06-22 17:29:19 INFO Jobs still running: ['4950682']. Waiting...
2026-06-22 17:29:34 INFO None 4950624: status FINISHED
2026-06-22 17:29:34 INFO None 4950627: status FINISHED
2026-06-22 17:29:34 INFO None 4950630: status FINISHED
2026-06-22 17:29:34 INFO None 4950634: status FINISHED
2026-06-22 17:29:34 INFO None 4950637: status FINISHED
2026-06-22 17:29:34 INFO None 4950641: status FINISHED
2026-06-22 17:29:34 INFO None 4950644: status FINISHED
2026-06-22 17:29:34 INFO None 4950649: status FINISHED
2026-06-22 17:29:34 INFO None 4950654: status FINISHED
2026-06-22 17:29:34 INFO None 4950660: status FINISHED
2026-06-22 17:29:34 INFO None 4950665: status FINISHED
2026-06-22 17:29:34 INFO None 4950669: status FINISHED
2026-06-22 17:29:34 INFO None 4950671: status FINISHED
2026-06-22 17:29:34 INFO None 4950678: status FINISHED
2026-06-22 17:29:34 INFO None 4950682: status RUNNING/PENDING
2026-06-22 17:29:34 INFO Jobs still running: ['4950682']. Waiting...
[2026-06-22T17:29:45.166] error: *** JOB 4950603 ON irene4869 CANCELLED AT 2026-06-22T17:29:45 DUE to SIGNAL Terminated ***
