+ SCRIPT_PID=3848280
+ /bin/bash -x /tmp/tmp.E5FRHSSlK9
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
2026-06-22 16:54:20 INFO 
███    ███ ██ ███    ███ ███████ ███████ ██
████  ████ ██ ████  ████ ██      ██      ██
██ ████ ██ ██ ██ ████ ██ █████   ███████ ██
██  ██  ██ ██ ██  ██  ██ ██           ██ ██
██      ██ ██ ██      ██ ███████ ███████ ██



2026-06-22 16:54:20 INFO [PIPELINE] =======================================
2026-06-22 16:54:20 INFO [PIPELINE] Starting chimere–DART orchestrator
2026-06-22 16:54:20 INFO [PIPELINE] Config file: config/config_irene_IM.yaml
2026-06-22 16:54:20 INFO [PIPELINE] Run dir: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart
2026-06-22 16:54:20 INFO [PIPELINE] Log file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/mimesi_orchestrator_logs/chimere_DART_20260622_165420.log
2026-06-22 16:54:20 INFO [PIPELINE] =======================================
2026-06-22 16:54:20 INFO Running assimilation with model_type=ModelType.CHIMERE
2026-06-22 16:54:20 INFO Using scheduler=Scheduler.SLURM, queue=rome
2026-06-22 16:54:20 INFO [STEP] ---- TIME LOOP START ----
2026-06-22 16:54:20 INFO [TIME] step_start current_time=2020-02-06 09:00:00 simulated_time=None dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-22 16:54:20 INFO [TIME] window start=2020-02-06 09:00:00 end=2020-02-06 11:00:00 run_hours=2 has_assimilation=True
2026-06-22 16:54:20 INFO ---------->>> Running CHIMERE model from 2020-02-06 09:00:00 to 2020-02-06 11:00:00
2026-06-22 16:54:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 16:54:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1
2026-06-22 16:54:20 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/end.2020020601_9_ENS1.nc
2026-06-22 16:54:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-06-22 16:54:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 16:54:20 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-06-22 16:54:21 INFO Queuing job for member 1...
2026-06-22 16:54:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 16:54:21 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-06-22 16:54:21 INFO Found: ['4949497']
2026-06-22 16:54:26 INFO [TGCC-IRENE] Submitted job with ID:['4949497']
2026-06-22 16:54:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 16:54:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2
2026-06-22 16:54:26 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/end.2020020601_9_ENS2.nc
2026-06-22 16:54:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-06-22 16:54:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 16:54:27 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-06-22 16:54:27 INFO Queuing job for member 2...
2026-06-22 16:54:27 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 16:54:27 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-06-22 16:54:27 INFO Found: ['4949500']
2026-06-22 16:54:32 INFO [TGCC-IRENE] Submitted job with ID:['4949500']
2026-06-22 16:54:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 16:54:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3
2026-06-22 16:54:32 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/end.2020020601_9_ENS3.nc
2026-06-22 16:54:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-06-22 16:54:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 16:54:32 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-06-22 16:54:32 INFO Queuing job for member 3...
2026-06-22 16:54:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 16:54:32 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-06-22 16:54:33 INFO Found: ['4949502']
2026-06-22 16:54:38 INFO [TGCC-IRENE] Submitted job with ID:['4949502']
2026-06-22 16:54:38 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 16:54:38 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4
2026-06-22 16:54:38 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/end.2020020601_9_ENS4.nc
2026-06-22 16:54:38 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-06-22 16:54:38 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 16:54:38 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-06-22 16:54:38 INFO Queuing job for member 4...
2026-06-22 16:54:38 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 16:54:38 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-06-22 16:54:39 INFO Found: ['4949504']
2026-06-22 16:54:44 INFO [TGCC-IRENE] Submitted job with ID:['4949504']
2026-06-22 16:54:44 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 16:54:44 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5
2026-06-22 16:54:44 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/end.2020020601_9_ENS5.nc
2026-06-22 16:54:44 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-06-22 16:54:44 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 16:54:44 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-06-22 16:54:44 INFO Queuing job for member 5...
2026-06-22 16:54:44 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 16:54:44 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-06-22 16:54:45 INFO Found: ['4949506']
2026-06-22 16:54:50 INFO [TGCC-IRENE] Submitted job with ID:['4949506']
2026-06-22 16:54:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 16:54:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6
2026-06-22 16:54:50 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/end.2020020601_9_ENS6.nc
2026-06-22 16:54:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-06-22 16:54:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 16:54:50 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-06-22 16:54:50 INFO Queuing job for member 6...
2026-06-22 16:54:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 16:54:50 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-06-22 16:54:50 INFO Found: ['4949509']
2026-06-22 16:54:55 INFO [TGCC-IRENE] Submitted job with ID:['4949509']
2026-06-22 16:54:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 16:54:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7
2026-06-22 16:54:55 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/end.2020020601_9_ENS7.nc
2026-06-22 16:54:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-06-22 16:54:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 16:54:55 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-06-22 16:54:55 INFO Queuing job for member 7...
2026-06-22 16:54:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 16:54:55 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-06-22 16:54:56 INFO Found: ['4949511']
2026-06-22 16:55:01 INFO [TGCC-IRENE] Submitted job with ID:['4949511']
2026-06-22 16:55:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 16:55:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8
2026-06-22 16:55:01 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/end.2020020601_9_ENS8.nc
2026-06-22 16:55:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-06-22 16:55:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 16:55:01 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-06-22 16:55:01 INFO Queuing job for member 8...
2026-06-22 16:55:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 16:55:01 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-06-22 16:55:02 INFO Found: ['4949514']
2026-06-22 16:55:07 INFO [TGCC-IRENE] Submitted job with ID:['4949514']
2026-06-22 16:55:07 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 16:55:07 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9
2026-06-22 16:55:07 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/end.2020020601_9_ENS9.nc
2026-06-22 16:55:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-06-22 16:55:07 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 16:55:07 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-06-22 16:55:07 INFO Queuing job for member 9...
2026-06-22 16:55:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 16:55:07 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-06-22 16:55:08 INFO Found: ['4949516']
2026-06-22 16:55:13 INFO [TGCC-IRENE] Submitted job with ID:['4949516']
2026-06-22 16:55:13 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 16:55:13 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10
2026-06-22 16:55:13 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/end.2020020601_9_ENS10.nc
2026-06-22 16:55:13 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-06-22 16:55:13 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 16:55:13 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-06-22 16:55:13 INFO Queuing job for member 10...
2026-06-22 16:55:13 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 16:55:13 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-06-22 16:55:13 INFO Found: ['4949518']
2026-06-22 16:55:18 INFO [TGCC-IRENE] Submitted job with ID:['4949518']
2026-06-22 16:55:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 16:55:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11
2026-06-22 16:55:18 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/end.2020020601_9_ENS11.nc
2026-06-22 16:55:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-06-22 16:55:18 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 16:55:18 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-06-22 16:55:18 INFO Queuing job for member 11...
2026-06-22 16:55:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 16:55:18 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-06-22 16:55:19 INFO Found: ['4949522']
2026-06-22 16:55:24 INFO [TGCC-IRENE] Submitted job with ID:['4949522']
2026-06-22 16:55:24 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 16:55:24 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12
2026-06-22 16:55:24 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/end.2020020601_9_ENS12.nc
2026-06-22 16:55:24 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-06-22 16:55:24 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 16:55:24 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-06-22 16:55:24 INFO Queuing job for member 12...
2026-06-22 16:55:24 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 16:55:24 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-06-22 16:55:25 INFO Found: ['4949527']
2026-06-22 16:55:30 INFO [TGCC-IRENE] Submitted job with ID:['4949527']
2026-06-22 16:55:30 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 16:55:30 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13
2026-06-22 16:55:30 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/end.2020020601_9_ENS13.nc
2026-06-22 16:55:30 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-06-22 16:55:30 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 16:55:30 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-06-22 16:55:30 INFO Queuing job for member 13...
2026-06-22 16:55:30 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 16:55:30 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-06-22 16:55:31 INFO Found: ['4949531']
2026-06-22 16:55:36 INFO [TGCC-IRENE] Submitted job with ID:['4949531']
2026-06-22 16:55:36 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 16:55:36 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14
2026-06-22 16:55:36 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/end.2020020601_9_ENS14.nc
2026-06-22 16:55:36 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-06-22 16:55:36 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 16:55:36 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-06-22 16:55:36 INFO Queuing job for member 14...
2026-06-22 16:55:36 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 16:55:36 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-06-22 16:55:36 INFO Found: ['4949534']
2026-06-22 16:55:41 INFO [TGCC-IRENE] Submitted job with ID:['4949534']
2026-06-22 16:55:41 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 16:55:41 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15
2026-06-22 16:55:41 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/end.2020020601_9_ENS15.nc
2026-06-22 16:55:41 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-06-22 16:55:41 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 16:55:41 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-06-22 16:55:41 INFO Queuing job for member 15...
2026-06-22 16:55:41 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 16:55:41 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-06-22 16:55:42 INFO Found: ['4949536']
2026-06-22 16:55:47 INFO [TGCC-IRENE] Submitted job with ID:['4949536']
2026-06-22 16:55:47 INFO Checking job status ...
2026-06-22 16:55:47 INFO None 4949497: status RUNNING/PENDING
2026-06-22 16:55:47 INFO None 4949500: status RUNNING/PENDING
2026-06-22 16:55:47 INFO None 4949502: status RUNNING/PENDING
2026-06-22 16:55:47 INFO None 4949504: status RUNNING/PENDING
2026-06-22 16:55:47 INFO None 4949506: status RUNNING/PENDING
2026-06-22 16:55:47 INFO None 4949509: status RUNNING/PENDING
2026-06-22 16:55:47 INFO None 4949511: status RUNNING/PENDING
2026-06-22 16:55:47 INFO None 4949514: status RUNNING/PENDING
2026-06-22 16:55:47 INFO None 4949516: status RUNNING/PENDING
2026-06-22 16:55:48 INFO None 4949518: status RUNNING/PENDING
2026-06-22 16:55:48 INFO None 4949522: status RUNNING/PENDING
2026-06-22 16:55:48 INFO None 4949527: status RUNNING/PENDING
2026-06-22 16:55:51 INFO None 4949531: status RUNNING/PENDING
2026-06-22 16:55:51 INFO None 4949534: status RUNNING/PENDING
2026-06-22 16:55:51 INFO None 4949536: status RUNNING/PENDING
2026-06-22 16:55:51 INFO Jobs still running: ['4949497', '4949500', '4949502', '4949504', '4949506', '4949509', '4949511', '4949514', '4949516', '4949518', '4949522', '4949527', '4949531', '4949534', '4949536']. Waiting...
2026-06-22 16:56:06 INFO None 4949497: status FINISHED
2026-06-22 16:56:06 INFO None 4949500: status RUNNING/PENDING
2026-06-22 16:56:07 INFO None 4949502: status RUNNING/PENDING
2026-06-22 16:56:07 INFO None 4949504: status RUNNING/PENDING
2026-06-22 16:56:07 INFO None 4949506: status RUNNING/PENDING
2026-06-22 16:56:07 INFO None 4949509: status FINISHED
2026-06-22 16:56:07 INFO None 4949511: status FINISHED
2026-06-22 16:56:07 INFO None 4949514: status RUNNING/PENDING
2026-06-22 16:56:07 INFO None 4949516: status FINISHED
2026-06-22 16:56:07 INFO None 4949518: status FINISHED
2026-06-22 16:56:07 INFO None 4949522: status FINISHED
2026-06-22 16:56:07 INFO None 4949527: status FINISHED
2026-06-22 16:56:07 INFO None 4949531: status RUNNING/PENDING
2026-06-22 16:56:07 INFO None 4949534: status RUNNING/PENDING
2026-06-22 16:56:07 INFO None 4949536: status RUNNING/PENDING
2026-06-22 16:56:07 INFO Jobs still running: ['4949500', '4949502', '4949504', '4949506', '4949514', '4949531', '4949534', '4949536']. Waiting...
2026-06-22 16:56:22 INFO None 4949497: status FINISHED
2026-06-22 16:56:22 INFO None 4949500: status RUNNING/PENDING
2026-06-22 16:56:22 INFO None 4949502: status RUNNING/PENDING
2026-06-22 16:56:22 INFO None 4949504: status RUNNING/PENDING
2026-06-22 16:56:22 INFO None 4949506: status RUNNING/PENDING
2026-06-22 16:56:22 INFO None 4949509: status FINISHED
2026-06-22 16:56:22 INFO None 4949511: status FINISHED
2026-06-22 16:56:22 INFO None 4949514: status RUNNING/PENDING
2026-06-22 16:56:22 INFO None 4949516: status FINISHED
2026-06-22 16:56:22 INFO None 4949518: status FINISHED
2026-06-22 16:56:22 INFO None 4949522: status FINISHED
2026-06-22 16:56:22 INFO None 4949527: status FINISHED
2026-06-22 16:56:22 INFO None 4949531: status RUNNING/PENDING
2026-06-22 16:56:22 INFO None 4949534: status RUNNING/PENDING
2026-06-22 16:56:22 INFO None 4949536: status RUNNING/PENDING
2026-06-22 16:56:22 INFO Jobs still running: ['4949500', '4949502', '4949504', '4949506', '4949514', '4949531', '4949534', '4949536']. Waiting...
2026-06-22 16:56:37 INFO None 4949497: status FINISHED
2026-06-22 16:56:37 INFO None 4949500: status FINISHED
2026-06-22 16:56:37 INFO None 4949502: status RUNNING/PENDING
2026-06-22 16:56:37 INFO None 4949504: status RUNNING/PENDING
2026-06-22 16:56:37 INFO None 4949506: status RUNNING/PENDING
2026-06-22 16:56:37 INFO None 4949509: status FINISHED
2026-06-22 16:56:37 INFO None 4949511: status FINISHED
2026-06-22 16:56:37 INFO None 4949514: status RUNNING/PENDING
2026-06-22 16:56:37 INFO None 4949516: status FINISHED
2026-06-22 16:56:37 INFO None 4949518: status FINISHED
2026-06-22 16:56:37 INFO None 4949522: status FINISHED
2026-06-22 16:56:37 INFO None 4949527: status FINISHED
2026-06-22 16:56:37 INFO None 4949531: status FINISHED
2026-06-22 16:56:37 INFO None 4949534: status RUNNING/PENDING
2026-06-22 16:56:37 INFO None 4949536: status RUNNING/PENDING
2026-06-22 16:56:37 INFO Jobs still running: ['4949502', '4949504', '4949506', '4949514', '4949534', '4949536']. Waiting...
2026-06-22 16:56:52 INFO None 4949497: status FINISHED
2026-06-22 16:56:52 INFO None 4949500: status FINISHED
2026-06-22 16:56:52 INFO None 4949502: status FINISHED
2026-06-22 16:56:52 INFO None 4949504: status FINISHED
2026-06-22 16:56:52 INFO None 4949506: status RUNNING/PENDING
2026-06-22 16:56:52 INFO None 4949509: status FINISHED
2026-06-22 16:56:52 INFO None 4949511: status FINISHED
2026-06-22 16:56:52 INFO None 4949514: status RUNNING/PENDING
2026-06-22 16:56:52 INFO None 4949516: status FINISHED
2026-06-22 16:56:52 INFO None 4949518: status FINISHED
2026-06-22 16:56:53 INFO None 4949522: status FINISHED
2026-06-22 16:56:53 INFO None 4949527: status FINISHED
2026-06-22 16:56:53 INFO None 4949531: status FINISHED
2026-06-22 16:56:53 INFO None 4949534: status RUNNING/PENDING
2026-06-22 16:56:53 INFO None 4949536: status RUNNING/PENDING
2026-06-22 16:56:53 INFO Jobs still running: ['4949506', '4949514', '4949534', '4949536']. Waiting...
2026-06-22 16:57:08 INFO None 4949497: status FINISHED
2026-06-22 16:57:08 INFO None 4949500: status FINISHED
2026-06-22 16:57:08 INFO None 4949502: status FINISHED
2026-06-22 16:57:08 INFO None 4949504: status FINISHED
2026-06-22 16:57:08 INFO None 4949506: status RUNNING/PENDING
2026-06-22 16:57:08 INFO None 4949509: status FINISHED
2026-06-22 16:57:08 INFO None 4949511: status FINISHED
2026-06-22 16:57:08 INFO None 4949514: status RUNNING/PENDING
2026-06-22 16:57:08 INFO None 4949516: status FINISHED
2026-06-22 16:57:08 INFO None 4949518: status FINISHED
2026-06-22 16:57:08 INFO None 4949522: status FINISHED
2026-06-22 16:57:08 INFO None 4949527: status FINISHED
2026-06-22 16:57:08 INFO None 4949531: status FINISHED
2026-06-22 16:57:08 INFO None 4949534: status FINISHED
2026-06-22 16:57:08 INFO None 4949536: status FINISHED
2026-06-22 16:57:08 INFO Jobs still running: ['4949506', '4949514']. Waiting...
2026-06-22 16:57:23 INFO None 4949497: status FINISHED
2026-06-22 16:57:23 INFO None 4949500: status FINISHED
2026-06-22 16:57:23 INFO None 4949502: status FINISHED
2026-06-22 16:57:23 INFO None 4949504: status FINISHED
2026-06-22 16:57:23 INFO None 4949506: status FINISHED
2026-06-22 16:57:23 INFO None 4949509: status FINISHED
2026-06-22 16:57:23 INFO None 4949511: status FINISHED
2026-06-22 16:57:23 INFO None 4949514: status RUNNING/PENDING
2026-06-22 16:57:23 INFO None 4949516: status FINISHED
2026-06-22 16:57:23 INFO None 4949518: status FINISHED
2026-06-22 16:57:23 INFO None 4949522: status FINISHED
2026-06-22 16:57:23 INFO None 4949527: status FINISHED
2026-06-22 16:57:23 INFO None 4949531: status FINISHED
2026-06-22 16:57:23 INFO None 4949534: status FINISHED
2026-06-22 16:57:23 INFO None 4949536: status FINISHED
2026-06-22 16:57:23 INFO Jobs still running: ['4949514']. Waiting...
2026-06-22 16:57:38 INFO None 4949497: status FINISHED
2026-06-22 16:57:38 INFO None 4949500: status FINISHED
2026-06-22 16:57:38 INFO None 4949502: status FINISHED
2026-06-22 16:57:38 INFO None 4949504: status FINISHED
2026-06-22 16:57:38 INFO None 4949506: status FINISHED
2026-06-22 16:57:38 INFO None 4949509: status FINISHED
2026-06-22 16:57:38 INFO None 4949511: status FINISHED
2026-06-22 16:57:38 INFO None 4949514: status RUNNING/PENDING
2026-06-22 16:57:38 INFO None 4949516: status FINISHED
2026-06-22 16:57:38 INFO None 4949518: status FINISHED
2026-06-22 16:57:38 INFO None 4949522: status FINISHED
2026-06-22 16:57:38 INFO None 4949527: status FINISHED
2026-06-22 16:57:38 INFO None 4949531: status FINISHED
2026-06-22 16:57:39 INFO None 4949534: status FINISHED
2026-06-22 16:57:39 INFO None 4949536: status FINISHED
2026-06-22 16:57:39 INFO Jobs still running: ['4949514']. Waiting...
2026-06-22 16:57:54 INFO None 4949497: status FINISHED
2026-06-22 16:57:54 INFO None 4949500: status FINISHED
2026-06-22 16:57:54 INFO None 4949502: status FINISHED
2026-06-22 16:57:54 INFO None 4949504: status FINISHED
2026-06-22 16:57:54 INFO None 4949506: status FINISHED
2026-06-22 16:57:54 INFO None 4949509: status FINISHED
2026-06-22 16:57:54 INFO None 4949511: status FINISHED
2026-06-22 16:57:54 INFO None 4949514: status FINISHED
2026-06-22 16:57:54 INFO None 4949516: status FINISHED
2026-06-22 16:57:54 INFO None 4949518: status FINISHED
2026-06-22 16:57:54 INFO None 4949522: status FINISHED
2026-06-22 16:57:54 INFO None 4949527: status FINISHED
2026-06-22 16:57:54 INFO None 4949531: status FINISHED
2026-06-22 16:57:54 INFO None 4949534: status FINISHED
2026-06-22 16:57:54 INFO None 4949536: status FINISHED
2026-06-22 16:57:54 INFO Jobs ['4949497', '4949500', '4949502', '4949504', '4949506', '4949509', '4949511', '4949514', '4949516', '4949518', '4949522', '4949527', '4949531', '4949534', '4949536'] have finished
2026-06-22 16:57:54 INFO Checking restart files were created ...
2026-06-22 16:57:54 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc(1002881383 bytes)
2026-06-22 16:57:54 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc(1002881383 bytes)
2026-06-22 16:57:54 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc(1002881383 bytes)
2026-06-22 16:57:54 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc(1002881383 bytes)
2026-06-22 16:57:54 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc(1002881383 bytes)
2026-06-22 16:57:54 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc(1002881383 bytes)
2026-06-22 16:57:54 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc(1002881383 bytes)
2026-06-22 16:57:54 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc(1002881383 bytes)
2026-06-22 16:57:54 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc(1002881383 bytes)
2026-06-22 16:57:54 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc(1002881383 bytes)
2026-06-22 16:57:54 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc(1002881383 bytes)
2026-06-22 16:57:54 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc(1002881383 bytes)
2026-06-22 16:57:54 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc(1002881383 bytes)
2026-06-22 16:57:54 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc(1002881383 bytes)
2026-06-22 16:57:54 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc(1002881383 bytes)
2026-06-22 16:57:54 INFO  Run_model() completed successfully.
2026-06-22 16:57:54 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-22 16:57:54 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 11:00:00 days=153072 seconds=39600
2026-06-22 16:57:54 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-06-22 16:57:54 INFO [TIME] increment current_time 2020-02-06 09:00:00 -> 2020-02-06 11:00:00
2026-06-22 16:57:54 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-22 16:57:54 INFO ---------->>> Running process_satellite_data()
2026-06-22 16:57:54 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12003.nc
2026-06-22 16:57:54 INFO ---------->>> Running run_obs_converter()
2026-06-22 16:57:54 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_39045_153072.out
2026-06-22 16:57:54 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_39045_153072.out
2026-06-22 16:57:54 INFO ---------->>> Running DART
2026-06-22 16:57:54 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-06-22 16:57:54 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/chim_ENS1_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/chim_ENS1_2020020611_1_out_toDART.nc
2026-06-22 16:58:03 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/chim_ENS2_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/chim_ENS2_2020020611_1_out_toDART.nc
2026-06-22 16:58:03 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/chim_ENS3_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/chim_ENS3_2020020611_1_out_toDART.nc
2026-06-22 16:58:04 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/chim_ENS4_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/chim_ENS4_2020020611_1_out_toDART.nc
2026-06-22 16:58:04 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/chim_ENS5_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/chim_ENS5_2020020611_1_out_toDART.nc
2026-06-22 16:58:04 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/chim_ENS6_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/chim_ENS6_2020020611_1_out_toDART.nc
2026-06-22 16:58:05 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/chim_ENS7_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/chim_ENS7_2020020611_1_out_toDART.nc
2026-06-22 16:58:05 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/chim_ENS8_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/chim_ENS8_2020020611_1_out_toDART.nc
2026-06-22 16:58:06 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/chim_ENS9_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/chim_ENS9_2020020611_1_out_toDART.nc
2026-06-22 16:58:06 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/chim_ENS10_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/chim_ENS10_2020020611_1_out_toDART.nc
2026-06-22 16:58:06 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/chim_ENS11_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/chim_ENS11_2020020611_1_out_toDART.nc
2026-06-22 16:58:07 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/chim_ENS12_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/chim_ENS12_2020020611_1_out_toDART.nc
2026-06-22 16:58:07 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/chim_ENS13_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/chim_ENS13_2020020611_1_out_toDART.nc
2026-06-22 16:58:07 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/chim_ENS14_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/chim_ENS14_2020020611_1_out_toDART.nc
2026-06-22 16:58:08 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/chim_ENS15_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/chim_ENS15_2020020611_1_out_toDART.nc
2026-06-22 16:58:08 INFO Replacement input_template.nml → input.nml completed successfully.
2026-06-22 16:58:08 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-06-22 16:58:08 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-06-22 16:58:08 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-06-22 16:58:08 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-06-22 16:58:08 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-06-22 16:58:22 INFO Found: []
2026-06-22 16:58:22 INFO No job id returned by command ./run_filter.bsh
2026-06-22 16:58:22 INFO No monitoring will be performed
2026-06-22 16:58:22 INFO Moving DART output files to analysis and preassim directories for date 2020020611 if present ...
2026-06-22 16:58:22 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020020611'
2026-06-22 16:58:22 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020020611'
2026-06-22 16:58:22 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020020611'
2026-06-22 16:58:22 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020020611'
2026-06-22 16:58:22 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020020611'
2026-06-22 16:58:22 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020020611'
2026-06-22 16:58:22 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020020611'
2026-06-22 16:58:22 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020020611'
2026-06-22 16:58:22 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020020611'
2026-06-22 16:58:22 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020020611'
2026-06-22 16:58:22 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020020611'
2026-06-22 16:58:22 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020020611'
2026-06-22 16:58:22 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020020611'
2026-06-22 16:58:22 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020020611'
2026-06-22 16:58:22 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020020611'
2026-06-22 16:58:22 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020020611'
2026-06-22 16:58:22 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020020611'
2026-06-22 16:58:22 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020020611'
2026-06-22 16:58:22 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020020611'
2026-06-22 16:58:22 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020020611'
2026-06-22 16:58:22 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020020611'
2026-06-22 16:58:22 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020020611'
2026-06-22 16:58:22 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020020611'
2026-06-22 16:58:22 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020020611'
2026-06-22 16:58:22 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020020611'
2026-06-22 16:58:22 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020020611'
2026-06-22 16:58:22 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020020611'
2026-06-22 16:58:22 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020020611'
2026-06-22 16:58:22 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020020611'
2026-06-22 16:58:22 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020020611'
2026-06-22 16:58:22 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020020611'
2026-06-22 16:58:22 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/analysis/2020020611'
2026-06-22 16:58:22 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020020611'
2026-06-22 16:58:23 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/preassim/2020020611'
2026-06-22 16:58:23 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-06-22 16:58:25 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-06-22 16:58:25 INFO run_dart() is DONE.
2026-06-22 16:58:25 INFO ---------->>> Running update_pollutant_in_end()
2026-06-22 16:58:26 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 16:58:26 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 16:58:27 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 16:58:27 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 16:58:27 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc
2026-06-22 16:58:30 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 16:58:31 INFO [Posterior Diff ENS1] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS1_2020020611.nc
2026-06-22 16:58:31 INFO [Posterior Diff ENS1] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS1_2020020611.relative.nc
2026-06-22 16:58:31 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 16:58:31 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 16:58:32 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 16:58:32 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 16:58:33 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc
2026-06-22 16:58:36 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 16:58:36 INFO [Posterior Diff ENS2] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS2_2020020611.nc
2026-06-22 16:58:36 INFO [Posterior Diff ENS2] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS2_2020020611.relative.nc
2026-06-22 16:58:37 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 16:58:37 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 16:58:37 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 16:58:37 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 16:58:38 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc
2026-06-22 16:58:41 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 16:58:41 INFO [Posterior Diff ENS3] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS3_2020020611.nc
2026-06-22 16:58:41 INFO [Posterior Diff ENS3] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS3_2020020611.relative.nc
2026-06-22 16:58:42 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 16:58:42 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 16:58:43 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 16:58:43 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 16:58:43 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc
2026-06-22 16:58:47 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 16:58:47 INFO [Posterior Diff ENS4] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS4_2020020611.nc
2026-06-22 16:58:47 INFO [Posterior Diff ENS4] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS4_2020020611.relative.nc
2026-06-22 16:58:47 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 16:58:47 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 16:58:48 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 16:58:48 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 16:58:49 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc
2026-06-22 16:58:52 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 16:58:52 INFO [Posterior Diff ENS5] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS5_2020020611.nc
2026-06-22 16:58:52 INFO [Posterior Diff ENS5] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS5_2020020611.relative.nc
2026-06-22 16:58:53 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 16:58:53 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 16:58:53 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 16:58:53 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 16:58:54 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc
2026-06-22 16:58:57 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 16:58:58 INFO [Posterior Diff ENS6] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS6_2020020611.nc
2026-06-22 16:58:58 INFO [Posterior Diff ENS6] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS6_2020020611.relative.nc
2026-06-22 16:58:58 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 16:58:58 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 16:58:59 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 16:58:59 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 16:59:00 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc
2026-06-22 16:59:03 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 16:59:03 INFO [Posterior Diff ENS7] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS7_2020020611.nc
2026-06-22 16:59:03 INFO [Posterior Diff ENS7] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS7_2020020611.relative.nc
2026-06-22 16:59:04 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 16:59:04 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 16:59:04 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 16:59:04 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 16:59:05 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc
2026-06-22 16:59:08 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 16:59:09 INFO [Posterior Diff ENS8] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS8_2020020611.nc
2026-06-22 16:59:09 INFO [Posterior Diff ENS8] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS8_2020020611.relative.nc
2026-06-22 16:59:09 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 16:59:09 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 16:59:10 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 16:59:10 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 16:59:11 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc
2026-06-22 16:59:14 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 16:59:15 INFO [Posterior Diff ENS9] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS9_2020020611.nc
2026-06-22 16:59:15 INFO [Posterior Diff ENS9] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS9_2020020611.relative.nc
2026-06-22 16:59:15 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 16:59:15 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 16:59:16 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 16:59:16 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 16:59:16 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc
2026-06-22 16:59:20 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 16:59:20 INFO [Posterior Diff ENS10] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS10_2020020611.nc
2026-06-22 16:59:20 INFO [Posterior Diff ENS10] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS10_2020020611.relative.nc
2026-06-22 16:59:21 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 16:59:21 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 16:59:21 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 16:59:21 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 16:59:22 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc
2026-06-22 16:59:25 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 16:59:26 INFO [Posterior Diff ENS11] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS11_2020020611.nc
2026-06-22 16:59:26 INFO [Posterior Diff ENS11] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS11_2020020611.relative.nc
2026-06-22 16:59:26 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 16:59:26 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 16:59:27 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 16:59:27 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 16:59:28 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc
2026-06-22 16:59:31 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 16:59:31 INFO [Posterior Diff ENS12] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS12_2020020611.nc
2026-06-22 16:59:31 INFO [Posterior Diff ENS12] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS12_2020020611.relative.nc
2026-06-22 16:59:31 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 16:59:31 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 16:59:32 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 16:59:32 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 16:59:33 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc
2026-06-22 16:59:36 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 16:59:37 INFO [Posterior Diff ENS13] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS13_2020020611.nc
2026-06-22 16:59:37 INFO [Posterior Diff ENS13] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS13_2020020611.relative.nc
2026-06-22 16:59:37 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 16:59:37 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 16:59:38 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 16:59:38 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 16:59:39 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc
2026-06-22 16:59:42 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 16:59:42 INFO [Posterior Diff ENS14] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS14_2020020611.nc
2026-06-22 16:59:42 INFO [Posterior Diff ENS14] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS14_2020020611.relative.nc
2026-06-22 16:59:43 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 16:59:43 INFO Scaled NO in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 16:59:44 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-06-22 16:59:44 INFO Scaled NO2 in /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc using posterior/prior ratio.
2026-06-22 16:59:44 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc
2026-06-22 16:59:48 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-06-22 16:59:48 INFO [Posterior Diff ENS15] Memory-optimized diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS15_2020020611.nc
2026-06-22 16:59:48 INFO [Posterior Diff ENS15] Memory-optimized relative diff saved to /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICEMIdmp_0615_15m_low_v2/posteriors/2020020611/diff_posterior_ENS15_2020020611.relative.nc
2026-06-22 16:59:48 INFO Next run starts from 2020-02-06 11:00:00
2026-06-22 16:59:48 INFO Cycle is DONE; starting a new loop!
2026-06-22 16:59:48 INFO [TIME] step_end current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-22 16:59:48 INFO [TIME] step_start current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-22 16:59:48 INFO [TIME] window start=2020-02-06 11:00:00 end=2020-02-06 13:00:00 run_hours=2 has_assimilation=True
2026-06-22 16:59:48 INFO ---------->>> Running CHIMERE model from 2020-02-06 11:00:00 to 2020-02-06 13:00:00
2026-06-22 16:59:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 16:59:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1
2026-06-22 16:59:48 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc
2026-06-22 16:59:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-06-22 16:59:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 16:59:48 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-06-22 16:59:49 INFO Queuing job for member 1...
2026-06-22 16:59:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 16:59:49 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-06-22 16:59:49 INFO Found: ['4949716']
2026-06-22 16:59:54 INFO [TGCC-IRENE] Submitted job with ID:['4949716']
2026-06-22 16:59:54 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 16:59:54 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2
2026-06-22 16:59:54 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc
2026-06-22 16:59:54 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-06-22 16:59:54 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 16:59:54 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-06-22 16:59:54 INFO Queuing job for member 2...
2026-06-22 16:59:54 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 16:59:54 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-06-22 16:59:55 INFO Found: ['4949719']
2026-06-22 17:00:00 INFO [TGCC-IRENE] Submitted job with ID:['4949719']
2026-06-22 17:00:00 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:00:00 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3
2026-06-22 17:00:00 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc
2026-06-22 17:00:00 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-06-22 17:00:00 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:00:00 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-06-22 17:00:00 INFO Queuing job for member 3...
2026-06-22 17:00:00 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:00:00 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-06-22 17:00:01 INFO Found: ['4949723']
2026-06-22 17:00:06 INFO [TGCC-IRENE] Submitted job with ID:['4949723']
2026-06-22 17:00:06 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:00:06 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4
2026-06-22 17:00:06 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc
2026-06-22 17:00:06 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-06-22 17:00:06 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:00:06 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-06-22 17:00:06 INFO Queuing job for member 4...
2026-06-22 17:00:06 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:00:06 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-06-22 17:00:06 INFO Found: ['4949726']
2026-06-22 17:00:11 INFO [TGCC-IRENE] Submitted job with ID:['4949726']
2026-06-22 17:00:11 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:00:11 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5
2026-06-22 17:00:11 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc
2026-06-22 17:00:11 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-06-22 17:00:11 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:00:11 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-06-22 17:00:12 INFO Queuing job for member 5...
2026-06-22 17:00:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:00:12 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-06-22 17:00:13 INFO Found: ['4949730']
2026-06-22 17:00:18 INFO [TGCC-IRENE] Submitted job with ID:['4949730']
2026-06-22 17:00:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:00:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6
2026-06-22 17:00:18 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc
2026-06-22 17:00:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-06-22 17:00:18 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:00:18 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-06-22 17:00:18 INFO Queuing job for member 6...
2026-06-22 17:00:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:00:18 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-06-22 17:00:18 INFO Found: ['4949734']
2026-06-22 17:00:23 INFO [TGCC-IRENE] Submitted job with ID:['4949734']
2026-06-22 17:00:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:00:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7
2026-06-22 17:00:23 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc
2026-06-22 17:00:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-06-22 17:00:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:00:23 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-06-22 17:00:23 INFO Queuing job for member 7...
2026-06-22 17:00:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:00:23 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-06-22 17:00:24 INFO Found: ['4949738']
2026-06-22 17:00:29 INFO [TGCC-IRENE] Submitted job with ID:['4949738']
2026-06-22 17:00:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:00:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8
2026-06-22 17:00:29 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc
2026-06-22 17:00:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-06-22 17:00:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:00:29 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-06-22 17:00:29 INFO Queuing job for member 8...
2026-06-22 17:00:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:00:29 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-06-22 17:00:30 INFO Found: ['4949740']
2026-06-22 17:00:35 INFO [TGCC-IRENE] Submitted job with ID:['4949740']
2026-06-22 17:00:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:00:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9
2026-06-22 17:00:35 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc
2026-06-22 17:00:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-06-22 17:00:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:00:35 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-06-22 17:00:35 INFO Queuing job for member 9...
2026-06-22 17:00:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:00:35 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-06-22 17:00:36 INFO Found: ['4949742']
2026-06-22 17:00:41 INFO [TGCC-IRENE] Submitted job with ID:['4949742']
2026-06-22 17:00:41 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:00:41 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10
2026-06-22 17:00:41 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc
2026-06-22 17:00:41 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-06-22 17:00:41 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:00:41 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-06-22 17:00:41 INFO Queuing job for member 10...
2026-06-22 17:00:41 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:00:41 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-06-22 17:00:41 INFO Found: ['4949745']
2026-06-22 17:00:46 INFO [TGCC-IRENE] Submitted job with ID:['4949745']
2026-06-22 17:00:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:00:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11
2026-06-22 17:00:46 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc
2026-06-22 17:00:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-06-22 17:00:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:00:46 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-06-22 17:00:46 INFO Queuing job for member 11...
2026-06-22 17:00:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:00:46 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-06-22 17:00:47 INFO Found: ['4949747']
2026-06-22 17:00:52 INFO [TGCC-IRENE] Submitted job with ID:['4949747']
2026-06-22 17:00:52 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:00:52 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12
2026-06-22 17:00:52 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc
2026-06-22 17:00:52 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-06-22 17:00:52 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:00:52 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-06-22 17:00:52 INFO Queuing job for member 12...
2026-06-22 17:00:52 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:00:52 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-06-22 17:00:53 INFO Found: ['4949751']
2026-06-22 17:00:58 INFO [TGCC-IRENE] Submitted job with ID:['4949751']
2026-06-22 17:00:58 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:00:58 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13
2026-06-22 17:00:58 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc
2026-06-22 17:00:58 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-06-22 17:00:58 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:00:58 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-06-22 17:00:58 INFO Queuing job for member 13...
2026-06-22 17:00:58 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:00:58 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-06-22 17:00:59 INFO Found: ['4949753']
2026-06-22 17:01:04 INFO [TGCC-IRENE] Submitted job with ID:['4949753']
2026-06-22 17:01:04 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:01:04 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14
2026-06-22 17:01:04 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc
2026-06-22 17:01:04 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-06-22 17:01:04 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:01:04 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-06-22 17:01:04 INFO Queuing job for member 14...
2026-06-22 17:01:04 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:01:04 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-06-22 17:01:05 INFO Found: ['4949756']
2026-06-22 17:01:10 INFO [TGCC-IRENE] Submitted job with ID:['4949756']
2026-06-22 17:01:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-22 17:01:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15
2026-06-22 17:01:10 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc
2026-06-22 17:01:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-06-22 17:01:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-22 17:01:10 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-06-22 17:01:10 INFO Queuing job for member 15...
2026-06-22 17:01:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-22 17:01:10 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-06-22 17:01:11 INFO Found: ['4949759']
2026-06-22 17:01:16 INFO [TGCC-IRENE] Submitted job with ID:['4949759']
2026-06-22 17:01:16 INFO Checking job status ...
2026-06-22 17:01:16 INFO None 4949716: status RUNNING/PENDING
2026-06-22 17:01:16 INFO None 4949719: status RUNNING/PENDING
2026-06-22 17:01:16 INFO None 4949723: status RUNNING/PENDING
2026-06-22 17:01:16 INFO None 4949726: status RUNNING/PENDING
2026-06-22 17:01:16 INFO None 4949730: status RUNNING/PENDING
2026-06-22 17:01:16 INFO None 4949734: status RUNNING/PENDING
2026-06-22 17:01:16 INFO None 4949738: status RUNNING/PENDING
2026-06-22 17:01:16 INFO None 4949740: status RUNNING/PENDING
2026-06-22 17:01:16 INFO None 4949742: status RUNNING/PENDING
2026-06-22 17:01:16 INFO None 4949745: status RUNNING/PENDING
2026-06-22 17:01:16 INFO None 4949747: status RUNNING/PENDING
2026-06-22 17:01:16 INFO None 4949751: status RUNNING/PENDING
2026-06-22 17:01:16 INFO None 4949753: status RUNNING/PENDING
2026-06-22 17:01:16 INFO None 4949756: status RUNNING/PENDING
2026-06-22 17:01:16 INFO None 4949759: status RUNNING/PENDING
2026-06-22 17:01:16 INFO Jobs still running: ['4949716', '4949719', '4949723', '4949726', '4949730', '4949734', '4949738', '4949740', '4949742', '4949745', '4949747', '4949751', '4949753', '4949756', '4949759']. Waiting...
2026-06-22 17:01:31 INFO None 4949716: status RUNNING/PENDING
2026-06-22 17:01:31 INFO None 4949719: status RUNNING/PENDING
2026-06-22 17:01:31 INFO None 4949723: status RUNNING/PENDING
2026-06-22 17:01:31 INFO None 4949726: status RUNNING/PENDING
2026-06-22 17:01:31 INFO None 4949730: status RUNNING/PENDING
2026-06-22 17:01:31 INFO None 4949734: status RUNNING/PENDING
2026-06-22 17:01:31 INFO None 4949738: status RUNNING/PENDING
2026-06-22 17:01:31 INFO None 4949740: status RUNNING/PENDING
2026-06-22 17:01:31 INFO None 4949742: status RUNNING/PENDING
2026-06-22 17:01:31 INFO None 4949745: status RUNNING/PENDING
2026-06-22 17:01:31 INFO None 4949747: status RUNNING/PENDING
2026-06-22 17:01:31 INFO None 4949751: status RUNNING/PENDING
2026-06-22 17:01:31 INFO None 4949753: status RUNNING/PENDING
2026-06-22 17:01:31 INFO None 4949756: status RUNNING/PENDING
2026-06-22 17:01:31 INFO None 4949759: status RUNNING/PENDING
2026-06-22 17:01:31 INFO Jobs still running: ['4949716', '4949719', '4949723', '4949726', '4949730', '4949734', '4949738', '4949740', '4949742', '4949745', '4949747', '4949751', '4949753', '4949756', '4949759']. Waiting...
2026-06-22 17:01:46 INFO None 4949716: status RUNNING/PENDING
2026-06-22 17:01:46 INFO None 4949719: status RUNNING/PENDING
2026-06-22 17:01:46 INFO None 4949723: status RUNNING/PENDING
2026-06-22 17:01:46 INFO None 4949726: status RUNNING/PENDING
2026-06-22 17:01:46 INFO None 4949730: status RUNNING/PENDING
2026-06-22 17:01:47 INFO None 4949734: status RUNNING/PENDING
2026-06-22 17:01:47 INFO None 4949738: status RUNNING/PENDING
2026-06-22 17:01:47 INFO None 4949740: status RUNNING/PENDING
2026-06-22 17:01:47 INFO None 4949742: status RUNNING/PENDING
2026-06-22 17:01:47 INFO None 4949745: status RUNNING/PENDING
2026-06-22 17:01:47 INFO None 4949747: status RUNNING/PENDING
2026-06-22 17:01:47 INFO None 4949751: status RUNNING/PENDING
2026-06-22 17:01:47 INFO None 4949753: status RUNNING/PENDING
2026-06-22 17:01:47 INFO None 4949756: status RUNNING/PENDING
2026-06-22 17:01:47 INFO None 4949759: status RUNNING/PENDING
2026-06-22 17:01:47 INFO Jobs still running: ['4949716', '4949719', '4949723', '4949726', '4949730', '4949734', '4949738', '4949740', '4949742', '4949745', '4949747', '4949751', '4949753', '4949756', '4949759']. Waiting...
2026-06-22 17:02:02 INFO None 4949716: status RUNNING/PENDING
2026-06-22 17:02:02 INFO None 4949719: status RUNNING/PENDING
2026-06-22 17:02:02 INFO None 4949723: status RUNNING/PENDING
2026-06-22 17:02:02 INFO None 4949726: status RUNNING/PENDING
2026-06-22 17:02:02 INFO None 4949730: status RUNNING/PENDING
2026-06-22 17:02:02 INFO None 4949734: status RUNNING/PENDING
2026-06-22 17:02:02 INFO None 4949738: status RUNNING/PENDING
2026-06-22 17:02:02 INFO None 4949740: status RUNNING/PENDING
2026-06-22 17:02:02 INFO None 4949742: status RUNNING/PENDING
2026-06-22 17:02:02 INFO None 4949745: status RUNNING/PENDING
2026-06-22 17:02:02 INFO None 4949747: status RUNNING/PENDING
2026-06-22 17:02:02 INFO None 4949751: status RUNNING/PENDING
2026-06-22 17:02:02 INFO None 4949753: status RUNNING/PENDING
2026-06-22 17:02:02 INFO None 4949756: status RUNNING/PENDING
2026-06-22 17:02:02 INFO None 4949759: status RUNNING/PENDING
2026-06-22 17:02:02 INFO Jobs still running: ['4949716', '4949719', '4949723', '4949726', '4949730', '4949734', '4949738', '4949740', '4949742', '4949745', '4949747', '4949751', '4949753', '4949756', '4949759']. Waiting...
[2026-06-22T17:02:09.050] error: *** JOB 4949478 ON irene5937 CANCELLED AT 2026-06-22T17:02:09 DUE to SIGNAL Terminated ***
