+ SCRIPT_PID=3937314
+ /bin/bash -x /tmp/tmp.dhheFrvX3p
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
2026-07-15 10:54:42 INFO 
███    ███ ██ ███    ███ ███████ ███████ ██
████  ████ ██ ████  ████ ██      ██      ██
██ ████ ██ ██ ██ ████ ██ █████   ███████ ██
██  ██  ██ ██ ██  ██  ██ ██           ██ ██
██      ██ ██ ██      ██ ███████ ███████ ██



2026-07-15 10:54:42 INFO [PIPELINE] =======================================
2026-07-15 10:54:42 INFO [PIPELINE] Starting chimere–DART orchestrator
2026-07-15 10:54:42 INFO [PIPELINE] Config file: config/config_irene_IM.yaml
2026-07-15 10:54:42 INFO [PIPELINE] Run dir: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart
2026-07-15 10:54:42 INFO [PIPELINE] Log file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/mimesi_orchestrator_logs/chimere_DART_20260715_105442.log
2026-07-15 10:54:42 INFO [PIPELINE] =======================================
2026-07-15 10:54:42 INFO Running assimilation with model_type=ModelType.CHIMERE
2026-07-15 10:54:42 INFO Using scheduler=Scheduler.SLURM, queue=rome
2026-07-15 10:54:42 INFO [STEP] ---- TIME LOOP START ----
2026-07-15 10:54:42 INFO [TIME] step_start current_time=2020-02-06 09:00:00 simulated_time=None dt=0 days 01:00:00 end_time=2020-02-06 23:00:00
2026-07-15 10:54:42 INFO [TIME] window start=2020-02-06 09:00:00 end=2020-02-06 11:00:00 run_hours=2 has_assimilation=True
2026-07-15 10:54:42 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-06 23:00:00
2026-07-15 10:54:42 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 11:00:00 days=153072 seconds=39600
2026-07-15 10:54:42 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-15 10:54:42 INFO [TIME] increment current_time 2020-02-06 09:00:00 -> 2020-02-06 11:00:00
2026-07-15 10:54:42 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-06 23:00:00
2026-07-15 10:54:42 INFO ---------->>> Running process_satellite_data()
2026-07-15 10:54:42 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12003.nc
2026-07-15 10:54:42 INFO ---------->>> Running run_obs_converter()
2026-07-15 10:54:42 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_39045_153072.out
2026-07-15 10:54:42 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_39045_153072.out
2026-07-15 10:54:42 INFO ---------->>> Running DART
2026-07-15 10:54:42 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-15 10:54:42 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/chim_ENS1_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/chim_ENS1_2020020611_1_out_toDART.nc
2026-07-15 10:54:50 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/chim_ENS2_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/chim_ENS2_2020020611_1_out_toDART.nc
2026-07-15 10:54:51 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/chim_ENS3_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/chim_ENS3_2020020611_1_out_toDART.nc
2026-07-15 10:54:51 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/chim_ENS4_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/chim_ENS4_2020020611_1_out_toDART.nc
2026-07-15 10:54:51 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/chim_ENS5_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/chim_ENS5_2020020611_1_out_toDART.nc
2026-07-15 10:54:52 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/chim_ENS6_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/chim_ENS6_2020020611_1_out_toDART.nc
2026-07-15 10:54:52 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/chim_ENS7_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/chim_ENS7_2020020611_1_out_toDART.nc
2026-07-15 10:54:53 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/chim_ENS8_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/chim_ENS8_2020020611_1_out_toDART.nc
2026-07-15 10:54:53 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/chim_ENS9_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/chim_ENS9_2020020611_1_out_toDART.nc
2026-07-15 10:54:53 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/chim_ENS10_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/chim_ENS10_2020020611_1_out_toDART.nc
2026-07-15 10:54:54 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/chim_ENS11_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/chim_ENS11_2020020611_1_out_toDART.nc
2026-07-15 10:54:54 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/chim_ENS12_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/chim_ENS12_2020020611_1_out_toDART.nc
2026-07-15 10:54:54 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/chim_ENS13_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/chim_ENS13_2020020611_1_out_toDART.nc
2026-07-15 10:54:55 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/chim_ENS14_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/chim_ENS14_2020020611_1_out_toDART.nc
2026-07-15 10:54:55 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/chim_ENS15_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/chim_ENS15_2020020611_1_out_toDART.nc
2026-07-15 10:54:55 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-15 10:54:56 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-15 10:54:56 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-15 10:54:56 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-15 10:54:56 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-15 10:54:56 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-15 10:55:10 INFO Found: []
2026-07-15 10:55:10 INFO No job id returned by command ./run_filter.bsh
2026-07-15 10:55:10 INFO No monitoring will be performed
2026-07-15 10:55:10 INFO Moving DART output files to analysis and preassim directories for date 2020020611 if present ...
2026-07-15 10:55:10 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-15 10:55:10 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-15 10:55:10 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-15 10:55:10 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-15 10:55:10 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-15 10:55:10 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-15 10:55:10 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-15 10:55:10 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-15 10:55:10 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-15 10:55:10 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-15 10:55:10 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-15 10:55:10 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-15 10:55:10 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-15 10:55:10 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-15 10:55:10 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-15 10:55:11 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-15 10:55:11 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-15 10:55:11 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-15 10:55:11 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-15 10:55:11 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-15 10:55:11 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-15 10:55:11 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-15 10:55:11 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-15 10:55:11 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-15 10:55:11 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-15 10:55:11 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-15 10:55:11 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-15 10:55:11 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-15 10:55:11 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-15 10:55:11 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-15 10:55:11 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-15 10:55:11 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020611'
2026-07-15 10:55:11 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-15 10:55:11 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020611'
2026-07-15 10:55:11 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-15 10:55:14 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-15 10:55:14 INFO run_dart() is DONE.
2026-07-15 10:55:14 INFO ---------->>> Running update_pollutant_in_end()
2026-07-15 10:55:14 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:55:14 INFO No previous orbit memory found.
2026-07-15 10:55:15 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:55:15 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS1.nc
2026-07-15 10:55:15 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:55:15 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:55:15 INFO No previous orbit memory found.
2026-07-15 10:55:15 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:55:15 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS1.nc
2026-07-15 10:55:16 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:55:16 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:55:16 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:55:16 INFO No previous orbit memory found.
2026-07-15 10:55:16 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:55:16 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS2.nc
2026-07-15 10:55:17 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:55:17 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:55:17 INFO No previous orbit memory found.
2026-07-15 10:55:17 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:55:17 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS2.nc
2026-07-15 10:55:18 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:55:18 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:55:18 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:55:18 INFO No previous orbit memory found.
2026-07-15 10:55:18 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:55:18 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS3.nc
2026-07-15 10:55:19 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:55:19 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:55:19 INFO No previous orbit memory found.
2026-07-15 10:55:19 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:55:19 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS3.nc
2026-07-15 10:55:19 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:55:19 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:55:20 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:55:20 INFO No previous orbit memory found.
2026-07-15 10:55:20 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:55:20 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS4.nc
2026-07-15 10:55:20 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:55:21 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:55:21 INFO No previous orbit memory found.
2026-07-15 10:55:21 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:55:21 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS4.nc
2026-07-15 10:55:21 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:55:21 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:55:22 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:55:22 INFO No previous orbit memory found.
2026-07-15 10:55:22 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:55:22 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS5.nc
2026-07-15 10:55:22 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:55:22 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:55:22 INFO No previous orbit memory found.
2026-07-15 10:55:22 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:55:22 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS5.nc
2026-07-15 10:55:23 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:55:23 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:55:23 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:55:23 INFO No previous orbit memory found.
2026-07-15 10:55:23 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:55:23 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS6.nc
2026-07-15 10:55:24 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:55:24 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:55:24 INFO No previous orbit memory found.
2026-07-15 10:55:24 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:55:24 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS6.nc
2026-07-15 10:55:24 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:55:24 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:55:25 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:55:25 INFO No previous orbit memory found.
2026-07-15 10:55:25 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:55:25 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS7.nc
2026-07-15 10:55:25 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:55:26 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:55:26 INFO No previous orbit memory found.
2026-07-15 10:55:26 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:55:26 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS7.nc
2026-07-15 10:55:26 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:55:26 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:55:27 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:55:27 INFO No previous orbit memory found.
2026-07-15 10:55:27 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:55:27 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS8.nc
2026-07-15 10:55:27 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:55:27 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:55:27 INFO No previous orbit memory found.
2026-07-15 10:55:27 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:55:27 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS8.nc
2026-07-15 10:55:28 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:55:28 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:55:28 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:55:28 INFO No previous orbit memory found.
2026-07-15 10:55:28 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:55:28 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS9.nc
2026-07-15 10:55:29 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:55:29 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:55:29 INFO No previous orbit memory found.
2026-07-15 10:55:29 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:55:29 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS9.nc
2026-07-15 10:55:30 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:55:30 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:55:30 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:55:30 INFO No previous orbit memory found.
2026-07-15 10:55:30 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:55:30 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS10.nc
2026-07-15 10:55:30 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:55:31 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:55:31 INFO No previous orbit memory found.
2026-07-15 10:55:31 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:55:31 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS10.nc
2026-07-15 10:55:31 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:55:31 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:55:32 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:55:32 INFO No previous orbit memory found.
2026-07-15 10:55:32 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:55:32 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS11.nc
2026-07-15 10:55:32 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:55:32 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:55:32 INFO No previous orbit memory found.
2026-07-15 10:55:32 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:55:32 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS11.nc
2026-07-15 10:55:33 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:55:33 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:55:33 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:55:33 INFO No previous orbit memory found.
2026-07-15 10:55:33 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:55:33 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS12.nc
2026-07-15 10:55:34 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:55:34 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:55:34 INFO No previous orbit memory found.
2026-07-15 10:55:34 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:55:34 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS12.nc
2026-07-15 10:55:35 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:55:35 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:55:35 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:55:35 INFO No previous orbit memory found.
2026-07-15 10:55:35 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:55:35 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS13.nc
2026-07-15 10:55:36 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:55:36 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:55:36 INFO No previous orbit memory found.
2026-07-15 10:55:36 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:55:36 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS13.nc
2026-07-15 10:55:36 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:55:36 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:55:37 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:55:37 INFO No previous orbit memory found.
2026-07-15 10:55:37 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:55:37 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS14.nc
2026-07-15 10:55:37 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:55:38 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:55:38 INFO No previous orbit memory found.
2026-07-15 10:55:38 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:55:38 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS14.nc
2026-07-15 10:55:38 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:55:38 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:55:38 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:55:38 INFO No previous orbit memory found.
2026-07-15 10:55:38 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:55:38 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS15.nc
2026-07-15 10:55:39 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:55:39 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:55:39 INFO No previous orbit memory found.
2026-07-15 10:55:39 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:55:39 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS15.nc
2026-07-15 10:55:40 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:55:40 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:55:40 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-15 10:55:40 INFO [TIME] step_end current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-06 23:00:00
2026-07-15 10:55:40 INFO [TIME] step_start current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-06 23:00:00
2026-07-15 10:55:40 INFO [TIME] window start=2020-02-06 11:00:00 end=2020-02-06 13:00:00 run_hours=2 has_assimilation=True
2026-07-15 10:55:40 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-06 23:00:00
2026-07-15 10:55:40 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 13:00:00 days=153072 seconds=46800
2026-07-15 10:55:40 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-15 10:55:40 INFO [TIME] increment current_time 2020-02-06 11:00:00 -> 2020-02-06 13:00:00
2026-07-15 10:55:40 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-06 23:00:00
2026-07-15 10:55:40 INFO ---------->>> Running process_satellite_data()
2026-07-15 10:55:40 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12004.nc
2026-07-15 10:55:40 INFO ---------->>> Running run_obs_converter()
2026-07-15 10:55:40 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_45135_153072.out
2026-07-15 10:55:40 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_45135_153072.out
2026-07-15 10:55:40 INFO ---------->>> Running DART
2026-07-15 10:55:40 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-15 10:55:40 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/chim_ENS1_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/chim_ENS1_2020020613_1_out_toDART.nc
2026-07-15 10:55:40 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/chim_ENS2_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/chim_ENS2_2020020613_1_out_toDART.nc
2026-07-15 10:55:40 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/chim_ENS3_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/chim_ENS3_2020020613_1_out_toDART.nc
2026-07-15 10:55:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/chim_ENS4_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/chim_ENS4_2020020613_1_out_toDART.nc
2026-07-15 10:55:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/chim_ENS5_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/chim_ENS5_2020020613_1_out_toDART.nc
2026-07-15 10:55:42 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/chim_ENS6_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/chim_ENS6_2020020613_1_out_toDART.nc
2026-07-15 10:55:42 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/chim_ENS7_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/chim_ENS7_2020020613_1_out_toDART.nc
2026-07-15 10:55:42 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/chim_ENS8_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/chim_ENS8_2020020613_1_out_toDART.nc
2026-07-15 10:55:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/chim_ENS9_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/chim_ENS9_2020020613_1_out_toDART.nc
2026-07-15 10:55:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/chim_ENS10_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/chim_ENS10_2020020613_1_out_toDART.nc
2026-07-15 10:55:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/chim_ENS11_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/chim_ENS11_2020020613_1_out_toDART.nc
2026-07-15 10:55:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/chim_ENS12_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/chim_ENS12_2020020613_1_out_toDART.nc
2026-07-15 10:55:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/chim_ENS13_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/chim_ENS13_2020020613_1_out_toDART.nc
2026-07-15 10:55:45 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/chim_ENS14_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/chim_ENS14_2020020613_1_out_toDART.nc
2026-07-15 10:55:45 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/chim_ENS15_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/chim_ENS15_2020020613_1_out_toDART.nc
2026-07-15 10:55:45 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-15 10:55:45 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-15 10:55:45 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-15 10:55:45 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-15 10:55:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-15 10:55:45 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-15 10:56:01 INFO Found: []
2026-07-15 10:56:01 INFO No job id returned by command ./run_filter.bsh
2026-07-15 10:56:01 INFO No monitoring will be performed
2026-07-15 10:56:01 INFO Moving DART output files to analysis and preassim directories for date 2020020613 if present ...
2026-07-15 10:56:01 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-15 10:56:01 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-15 10:56:01 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-15 10:56:01 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-15 10:56:01 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-15 10:56:01 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-15 10:56:01 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-15 10:56:01 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-15 10:56:01 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-15 10:56:01 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-15 10:56:01 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-15 10:56:01 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-15 10:56:01 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-15 10:56:01 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-15 10:56:01 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-15 10:56:01 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-15 10:56:01 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-15 10:56:01 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-15 10:56:01 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-15 10:56:01 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-15 10:56:01 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-15 10:56:01 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-15 10:56:01 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-15 10:56:01 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-15 10:56:01 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-15 10:56:01 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-15 10:56:01 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-15 10:56:01 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-15 10:56:01 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-15 10:56:01 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-15 10:56:01 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-15 10:56:01 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020613'
2026-07-15 10:56:01 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-15 10:56:01 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020613'
2026-07-15 10:56:01 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-15 10:56:01 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-15 10:56:01 INFO run_dart() is DONE.
2026-07-15 10:56:01 INFO ---------->>> Running update_pollutant_in_end()
2026-07-15 10:56:02 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:02 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS1.nc
2026-07-15 10:56:02 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:02 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:02 INFO ['2020-02-06T11:00:00.000000000']
2026-07-15 10:56:02 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:02 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS1.nc
2026-07-15 10:56:02 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:03 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:03 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS1.nc
2026-07-15 10:56:03 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:03 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:03 INFO ['2020-02-06T11:00:00.000000000']
2026-07-15 10:56:03 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:03 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS1.nc
2026-07-15 10:56:03 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:03 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:56:03 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:03 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS2.nc
2026-07-15 10:56:03 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:03 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:03 INFO ['2020-02-06T11:00:00.000000000']
2026-07-15 10:56:03 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:03 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS2.nc
2026-07-15 10:56:04 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:04 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:04 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS2.nc
2026-07-15 10:56:04 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:04 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:04 INFO ['2020-02-06T11:00:00.000000000']
2026-07-15 10:56:04 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:04 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS2.nc
2026-07-15 10:56:05 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:05 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:56:05 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:05 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS3.nc
2026-07-15 10:56:05 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:05 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:05 INFO ['2020-02-06T11:00:00.000000000']
2026-07-15 10:56:05 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:05 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS3.nc
2026-07-15 10:56:06 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:06 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:06 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS3.nc
2026-07-15 10:56:06 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:06 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:06 INFO ['2020-02-06T11:00:00.000000000']
2026-07-15 10:56:06 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:06 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS3.nc
2026-07-15 10:56:06 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:06 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:56:07 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:07 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS4.nc
2026-07-15 10:56:07 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:07 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:07 INFO ['2020-02-06T11:00:00.000000000']
2026-07-15 10:56:07 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:07 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS4.nc
2026-07-15 10:56:07 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:08 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:08 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS4.nc
2026-07-15 10:56:08 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:08 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:08 INFO ['2020-02-06T11:00:00.000000000']
2026-07-15 10:56:08 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:08 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS4.nc
2026-07-15 10:56:08 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:08 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:56:09 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:09 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS5.nc
2026-07-15 10:56:09 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:09 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:09 INFO ['2020-02-06T11:00:00.000000000']
2026-07-15 10:56:09 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:09 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS5.nc
2026-07-15 10:56:09 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:10 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:10 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS5.nc
2026-07-15 10:56:10 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:10 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:10 INFO ['2020-02-06T11:00:00.000000000']
2026-07-15 10:56:10 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:10 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS5.nc
2026-07-15 10:56:10 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:10 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:56:11 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:11 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS6.nc
2026-07-15 10:56:11 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:11 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:11 INFO ['2020-02-06T11:00:00.000000000']
2026-07-15 10:56:11 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:11 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS6.nc
2026-07-15 10:56:11 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:12 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:12 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS6.nc
2026-07-15 10:56:12 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:12 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:12 INFO ['2020-02-06T11:00:00.000000000']
2026-07-15 10:56:12 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:12 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS6.nc
2026-07-15 10:56:12 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:12 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:56:12 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:12 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS7.nc
2026-07-15 10:56:12 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:12 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:12 INFO ['2020-02-06T11:00:00.000000000']
2026-07-15 10:56:12 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:12 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS7.nc
2026-07-15 10:56:13 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:13 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:13 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS7.nc
2026-07-15 10:56:13 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:13 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:13 INFO ['2020-02-06T11:00:00.000000000']
2026-07-15 10:56:13 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:13 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS7.nc
2026-07-15 10:56:14 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:14 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:56:14 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:14 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS8.nc
2026-07-15 10:56:14 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:14 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:14 INFO ['2020-02-06T11:00:00.000000000']
2026-07-15 10:56:14 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:14 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS8.nc
2026-07-15 10:56:15 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:15 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:15 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS8.nc
2026-07-15 10:56:15 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:15 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:15 INFO ['2020-02-06T11:00:00.000000000']
2026-07-15 10:56:15 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:15 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS8.nc
2026-07-15 10:56:16 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:16 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:56:16 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:16 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS9.nc
2026-07-15 10:56:16 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:16 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:16 INFO ['2020-02-06T11:00:00.000000000']
2026-07-15 10:56:16 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:16 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS9.nc
2026-07-15 10:56:17 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:17 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:17 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS9.nc
2026-07-15 10:56:17 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:17 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:17 INFO ['2020-02-06T11:00:00.000000000']
2026-07-15 10:56:17 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:17 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS9.nc
2026-07-15 10:56:18 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:18 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:56:18 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:18 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS10.nc
2026-07-15 10:56:18 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:18 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:18 INFO ['2020-02-06T11:00:00.000000000']
2026-07-15 10:56:18 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:18 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS10.nc
2026-07-15 10:56:19 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:19 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:19 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS10.nc
2026-07-15 10:56:19 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:19 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:19 INFO ['2020-02-06T11:00:00.000000000']
2026-07-15 10:56:19 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:19 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS10.nc
2026-07-15 10:56:19 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:19 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:56:20 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:20 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS11.nc
2026-07-15 10:56:20 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:20 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:20 INFO ['2020-02-06T11:00:00.000000000']
2026-07-15 10:56:20 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:20 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS11.nc
2026-07-15 10:56:20 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:21 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:21 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS11.nc
2026-07-15 10:56:21 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:21 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:21 INFO ['2020-02-06T11:00:00.000000000']
2026-07-15 10:56:21 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:21 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS11.nc
2026-07-15 10:56:21 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:21 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:56:21 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:21 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS12.nc
2026-07-15 10:56:21 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:21 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:21 INFO ['2020-02-06T11:00:00.000000000']
2026-07-15 10:56:21 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:21 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS12.nc
2026-07-15 10:56:22 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:22 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:22 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS12.nc
2026-07-15 10:56:22 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:22 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:22 INFO ['2020-02-06T11:00:00.000000000']
2026-07-15 10:56:22 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:22 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS12.nc
2026-07-15 10:56:23 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:23 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:56:23 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:23 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS13.nc
2026-07-15 10:56:23 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:23 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:23 INFO ['2020-02-06T11:00:00.000000000']
2026-07-15 10:56:23 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:23 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS13.nc
2026-07-15 10:56:24 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:24 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:24 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS13.nc
2026-07-15 10:56:24 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:24 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:24 INFO ['2020-02-06T11:00:00.000000000']
2026-07-15 10:56:24 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:24 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS13.nc
2026-07-15 10:56:25 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:25 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:56:25 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:25 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS14.nc
2026-07-15 10:56:25 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:25 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:25 INFO ['2020-02-06T11:00:00.000000000']
2026-07-15 10:56:25 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:25 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS14.nc
2026-07-15 10:56:26 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:26 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:26 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS14.nc
2026-07-15 10:56:26 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:26 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:26 INFO ['2020-02-06T11:00:00.000000000']
2026-07-15 10:56:26 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:26 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS14.nc
2026-07-15 10:56:26 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:26 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:56:27 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:27 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS15.nc
2026-07-15 10:56:27 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:27 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:27 INFO ['2020-02-06T11:00:00.000000000']
2026-07-15 10:56:27 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:27 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS15.nc
2026-07-15 10:56:27 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:28 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:28 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS15.nc
2026-07-15 10:56:28 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:28 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:28 INFO ['2020-02-06T11:00:00.000000000']
2026-07-15 10:56:28 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:28 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS15.nc
2026-07-15 10:56:28 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:28 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:56:28 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-15 10:56:28 INFO [TIME] step_end current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-06 23:00:00
2026-07-15 10:56:28 INFO [TIME] step_start current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-06 23:00:00
2026-07-15 10:56:28 INFO [TIME] window start=2020-02-06 13:00:00 end=2020-02-06 14:00:00 run_hours=1 has_assimilation=True
2026-07-15 10:56:28 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-06 23:00:00
2026-07-15 10:56:28 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 14:00:00 days=153072 seconds=50400
2026-07-15 10:56:28 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-15 10:56:28 INFO [TIME] increment current_time 2020-02-06 13:00:00 -> 2020-02-06 14:00:00
2026-07-15 10:56:28 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-06 23:00:00
2026-07-15 10:56:28 INFO ---------->>> Running process_satellite_data()
2026-07-15 10:56:28 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12005.nc
2026-07-15 10:56:28 INFO ---------->>> Running run_obs_converter()
2026-07-15 10:56:28 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_51225_153072.out
2026-07-15 10:56:28 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_51225_153072.out
2026-07-15 10:56:28 INFO ---------->>> Running DART
2026-07-15 10:56:28 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-15 10:56:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/chim_ENS1_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/chim_ENS1_2020020614_1_out_toDART.nc
2026-07-15 10:56:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/chim_ENS2_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/chim_ENS2_2020020614_1_out_toDART.nc
2026-07-15 10:56:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/chim_ENS3_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/chim_ENS3_2020020614_1_out_toDART.nc
2026-07-15 10:56:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/chim_ENS4_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/chim_ENS4_2020020614_1_out_toDART.nc
2026-07-15 10:56:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/chim_ENS5_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/chim_ENS5_2020020614_1_out_toDART.nc
2026-07-15 10:56:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/chim_ENS6_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/chim_ENS6_2020020614_1_out_toDART.nc
2026-07-15 10:56:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/chim_ENS7_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/chim_ENS7_2020020614_1_out_toDART.nc
2026-07-15 10:56:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/chim_ENS8_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/chim_ENS8_2020020614_1_out_toDART.nc
2026-07-15 10:56:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/chim_ENS9_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/chim_ENS9_2020020614_1_out_toDART.nc
2026-07-15 10:56:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/chim_ENS10_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/chim_ENS10_2020020614_1_out_toDART.nc
2026-07-15 10:56:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/chim_ENS11_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/chim_ENS11_2020020614_1_out_toDART.nc
2026-07-15 10:56:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/chim_ENS12_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/chim_ENS12_2020020614_1_out_toDART.nc
2026-07-15 10:56:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/chim_ENS13_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/chim_ENS13_2020020614_1_out_toDART.nc
2026-07-15 10:56:33 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/chim_ENS14_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/chim_ENS14_2020020614_1_out_toDART.nc
2026-07-15 10:56:33 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/chim_ENS15_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/chim_ENS15_2020020614_1_out_toDART.nc
2026-07-15 10:56:33 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-15 10:56:33 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-15 10:56:33 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-15 10:56:33 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-15 10:56:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-15 10:56:33 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-15 10:56:43 INFO Found: []
2026-07-15 10:56:43 INFO No job id returned by command ./run_filter.bsh
2026-07-15 10:56:43 INFO No monitoring will be performed
2026-07-15 10:56:43 INFO Moving DART output files to analysis and preassim directories for date 2020020614 if present ...
2026-07-15 10:56:43 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-15 10:56:43 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-15 10:56:43 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-15 10:56:43 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-15 10:56:43 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-15 10:56:43 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-15 10:56:43 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-15 10:56:43 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-15 10:56:43 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-15 10:56:43 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-15 10:56:43 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-15 10:56:43 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-15 10:56:43 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-15 10:56:43 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-15 10:56:43 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-15 10:56:43 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-15 10:56:43 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-15 10:56:43 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-15 10:56:43 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-15 10:56:43 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-15 10:56:43 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-15 10:56:43 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-15 10:56:43 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-15 10:56:43 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-15 10:56:43 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-15 10:56:43 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-15 10:56:44 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-15 10:56:44 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-15 10:56:44 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-15 10:56:44 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-15 10:56:44 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-15 10:56:44 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020614'
2026-07-15 10:56:44 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-15 10:56:44 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020614'
2026-07-15 10:56:44 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-15 10:56:44 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-15 10:56:44 INFO run_dart() is DONE.
2026-07-15 10:56:44 INFO ---------->>> Running update_pollutant_in_end()
2026-07-15 10:56:44 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:44 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS1.nc
2026-07-15 10:56:44 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:44 INFO 2020-02-06T13:00:00.000000000
2026-07-15 10:56:44 INFO ['2020-02-06T11:00:00.000000000' '2020-02-06T13:00:00.000000000']
2026-07-15 10:56:44 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:44 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS1.nc
2026-07-15 10:56:45 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:45 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:45 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS1.nc
2026-07-15 10:56:45 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:45 INFO 2020-02-06T13:00:00.000000000
2026-07-15 10:56:45 INFO ['2020-02-06T11:00:00.000000000' '2020-02-06T13:00:00.000000000']
2026-07-15 10:56:45 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:45 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS1.nc
2026-07-15 10:56:45 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:45 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:56:46 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:46 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS2.nc
2026-07-15 10:56:46 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:46 INFO 2020-02-06T13:00:00.000000000
2026-07-15 10:56:46 INFO ['2020-02-06T11:00:00.000000000' '2020-02-06T13:00:00.000000000']
2026-07-15 10:56:46 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:46 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS2.nc
2026-07-15 10:56:46 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:47 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:47 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS2.nc
2026-07-15 10:56:47 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:47 INFO 2020-02-06T13:00:00.000000000
2026-07-15 10:56:47 INFO ['2020-02-06T11:00:00.000000000' '2020-02-06T13:00:00.000000000']
2026-07-15 10:56:47 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:47 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS2.nc
2026-07-15 10:56:47 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:47 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:56:48 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:48 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS3.nc
2026-07-15 10:56:48 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:48 INFO 2020-02-06T13:00:00.000000000
2026-07-15 10:56:48 INFO ['2020-02-06T11:00:00.000000000' '2020-02-06T13:00:00.000000000']
2026-07-15 10:56:48 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:48 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS3.nc
2026-07-15 10:56:48 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:49 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:49 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS3.nc
2026-07-15 10:56:49 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:49 INFO 2020-02-06T13:00:00.000000000
2026-07-15 10:56:49 INFO ['2020-02-06T11:00:00.000000000' '2020-02-06T13:00:00.000000000']
2026-07-15 10:56:49 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:49 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS3.nc
2026-07-15 10:56:49 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:49 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:56:50 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:50 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS4.nc
2026-07-15 10:56:50 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:50 INFO 2020-02-06T13:00:00.000000000
2026-07-15 10:56:50 INFO ['2020-02-06T11:00:00.000000000' '2020-02-06T13:00:00.000000000']
2026-07-15 10:56:50 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:50 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS4.nc
2026-07-15 10:56:50 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:51 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:51 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS4.nc
2026-07-15 10:56:51 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:51 INFO 2020-02-06T13:00:00.000000000
2026-07-15 10:56:51 INFO ['2020-02-06T11:00:00.000000000' '2020-02-06T13:00:00.000000000']
2026-07-15 10:56:51 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:51 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS4.nc
2026-07-15 10:56:51 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:51 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:56:51 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:51 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS5.nc
2026-07-15 10:56:51 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:51 INFO 2020-02-06T13:00:00.000000000
2026-07-15 10:56:51 INFO ['2020-02-06T11:00:00.000000000' '2020-02-06T13:00:00.000000000']
2026-07-15 10:56:51 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:51 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS5.nc
2026-07-15 10:56:52 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:52 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:52 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS5.nc
2026-07-15 10:56:52 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:52 INFO 2020-02-06T13:00:00.000000000
2026-07-15 10:56:52 INFO ['2020-02-06T11:00:00.000000000' '2020-02-06T13:00:00.000000000']
2026-07-15 10:56:52 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:52 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS5.nc
2026-07-15 10:56:53 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:53 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:56:53 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:53 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS6.nc
2026-07-15 10:56:53 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:53 INFO 2020-02-06T13:00:00.000000000
2026-07-15 10:56:53 INFO ['2020-02-06T11:00:00.000000000' '2020-02-06T13:00:00.000000000']
2026-07-15 10:56:53 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:53 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS6.nc
2026-07-15 10:56:54 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:55 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:55 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS6.nc
2026-07-15 10:56:55 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:55 INFO 2020-02-06T13:00:00.000000000
2026-07-15 10:56:55 INFO ['2020-02-06T11:00:00.000000000' '2020-02-06T13:00:00.000000000']
2026-07-15 10:56:55 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:55 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS6.nc
2026-07-15 10:56:55 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:55 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:56:55 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:55 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS7.nc
2026-07-15 10:56:56 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:56 INFO 2020-02-06T13:00:00.000000000
2026-07-15 10:56:56 INFO ['2020-02-06T11:00:00.000000000' '2020-02-06T13:00:00.000000000']
2026-07-15 10:56:56 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:56 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS7.nc
2026-07-15 10:56:56 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:56 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:56 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS7.nc
2026-07-15 10:56:56 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:56 INFO 2020-02-06T13:00:00.000000000
2026-07-15 10:56:56 INFO ['2020-02-06T11:00:00.000000000' '2020-02-06T13:00:00.000000000']
2026-07-15 10:56:56 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:56 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS7.nc
2026-07-15 10:56:57 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:57 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:56:57 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:57 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS8.nc
2026-07-15 10:56:57 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:57 INFO 2020-02-06T13:00:00.000000000
2026-07-15 10:56:57 INFO ['2020-02-06T11:00:00.000000000' '2020-02-06T13:00:00.000000000']
2026-07-15 10:56:57 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:57 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS8.nc
2026-07-15 10:56:58 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:58 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:58 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS8.nc
2026-07-15 10:56:58 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:58 INFO 2020-02-06T13:00:00.000000000
2026-07-15 10:56:58 INFO ['2020-02-06T11:00:00.000000000' '2020-02-06T13:00:00.000000000']
2026-07-15 10:56:58 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:58 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS8.nc
2026-07-15 10:56:59 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:56:59 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:56:59 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:56:59 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS9.nc
2026-07-15 10:56:59 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:56:59 INFO 2020-02-06T13:00:00.000000000
2026-07-15 10:56:59 INFO ['2020-02-06T11:00:00.000000000' '2020-02-06T13:00:00.000000000']
2026-07-15 10:56:59 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:56:59 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS9.nc
2026-07-15 10:57:00 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:57:00 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:57:00 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS9.nc
2026-07-15 10:57:00 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:57:00 INFO 2020-02-06T13:00:00.000000000
2026-07-15 10:57:00 INFO ['2020-02-06T11:00:00.000000000' '2020-02-06T13:00:00.000000000']
2026-07-15 10:57:00 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:57:00 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS9.nc
2026-07-15 10:57:01 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:57:01 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:57:01 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:57:01 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS10.nc
2026-07-15 10:57:01 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:57:01 INFO 2020-02-06T13:00:00.000000000
2026-07-15 10:57:01 INFO ['2020-02-06T11:00:00.000000000' '2020-02-06T13:00:00.000000000']
2026-07-15 10:57:01 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:57:01 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS10.nc
2026-07-15 10:57:02 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:57:02 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:57:02 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS10.nc
2026-07-15 10:57:02 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:57:02 INFO 2020-02-06T13:00:00.000000000
2026-07-15 10:57:02 INFO ['2020-02-06T11:00:00.000000000' '2020-02-06T13:00:00.000000000']
2026-07-15 10:57:02 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:57:02 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS10.nc
2026-07-15 10:57:03 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:57:03 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:57:03 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:57:03 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS11.nc
2026-07-15 10:57:03 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:57:03 INFO 2020-02-06T13:00:00.000000000
2026-07-15 10:57:03 INFO ['2020-02-06T11:00:00.000000000' '2020-02-06T13:00:00.000000000']
2026-07-15 10:57:03 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:57:03 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS11.nc
2026-07-15 10:57:04 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:57:04 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:57:04 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS11.nc
2026-07-15 10:57:04 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:57:04 INFO 2020-02-06T13:00:00.000000000
2026-07-15 10:57:04 INFO ['2020-02-06T11:00:00.000000000' '2020-02-06T13:00:00.000000000']
2026-07-15 10:57:04 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:57:04 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS11.nc
2026-07-15 10:57:04 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:57:04 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:57:05 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:57:05 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS12.nc
2026-07-15 10:57:05 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:57:05 INFO 2020-02-06T13:00:00.000000000
2026-07-15 10:57:05 INFO ['2020-02-06T11:00:00.000000000' '2020-02-06T13:00:00.000000000']
2026-07-15 10:57:05 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:57:05 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS12.nc
2026-07-15 10:57:06 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:57:06 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:57:06 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS12.nc
2026-07-15 10:57:06 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:57:06 INFO 2020-02-06T13:00:00.000000000
2026-07-15 10:57:06 INFO ['2020-02-06T11:00:00.000000000' '2020-02-06T13:00:00.000000000']
2026-07-15 10:57:06 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:57:06 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS12.nc
2026-07-15 10:57:07 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:57:07 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:57:07 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:57:07 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS13.nc
2026-07-15 10:57:07 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:57:07 INFO 2020-02-06T13:00:00.000000000
2026-07-15 10:57:07 INFO ['2020-02-06T11:00:00.000000000' '2020-02-06T13:00:00.000000000']
2026-07-15 10:57:07 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:57:08 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS13.nc
2026-07-15 10:57:08 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:57:08 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:57:08 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS13.nc
2026-07-15 10:57:08 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:57:08 INFO 2020-02-06T13:00:00.000000000
2026-07-15 10:57:08 INFO ['2020-02-06T11:00:00.000000000' '2020-02-06T13:00:00.000000000']
2026-07-15 10:57:08 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:57:08 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS13.nc
2026-07-15 10:57:09 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:57:09 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:57:09 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:57:09 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS14.nc
2026-07-15 10:57:09 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:57:09 INFO 2020-02-06T13:00:00.000000000
2026-07-15 10:57:09 INFO ['2020-02-06T11:00:00.000000000' '2020-02-06T13:00:00.000000000']
2026-07-15 10:57:09 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:57:10 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS14.nc
2026-07-15 10:57:10 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:57:10 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:57:10 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS14.nc
2026-07-15 10:57:10 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:57:10 INFO 2020-02-06T13:00:00.000000000
2026-07-15 10:57:10 INFO ['2020-02-06T11:00:00.000000000' '2020-02-06T13:00:00.000000000']
2026-07-15 10:57:10 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:57:10 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS14.nc
2026-07-15 10:57:11 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:57:11 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:57:12 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:57:12 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS15.nc
2026-07-15 10:57:12 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:57:12 INFO 2020-02-06T13:00:00.000000000
2026-07-15 10:57:12 INFO ['2020-02-06T11:00:00.000000000' '2020-02-06T13:00:00.000000000']
2026-07-15 10:57:12 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:57:12 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISA_ENS15.nc
2026-07-15 10:57:12 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:57:13 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-15 10:57:13 INFO Loading orbit memory: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS15.nc
2026-07-15 10:57:13 INFO 2020-02-06T11:00:00.000000000
2026-07-15 10:57:13 INFO 2020-02-06T13:00:00.000000000
2026-07-15 10:57:13 INFO ['2020-02-06T11:00:00.000000000' '2020-02-06T13:00:00.000000000']
2026-07-15 10:57:13 INFO Emission correction applied with pixel-based damping.
2026-07-15 10:57:13 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory/ratio_memory_file_EMISB_ENS15.nc
2026-07-15 10:57:13 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-15 10:57:13 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 10:57:13 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-15 10:57:13 INFO [TIME] step_end current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-06 23:00:00
2026-07-15 10:57:13 INFO [TIME] step_start current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-06 23:00:00
2026-07-15 10:57:13 INFO [TIME] window start=2020-02-06 14:00:00 end=2020-02-07 00:00:00 run_hours=10 has_assimilation=False
2026-07-15 10:57:13 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 14:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-06 23:00:00
2026-07-15 10:57:13 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 00:00:00 days=153073 seconds=0
2026-07-15 10:57:13 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-15 10:57:13 INFO [TIME] increment current_time 2020-02-06 14:00:00 -> 2020-02-07 00:00:00
2026-07-15 10:57:13 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 00:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-06 23:00:00
2026-07-15 10:57:13 INFO ---------->>> Running process_satellite_data()
2026-07-15 10:57:13 INFO [DART] No satellite data found, skipping assimilation
2026-07-15 10:57:13 INFO after_assimilation() skipped
2026-07-15 10:57:13 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-15 10:57:13 INFO [TIME] step_end current_time=2020-02-07 00:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-06 23:00:00
2026-07-15 10:57:13 INFO [PIPELINE] ---- TIME LOOP END ----
+ exit 0
