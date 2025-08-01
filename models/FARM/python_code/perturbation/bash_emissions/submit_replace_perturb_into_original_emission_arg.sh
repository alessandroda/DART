#!/bin/bash
source /gporq3/minni/FARM-DART/miniconda3/etc/profile.d/conda.sh
date_start=@date_start
date_end=@date_end
case_dir=@case_dir
var_replace=@emission_var_to_replace
sub_dir_name=@sub_dir_name
cresco_queue=@cresco_queue
for mem in {0..19}
do
    bsub -q  $cresco_queue -n 24 -e err_${mem}_${date_start}.log -o out_${mem}_${date_end}.log /gporq3/minni/FARM-DART/miniconda3/bin/python3.10 /gporq3/minni/FARM-DART/DART/models/FARM/python_code/perturbation/replace_perturb_into_original_emission_arg.py  $date_start $date_end $mem $case_dir $sub_dir_name $var_replace
done

