#!/bin/bash
source /gporq3/minni/FARM-DART/miniconda3/etc/profile.d/conda.sh
date_start='2023081100'
date_end='2023081300'
data_dir='data_3days_3000'
var_replace='veSO2'
sub_dir_name='veSO2_3000'
for mem in {0..19}
#mem=11
do
    bsub -q cameo_h144 -n 24 -e err_${mem}_${date_start}.log -o out_${mem}_${date_end}.log /gporq3/minni/FARM-DART/miniconda3/bin/python3.10 /gporq3/minni/FARM-DART/DART/models/FARM/python_code/perturbation/replace_perturb_into_original_emission_arg.py  $date_start $date_end $mem $data_dir $sub_dir_name $var_replace
done

