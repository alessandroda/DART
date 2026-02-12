import xarray as xr
from pathlib import Path
import pandas as pd
from tqdm import tqdm
from datetime import datetime, timedelta
import os

file_path_original = Path('/g100_work/ARPAE_AQM/MIMESI/runs_chimere/ITA7/basecase/data_ITA7_ITA7/')
file_path_to_be_replaced = Path(f'/g100_work/ARPAE_AQM/MIMESI/perturbation_fields/emi_mems')
date_start = '2026012500'
date_end = '2026012600'
data_dir = 'ITA7'
start_dt = datetime.strptime(date_start, '%Y%m%d%H')
end_dt = datetime.strptime(date_end, '%Y%m%d%H')
var_replace = 'NO2'
sub_dir_name = 'NO2_1000'
current_dt = start_dt
no_mems = 3
while current_dt < end_dt:
    date = current_dt.strftime('%Y%m%d%H')
    datep1 = (current_dt + timedelta(days=1)).strftime('%Y%m%d%H')
    file_original_name = f'AEMISSIONS.{date}_{datep1}_ITA7.nc'
    ds = xr.open_dataset(file_path_original / file_original_name)
    for mem in tqdm(range(no_mems)):
        file_dest = Path(f'/g100_work/ARPAE_AQM/MIMESI/runs_chimere/{data_dir}/RUN_{mem}/EMISSION_{mem}/')
        if not os.path.exists(file_dest):
            print(f'dir_path: {file_dest}')
            os.makedirs(file_dest)
        ds_perturbed = xr.open_dataset(file_path_to_be_replaced / f'emi_{mem}'/ sub_dir_name / f'AEMISSIONS.{date}_{datep1}_ITA7_{mem}.nc')
        ds[var_replace].values = ds_perturbed[var_replace].values
        ds.to_netcdf( file_dest / file_original_name)   
    ds.close()
    current_dt += timedelta(days=1)

