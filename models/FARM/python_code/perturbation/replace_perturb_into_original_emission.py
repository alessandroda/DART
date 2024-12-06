import xarray as xr
from pathlib import Path
import pandas as pd
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta
file_path_original = Path('/gporq3/minni/FARM-DART/perturbation_fields/data/emission_base_test/2023/emi/08/')
file_path_to_be_replaced = Path(f'/gporq3/minni/FARM-DART/perturbation_fields/emi_mems')
date_start = '2023081100'
date_end = '2023082000'

start_dt = datetime.strptime(date_start, '%Y%m%d%H')
end_dt = datetime.strptime(date_end, '%Y%m%d%H')

current_dt = start_dt
no_mems = 20
while current_dt <= end_dt:
    date = current_dt.strftime('%Y%m%d%H')
    var_replace = 'veSO2'
    file_original_name = f'HERMESv3_{date}.nc'
    ds = xr.open_dataset(file_path_original / file_original_name)
    for mem in tqdm(range(no_mems)):
        file_dest = Path(f'/gporq3/minni/FARM-DART/RUN/data/INPUT/HERMES/emi_{mem}')
        ds_perturbed = xr.open_dataset(file_path_to_be_replaced / f'emi_{mem}'/  f'HERMESv3_{date}_{mem}.nc')
        ds[var_replace].values = ds_perturbed[var_replace].values
        ds.to_netcdf( file_dest / file_original_name)   
    ds.close()
    current_dt += timedelta(days=1)

