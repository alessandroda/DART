import xarray as xr

# Load the original dataset
ds = xr.open_dataset(
    "/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/work/test_24h/original_dataset.nc"
)

# Create a new dataset to store shifted data
new_ds = ds.copy(deep=True)

# Loop through all time steps except the first one
for i in range(1, len(ds["time"])):
    # Copy data from the current time step to the previous time step
    new_ds["c_NO2"][24] = ds["c_NO2"][i]
    # Save the new dataset to a new file
    new_ds.to_netcdf(
        f"/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/work/test_24h/{i}.nc"
    )
