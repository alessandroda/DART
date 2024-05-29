import numpy as np
import xarray as xr


def find_closest_point(x, y, x_series, y_series):
    """
    Find the closest point in x_series and y_series
    to the given x and y coordinates.
    """
    idx_lon = np.abs(np.unique(x_series.compute().values) - x).argmin()
    idx_lat = np.abs(np.unique(y_series.compute().values) - y).argmin()

    x_closest = x_series.compute().values[idx_lon]
    y_closest = y_series.compute().values[idx_lat]
    return x_closest, y_closest


analysis = xr.open_dataset(
    "/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/work/rf_ens_experiment/analysis.nc"
)
preassim = xr.open_dataset(
    "/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/work/rf_ens_experiment/preassim.nc"
)
x = 17.2536
y = 40.5211
x_closest, y_closest = find_closest_point(x, y, analysis["x"], analysis["y"])
print(x_closest, y_closest)
an_no2 = analysis.sel(x=x_closest, y=y_closest)
pr_no2 = preassim.sel(x=x_closest, y=y_closest)
print(an_no2["c_NO2"][:, 0, 0].values)

print(pr_no2["c_NO2"][:, 0, 0].values)
