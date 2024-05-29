import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
import pandas as pd
from pydantic import BaseModel
import matplotlib.colors as mcolors
from datetime import datetime


class Settings(BaseModel):
    lonmin: float = 5.325  # 5.325
    lonmax: float = 19
    latmin: float = 35.380
    latmax: float = 47.8

    range_conc = {
        "obs_s5p_tropomi": (0, 0.0002),
        "prior_ensemble_mean": (0, 0.0002),
        "posterior_ensemble_mean": (0, 0.0002),
        "prior_ensemble_spread": (0, 0.0002),
        "posterior_ensemble_spread": (0, 0.0002),
        "prior_ensemble_member_1": (0, 0.0002),
        "posterior_ensemble_member_1": (0, 0.0002),
    }


colors = [
    (0, "blue"),
    (0.15, "cyan"),
    (0.3, "lightgreen"),
    (0.45, "yellow"),
    (0.6, "orange"),
    (0.75, "red"),
    (1, "purple"),
]

cmap = mcolors.LinearSegmentedColormap.from_list("custom_map", colors)

settings = Settings()
workdir = "/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/IRIDE_ens/posteriors"
out_dir = "/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/IRIDE_ens/plots"
os.makedirs(out_dir, exist_ok=True)
subdirectories = [
    d
    for d in os.listdir(workdir)
    if os.path.isdir(os.path.join(workdir, d)) and d != "plots"
]
values = {
    "obs_s5p_tropomi": 1,
    "prior_ensemble_mean": 2,
    "posterior_ensemble_mean": 3,
    "prior_ensemble_spread": 4,
    "posterior_ensemble_spread": 5,
    "prior_ensemble_member_1": 6,
    "posterior_ensemble_member_1": 7,
}

# Create a figure with subplots

key = "subplots"

# for key, ax in zip(
#     ["obs_s5p_tropomi", "prior_ensemble_member_1", "posterior_ensemble_member_1"], axs
# ):
list_obs_vcd = [
    "obs_s5p_tropomi",
    "prior_ensemble_member_1",
    "posterior_ensemble_member_1",
]
out_dir_key = os.path.join(out_dir, key)
os.makedirs(out_dir_key, exist_ok=True)
time_values = []
# only subdirectories from 0 to 10
obs_total = {key: [] for key in list_obs_vcd}
for subdir in subdirectories:

    for key in list_obs_vcd:
        obs_values = []
        try:
            with open(
                os.path.join(workdir, subdir, f"obs_seq_{subdir}.final"), "r"
            ) as file:
                lines = file.readlines()
                for i in range(len(lines)):
                    line = lines[i]
                    if line.startswith("loc3d"):
                        next_line = lines[i + 1]
                        data = next_line.split()
                    elif line.startswith(" OBS"):
                        obs_values.extend(map(float, lines[i + values[key]].split()))
                obs_total[key].append(np.array(obs_values).mean())
                time_values.append(
                    pd.to_datetime(subdir, format="%Y%m%d_%H%M%S")
                )  # Store time values
        except FileNotFoundError:
            print(f"File {subdir}/obs_seq_{str(subdir)}.final not found")
        # plot timeseries plot
        # Plot timeseries plot
    # Plot timeseries plot for each observation
time_values = list(set(time_values))


# Sort the time_values list
time_values.sort()

# Convert time_values to datetime objects
time_values_dt = [datetime.strptime(str(t), "%Y-%m-%d %H:%M:%S") for t in time_values]

jan_2021_mask = [(dt.year == 2021 and dt.month == 2) for dt in time_values_dt]
time_values = [t for t, mask in zip(time_values, jan_2021_mask) if mask]
for key in list_obs_vcd:
    obs_total[key] = [obs for obs, mask in zip(obs_total[key], jan_2021_mask) if mask]

# Convert time_values to just dates
dates = [dt.date() for dt in time_values]

# Initialize dictionaries to store daily mean values for each observation
daily_means = {key: [] for key in list_obs_vcd}

# Compute daily means
for key in list_obs_vcd:
    daily_values = {}
    for dt, value in zip(dates, obs_total[key]):
        if dt not in daily_values:
            daily_values[dt] = []
        daily_values[dt].append(value)
    daily_means[key] = [np.mean(vals) for vals in daily_values.values()]


# Plot daily mean for each observation
plt.figure(figsize=(20, 6))  # Adjust figure size as needed
for key in list_obs_vcd:
    plt.scatter(
        list(daily_values.keys()),
        daily_means[key],
        label=key,
        marker="o",
        linestyle="--",
    )
# tick interval 1 day
plt.xlabel("Date")
plt.ylabel("vcd [mol/m2]")
plt.legend()
plt.xticks(list(daily_values.keys()), rotation="vertical")
plt.tight_layout()
plt.show()
plt.savefig("test.png")
print("did")
