import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
from pydantic import BaseModel
import matplotlib.colors as mcolors


class Settings(BaseModel):

    lonmin: float = 5.325  # 5.325
    lonmax: float = 19
    latmin: float = 35.380
    latmax: float = 48
    # Domain boundary Spott UTM32
    range_conc = {
        "obs_s5p_tropomi": (0, 0.0002),
        "prior_ensemble_mean": (0, 0.0002),
        "posterior_ensemble_mean": (0, 0.0002),
        "prior_ensemble_spread": (0, 0.0002),
        "posterior_ensemble_spread": (0, 0.0002),
        "prior_ensemble_member_1": (0, 0.0002),
        "posterior_ensemble_member_1": (0, 0.0002),
    }


# Generate a similar colormap
# Define colors and corresponding values
colors = [
    (0, "blue"),
    (0.15, "cyan"),
    (0.3, "lightgreen"),
    (0.45, "yellow"),
    (0.6, "orange"),
    (0.75, "red"),
    (1, "purple"),
]

# Create a colormap using linear interpolation
cmap = mcolors.LinearSegmentedColormap.from_list("custom_map", colors)

settings = Settings()
workdir = "/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/IRIDE_ens/posteriors"
out_dir = workdir + "/plots"
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
for key, value in values.items():
    out_dir_key = out_dir + "/" + key
    os.makedirs(out_dir_key, exist_ok=True)
    for subdir in subdirectories:
        # Extracting data from the file
        x_values = []
        y_values = []
        obs_values = []
        try:
            with open(
                f"/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/IRIDE_ens/posteriors/{subdir}/obs_seq_{str(subdir)}.final",
                "r",
            ) as file:
                lines = file.readlines()
                for i in range(len(lines)):
                    line = lines[i]
                    if line.startswith("loc3d"):
                        next_line = lines[i + 1]
                        data = next_line.split()
                        x_values.append(float(data[0]))
                        y_values.append(float(data[1]))
                    elif line.startswith(" OBS"):
                        obs_values.extend(map(float, lines[i + value].split()))

            # Plotting the map
            plt.figure(figsize=(10, 8))

            # Map
            ax = plt.axes(projection=ccrs.PlateCarree())
            ax.set_title(f"{key}").set_fontsize(20)
            ax.set_extent(
                [
                    settings.lonmin,
                    settings.lonmax,
                    settings.latmin,
                    settings.latmax,
                ],
                ccrs.PlateCarree(),
            )
            # Adding coastlines and features
            ax.coastlines(resolution="10m")
            ax.add_feature(cfeature.BORDERS, linewidth=0.5)
            ax.add_feature(cfeature.LAND, facecolor="white")
            ax.add_feature(cfeature.OCEAN, facecolor="lightblue")
            ax.add_feature(cfeature.LAKES, facecolor="lightblue")
            ax.add_feature(cfeature.RIVERS, edgecolor="lightblue")

            # Convert longitude and latitude from radians to degrees
            x_values_deg = np.degrees(x_values)
            y_values_deg = np.degrees(y_values)
            norm = mcolors.PowerNorm(
                gamma=0.5,
                vmin=settings.range_conc[key][0],
                vmax=settings.range_conc[key][1],
            )
            # Scatter plot on the map
            sc = ax.scatter(
                x_values_deg,
                y_values_deg,
                c=obs_values,
                marker="s",
                s=4,
                vmin=settings.range_conc[key][0],
                vmax=settings.range_conc[key][1],
                cmap=cmap,
                transform=ccrs.PlateCarree(),
            )

            # Colorbar
            cbar = plt.colorbar(sc, orientation="horizontal", shrink=0.8, pad=0.05)
            cbar.set_label(key, fontsize=16)

            cbar.set_label(f"NO2 [mol/m2]", fontsize=16)
            # Adjust labels
            ax.set_xlabel("Longitude")
            ax.set_ylabel("Latitude")
            plt.savefig(out_dir_key + f"/{subdir}_{key}.png")
        except FileNotFoundError:
            print(f"File {subdir}/obs_seq_{str(subdir)}.final not found")
