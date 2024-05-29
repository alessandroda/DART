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
out_dir_key = os.path.join(out_dir, key)
os.makedirs(out_dir_key, exist_ok=True)
for subdir in subdirectories:
    fig, axs = plt.subplots(
        1, 3, figsize=(18, 8), subplot_kw={"projection": ccrs.PlateCarree()}
    )
    for key, ax in zip(
        ["obs_s5p_tropomi", "prior_ensemble_member_1", "posterior_ensemble_member_1"],
        axs,
    ):
        x_values = []
        y_values = []
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
                        x_values.append(float(data[0]))
                        y_values.append(float(data[1]))
                    elif line.startswith(" OBS"):
                        obs_values.extend(map(float, lines[i + values[key]].split()))

            ax.set_title(f"{key} time: {subdir}")
            ax.set_extent(
                [settings.lonmin, settings.lonmax, settings.latmin, settings.latmax],
                ccrs.PlateCarree(),
            )
            ax.coastlines(resolution="10m")
            ax.add_feature(cfeature.BORDERS, linewidth=0.5)
            ax.add_feature(cfeature.LAND, facecolor="white")
            ax.add_feature(cfeature.OCEAN, facecolor="lightblue")
            ax.add_feature(cfeature.LAKES, facecolor="lightblue")
            ax.add_feature(cfeature.RIVERS, edgecolor="lightblue")

            x_values_deg = np.degrees(x_values)
            y_values_deg = np.degrees(y_values)
            norm = mcolors.PowerNorm(
                gamma=0.5,
                vmin=settings.range_conc[key][0],
                vmax=settings.range_conc[key][1],
            )

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

            x_ticks = np.linspace(settings.lonmin, settings.lonmax, 10)
            y_ticks = np.linspace(settings.latmin, settings.latmax, 10)
            ax.set_xticks(np.floor(x_ticks))
            ax.set_yticks(np.floor(y_ticks))
            ax.set_xticklabels(np.floor(x_ticks), rotation=90, fontsize=8)
            ax.set_yticklabels(np.floor(y_ticks), fontsize=8)

            ax.set_xlabel("Longitude", fontsize=8)
            ax.set_ylabel("Latitude", fontsize=8)
            # ax.ticklabel_format(style="plain")

            cbar = plt.colorbar(sc, orientation="horizontal", shrink=0.5, ax=ax)
            # change fontsize to 8 cbar
            cbar.ax.tick_params(labelsize=8)
            cbar.set_label(f"NO2 [mol/m2]", fontsize=8)

        except FileNotFoundError:
            print(f"File {subdir}/obs_seq_{str(subdir)}.final not found")
    plt.savefig(os.path.join(out_dir_key, f"{subdir}_{key}.png"))
    plt.close()
