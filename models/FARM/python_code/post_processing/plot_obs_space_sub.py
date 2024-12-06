from matplotlib.colorbar import ColorbarBase
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
from pydantic import BaseModel
import matplotlib.colors as mcolors

def generate_file_name(input_code):
    last_two_digits=int(input_code[-2:])
    if 0 <= last_two_digits <= 9:
        case_suffix = f"0{last_two_digits}0000"
    else:
        case_suffix = f"{last_two_digits}0000"

    file_name = f"{input_code[:-2]}_{case_suffix}"
    return file_name



class Settings():
    lonmin: float = -25  # 5.325
    lonmax: float = 45
    latmin: float = 30
    latmax: float = 72

    range_conc = {
        "obs_s5p_tropomi": (0, 0.0008),
        "prior_ensemble_mean": (0, 0.0008),
        "posterior_ensemble_mean": (0, 0.0008),
        "prior_ensemble_spread": (0, 0.0008),
        "posterior_ensemble_spread": (0, 0.0008),
        "prior_ensemble_member_1": (0, 0.0008),
        "posterior_ensemble_member_1": (0, 0.0008),
    }


colors = [
    (0, "blue"),
    (0.15, "cyan"),
    (0.3, "darkgreen"),
    (0.45, "yellow"),
    (0.6, "orange"),
    (0.75, "red"),
    (1, "purple"),
]

cmap = mcolors.LinearSegmentedColormap.from_list("custom_map", colors)

settings = Settings()
workdir = "/gporq3/minni/FARM-DART/RUN/data/posteriors"
out_dir = "/gporq3/minni/FARM-DART/post_proc/report_plots"
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
norm = None
key = "subplots"

# for key, ax in zip(
#     ["obs_s5p_tropomi", "prior_ensemble_member_1", "posterior_ensemble_member_1"], axs
# ):
out_dir_key = os.path.join(out_dir, key)
os.makedirs(out_dir_key, exist_ok=True)
for subdir in subdirectories:
    fig, axs = plt.subplots(
        1, 3, figsize=(22, 10), subplot_kw={"projection": ccrs.PlateCarree()}
    )
    for key, ax in zip(
        ["obs_s5p_tropomi", "prior_ensemble_mean", "posterior_ensemble_mean"],
        axs,
    ):
        x_values = []
        y_values = []
        obs_values = []
        try:
            date_full = generate_file_name(str(subdir))
            with open(
                os.path.join(workdir, subdir, f"obs_seq_{date_full}.final"), "r"
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
            ax.add_feature(cfeature.LAND, facecolor="lightgrey")
            ax.add_feature(cfeature.OCEAN, facecolor="lightgrey")
            ax.add_feature(cfeature.LAKES, facecolor="lightgrey")
            ax.add_feature(cfeature.RIVERS, edgecolor="lightgrey")

            x_values_deg = np.degrees(x_values)
            y_values_deg = np.degrees(y_values)
            norm = mcolors.Normalize(
                vmin=settings.range_conc[key][0],
                vmax=settings.range_conc[key][1],
            )

            sc = ax.scatter(
                x_values_deg,
                y_values_deg,
                # multiply array by 100
                c=obs_values,
                marker="s",
                s=1,
                vmin=settings.range_conc[key][0],
                vmax=settings.range_conc[key][1],
                cmap=cmap,
                transform=ccrs.PlateCarree(),
            )

            x_ticks = np.linspace(settings.lonmin, settings.lonmax, 15)
            y_ticks = np.linspace(settings.latmin, settings.latmax, 15)
            ax.set_xticks(np.floor(x_ticks))
            ax.set_yticks(np.floor(y_ticks))
            ax.set_xticklabels(np.floor(x_ticks), rotation=90, fontsize=12)
            ax.set_yticklabels(np.floor(y_ticks), fontsize=12)

            ax.set_xlabel("Longitude", fontsize=12)
            ax.set_ylabel("Latitude", fontsize=12)
            # ax.ticklabel_format(style="plain")

        except FileNotFoundError:
            print(f"File {subdir}/obs_seq_{str(subdir)}.final not found")
    # Add a single colorbar for the entire figure
    fig.tight_layout(pad=1.0)
    cbar_ax = fig.add_axes([0.25, 0.05, 0.5, 0.03])  # [left, bottom, width, height]
    # change size of ColorBase

    cbar = ColorbarBase(cbar_ax, cmap=cmap, norm=norm, orientation="horizontal")
    # set cbar font size
    cbar.ax.tick_params(labelsize=12)
    #cbar.set_ticks(np.linspace(0, 0.0002, 4))
    cbar.set_label(f"SO2 [mol/m2]", fontsize=12)
    plt.savefig(os.path.join(out_dir_key, f"{subdir}_{key}.png"))
    plt.close()
