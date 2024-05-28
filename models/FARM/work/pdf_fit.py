from scipy.stats import norm
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
out_dir = "/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/IRIDE_ens/loc_plots"
os.makedirs(out_dir, exist_ok=True)
subdirectories = [
    d
    for d in os.listdir(workdir)
    if os.path.isdir(os.path.join(workdir, d)) and d != "plots"
]
values_prior = {f"prior_ensemble_member_{i}": i * 2 + 4 for i in range(1, 17)}
values_posterior = {f"prior_ensemble_member_{i}": i * 2 + 5 for i in range(1, 17)}
# Create a figure with subplots

key = "subplots"

# for key, ax in zip(
#     ["obs_s5p_tropomi", "prior_ensemble_member_1", "posterior_ensemble_member_1"], axs
# ):
out_dir_key = os.path.join(out_dir, key)
os.makedirs(out_dir_key, exist_ok=True)

for subdir in subdirectories:
    if subdir != "20210108_121500":
        continue
    # fig, axs = plt.subplots(
    #     3, 1, figsize=(18, 6), subplot_kw={"projection": ccrs.PlateCarree()}
    # )
    for values_dict in [values_posterior]:
        location_values = {}
        values = []
        try:
            with open(
                os.path.join(workdir, subdir, f"obs_seq_{subdir}.final"), "r"
            ) as file:
                lines = file.readlines()
                for i in range(len(lines)):
                    line = lines[i]
                    if line.startswith(" OBS"):
                        loc_line = lines[i + 43]
                        loc_data = loc_line.split()
                        # location_values[float(loc_data[0]), float(loc_data[1])] = [
                        #     float(lines[i + value]) for value in values_dict.values()
                        # ]
                        data = [
                            float(lines[i + value]) for value in values_dict.values()
                        ]
                        mu, std = norm.fit(data)
                        plt.hist(data, bins=25, density=False, alpha=0.6, color="g")
                        xmin, xmax = plt.xlim()
                        x = np.linspace(xmin, xmax, 100)
                        p = norm.pdf(x, mu, std)
                        plt.plot(x, p, "k", linewidth=2)
                        for i, (x, y) in enumerate(zip(data, np.zeros_like(data))):
                            plt.scatter(x, y, color="red", label="Data points")
                            plt.annotate(
                                f"{i+1}",
                                (x, y),
                                textcoords="offset points",
                                xytext=(0, 10),
                                ha="center",
                            )
                        plt.xlabel("Value")
                        plt.title(
                            "Fit result: mean = {:.5e},  std = {:.5e}".format(mu, std)
                        )
                        plt.grid(True)
                        plt.savefig(
                            os.path.join(
                                out_dir_key,
                                f"posterior_{loc_data[0]}_{loc_data[1]}.png",
                            )
                        )
                        plt.close()
        except FileNotFoundError:
            print(f"File {subdir}/obs_seq_{str(subdir)}.final not found")
