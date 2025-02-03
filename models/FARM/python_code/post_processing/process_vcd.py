
from matplotlib.colorbar import ColorbarBase
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
from pydantic import BaseModel
import matplotlib.colors as mcolors
import pandas as pd
import xarray as xr
from tqdm import tqdm

def total_vcd_plot(all_dfs: dict):
    '''
    Create a dataframe of the structure
    time obs_so2, prior_so2, posterior_so2
    '''
    df = pd.DataFrame(columns = ['time', 'obs_so2', 'prior_so2',
        'posterior_so2'])
    for key, ax in zip(all_dfs.keys(), axs):
        print(f'processing time: {all_dfs[key]["time"].unique()}')





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
    min_scale: float = 0
    max_scale: float = 0.001
    start_time = '2023080100'
    end_time = '2023083123'
    range_conc = {
        "obs_s5p_tropomi": (0, 0.0008),
        "prior_ensemble_mean": (0, 0.0008),
        "posterior_ensemble_mean": (0, 0.0008),
        "prior_ensemble_spread": (0, 0.0008),
        "posterior_ensemble_spread": (0, 0.0008),
        "prior_ensemble_member_1": (0, 0.0008),
        "posterior_ensemble_member_1": (0, 0.0008),
    }

    values = {
        "obs_s5p_tropomi": 1,
        "prior_ensemble_mean": 2,
        "posterior_ensemble_mean": 3,
        "prior_ensemble_spread": 4,
        "posterior_ensemble_spread": 5,
        "prior_ensemble_member_1": 6,
        "posterior_ensemble_member_1": 7,
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

#cmap = mcolors.LinearSegmentedColormap.from_list("custom_map", colors)
cmap = plt.cm.jet
settings = Settings()
workdir = "/gporq3/minni/FARM-DART/RUN/data/posteriors"
out_dir = "/gporq3/minni/FARM-DART/post_proc/report_plots"
os.makedirs(out_dir, exist_ok=True)
subdirectories = [
    d
    for d in os.listdir(workdir)
    if os.path.isdir(os.path.join(workdir, d)) and d != "plots"
]
key = "subplots"

out_dir_key = os.path.join(out_dir, key)
os.makedirs(out_dir_key, exist_ok=True)

def regular_grid():
    # Define grid boundaries
    lon_bins = np.linspace(settings.lonmin, settings.lonmax, 461)  # 100 grid cells for longitude
    lat_bins = np.linspace(settings.latmin, settings.latmax, 421)  # 100 grid cells for latitude

    # # Create a meshgrid for plotting
    # lon_centers = 0.5 * (lon_bins[:-1] + lon_bins[1:])
    # lat_centers = 0.5 * (lat_bins[:-1] + lat_bins[1:])
    # return np.meshgrid(lon_centers, lat_centers)
    return lon_bins, lat_bins


def transform_obs_seq_out(subdirectories, lon_grid, lat_grid):
    start_time = pd.to_datetime(settings.start_time, format="%Y%m%d%H")
    end_time = pd.to_datetime(settings.end_time, format="%Y%m%d%H")
    index = pd.date_range(start=start_time, end=end_time, freq="H")
    result =  pd.DataFrame(columns = ['obs_s5p_tropomi', 'prior_ensemble_mean',
        'posterior_ensemble_mean'], index=index)
    paths = {
        "obs_s5p_tropomi" : [],
        "prior_ensemble_mean" : [],
        "posterior_ensemble_mean" : []
    }
    for subdir in tqdm(subdirectories):
        for key in ["obs_s5p_tropomi", "prior_ensemble_mean", "posterior_ensemble_mean"]:
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
                            obs_values.extend(map(float, lines[i + settings.values[key]].split()))
                if x_values and y_values and obs_values:
                    from scipy.stats import binned_statistic_2d

                    # Compute grid-averaged values
                    x_values_deg = np.degrees(x_values)
                    y_values_deg = np.degrees(y_values)
                    stat, x_bin,y_bin, _ = binned_statistic_2d(x_values_deg,
                            y_values_deg, obs_values,  statistic='mean',
                            bins=[lon_grid, lat_grid])

                    import xarray as xr
                    # Convert to an xarray DataArray for better handling of coordinates and rename var with key
                    grid_data = xr.DataArray(
                        stat.T,  # Transposed to match (lat, lon) convention
                        name=key,
                        coords={'latitude': lat_grid[:-1], 'longitude':
                            lon_grid[:-1], 'time': pd.to_datetime(subdir, format="%Y%m%d%H")},
                        dims=['latitude', 'longitude']
                    )
                    #average variable on coords
                    mean_on_domain = grid_data.mean(dim=['latitude', 'longitude'])
                    #at corresponding time index set the mean_on_domain
                    result.loc[pd.to_datetime(subdir, format="%Y%m%d%H"), key] = mean_on_domain.values.item()
                    # Save to NetCDF
                    paths[key].append(os.path.join(out_dir_key, f"{key}_{subdir}.nc"))
                    grid_data.to_netcdf(os.path.join(out_dir_key, f"{key}_{subdir}.nc"))

                else:
                    print(f"No data found in subdir {subdir}")
            except FileNotFoundError:
                print(f"File {subdir}/obs_seq_{str(subdir)}.final not found")
    return result, paths

#def subplot_observation_space(paths: dict):
#    fig, axs = plt.subplots(
#        1, 3, figsize=(22, 10), subplot_kw={"projection": ccrs.PlateCarree()}
#    )
#    for name, list_paths, ax in zip(paths.keys(), paths.values(), axs):
#        for path_nc in list_paths:
#            dataset = xr.open_dataset(path_nc)
#            breakpoint()
#            ax.set_title(f"{key} time: {dataset.time.values}")
#            ax.set_extent(
#                [settings.lonmin, settings.lonmax, settings.latmin, settings.latmax],
#                ccrs.PlateCarree(),
#            )
#            ax.coastlines(resolution="10m")
#            ax.add_feature(cfeature.BORDERS, linewidth=0.5)
#            ax.add_feature(cfeature.LAND, facecolor="lightgrey")
#            ax.add_feature(cfeature.OCEAN, facecolor="lightgrey")
#            ax.add_feature(cfeature.LAKES, facecolor="lightgrey")
#            ax.add_feature(cfeature.RIVERS, edgecolor="lightgrey")
#            breakpoint()
#
#            norm = mcolors.Normalize(
#                vmin=settings.range_conc[key][0],
#                vmax=settings.range_conc[key][1],
#            )
#
#            # Create the plot
#            img = plt.pcolormesh(
#                dataset.variables["x"],
#                dataset.variables["y"],
#                dataset.variables[name],
#                alpha=0.8,
#                cmap=cmap,
#                vmin=settings.min_scale,
#                vmax=settings.max_scale,
#                shading="auto",
#                zorder=1,
#            )
#
#            # Add a colorbar
#            cbar = plt.colorbar(img, ax=plt.gca())
#
#
#            x_ticks = np.linspace(settings.lonmin, settings.lonmax, 15)
#            y_ticks = np.linspace(settings.latmin, settings.latmax, 15)
#            ax.set_xticks(np.floor(x_ticks))
#            ax.set_yticks(np.floor(y_ticks))
#            ax.set_xticklabels(np.floor(x_ticks), rotation=90, fontsize=12)
#            ax.set_yticklabels(np.floor(y_ticks), fontsize=12)
#
#            ax.set_xlabel("Longitude", fontsize=12)
#            ax.set_ylabel("Latitude", fontsize=12)
#
#            fig.tight_layout(pad=1.0)
#            cbar_ax = fig.add_axes([0.25, 0.05, 0.5, 0.03])  # [left, bottom, width, height]
#            # change size of ColorBase
#
#            cbar = ColorbarBase(cbar_ax, cmap=cmap, norm=norm, orientation="horizontal")
#            # set cbar font size
#            cbar.ax.tick_params(labelsize=12)
#            #cbar.set_ticks(np.linspace(0, 0.0002, 4))
#            cbar.set_label(f"SO2 [mol/m2]", fontsize=12)
#            # replace .nc in path_name with .png
#            path_name = path_nc.replace(".nc", ".png")
#            plt.savefig(path_name)
#            plt.close()

def plot_observation_space(paths: dict):
    for name, list_paths in paths.items():
        for path_nc in list_paths:
            fig, ax = plt.subplots(
                1, 1, figsize=(10, 8), subplot_kw={"projection": ccrs.PlateCarree()}
            )

            dataset = xr.open_dataset(path_nc)
            ax.set_title(f"{name} time: {dataset.time.values}")
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

            norm = mcolors.Normalize(
                vmin=settings.range_conc[name][0],
                vmax=settings.range_conc[name][1],
            )

            # Create the plot
            img = ax.pcolormesh(
                dataset.variables["longitude"],
                dataset.variables["latitude"],
                dataset.variables[name],
                alpha=0.8,
                cmap=cmap,
                vmin=settings.min_scale,
                vmax=settings.max_scale,
                shading="auto",
                zorder=1,
            )

            # Add a colorbar
            cbar = plt.colorbar(img, ax=ax, orientation="horizontal", pad=0.05)

            x_ticks = np.linspace(settings.lonmin, settings.lonmax, 15)
            y_ticks = np.linspace(settings.latmin, settings.latmax, 15)
            ax.set_xticks(np.floor(x_ticks))
            ax.set_yticks(np.floor(y_ticks))
            ax.set_xticklabels(np.floor(x_ticks), rotation=90, fontsize=12)
            ax.set_yticklabels(np.floor(y_ticks), fontsize=12)

            ax.set_xlabel("Longitude", fontsize=12)
            ax.set_ylabel("Latitude", fontsize=12)

            fig.tight_layout(pad=1.0)
            cbar.set_label(f"SO2 [mol/m2]", fontsize=12)

            # Replace .nc in path_name with .png
            path_name = path_nc.replace(".nc", ".png")
            plt.savefig(path_name)
            plt.close()

lon_grid, lat_grid = regular_grid()
domain_averaged, paths = transform_obs_seq_out(subdirectories, lon_grid, lat_grid)
domain_averaged.to_csv('test.csv')
plot_observation_space(paths)
#plot_ts_domain_averaged(domain_averaged)

