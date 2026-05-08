from matplotlib.colorbar import ColorbarBase
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
from pydantic import BaseModel
import matplotlib.colors as mcolors
import pandas as pd
import xarray as xr
from tqdm import tqdm
from scipy.interpolate import griddata
import os
import cartopy.feature as cfeature
Sicily = {
    "name": "Sicily",
    "lonmin_grid": 12.0,
    "lonmax_grid": 19.0,
    "latmin_grid": 33.0,
    "latmax_grid": 38.5,
}

Domain_minni = {
    "name": "Minni",
    "lonmin_grid": -25.0,
    "lonmax_grid": 45.0,
    "latmin_grid": 30.0,
    "latmax_grid": 72.0,
}

Domain_italy = {
    "name": "Italy",
    "lonmin_grid": 6.6,
    "lonmax_grid": 18.8,
    "latmin_grid": 36.6,
    "latmax_grid": 47.1,
}

Domain_kairos = {

    "name" : "kairos",
    "lonmin_grid" : 6.6,
    "lonmax_grid" : 18.8,
    "latmin_grid" : 36.6,
    "latmax_grid" : 47.1,
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

def plot_monthly_mean_from_nc(out_dir, month="2025-11"):
    import glob

    # ==========================================
    # 1. LOAD ALL DAILY NETCDF
    # ==========================================
    files = sorted(glob.glob(os.path.join(out_dir, "*.nc")))

    if not files:
        raise ValueError("No NetCDF files found!")

    ds = xr.open_mfdataset(files, combine="nested", concat_dim="time")

    print(ds)

    # ==========================================
    # 2. SELECT MONTH
    # ==========================================
    ds_month = ds.sel(time=month)

    print(f"Selected {ds_month.time.size} timesteps")

    # ==========================================
    # 3. FILTER LOW COVERAGE (IMPORTANT)
    # ==========================================
    if "nobs" in ds_month:
        mask = ds_month["nobs"].sum(dim="time") > 20
        ds_month = ds_month.where(mask)

    # ==========================================
    # 4. MONTHLY MEAN
    # ==========================================
    ds_mean = ds_month.mean(dim="time", skipna=True)

    # ==========================================
    # 5. VARIABLES
    # ==========================================
    obs = ds_mean["obs_s5p_tropomi"]
    prior = ds_mean["prior_ensemble_mean"]
    post = ds_mean["posterior_ensemble_mean"]

    diff = post - prior  # safer than %

    # ==========================================
    # 6. SMOOTHING (optional but recommended)
    # ==========================================
    obs = obs.rolling(latitude=3, longitude=3, center=True).mean()
    prior = prior.rolling(latitude=3, longitude=3, center=True).mean()
    post = post.rolling(latitude=3, longitude=3, center=True).mean()
    diff = diff.rolling(latitude=3, longitude=3, center=True).mean()

    # ==========================================
    # 7. PLOT
    # ==========================================
    fig, axs = plt.subplots(
        1, 4, figsize=(18, 5),
        subplot_kw={"projection": ccrs.PlateCarree()}
    )

    titles = [
        "Satellite observations",
        "Model forecast",
        "Assimilated analysis",
        "Analysis - Background"
    ]

    data_list = [obs, prior, post, diff]
    cmaps = ["viridis", "viridis", "viridis", "RdBu_r"]

    for ax, data, title, cmap in zip(axs, data_list, titles, cmaps):

        ax.set_extent([
            settings.grid_bounds["lonmin_grid"],
            settings.grid_bounds["lonmax_grid"],
            settings.grid_bounds["latmin_grid"],
            settings.grid_bounds["latmax_grid"]
        ], crs=ccrs.PlateCarree())

        ax.coastlines(resolution="10m", linewidth=0.5)
        ax.add_feature(cfeature.BORDERS, linewidth=0.3)

        if title == "Analysis - Background":
            img = ax.pcolormesh(
                data.longitude, data.latitude, data,
                cmap=cmap, vmin=-2e-5, vmax=2e-5,
                shading="auto"
            )
        else:
            img = ax.pcolormesh(
                data.longitude, data.latitude, data,
                cmap=cmap, vmin=0, vmax=5e-5,
                shading="auto"
            )

        ax.set_title(title)

        cbar = plt.colorbar(img, ax=ax, orientation="horizontal", pad=0.05)
        if title == "Analysis - Background":
            cbar.set_label("Δ NO₂ [mol/m²]")
        else:
            cbar.set_label("NO₂ [mol/m²]")

    fig.suptitle(f"Monthly mean NO₂ — {month}")

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"monthly_mean_{month}.png"), dpi=300)
    plt.close()

    print("Monthly plot saved")

def centers_to_edges(centers):
    centers = np.asarray(centers)
    # midpoints between centers
    mids = (centers[1:] + centers[:-1]) / 2.0
    # extrapolate edges at the ends
    first = centers[0] - (centers[1] - centers[0]) / 2.0
    last = centers[-1] + (centers[-1] - centers[-2]) / 2.0
    edges = np.concatenate(([first], mids, [last]))
    return edges


def generate_file_name(input_code):
    last_two_digits = int(input_code[-2:])
    if 0 <= last_two_digits <= 9:
        case_suffix = f"0{last_two_digits}0000"
    else:
        case_suffix = f"{last_two_digits}0000"

    file_name = f"{input_code[:-2]}_{case_suffix}"
    return file_name


class Settings:
    lonmin: float = -25.0  # 5.325
    lonmax: float = 45.0
    latmin: float = 30.0
    latmax: float = 72.0

    grid_bounds: dict = Domain_kairos

    min_scale: float = 0
    max_scale: float = 0.00001
    start_time = "2025110200"
    end_time = "2025113000"
    range_conc = {
        "obs_s5p_tropomi": (0, 0.0001),
        "prior_ensemble_mean": (0, 0.0001),
        "posterior_ensemble_mean": (0, 0.0001),
        "prior_ensemble_spread": (0, 0.0001),
        "posterior_ensemble_spread": (0, 0.0001),
        "prior_ensemble_member_1": (0, 0.0001),
        "posterior_ensemble_member_1": (0, 0.0001),
    }

    values = {
        "obs_s5p_tropomi": 1,
        "prior_ensemble_mean": 2,
        "posterior_ensemble_mean": 3,
        "prior_ensemble_spread": 4,
        "posterior_ensemble_spread": 5,
        "prior_ensemble_member_1": 6,
        "posterior_ensemble_member_1": 7,
        "prior_ensemble_member_2": 8,
        "posterior_ensemble_member_2": 9,
        "prior_ensemble_member_3": 10,
        "posterior_ensemble_member_3": 11,
        "prior_ensemble_member_4": 12,
        "posterior_ensemble_member_4": 13,
        "prior_ensemble_member_5": 14,
        "posterior_ensemble_member_5": 15,
        "prior_ensemble_member_6": 16,
        "posterior_ensemble_member_6": 17,
        "prior_ensemble_member_7": 18,
        "posterior_ensemble_member_7": 19,
        "prior_ensemble_member_8": 20,
        "posterior_ensemble_member_8": 21,
        "prior_ensemble_member_9": 22,
        "posterior_ensemble_member_9": 23,
        "prior_ensemble_member_10": 24,
        "posterior_ensemble_member_10": 25,
        "prior_ensemble_member_11": 26,
        "posterior_ensemble_member_11": 27,
        "prior_ensemble_member_12": 28,
        "posterior_ensemble_member_12": 29,
    }


def get_regular_from_nc(nc_template):
    ds = xr.open_dataset(nc_template)

    lon = ds["lon"][0,:].values
    lat = ds["lat"][:,0].values

    if not (
        np.allclose(np.diff(lon), np.diff(lon)[0])
        and np.allclose(np.diff(lat), np.diff(lat)[0])
    ):
        raise ValueError("Template grid is not regular")

    lon_bins = np.linspace(lon.min(), lon.max(), lon.size + 1)
    lat_bins = np.linspace(lat.min(), lat.max(), lat.size + 1)

    return lon_bins, lat_bins


def regular_grid():
    path = (
        "/gporq3/minni/FARM-DART/RUN/data_202308_so2/OUTPUT_0/OUT/conc_g1_2023080311.nc"
    )
    ds_minni = xr.open_dataset(path)
    lon_bins = ds_minni.lon.values
    lat_bins = ds_minni.lat.values

    # dlon, dlat = 0.05, 0.05
    # breakpoint()
    # nlon = int(round((settings.lonmax - settings.lonmin) / dlon))
    # nlat  = int(round((settings.latmax - settings.latmin) / dlat))

    # print(nlon, nlat)
    # Define grid boundaries
    # lon_bins = np.linspace(settings.lonmin, settings.lonmax, nlon +1)
    # lat_bins = np.linspace(settings.latmin, settings.latmax, nlat +1 )
    # breakpoint()
    # lon_grid = ds_PRODUCT['longitude'][0,:,:].values[0,:].tolist()
    # lat_grid = ds_PRODUCT['latitude'][0,:,:].values[:,0].tolist()

    # breakpoint()
    # lon_bins = centers_to_edges(lon_grid)
    # lat_bins = centers_to_edges(lat_grid)
    # # Create a meshgrid for plotting
    # lon_centers = 0.5 * (lon_bins[:-1] + lon_bins[1:])
    # lat_centers = 0.5 * (lat_bins[:-1] + lat_bins[1:])
    # return np.meshgrid(lon_centers, lat_centers)
    return lon_bins, lat_bins


def transform_obs_seq_out(subdirectories, lon_grid, lat_grid, **kwargs):
    specific_time = kwargs.get("specific_time", None)
    print(specific_time)
    start_time = pd.to_datetime(settings.start_time, format="%Y%m%d%H")
    end_time = pd.to_datetime(settings.end_time, format="%Y%m%d%H")
    index = pd.date_range(start=start_time, end=end_time, freq="h")
    result = pd.DataFrame(
        columns=[
            "obs_s5p_tropomi",
            "prior_ensemble_mean",
            "posterior_ensemble_mean",
            "prior_ensemble_spread",
            "posterior_ensemble_spread",
            "prior_ensemble_member_1",
            "posterior_ensemble_member_1",
            "prior_ensemble_member_2",
            "posterior_ensemble_member_2",
            "prior_ensemble_member_3",
            "posterior_ensemble_member_3",
            "prior_ensemble_member_4",
            "posterior_ensemble_member_4",
            "prior_ensemble_member_5",
            "posterior_ensemble_member_5",
            "prior_ensemble_member_6",
            "posterior_ensemble_member_6",
            "prior_ensemble_member_7",
            "posterior_ensemble_member_7",
            "prior_ensemble_member_8",
            "posterior_ensemble_member_8",
            "prior_ensemble_member_9",
            "posterior_ensemble_member_9",
            "prior_ensemble_member_10",
            "posterior_ensemble_member_10",
            "prior_ensemble_member_11",
            "posterior_ensemble_member_11",
            "prior_ensemble_member_12",
            "posterior_ensemble_member_12",
        ],
        index=index,
    )
    paths = {
        "obs_s5p_tropomi": [],
        "prior_ensemble_mean": [],
        "posterior_ensemble_mean": [],
        "prior_ensemble_spread": [],
        "posterior_ensemble_spread": [],
        "prior_ensemble_member_1": [],
        "posterior_ensemble_member_1": [],
        "prior_ensemble_member_2": [],
        "posterior_ensemble_member_2": [],
        "prior_ensemble_member_3": [],
        "posterior_ensemble_member_3": [],
        "prior_ensemble_member_4": [],
        "posterior_ensemble_member_4": [],
        "prior_ensemble_member_5": [],
        "posterior_ensemble_member_5": [],
        "prior_ensemble_member_6": [],
        "posterior_ensemble_member_6": [],
        "prior_ensemble_member_7": [],
        "posterior_ensemble_member_7": [],
        "prior_ensemble_member_8": [],
        "posterior_ensemble_member_8": [],
        "prior_ensemble_member_9": [],
        "posterior_ensemble_member_9": [],
        "prior_ensemble_member_10": [],
        "posterior_ensemble_member_10": [],
        "prior_ensemble_member_11": [],
        "posterior_ensemble_member_11": [],
        "prior_ensemble_member_12": [],
        "posterior_ensemble_member_12": [],
    }

    for subdir in tqdm(subdirectories):
        if isinstance(specific_time, str):
            if subdir != specific_time:
                continue
        print(subdir)
        for key in [
            "obs_s5p_tropomi",
            "prior_ensemble_mean",
            "posterior_ensemble_mean",
            "prior_ensemble_spread",
            "posterior_ensemble_spread",
            "prior_ensemble_member_1",
            "posterior_ensemble_member_1",
            "prior_ensemble_member_2",
            "posterior_ensemble_member_2",
            "prior_ensemble_member_3",
            "posterior_ensemble_member_3",
            "prior_ensemble_member_4",
            "posterior_ensemble_member_4",
            "prior_ensemble_member_5",
            "posterior_ensemble_member_5",
            "prior_ensemble_member_6",
            "posterior_ensemble_member_6",
            "prior_ensemble_member_7",
            "posterior_ensemble_member_7",
            "prior_ensemble_member_8",
            "posterior_ensemble_member_8",
            "prior_ensemble_member_9",
            "posterior_ensemble_member_9",
            "prior_ensemble_member_10",
            "posterior_ensemble_member_10",
            "prior_ensemble_member_11",
            "posterior_ensemble_member_11",
            "prior_ensemble_member_12",
            "posterior_ensemble_member_12",
        ]:
            # if os.path.exists(os.path.join(out_dir_key1, f"{key}_{subdir}_{settings.grid_bounds['name']}.nc")):
            #    continue
            x_values = []
            y_values = []
            obs_values = []
            x_values_deg = []
            try:
                date_full = generate_file_name(str(subdir))
                print(f"obs_seq_{date_full}.final")

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
                            obs_values.extend(
                                map(float, lines[i + settings.values[key]].split())
                            )

                if x_values and y_values and obs_values:
                    from scipy.stats import binned_statistic_2d
                    from scipy.interpolate import interp1d

                    # Compute grid-averaged values
                    # breakpoint()
                    # for val in x_values:
                    #    if val > np.pi:
                    #        x_values_deg.append(np.degrees(val)-360)
                    #    else:
                    #        x_values_deg.append(np.degrees(val))
                    # x_values_deg = np.degrees(x_values) - 180

                    # Filter negatives and zeros before regridding
                    obs_values = np.array(obs_values, dtype=float)
                    x_values_deg = np.degrees(x_values)
                    x_values_deg = (
                        x_values_deg + 180
                    ) % 360 - 180  # Convert to [-180, 180] consistently
                    y_values_deg = np.degrees(y_values)

                    print(f"{key}: {np.nanmin(obs_values)}")
                    valid_mask = obs_values > 0

                    x_values_deg = x_values_deg[valid_mask]
                    y_values_deg = y_values_deg[valid_mask]
                    obs_values = obs_values[valid_mask]

                    print(f"min after mask negatives {key}: {np.nanmin(obs_values)}")
                    print(f"max after mask negatives {key}: {np.nanmax(obs_values)}")
                    print(f"{key}, {len(obs_values)}")
                    stat, x_bin, y_bin, _ = binned_statistic_2d(
                        x_values_deg,
                        y_values_deg,
                        obs_values,
                        statistic="mean",
                        bins=[lon_grid, lat_grid],
                    )

                    count, x_bin, y_bin, _ = binned_statistic_2d(
                        x_values_deg,
                        y_values_deg,
                        obs_values,
                        statistic="count",
                        bins=[lon_grid, lat_grid],
                    )
                    print(f"Number of bins with data for {key}: {np.sum(count>0)}")
                    print(f"each bin count max {key}: {np.nanmax(count)}")
                    # x_centers = (x_bin[:-1] + x_bin[1:] / 2)
                    # y_centers = (y_bin[:-1] + y_bin[1:] / 2)
                    # X, Y = np.meshgrid(x_centers, y_centers, indexing='ij')
                    # valid_mask = stat != 0

                    # if np.any(~valid_mask): # Only interpolate if there are missing values
                    #    print(key)
                    #    points = np.array([X[valid_mask], Y[valid_mask]]).T

                    #    values = stat[valid_mask] # Convert to an xarray DataArray for better handling of coordinates and rename var with key
                    # Interpolate missing values using griddata stat.T[...,np.newaxis], # Transposed to match (lat, lon) convention
                    #    stat_filled = griddata(points, values, (Y, X), method='cubic')
                    #    if np.any(np.isnan(stat_filled)):
                    #        stat_filled = griddata(points, values, (Y, X), method='nearest')
                    # else: #average variable on coords
                    #    stat_filled = stat  # No missing values, use original

                    grid_data = xr.DataArray(
                        stat.T[
                            ..., np.newaxis
                        ],  # Transposed to match (lat, lon) convention
                        name=key,
                        coords={
                            "latitude": lat_grid[:-1],
                            "longitude": lon_grid[:-1],
                            "time": pd.to_datetime([subdir], format="%Y%m%d%H"),
                        },
                        dims=["latitude", "longitude", "time"],
                    )

                    subset = grid_data.sel(
                        longitude=slice(
                            settings.grid_bounds["lonmin_grid"],
                            settings.grid_bounds["lonmax_grid"],
                        ),
                        latitude=slice(
                            settings.grid_bounds["latmin_grid"],
                            settings.grid_bounds["latmax_grid"],
                        ),
                    )

                    subset = subset.where(subset > 0)
                    mean_on_domain = subset.mean(dim=["latitude", "longitude"])
                    print(f"{key}: {mean_on_domain}")
                    # at corresponding time index set the mean_on_domain
                    result.loc[pd.to_datetime(subdir, format="%Y%m%d%H"), key] = (
                        mean_on_domain.values.item()
                    )
                    # Save to NetCDF
                    paths[key].append(
                        os.path.join(
                            out_dir, f"{key}_{subdir}_{settings.grid_bounds['name']}.nc"
                        )
                    )
                    subset.to_netcdf(
                        os.path.join(
                            out_dir, f"{key}_{subdir}_{settings.grid_bounds['name']}.nc"
                        ),
                        engine="netcdf4",
                    )
                else:
                    print(f"No data found in subdir {subdir}")
            except FileNotFoundError:
                print(f"File {subdir}/obs_seq_{str(subdir)}.final not found")
    return result, paths


def plot_observation_space(paths: dict):
    for name, list_paths in paths.items():
        for path_nc in list_paths:

            fig = plt.figure(figsize=(15, 11))

            ax = plt.axes(projection=ccrs.PlateCarree())
            plt.rcParams.update({"font.size": 24})
            dataset = xr.open_dataset(path_nc)
            # ax.set_title(f"{name} time: {dataset.time.values}", fontsize =20)
            ax.set_extent(
                [
                    settings.grid_bounds["lonmin_grid"],
                    settings.grid_bounds["lonmax_grid"],
                    settings.grid_bounds["latmin_grid"],
                    settings.grid_bounds["latmax_grid"],
                ],
                ccrs.PlateCarree(),
            )
            ax.coastlines(resolution="10m")
            ax.add_feature(cfeature.BORDERS, linewidth=0.5)
            ax.add_feature(cfeature.LAND, facecolor="white")
            # ax.add_feature(cfeature.OCEAN, facecolor="lightgrey")
            # ax.add_feature(cfeature.LAKES, facecolor="lightgrey")
            # ax.add_feature(cfeature.RIVERS, edgecolor="lightgrey")

            norm = mcolors.Normalize(
                vmin=settings.range_conc[name][0],
                vmax=settings.range_conc[name][1],
            )
            var_data = dataset.variables[name][:, :, 0].data
            # Create the plot
            img = ax.pcolormesh(
                dataset.variables["longitude"],
                dataset.variables["latitude"],
                var_data,
                alpha=0.8,
                cmap=cmap,
                vmin=settings.min_scale,
                vmax=settings.max_scale,
                shading="nearest",
                zorder=1,
            )

            # Add a colorbar
            cbar = plt.colorbar(
                img, ax=ax, orientation="horizontal", shrink=0.5, pad=0.2
            )

            x_ticks = np.arange(
                settings.grid_bounds["lonmin_grid"],
                settings.grid_bounds["lonmax_grid"],
                5,
            )
            y_ticks = np.arange(
                settings.grid_bounds["latmin_grid"],
                settings.grid_bounds["latmax_grid"],
                5,
            )
            ax.set_xticks(np.floor(x_ticks))
            ax.set_yticks(np.floor(y_ticks))
            ax.set_xticklabels(np.floor(x_ticks), rotation=90, fontsize=24)
            ax.set_yticklabels(np.floor(y_ticks), fontsize=24)
            ax.gridlines(xlocs=x_ticks, ylocs=y_ticks, linestyle="--", color="grey")
            ax.set_xlabel("Lon", fontsize=24)
            ax.set_ylabel("Lat", fontsize=24)

            # fig.tight_layout(pad=1.0)
            cbar.set_label(f"NO2 [mol/m2]", fontsize=24)
            cbar.ax.tick_params(labelsize=24)
            # Replace .nc in path_name with .png
            filename = os.path.basename(path_nc)
            path_name = filename.replace(".nc", f".png")
            plt.tight_layout()
            plt.savefig(os.path.join(out_dir_key2, filename))
            plt.close()


# cmap = mcolors.LinearSegmentedColormap.from_list("custom_map", colors)
cmap = plt.cm.jet
settings = Settings()

# for i in ['1000', '0100', '0010', '3000', '3000_noinfl', '1000_noinfl',
#        '3000_restart_infl', '3000_vert_loc']:
#for i in ["0010"]:
#    case_name = f"ITA7"
#    workdir = f"/g100_scratch/userexternal/adausili/ITA7/posteriors"
out_dir = f"/g100_scratch/userexternal/adausili/ITA7/pp_old"
#    nc_template = (
#        "/g100_scratch/userexternal/adausili/ITA7/RUN_0/out.2025113000_2025113001_ITA7.nc"
#    )
#    os.makedirs(out_dir, exist_ok=True)
#
#    subdirectories = [
#        d
#        for d in os.listdir(workdir)
#        if os.path.isdir(os.path.join(workdir, d)) and d != "plots"
#    ]
#
#    lon_grid, lat_grid = get_regular_from_nc(nc_template)
#    print(lon_grid, lat_grid)
#    try:
#        domain_averaged, paths = transform_obs_seq_out(
#            subdirectories, lon_grid, lat_grid
#        )  # kwargs : specific_time
#        domain_averaged.to_csv(
#            out_dir + f'/{case_name}_domain_avg_{settings.grid_bounds["name"]}.csv'
#        )
#
#    except:
#        print("error")

plot_monthly_mean_from_nc(out_dir, month="2025-11")
# plot_ts_domain_averaged(domain_averaged)
