from pathlib import Path
from matplotlib import animation, colors
from matplotlib import pyplot as plt
import pandas as pd
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
import warnings
import contextily as ctx
import os
from pydantic import BaseModel
import matplotlib.cm as cm
import imageio.v2 as imageio

from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap

warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=FutureWarning)


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
cmap = LinearSegmentedColormap.from_list("custom_map", colors)


blue_to_purple_20 = [
    (0.0, 0.6, 0.7607843137254902, 1),
    (0.0, 0.6, 0.7607843137254902, 1),
    (0.01568627450980392, 0.8666666666666667, 0.9058823529411765, 1),
    (0.13333333333333333, 0.9529411764705882, 0.8196078431372549, 1),
    (0.34509803921568627, 0.8862745098039215, 0.5372549019607843, 1),
    (0.5411764705882353, 0.8313725490196079, 0.2901960784313726, 1),
    (0.6980392156862745, 0.8509803921568627, 0.14901960784313725, 1),
    (0.8235294117647058, 0.9137254901960784, 0.08627450980392157, 1),
    (0.9254901960784314, 0.9176470588235294, 0.03529411764705882, 1),
    (0.9725490196078431, 0.792156862745098, 0.011764705882352941, 1),
    (0.9882352941176471, 0.5882352941176471, 0.00784313725490196, 1),
    (0.996078431372549, 0.403921568627451, 0.0, 1),
    (1.0, 0.25098039215686274, 0.0, 1),
    (1.0, 0.10980392156862745, 0.0, 1),
    (1.0, 0.0196078431372549, 0.12549019607843137, 1),
    (1.0, 0.0, 0.4235294117647059, 1),
]
dict_layers = {
    0: 20,
    1: 65,
    2: 125,
    3: 210,
    4: 325,
    5: 480,
    6: 690,
    7: 975,
    8: 1360,
    9: 1880,
    10: 2580,
    11: 3525,
    12: 4805,
    13: 6290,
    14: 8040,
    15: 9790,
    16: 11790,
}

n_colors = len(blue_to_purple_20)
jet_colors_with_alpha = plt.cm.jet(np.linspace(0, 1, n_colors))

# Add transparency (alpha) values to the colors
alpha_values = [color[3] for color in blue_to_purple_20]
jet_colors_with_alpha[:, 3] = alpha_values

# Convert RGBA values to tuple format and create a dictionary
jet_colors_dict_with_alpha = [color for color in jet_colors_with_alpha]
jet_colors_dict_with_alpha[0] = (1.0, 1.0, 1.0, 1.0)


def find_grid_cell(dataset, x, y):
    # Find the nearest grid cell
    nearest_x = dataset["x"].sel(x=x, method="nearest").values
    nearest_y = dataset["y"].sel(y=y, method="nearest").values

    # Find the indices of the nearest grid cell
    x_index = np.where(dataset["x"].values == nearest_x)[0][0]
    y_index = np.where(dataset["y"].values == nearest_y)[0][0]

    return x_index, y_index


import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score


def r2(observed, predicted):
    predicted = predicted[~np.isnan(observed)]
    observed = observed[~np.isnan(observed)]
    # observed = observed[~np.isnan(predicted)]
    # predicted = predicted[~np.isnan(predicted)]
    cov = np.cov(observed, predicted)
    return (cov[0, 1] / np.sqrt((cov[0, 0] * cov[1, 1]))) ** 2


def make_obs_vs_data(settings, dataset_variable, dataset, station_points, **kwargs):

    name_process = kwargs.get("name", "")
    # 1. Query cells where the station belongs to
    import pandas as pd

    # Filter out negative values, NaNs, and values higher than 200

    values_model = pd.DataFrame(columns=["c_NO2"])
    # Format output string date
    formatted_str = str(dataset.time.values)
    print("processing: " + formatted_str)
    filtered_time_points = station_points[
        station_points.index == np.datetime64(dataset.time.values[0])
    ]
    filtered_time_points["model"] = 0
    for lon, lat in zip(filtered_time_points["lon"], filtered_time_points["lat"]):
        x_index, y_index = find_grid_cell(dataset, lon, lat)
        model_value = dataset_variable[0, y_index, x_index].data
        filtered_time_points.loc[
            (filtered_time_points["lat"] == lat) & (filtered_time_points["lon"] == lon),
            "model",
        ] = model_value
        # Append model_value to values_model DataFrame
        values_model = pd.concat(
            [values_model, pd.DataFrame({"c_NO2": [model_value]})],
            ignore_index=True,
        )
    # 2. Plot obs vs data
    # Remove NaN values

    values_obs = filtered_time_points["value"].values

    values_obs_filtered = values_obs[~(values_obs == -9999.0)]
    values_model_filtered = values_model[~(values_obs == -9999.0)]
    plt.figure(figsize=(8, 6))
    plt.scatter(values_obs, values_model, color="blue", label="Observed vs. Model Data")
    plt.xlabel("Observed NO2")
    plt.ylabel("Model NO2")
    plt.xlim(right=50, left=0)
    plt.ylim(bottom=0, top=50)
    plt.title("Observed vs. Model NO2")
    plt.legend()
    plt.grid(True)
    plt.savefig(
        Path(settings.output_directory)
        / f"{name_process}_{str(dataset.time.values[0])}.png"
    )
    # r2 coefficient without NaN
    print(r2_score(values_obs_filtered, values_model_filtered))
    print(
        r2(
            values_obs_filtered,
            np.array(values_model_filtered["c_NO2"].values.tolist()),
        )
    )


class Settings(BaseModel):
    # NETCDF
    # ds_year = ds.sortby('time').resample(time="Y").mean()
    word = "preassim"
    path_netcdf: str = (
        f"/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/IRIDE_ens/{word}/"
    )
    file_name = f"diff_member_0001.nc"
    output_directory: str = (
        "/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/IRIDE_ens/3D_plots/"
    )
    type_statistics = ""  # options
    model: str = f"{word}"
    engine: str = "netcdf4"
    # list 16 integers from 1 to 16
    layers = list(range(0, 1))
    conv_factor: int = 1
    dest_crs: int = 4326
    source_crs: int = 4326
    name_variable: str = "c_NO2"
    obs_variable: str = "value"
    start_time: str = "2021-01-06T12:00:00"
    end_time: str = "2021-01-08T12:00:00"
    # Domain Italy
    lonmin: float = 5.4  # 5.325
    lonmax: float = 19
    latmin: float = 35.380
    latmax: float = 47.8
    # # Domain North
    # lonmin: float = 8  # 5.325
    # lonmax: float = 14
    # latmin: float = 43
    # latmax: float = 47.8
    # Map style
    label_size: int = 18
    title_font_size: int = 20
    size_figure: tuple = (15, 12)
    min_scale: int = -20
    max_scale: int = 20
    title: str = "diff preassim/analysis"
    var_units: str = "NO2 ugm3"  # "NO. exceed."

    make_animation: bool = False
    clean_pngs: bool = False
    frame_per_second: int = 1
    type_colorbar: str = "YOB_transparent"  # YOB_transparent, Viridis_transparent
    power_exp: int = 4  # 0 -> no transparency inf -> full transparent
    basemap_provider: str = (
        "CartoDB.Positron"  # OpenStreetMap.HOT, MapTilesAPI.OSMEnglish, OpenTopoMap, OpenRailwayMap, JusticeMap.black others from https://github.com/geopandas/contextily/blob/main/notebooks/providers_deepdive.ipynb
    )
    logo_path: str = "/mnt/mumbai_n4r5/dausilio/ncmakemapspy/arianet_logo.png"


def get_basemap_provider(basemap_provider_name: str):
    if basemap_provider_name == "CartoDB.Positron":
        return ctx.providers.CartoDB.Positron
    elif basemap_provider_name == "OpenStreetMap.HOT":
        return ctx.providers.OpenStreetMap.HOT
    elif basemap_provider_name == "MapTilesAPI.OSMEnglish":
        return ctx.providers.MapTilesAPI.OSMEnglish
    elif basemap_provider_name == "OpenTopoMap":
        return ctx.providers.MapTilesAPI.OpenTopoMap
    elif basemap_provider_name == "OpenStreetMap.BZH":
        return ctx.providers.MapTilesAPI.OpenStreetMap.BZH


def get_my_cmap(type_colorbar: str):
    if type_colorbar is None:
        return
    if type_colorbar == "YOB_transparent":
        cm_values = np.linspace(0, 1, 16536)
        # use the Yellow Orange Brown color map as reference
        alpha_cm = plt.cm.jet(cm_values)
        # Change alpha value to follow a square root law
        alpha_cm[:, -1] = np.power(cm_values, 4)
        # build new color map
        return colors.ListedColormap(alpha_cm)
    elif type_colorbar == "Viridis_transparent":
        # Fill this branch using the "viridis" colormap
        cm_values = np.linspace(0, 1, 16536)
        # Use the "viridis" colormap as reference
        alpha_cm = plt.cm.viridis(cm_values)
        # Change alpha value to follow a square root law
        alpha_cm[:, -1] = np.power(cm_values, 0.3)
        # Build a new color map
        return colors.ListedColormap(alpha_cm)
    else:
        return plt.cm.jet


def make_single_pngs(settings: Settings, **kwargs):
    station_points = kwargs.get("station_points", None)
    shape_feature = kwargs.get("state_feature", None)
    filenames = []
    dirlist = []
    for dir in os.listdir(settings.path_netcdf):
        if dir.endswith("00"):
            dirlist.append(dir)
    for dirname in dirlist:
        # check folder not empty
        if not os.listdir(settings.path_netcdf + dirname):
            continue
        time = pd.to_datetime(dirname, format="%Y%m%d_%H%M%S")
        # continue the loop if dataset empty
        if not os.path.exists(
            settings.path_netcdf + dirname + "/" + settings.file_name
        ):
            continue
        dataset = open_dataset(
            settings.path_netcdf + dirname + "/" + settings.file_name, settings
        )
        dataset_variable = dataset[settings.name_variable]
        # make_obs_vs_data(
        #     settings=settings,
        #     dataset=dataset,
        #     dataset_variable=dataset_variable,
        #     station_points=station_points,
        #     name=settings.word,
        # )
        dataset = dataset.sel(time=time, method="nearest")
        for layer_i in settings.layers:
            # if dataset.time[i].values >= np.datetime64(
            #     settings.start_time
            # ) and dataset.time[i].values <= np.datetime64(settings.end_time):
            # Extract the data for the current time step
            var_data = dataset_variable[layer_i, :, :].data

            # Format output string date
            formatted_str = str(dataset.time.dt.strftime("%Y-%m-%d_%H:%M:%S").data)
            print("processing: " + formatted_str)

            fig = plt.figure(figsize=settings.size_figure)
            ax = plt.axes(projection=ccrs.PlateCarree())
            ax.set_extent(
                [
                    settings.lonmin,
                    settings.lonmax,
                    settings.latmin,
                    settings.latmax,
                ],
                ccrs.PlateCarree(),
            )
            x_ticks = np.linspace(
                settings.lonmin,
                settings.lonmax,
                10,
            )
            y_ticks = np.linspace(
                settings.latmin,
                settings.latmax,
                10,
            )
            ax.set_xticks(x_ticks)
            ax.set_yticks(y_ticks)

            plt.xticks(rotation=90, fontsize=18)
            plt.yticks(fontsize=18)
            plt.xlabel("Lon [deg]", fontsize=18)
            plt.ylabel("Lon [deg]", fontsize=18)
            ax.ticklabel_format(style="plain")
            # plot pcolormesh using dataset['x'], dataset['y'] and var_data
            # cmap = plt.cm.jet
            # cmap = get_my_cmap(settings.type_colorbar)
            # cmap = plt.cm.YlOrBr
            plt.title(
                f"{settings.word} {str(dict_layers[layer_i])} m {formatted_str}",
                fontdict={
                    "fontsize": settings.title_font_size,
                    "fontweight": "bold",
                },
            )

            img = plt.pcolormesh(
                dataset["x"],
                dataset["y"],
                var_data,
                alpha=0.75,
                vmin=settings.min_scale,
                vmax=settings.max_scale,
                cmap=cmap,
                shading="auto",
                zorder=1,
            )
            # img = plt.pcolormesh(
            #     dataset["x"],
            #     dataset["y"],
            #     var_data.transpose(),
            #     alpha=0.75,
            #     cmap=,
            #     transform=ccrs.PlateCarree(),
            #     vmin=settings.min_scale,
            #     vmax=settings.max_scale,
            #     shading="auto",
            #     zorder=1,
            # )
            # ax.gridlines(xlocs=x_ticks, ylocs=y_ticks, linestyle="--", color="black")
            plt.grid(zorder=7)
            # ctx.add_basemap(
            #     ax, crs=ccrs.UTM(32), zoom=1, source=get_basemap_provider(settings.basemap_provider), alpha=1
            # )
            # graph_proj = ox.project_graph(graph,to_crs = 'EPSG:4326', to_latlong=True)
            # ox.plot_graph(graph_proj, ax=ax, node_size=0, edge_color='black', edge_linewidth=0.1)
            # logo = plt.imread(settings.logo_path)
            # logo_height = 0.5 # Adjust the height of the logo as needed
            # logo_width = logo_height * logo.shape[1] / logo.shape[0]
            # ax.imshow(logo, extent=(settings.latmax - logo_width,settings.latmax, settings.lonmax - logo_height,settings.lonmax), transform=ccrs.PlateCarree())
            if station_points is not None:
                filtered_time_points = station_points[
                    station_points.index == np.datetime64(dataset.time.data)
                ]
                plt.scatter(
                    filtered_time_points["lon"],
                    filtered_time_points["lat"],
                    c=filtered_time_points[settings.obs_variable],
                    s=150,
                    transform=ccrs.PlateCarree(),
                    cmap=cmap,
                    vmin=settings.min_scale,
                    vmax=settings.max_scale,
                    zorder=7,
                    edgecolors="black",
                )

                # for index, row in filtered_time_points.iterrows():
                #     txt = int(row["value"])
                #     offset = 0  # Adjust this value as needed
                #     plt.text(
                #         row["lon"] - offset,  # Adjust x position
                #         row["lat"] - offset,  # Adjust y position
                #         txt,
                #         ha="center",  # Adjust horizontal alignment
                #         va="center",  # Adjust vertical alignment
                #         fontsize=12,
                #         zorder=7,
                #         color="black",
                #     )
            if shape_feature is not None:
                ax.add_feature(
                    shape_feature, facecolor=(1, 1, 1, 0), linewidth=0.5, zorder=4
                )

            # Format colorbar
            cbar = fig.colorbar(
                img, ax=ax, orientation="vertical", fraction=0.04, pad=0.02
            )
            cbar.set_label(
                f"{settings.title}  [{settings.var_units}]",
                fontsize=settings.label_size,
            )
            # number ticks to 10
            # cbar.set_ticks(np.linspace(0, 50, 10))
            cbar.ax.tick_params(labelsize=settings.label_size)

            # Save the image
            image_filename = os.path.join(
                settings.output_directory,
                f"{settings.title}_{formatted_str}_{layer_i}_{settings.model}_north.png",
            )
            filenames.append(image_filename)
            # add cfeature roads
            # ax.add_feature(cfeature.OCEAN.with_scale("10m"), zorder=2)
            ax.add_feature(cfeature.LAKES.with_scale("10m"), zorder=2)
            ax.add_feature(cfeature.BORDERS.with_scale("10m"), zorder=2)
            ax.add_feature(cfeature.RIVERS.with_scale("10m"), zorder=2)
            # ax.add_feature(cfeature.LAND.with_scale("10m"), zorder=2)
            ax.add_feature(cfeature.COASTLINE.with_scale("10m"), zorder=2)

            plt.savefig(image_filename)
            plt.close()
    return filenames


def open_dataset(path_netcdf, settings):
    dataset = xr.open_dataset(path_netcdf, engine=settings.engine)
    dataset = dataset.rio.write_crs(settings.source_crs)
    dataset = dataset.assign_coords(
        x=(dataset.x.values * 0.15) + 5.4,
        y=(dataset.y.values * 0.104) + 35.43,
        # x=dataset.x.values,
        # y=dataset.y.values,
    )
    if settings.dest_crs != settings.source_crs:
        dataset = dataset.rio.reproject(settings.dest_crs)
    return dataset


def main():
    settings = Settings()
    # subset_variable = dataset[settings.name_variable]
    output_directory = settings.output_directory
    os.makedirs(output_directory, exist_ok=True)
    import pandas as pd

    # station_points = pd.read_csv(
    #     "/mnt/mumbai_n4r5/dausilio/projects/eea_app/csv_table/stations_eea_2021-01-01000000_2021-01-01000000_IT_HR_CH_SI_AT_SK_NO2.csv"
    # )
    # station_points["datetime"] = pd.to_datetime(station_points["datetime"])
    # station_points = station_points[
    #     (station_points["type"] == "background")
    #     & (
    #         (station_points["type_area"] == "rural")
    #         | (station_points["type_area"] == "suburban")
    #     )
    # ]
    # station_points.set_index("datetime", inplace=True)

    filenames = make_single_pngs(
        settings=settings,
        # station_points=station_points,
    )

    if settings.make_animation:
        # use filnames to make pngs
        for filename in filenames:
            image = imageio.imread(filename)
            imageio.imwrite(filename, image)
        gif_filename = os.path.join(
            settings.output_directory,
            f"{settings.title}_{settings.model}_animation.gif",
        )
        with imageio.get_writer(gif_filename, fps=settings.frame_per_second) as writer:
            for filename in filenames:
                image = imageio.imread(filename)
                writer.append_data(image)
    if settings.clean_pngs:
        filelist = [f for f in os.listdir(output_directory) if f.endswith(".png")]
        for f in filelist:
            os.remove(os.path.join(output_directory, f))


if __name__ == "__main__":
    main()
