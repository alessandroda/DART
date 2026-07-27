import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.transforms import Bbox
from matplotlib.patches import Rectangle
try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
except ImportError:  # cartopy is optional for non-map outputs
    ccrs = None
    cfeature = None
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional
import argparse
import re


# ==========================================
# 1. CONFIGURATION BLOCK
# ==========================================
@dataclass
class Config:
    # -- Paths and File Names --
    workdir: Path = Path("/g100_scratch/userexternal/adausili/ITA7/posteriors")
    sub_dir_name: str = "process_obs_seq_file"
    output_dir: Optional[Path] = None
    # Reference model grid (used for monthly-mean CAMS plots).
    # Extracted from: /home/dausilia/forGiorgia/out.2025110200_2025120100_ITA7.nc (lat/lon variables).
    # Kept hardcoded to avoid dependency on external files at runtime.
    reference_lon_start: float = 4.36
    reference_lon_step: float = 0.09
    reference_lon_count: int = 174
    reference_lat_start: float = 35.20
    reference_lat_step: float = 0.06
    reference_lat_count: int = 229

    reference_longitude: xr.DataArray = field(init=False, repr=False)
    reference_latitude: xr.DataArray = field(init=False, repr=False)

    experiment_name: str = field(init=False)

    # -- Spatial Domain --
    lonmin: float = 6.6
    lonmax: float = 18.8
    latmin: float = 36.6
    latmax: float = 47.1
    round_decimals: int = 4

    # -- Unit Conversion --
    # mol/m2 -> molecules/cm2
    conversion_factor: float = 6.02214e19
    csv_and_plots_units: str = "mol/m2"  # "mol/m2" or "molecules/cm2"

    # -- Ensemble Members in obs_seq_*.final --
    num_members: int = 12
    members_start_offset: int = 6  # first member line relative to " OBS" line
    debug_timestamp: Optional[str] = None  # YYYYMMDDHH
    debug_obs_index: int = 0
    debug_dump_obs: bool = False
    debug_break: bool = False
    only_timestamp: Optional[str] = None  # YYYYMMDDHH (process a single folder)

    # -- Variables to Extract from Text --
    variables_of_interest: List[str] = field(default_factory=lambda: [
        "obs_s5p_tropomi",
        "prior_ensemble_mean",
        "posterior_ensemble_mean",
        "obs_s5p_tropomi_errstd",
        "prior_ensemble_spread",
        "posterior_ensemble_spread",
    ])

    # -- Variables to Plot --
    plot_variables: List[str] = field(default_factory=lambda: [
        "obs_s5p_tropomi",
        "prior_ensemble_mean",
        "posterior_ensemble_mean",
        "mean_diff_pct",

        "obs_s5p_tropomi_errstd",
        "prior_ensemble_spread",
        "posterior_ensemble_spread",
        "spread_diff_pct",

        "prior_rmse",
        "posterior_rmse",
        "prior_totalspread",
        "posterior_totalspread",

        "prior_consistency_ratio",
        "posterior_consistency_ratio",
    ])

    # fixed offsets from line " OBS"
    values_mapping: Dict[str, int] = field(default_factory=lambda: {
        "obs_s5p_tropomi": 1,
        "prior_ensemble_mean": 2,
        "posterior_ensemble_mean": 3,
        "prior_ensemble_spread": 4,
        "posterior_ensemble_spread": 5,
    })

    # -- Plot Settings --
    cmap_name: str = "jet"
    cbar_label: str = "NO₂ [molecules/cm²]"
    map_resolution: str = "10m"

    range_conc: Dict[str, Tuple[float, float]] = field(default_factory=lambda: {
        "obs_s5p_tropomi": (0, 5e15),
        "prior_ensemble_mean": (0, 5e15),
        "posterior_ensemble_mean": (0, 5e15),

        "obs_s5p_tropomi_errstd": (5e13, 5e15),
        "prior_ensemble_spread": (5e13, 5e15),
        "posterior_ensemble_spread": (5e13, 5e15),

        "prior_rmse": (5e13, 5e15),
        "posterior_rmse": (5e13, 5e15),
        "prior_totalspread": (5e13, 5e15),
        "posterior_totalspread": (5e13, 5e15),

        "mean_diff_pct": (-50, 50),
        "spread_diff_pct": (-50, 50),

        "prior_consistency_ratio": (0.5, 2.0),
        "posterior_consistency_ratio": (0.5, 2.0),
    })

    def __post_init__(self):
        self.experiment_name = self.workdir.parent.name
        if self.output_dir is None:
            self.output_dir = self.workdir
        # Precompute reference 1D coordinate vectors (model grid)
        lon = np.round(self.reference_lon_start + self.reference_lon_step * np.arange(self.reference_lon_count), 2)
        lat = np.round(self.reference_lat_start + self.reference_lat_step * np.arange(self.reference_lat_count), 2)
        self.reference_longitude = xr.DataArray(lon, dims=("longitude",), name="longitude")
        self.reference_latitude = xr.DataArray(lat, dims=("latitude",), name="latitude")


# ==========================================
# 2. DATA PROCESSING CORE
# ==========================================
class DataProcessor:
    def __init__(self, config: Config):
        self.cfg = config

    @staticmethod
    def _generate_file_name(input_code: str) -> str:
        """Formats the date string to match the raw text file name."""
        last_two = int(input_code[-2:])
        suffix = f"0{last_two}0000" if 0 <= last_two <= 9 else f"{last_two}0000"
        return f"{input_code[:-2]}_{suffix}"

    def process_all_directories(self) -> List[Path]:
        """Scans workdir for timestamp folders, extracts data, computes diagnostics, and saves NetCDFs."""
        nc_paths = []

        subdirectories = [
            d for d in self.cfg.workdir.iterdir()
            if d.is_dir() and d.name.isdigit() and len(d.name) == 10
        ]
        if self.cfg.only_timestamp is not None:
            subdirectories = [d for d in subdirectories if d.name == self.cfg.only_timestamp]

        if not subdirectories:
            print(f"No timestamp directories found in {self.cfg.workdir}")
            return []

        for subdir in sorted(subdirectories):
            timestamp = subdir.name
            print(f"\n--- Processing Timestamp: {timestamp} ---")

            current_time = pd.to_datetime(timestamp, format="%Y%m%d%H")
            date_full = self._generate_file_name(timestamp)

            file_path = subdir / f"obs_seq_{date_full}.final"
            print(f"Looking for file: {file_path}")

            if not file_path.exists():
                print(f"Warning: The file {file_path} does not exist. Skipping this folder.")
                continue

            try:
                with open(file_path, "r") as file:
                    lines = file.readlines()
            except IOError as e:
                print(f"Error reading {file_path}: {e}")
                continue

            dataset = self._parse_file_to_dataset(lines, current_time, timestamp, file_path)

            if dataset is not None:
                out_dir_key = subdir / self.cfg.sub_dir_name
                out_dir_key.mkdir(parents=True, exist_ok=True)

                nc_file_path = out_dir_key / f"{self.cfg.experiment_name}_{timestamp}.nc"
                dataset.to_netcdf(nc_file_path)
                nc_paths.append(nc_file_path)
                print(f"Data successfully extracted and saved to: {nc_file_path}")

                self._print_domain_summary(dataset, timestamp)
            else:
                print(f"No valid data found in {file_path}.")

        return nc_paths

    def _parse_file_to_dataset(
        self,
        lines: List[str],
        current_time: pd.Timestamp,
        timestamp: str,
        file_path: Path,
    ) -> Optional[xr.Dataset]:
        """Extract observations and compute DART-style diagnostics."""
        obs_indices = [i for i, line in enumerate(lines) if line.startswith(" OBS")]
        if not obs_indices:
            return None

        x_vals, y_vals = [], []
        obs_dict = {var: [] for var in self.cfg.variables_of_interest}
        prior_member_cols = [f"prior_ens_member_{k:02d}" for k in range(1, self.cfg.num_members + 1)]
        post_member_cols = [f"posterior_ens_member_{k:02d}" for k in range(1, self.cfg.num_members + 1)]
        for col in prior_member_cols + post_member_cols:
            obs_dict[col] = []

        for idx, start_idx in enumerate(obs_indices):
            end_idx = obs_indices[idx + 1] if idx + 1 < len(obs_indices) else len(lines)

            # 1. Extract Coordinates (loc3d)
            loc_found = False
            for j in range(start_idx, end_idx):
                if lines[j].startswith("loc3d"):
                    try:
                        coords = lines[j + 1].split()
                        x_vals.append(float(coords[0]))
                        y_vals.append(float(coords[1]))
                        loc_found = True
                    except (IndexError, ValueError):
                        pass
                    break

            if not loc_found:
                continue

            # 2. Extract fixed-offset variables
            for var, offset in self.cfg.values_mapping.items():
                if start_idx + offset < end_idx:
                    try:
                        val = float(lines[start_idx + offset].split()[0]) * self.cfg.conversion_factor
                        obs_dict[var].append(val)
                    except ValueError:
                        obs_dict[var].append(np.nan)
                else:
                    obs_dict[var].append(np.nan)

            # 3. Extract observation error variance from last non-empty line
            errstd_idx = -1
            if "obs_s5p_tropomi_errstd" in self.cfg.variables_of_interest:
                errstd_idx = end_idx - 1
                while errstd_idx > start_idx and not lines[errstd_idx].strip():
                    errstd_idx -= 1

                try:
                    # obs_seq.final stores observation error variance
                    err_var = float(lines[errstd_idx].split()[0])
                    err_std = np.sqrt(err_var) * self.cfg.conversion_factor
                    obs_dict["obs_s5p_tropomi_errstd"].append(err_std)
                except (ValueError, IndexError):
                    obs_dict["obs_s5p_tropomi_errstd"].append(np.nan)

            # 4. Extract ensemble members (prior/post) starting from members_start_offset.
            # Supports two common layouts:
            # - Layout A: each member line has "prior posterior" (two floats)
            # - Layout B: first N floats are prior members, next N floats are posterior members
            prior_members: List[float] = []
            post_members: List[float] = []
            member_start = start_idx + self.cfg.members_start_offset
            if member_start < end_idx:
                floats_by_line: List[List[float]] = []
                first_floats_by_line: List[float] = []
                for j in range(member_start, end_idx):
                    parts = lines[j].split()
                    if not parts:
                        continue
                    vals: List[float] = []
                    for p in parts:
                        try:
                            vals.append(float(p))
                        except ValueError:
                            continue
                    if not vals:
                        continue
                    floats_by_line.append(vals)
                    first_floats_by_line.append(vals[0])

                    # stop early if we already have enough numbers for layout B
                    if len(first_floats_by_line) >= 2 * self.cfg.num_members:
                        break

                # Detect layout A if we have at least N lines each with >=2 floats
                layout_a = False
                if len(floats_by_line) >= self.cfg.num_members:
                    first_n = floats_by_line[: self.cfg.num_members]
                    if all(len(v) >= 2 for v in first_n):
                        layout_a = True

                if layout_a:
                    for v in floats_by_line[: self.cfg.num_members]:
                        prior_members.append(v[0] * self.cfg.conversion_factor)
                        post_members.append(v[1] * self.cfg.conversion_factor)
                else:
                    # Layout B: sequential member lines (prior block then posterior block).
                    # Use the first float from each line to avoid capturing extra numeric fields.
                    for v in first_floats_by_line[: self.cfg.num_members]:
                        prior_members.append(v * self.cfg.conversion_factor)
                    post_block = first_floats_by_line[self.cfg.num_members : 2 * self.cfg.num_members]
                    for v in post_block:
                        post_members.append(v * self.cfg.conversion_factor)

                if (
                    self.cfg.debug_timestamp is not None
                    and timestamp == self.cfg.debug_timestamp
                    and idx == self.cfg.debug_obs_index
                ):
                    if self.cfg.debug_dump_obs:
                        out_dir = Path(self.cfg.output_dir)
                        out_dir.mkdir(parents=True, exist_ok=True)
                        out_txt = out_dir / f"{self.cfg.experiment_name}_{timestamp}_obs{idx:05d}_debug.txt"
                        with open(out_txt, "w") as f:
                            f.write(f"file_path: {file_path}\n")
                            f.write(f"timestamp: {timestamp}\n")
                            f.write(f"obs_index: {idx}\n")
                            f.write(f"start_idx: {start_idx}\n")
                            f.write(f"end_idx  : {end_idx}\n")
                            f.write(f"members_start_offset: {self.cfg.members_start_offset}\n")
                            f.write(f"member_start_line   : {member_start}\n")
                            f.write(f"num_members         : {self.cfg.num_members}\n")
                            f.write(f"layout_detected     : {'A (two-col)' if layout_a else 'B (blocks)'}\n")
                            f.write("\n--- Raw OBS block lines ---\n")
                            for k in range(start_idx, end_idx):
                                f.write(f"{k:8d}: {lines[k]}")
                            f.write("\n--- Parsed members (mol/m2) ---\n")
                            pm = [v / self.cfg.conversion_factor if pd.notna(v) else np.nan for v in prior_members]
                            qm = [v / self.cfg.conversion_factor if pd.notna(v) else np.nan for v in post_members]
                            f.write(
                                "prior_members: "
                                + " ".join(f"{v:.6e}" if np.isfinite(v) else "nan" for v in pm)
                                + "\n"
                            )
                            f.write(
                                "post_members : "
                                + " ".join(f"{v:.6e}" if np.isfinite(v) else "nan" for v in qm)
                                + "\n"
                            )
                            prior_mean_file = obs_dict["prior_ensemble_mean"][-1]
                            post_mean_file = obs_dict["posterior_ensemble_mean"][-1]
                            f.write("\n--- Consistency checks (mol/m2) ---\n")
                            f.write(
                                f"prior_mean_file     : {prior_mean_file / self.cfg.conversion_factor:.6e}\n"
                            )
                            f.write(
                                f"prior_mean_members  : {np.nanmean(pm):.6e}\n"
                            )
                            f.write(
                                f"posterior_mean_file : {post_mean_file / self.cfg.conversion_factor:.6e}\n"
                            )
                            f.write(
                                f"posterior_mean_members: {np.nanmean(qm):.6e}\n"
                            )
                        print(f"Debug OBS dump written to: {out_txt}")

                    if self.cfg.debug_break:
                        print("Entering breakpoint() for member parsing inspection...")
                        breakpoint()

            if len(prior_members) < self.cfg.num_members:
                prior_members.extend([np.nan] * (self.cfg.num_members - len(prior_members)))
            if len(post_members) < self.cfg.num_members:
                post_members.extend([np.nan] * (self.cfg.num_members - len(post_members)))

            for k, col in enumerate(prior_member_cols):
                obs_dict[col].append(prior_members[k])
            for k, col in enumerate(post_member_cols):
                obs_dict[col].append(post_members[k])

            # Optional debug on first obs
            if idx in [0, 1]:
                print(f"\n  > SAFETY CHECK: OBS {idx + 1}")
                print(f"  > Raw lon (rad): {x_vals[-1]}, Raw lat (rad): {y_vals[-1]}")
                print(f"  > Converted lon (deg): {np.round(np.degrees(x_vals[-1]), self.cfg.round_decimals)}")
                print(f"  > Converted lat (deg): {np.round(np.degrees(y_vals[-1]), self.cfg.round_decimals)}")
                print("  " + "-" * 28)
                for var in self.cfg.variables_of_interest:
                    val = obs_dict[var][-1]
                    val_molm2 = val / self.cfg.conversion_factor if pd.notna(val) else np.nan
                    print(f"  > {var}: {val_molm2} [mol/m2]")
                if self.cfg.num_members > 0:
                    pm = prior_members[0] / self.cfg.conversion_factor if pd.notna(prior_members[0]) else np.nan
                    qm = post_members[0] / self.cfg.conversion_factor if pd.notna(post_members[0]) else np.nan
                    print(f"  > prior member 01: {pm} [mol/m2]")
                    print(f"  > post  member 01: {qm} [mol/m2]")
                print("  " + "-" * 28)

        if not x_vals:
            return None

        # Build dataframe in observation space
        df = pd.DataFrame({
            "longitude": np.round(np.degrees(x_vals), self.cfg.round_decimals),
            "latitude": np.round(np.degrees(y_vals), self.cfg.round_decimals),
            **obs_dict
        })

        # Remove rows with missing key values
        required_cols = [
            "obs_s5p_tropomi",
            "prior_ensemble_mean",
            "posterior_ensemble_mean",
            "obs_s5p_tropomi_errstd",
            "prior_ensemble_spread",
            "posterior_ensemble_spread",
        ]
        df = df.dropna(subset=required_cols)

        if df.empty:
            return None

        # ------------------------------------------
        # DART-style per-observation diagnostics
        # ------------------------------------------
        df["prior_sqerr"] = (df["prior_ensemble_mean"] - df["obs_s5p_tropomi"]) ** 2
        df["posterior_sqerr"] = (df["posterior_ensemble_mean"] - df["obs_s5p_tropomi"]) ** 2

        df["prior_variance"] = df["prior_ensemble_spread"] ** 2
        df["posterior_variance"] = df["posterior_ensemble_spread"] ** 2

        df["obs_error_variance"] = df["obs_s5p_tropomi_errstd"] ** 2

        df["prior_total_variance"] = df["prior_variance"] + df["obs_error_variance"]
        df["posterior_total_variance"] = df["posterior_variance"] + df["obs_error_variance"]

        # ------------------------------------------
        # Aggregate by grid cell
        # ------------------------------------------
        grouped = df.groupby(["latitude", "longitude"])

        mean_cols = [
            "obs_s5p_tropomi",
            "prior_ensemble_mean",
            "posterior_ensemble_mean",
            "obs_s5p_tropomi_errstd",
            "prior_ensemble_spread",
            "posterior_ensemble_spread",
            "prior_sqerr",
            "posterior_sqerr",
            "prior_variance",
            "posterior_variance",
            "obs_error_variance",
            "prior_total_variance",
            "posterior_total_variance",
        ]
        mean_cols.extend([c for c in prior_member_cols + post_member_cols if c in df.columns])

        grid_dataset = grouped[mean_cols].mean().to_xarray()
        grid_dataset["nobs"] = grouped.size().to_xarray()

        # ------------------------------------------
        # Final diagnostics per grid cell
        # ------------------------------------------
        grid_dataset["prior_rmse"] = np.sqrt(grid_dataset["prior_sqerr"])
        grid_dataset["posterior_rmse"] = np.sqrt(grid_dataset["posterior_sqerr"])

        grid_dataset["prior_spread_diag"] = np.sqrt(grid_dataset["prior_variance"])
        grid_dataset["posterior_spread_diag"] = np.sqrt(grid_dataset["posterior_variance"])

        grid_dataset["prior_totalspread"] = np.sqrt(grid_dataset["prior_total_variance"])
        grid_dataset["posterior_totalspread"] = np.sqrt(grid_dataset["posterior_total_variance"])

        grid_dataset["prior_consistency_ratio"] = (
            grid_dataset["prior_rmse"] / grid_dataset["prior_totalspread"]
        )
        grid_dataset["posterior_consistency_ratio"] = (
            grid_dataset["posterior_rmse"] / grid_dataset["posterior_totalspread"]
        )

        # ------------------------------------------
        # Relative changes
        # ------------------------------------------
        prior_mean = grid_dataset["prior_ensemble_mean"]
        post_mean = grid_dataset["posterior_ensemble_mean"]
        grid_dataset["mean_diff_pct"] = (
            (post_mean - prior_mean) / prior_mean.where(prior_mean != 0)
        ) * 100.0

        prior_spread = grid_dataset["prior_ensemble_spread"]
        post_spread = grid_dataset["posterior_ensemble_spread"]
        grid_dataset["spread_diff_pct"] = (
            (post_spread - prior_spread) / prior_spread.where(prior_spread != 0)
        ) * 100.0

        # Keep a time dimension (length=1) for easier concatenation downstream.
        grid_dataset = grid_dataset.expand_dims(time=[current_time])

        return grid_dataset

    def _print_domain_summary(self, dataset: xr.Dataset, timestamp: str) -> None:
        """Print domain-wide summary diagnostics."""
        def safe_mean(varname: str) -> float:
            if varname not in dataset:
                return np.nan
            return float(dataset[varname].mean(skipna=True).values)

        prior_rmse = safe_mean("prior_rmse")
        posterior_rmse = safe_mean("posterior_rmse")
        prior_totalspread = safe_mean("prior_totalspread")
        posterior_totalspread = safe_mean("posterior_totalspread")
        prior_ratio = safe_mean("prior_consistency_ratio")
        posterior_ratio = safe_mean("posterior_consistency_ratio")
        nobs = safe_mean("nobs")

        print(f"\n--- DOMAIN SUMMARY {timestamp} ---")
        print(f"mean nobs/cell               : {nobs:.2f}")
        print(f"prior rmse                   : {prior_rmse:.3e}")
        print(f"posterior rmse               : {posterior_rmse:.3e}")
        print(f"prior totalspread            : {prior_totalspread:.3e}")
        print(f"posterior totalspread        : {posterior_totalspread:.3e}")
        print(f"prior rmse/totalspread       : {prior_ratio:.3f}")
        print(f"posterior rmse/totalspread   : {posterior_ratio:.3f}")
        print("-" * 35)


# ==========================================
# 3. VISUALIZATION CORE
# ==========================================
class Visualizer:
    def __init__(self, config: Config):
        self.cfg = config

    def export_domain_timeseries_csv(self, nc_paths: List[Path]) -> Optional[Path]:
        """
        Export a domain-aggregated time series CSV.

        This keeps the existing processing logic intact by reusing `_build_domain_timeseries()`,
        which also handles unit conversion according to `cfg.csv_and_plots_units`.
        """
        df = self._build_domain_timeseries(nc_paths)
        if df is None or df.empty:
            return None

        out_path = Path(self.cfg.output_dir) / f"{self.cfg.experiment_name}_domain_timeseries.csv"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_path, index=True)
        print(f"Domain time series CSV saved to: {out_path}")
        return out_path

    def plot_cams_style_maps(
        self,
        nc_paths: List[Path],
        month: str = "2025-12",
        add_difference_panel: bool = True,
        add_monthly_mean: bool = False,
        separate_figures: bool = False,
    ) -> Optional[Path]:
        """
        Create CAMS/ECMWF-style publication-quality NO2 diagnostics maps.

        Notes:
        - Produces one figure per orbit/time (per NetCDF) within the requested month.
        - Additionally writes/overwrites a month-level filename
          "{experiment}_CAMS_style_{month}.png" pointing to the latest processed orbit,
          to preserve backward compatibility with earlier expectations.
        """
        if ccrs is None or cfeature is None:
            print("Cartopy not available; skipping CAMS-style maps.")
            return None
        if not nc_paths:
            return None

        try:
            month_period = pd.Period(month, freq="M")
        except Exception:
            print(f"Invalid month '{month}' (expected YYYY-MM).")
            return None

        def _get_cmap(preferred: str, fallback: str):
            try:
                return plt.get_cmap(preferred)
            except Exception:
                return plt.get_cmap(fallback)

        # CAMS-like palette (light cyan -> yellow -> orange -> brown), discretized for a
        # clean publication-style colorbar reminiscent of CAMS/ECMWF products.
        cams_colors = [
            "#cfeff7",  # very light cyan
            "#b9e4f0",
            "#9fd6e8",
            "#7fc1d6",
            "#bfe3d3",  # pale greenish
            "#eef3b6",  # light yellow-green
            "#fee68c",  # light yellow
            "#f7cd5a",  # yellow/orange
            "#f0ad34",  # orange
            "#e07a1c",  # deep orange
            "#8c2d16",  # brown
        ]
        cmap_no2 = mcolors.ListedColormap(cams_colors, name="cams_no2")
        cmap_diff = plt.get_cmap("RdBu_r")

        no2_vmin, no2_vmax = 0.0, 10e15
        increment_vmax = 1.5e15
        increment_vmin, increment_vmax = -increment_vmax, increment_vmax
        no2_bounds = np.linspace(no2_vmin, no2_vmax, len(cams_colors) + 1)
        no2_norm = mcolors.BoundaryNorm(no2_bounds, cmap_no2.N, clip=False)
        increment_min_abs = 2e14

        titles_no2 = [
            "Satellite observations",
            "Model forecast",
            "Assimilated analysis",
        ]

        def _plot_single_panel(
            *,
            lon: xr.DataArray,
            lat: xr.DataArray,
            field: xr.DataArray,
            title: str,
            out_path: Path,
            cmap,
            norm=None,
            vmin=None,
            vmax=None,
            cbar_label: Optional[str] = None,
            extend: str = "both",
            ticks: Optional[List[float]] = None,
            figsize: Tuple[float, float] = (6.0, 5.0),
        ) -> Path:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
            fig.patch.set_facecolor("white")
            ax.set_facecolor("white")
            ax.set_extent([self.cfg.lonmin, self.cfg.lonmax, self.cfg.latmin, self.cfg.latmax], crs=ccrs.PlateCarree())
            ax.add_feature(cfeature.LAND, facecolor="#f2f2f2", edgecolor="none", zorder=0)
            ax.add_feature(cfeature.OCEAN, facecolor="white", edgecolor="none", zorder=0)
            ax.coastlines(resolution="10m", linewidth=0.8, color="black")
            ax.add_feature(cfeature.BORDERS, linewidth=0.6, edgecolor="black")
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
            if hasattr(ax, "outline_patch"):
                ax.outline_patch.set_visible(False)

            mesh = ax.pcolormesh(
                lon,
                lat,
                field,
                cmap=cmap,
                norm=norm,
                vmin=vmin,
                vmax=vmax,
                shading="auto",
            )
            ax.set_title(title, fontsize=12, pad=8)

            # Compact centered colorbar
            cax = fig.add_axes([0.18, 0.10, 0.64, 0.045])
            cbar = fig.colorbar(mesh, cax=cax, orientation="horizontal", extend=extend, ticks=ticks)
            cbar.outline.set_edgecolor("black")
            cbar.outline.set_linewidth(0.6)
            cbar.ax.tick_params(labelsize=9, length=2, width=0.6, colors="black")
            if cbar_label:
                cbar.set_label(cbar_label, fontsize=11, color="black")

            out_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
            plt.close(fig)
            return out_path

        last_orbit_path: Optional[Path] = None
        month_obs: List[xr.DataArray] = []
        month_prior: List[xr.DataArray] = []
        month_post: List[xr.DataArray] = []
        ref_lon: Optional[xr.DataArray] = self.cfg.reference_longitude if add_monthly_mean else None
        ref_lat: Optional[xr.DataArray] = self.cfg.reference_latitude if add_monthly_mean else None

        for path_nc in sorted(nc_paths):
            try:
                with xr.open_dataset(path_nc) as ds:
                    time_val = self._extract_time(ds, path_nc)
                    if time_val is None:
                        continue
                    ts = pd.Timestamp(time_val)
                    if ts.to_period("M") != month_period:
                        continue

                    required = ["obs_s5p_tropomi", "prior_ensemble_mean", "posterior_ensemble_mean"]
                    if any(v not in ds.variables for v in required):
                        print(f"Missing required variables in {path_nc.name}; skipping CAMS-style map.")
                        continue

                    def _as_2d(da: xr.DataArray) -> xr.DataArray:
                        return da.isel(time=0) if "time" in da.dims else da

                    obs = _as_2d(ds["obs_s5p_tropomi"])
                    prior = _as_2d(ds["prior_ensemble_mean"])
                    post = _as_2d(ds["posterior_ensemble_mean"])

                    if add_monthly_mean:
                        if ref_lon is None and "longitude" in ds:
                            ref_lon = ds["longitude"]
                        if ref_lat is None and "latitude" in ds:
                            ref_lat = ds["latitude"]

                        def _to_ref_grid(da: xr.DataArray) -> xr.DataArray:
                            if ref_lon is None or ref_lat is None:
                                return da
                            if "longitude" not in da.coords or "latitude" not in da.coords:
                                return da
                            try:
                                lon = da["longitude"]
                                lat = da["latitude"]
                                # If the grid already matches, avoid any resampling.
                                if lon.shape == ref_lon.shape and lat.shape == ref_lat.shape:
                                    if np.array_equal(lon.values, ref_lon.values) and np.array_equal(lat.values, ref_lat.values):
                                        return da
                                # Most common case here is 1D lat/lon coordinates.
                                if lon.ndim == 1 and lat.ndim == 1 and ref_lon.ndim == 1 and ref_lat.ndim == 1:
                                    return da.interp(longitude=ref_lon, latitude=ref_lat, method="nearest")
                                # Fallback: nearest reindex (works for 1D coords).
                                return da.reindex({"longitude": ref_lon, "latitude": ref_lat}, method="nearest")
                            except Exception:
                                # Last resort: force coords to reference if shapes match (prevents union-grid concat).
                                try:
                                    if da["longitude"].shape == ref_lon.shape and da["latitude"].shape == ref_lat.shape:
                                        return da.assign_coords(longitude=ref_lon, latitude=ref_lat)
                                except Exception:
                                    pass
                                return da

                        month_obs.append(_to_ref_grid(obs))
                        month_prior.append(_to_ref_grid(prior))
                        month_post.append(_to_ref_grid(post))

                    increment = (post - prior).where(np.abs(post - prior) > increment_min_abs)

                    if separate_figures:
                        orbit_tag = ts.strftime("%Y%m%d%H")
                        base = Path(self.cfg.output_dir) / f"{self.cfg.experiment_name}_CAMS_style_{month}_{orbit_tag}"
                        out_obs = _plot_single_panel(
                            lon=ds["longitude"],
                            lat=ds["latitude"],
                            field=obs,
                            title="Observations",
                            out_path=base.with_name(base.name + "_obs.png"),
                            cmap=cmap_no2,
                            norm=no2_norm,
                            cbar_label="NO$_2$ column [molecules/cm$^2$]  (×10$^{15}$)",
                            ticks=[0, 1e15, 2e15, 3e15, 4e15, 5e15],
                        )
                        _plot_single_panel(
                            lon=ds["longitude"],
                            lat=ds["latitude"],
                            field=prior,
                            title="Forecast",
                            out_path=base.with_name(base.name + "_forecast.png"),
                            cmap=cmap_no2,
                            norm=no2_norm,
                            cbar_label="NO$_2$ column [molecules/cm$^2$]  (×10$^{15}$)",
                            ticks=[0, 1e15, 2e15, 3e15, 4e15, 5e15],
                        )
                        _plot_single_panel(
                            lon=ds["longitude"],
                            lat=ds["latitude"],
                            field=post,
                            title="Analysis",
                            out_path=base.with_name(base.name + "_analysis.png"),
                            cmap=cmap_no2,
                            norm=no2_norm,
                            cbar_label="NO$_2$ column [molecules/cm$^2$]  (×10$^{15}$)",
                            ticks=[0, 1e15, 2e15, 3e15, 4e15, 5e15],
                        )
                        if add_difference_panel:
                            _plot_single_panel(
                                lon=ds["longitude"],
                                lat=ds["latitude"],
                                field=increment,
                                title="Assimilation increment",
                                out_path=base.with_name(base.name + "_rel_diff.png"),
                                cmap=cmap_diff,
                                vmin=increment_vmin,
                                vmax=increment_vmax,
                                cbar_label="Analysis - Forecast [molecules/cm$^2$]",
                                ticks=[-1.5e15, -1.0e15, -5e14, 0.0, 5e14, 1.0e15, 1.5e15],
                            )

                        # Keep backward compatible month-level filename pointing to latest "obs" panel.
                        month_level = Path(self.cfg.output_dir) / f"{self.cfg.experiment_name}_CAMS_style_{month}.png"
                        try:
                            import shutil
                            shutil.copyfile(out_obs, month_level)
                        except Exception:
                            pass

                        last_orbit_path = out_obs
                        continue

                    if add_difference_panel:
                        ncols = 4
                        figsize = (20, 5)
                        titles = titles_no2 + ["Assimilation increment"]
                        data_list = [obs, prior, post, increment]
                    else:
                        ncols = 3
                        figsize = (16, 5)
                        titles = titles_no2
                        data_list = [obs, prior, post]

                    # GridSpec gives better control over map/cbar spacing than tight_layout.
                    fig = plt.figure(figsize=figsize)
                    gs = fig.add_gridspec(
                        2,
                        ncols,
                        height_ratios=[1.0, 0.10],
                        hspace=0.05,
                        wspace=0.02,
                    )
                    axs = np.array([fig.add_subplot(gs[0, i], projection=ccrs.PlateCarree()) for i in range(ncols)])
                    fig.patch.set_facecolor("white")

                    mesh_no2 = None
                    mesh_diff = None

                    for ax, data, title in zip(axs, data_list, titles):
                        ax.set_facecolor("white")
                        ax.set_extent(
                            [self.cfg.lonmin, self.cfg.lonmax, self.cfg.latmin, self.cfg.latmax],
                            crs=ccrs.PlateCarree(),
                        )
                        ax.add_feature(cfeature.LAND, facecolor="#f2f2f2", edgecolor="none", zorder=0)
                        ax.add_feature(cfeature.OCEAN, facecolor="white", edgecolor="none", zorder=0)
                        ax.coastlines(resolution="10m", linewidth=0.8, color="black")
                        ax.add_feature(cfeature.BORDERS, linewidth=0.6, edgecolor="black")
                        ax.set_xticks([])
                        ax.set_yticks([])
                        for spine in ax.spines.values():
                            spine.set_visible(False)
                        if hasattr(ax, "outline_patch"):
                            ax.outline_patch.set_visible(False)

                        if title == "Assimilation increment":
                            mesh_diff = ax.pcolormesh(
                                ds["longitude"],
                                ds["latitude"],
                                data,
                                cmap=cmap_diff,
                                vmin=increment_vmin,
                                vmax=increment_vmax,
                                shading="auto",
                            )
                        else:
                            mesh_no2 = ax.pcolormesh(
                                ds["longitude"],
                                ds["latitude"],
                                data,
                                cmap=cmap_no2,
                                norm=no2_norm,
                                shading="auto",
                            )

                        ax.set_title(title, fontsize=12, pad=6)

                    if mesh_no2 is None:
                        plt.close(fig)
                        continue

                    # Centered, compact colorbars, matching size across bars.
                    fig.canvas.draw()
                    if add_difference_panel:
                        bbox_no2 = Bbox.union([axs[0].get_position(), axs[1].get_position(), axs[2].get_position()])
                        bbox_diff = axs[3].get_position()
                        cbar_width = min(bbox_no2.width, bbox_diff.width)
                        cbar_height = 0.035
                        cbar_pad = 0.020

                        no2_left = (bbox_no2.x0 + bbox_no2.x1) / 2.0 - cbar_width / 2.0
                        diff_left = (bbox_diff.x0 + bbox_diff.x1) / 2.0 - cbar_width / 2.0
                        cbar_bottom = min(bbox_no2.y0, bbox_diff.y0) - cbar_pad - cbar_height

                        no2_cax = fig.add_axes([no2_left, cbar_bottom, cbar_width, cbar_height])
                        no2_ticks = np.linspace(no2_vmin, no2_vmax, 6)
                        cbar = fig.colorbar(
                            mesh_no2,
                            cax=no2_cax,
                            orientation="horizontal",
                            extend="both",
                            boundaries=no2_bounds,
                            ticks=no2_ticks,
                            spacing="proportional",
                        )
                        cbar.outline.set_edgecolor("black")
                        cbar.outline.set_linewidth(0.6)
                        cbar.ax.tick_params(labelsize=9, length=2, width=0.6, colors="black")
                        cbar.ax.set_xticklabels([f"{v/1e15:.0f}" for v in no2_ticks])
                        cbar.set_label("NO$_2$ column [molecules/cm$^2$]  (×10$^{15}$)", fontsize=11, color="black")

                        if mesh_diff is not None:
                            diff_cax = fig.add_axes([diff_left, cbar_bottom, cbar_width, cbar_height])
                            cbar_diff = fig.colorbar(
                                mesh_diff,
                                cax=diff_cax,
                                orientation="horizontal",
                                extend="both",
                                ticks=[-40, -20, 0, 20, 40],
                            )
                            cbar_diff.outline.set_edgecolor("black")
                            cbar_diff.outline.set_linewidth(0.6)
                            cbar_diff.ax.tick_params(labelsize=9, length=2, width=0.6, colors="black")
                            cbar_diff.set_label("Analysis - Forecast [molecules/cm$^2$]", fontsize=11, color="black")
                    else:
                        bbox_no2 = Bbox.union([ax.get_position() for ax in axs])
                        cbar_width = bbox_no2.width * 0.55
                        cbar_height = 0.035
                        cbar_pad = 0.020
                        no2_left = (bbox_no2.x0 + bbox_no2.x1) / 2.0 - cbar_width / 2.0
                        cbar_bottom = bbox_no2.y0 - cbar_pad - cbar_height
                        no2_cax = fig.add_axes([no2_left, cbar_bottom, cbar_width, cbar_height])
                        no2_ticks = np.linspace(no2_vmin, no2_vmax, 6)
                        cbar = fig.colorbar(
                            mesh_no2,
                            cax=no2_cax,
                            orientation="horizontal",
                            extend="both",
                            boundaries=no2_bounds,
                            ticks=no2_ticks,
                            spacing="proportional",
                        )
                        cbar.outline.set_edgecolor("black")
                        cbar.outline.set_linewidth(0.6)
                        cbar.ax.tick_params(labelsize=9, length=2, width=0.6, colors="black")
                        cbar.ax.set_xticklabels([f"{v/1e15:.0f}" for v in no2_ticks])
                        cbar.set_label("NO$_2$ column [molecules/cm$^2$]  (×10$^{15}$)", fontsize=11, color="black")

                    # Title format requested
                    fig.suptitle(
                        f"Time: {ts:%a %d %b %Y %H} UTC. Variable: NO$_2$ Column. "
                        f"Area: Italy, Experiment {self.cfg.experiment_name}",
                        fontsize=12,
                        y=0.98,
                    )

                    # Clean margins
                    fig.subplots_adjust(left=0.01, right=0.99, top=0.90, bottom=0.22)

                    orbit_tag = ts.strftime("%Y%m%d%H")
                    out_orbit = Path(self.cfg.output_dir) / f"{self.cfg.experiment_name}_CAMS_style_{month}_{orbit_tag}.png"
                    fig.savefig(out_orbit, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
                    plt.close(fig)

                    month_level = Path(self.cfg.output_dir) / f"{self.cfg.experiment_name}_CAMS_style_{month}.png"
                    try:
                        import shutil
                        shutil.copyfile(out_orbit, month_level)
                    except Exception:
                        pass

                    last_orbit_path = out_orbit
            except Exception as e:
                print(f"Failed CAMS-style map for {path_nc.name}: {e}")

        if last_orbit_path is None:
            print(f"No CAMS-style maps produced for month {month}.")
            return None

        if add_monthly_mean and month_obs and month_prior and month_post:
            try:
                obs_mean = xr.concat(month_obs, dim="orbit").mean(dim="orbit", skipna=True)
                prior_mean = xr.concat(month_prior, dim="orbit").mean(dim="orbit", skipna=True)
                post_mean = xr.concat(month_post, dim="orbit").mean(dim="orbit", skipna=True)
                increment_mean = (post_mean - prior_mean).where(np.abs(post_mean - prior_mean) > increment_min_abs)

                if separate_figures:
                    base = Path(self.cfg.output_dir) / f"{self.cfg.experiment_name}_CAMS_style_{month}_monthly_mean"
                    _plot_single_panel(
                        lon=obs_mean["longitude"],
                        lat=obs_mean["latitude"],
                        field=obs_mean,
                        title="Monthly mean Observations",
                        out_path=base.with_name(base.name + "_obs.png"),
                        cmap=cmap_no2,
                        norm=no2_norm,
                        cbar_label="NO$_2$ column [molecules/cm$^2$]  (×10$^{15}$)",
                        ticks=[0, 1e15, 2e15, 3e15, 4e15, 5e15],
                        figsize=(6.5, 5.5),
                    )
                    _plot_single_panel(
                        lon=prior_mean["longitude"],
                        lat=prior_mean["latitude"],
                        field=prior_mean,
                        title="Monthly mean Forecast",
                        out_path=base.with_name(base.name + "_forecast.png"),
                        cmap=cmap_no2,
                        norm=no2_norm,
                        cbar_label="NO$_2$ column [molecules/cm$^2$]  (×10$^{15}$)",
                        ticks=[0, 1e15, 2e15, 3e15, 4e15, 5e15],
                        figsize=(6.5, 5.5),
                    )
                    _plot_single_panel(
                        lon=post_mean["longitude"],
                        lat=post_mean["latitude"],
                        field=post_mean,
                        title="Monthly mean Analysis",
                        out_path=base.with_name(base.name + "_analysis.png"),
                        cmap=cmap_no2,
                        norm=no2_norm,
                        cbar_label="NO$_2$ column [molecules/cm$^2$]  (×10$^{15}$)",
                        ticks=[0, 1e15, 2e15, 3e15, 4e15, 5e15],
                        figsize=(6.5, 5.5),
                    )
                    if add_difference_panel:
                        _plot_single_panel(
                            lon=increment_mean["longitude"],
                            lat=increment_mean["latitude"],
                            field=increment_mean,
                            title="Monthly mean Assimilation increment",
                            out_path=base.with_name(base.name + "_rel_diff.png"),
                            cmap=cmap_diff,
                            vmin=increment_vmin,
                            vmax=increment_vmax,
                            cbar_label="Analysis - Forecast [molecules/cm$^2$]",
                            ticks=[-1.5e15, -1.0e15, -5e14, 0.0, 5e14, 1.0e15, 1.5e15],
                            figsize=(6.5, 5.5),
                        )
                    return last_orbit_path

                # Monthly mean: use 2x2 layout for better use of space
                if add_difference_panel:
                    nrows, ncols = 2, 2
                    figsize = (14, 7)
                    titles = [
                        "Observations",
                        "Forecast",
                        "Analysis",
                        "Assimilation increment",
                    ]
                    data_list = [obs_mean, prior_mean, post_mean, increment_mean]
                else:
                    nrows, ncols = 2, 2
                    figsize = (14, 7)
                    titles = [
                        "Observations",
                        "Forecast",
                        "Analysis",
                        "",
                    ]
                    data_list = [obs_mean, prior_mean, post_mean, None]

                fig = plt.figure(figsize=figsize)
                gs = fig.add_gridspec(nrows, ncols, hspace=0.08, wspace=0.03)
                axs = np.empty((nrows, ncols), dtype=object)
                for r in range(nrows):
                    for c in range(ncols):
                        axs[r, c] = fig.add_subplot(gs[r, c], projection=ccrs.PlateCarree())
                fig.patch.set_facecolor("white")

                mesh_no2 = None
                mesh_diff = None
                flat_axes = [axs[0, 0], axs[0, 1], axs[1, 0], axs[1, 1]]
                for ax, data, title in zip(flat_axes, data_list, titles):
                    if data is None:
                        ax.set_visible(False)
                        continue
                    ax.set_facecolor("white")
                    ax.set_extent(
                        [self.cfg.lonmin, self.cfg.lonmax, self.cfg.latmin, self.cfg.latmax],
                        crs=ccrs.PlateCarree(),
                    )
                    ax.add_feature(cfeature.LAND, facecolor="#f2f2f2", edgecolor="none", zorder=0)
                    ax.add_feature(cfeature.OCEAN, facecolor="white", edgecolor="none", zorder=0)
                    ax.coastlines(resolution="10m", linewidth=0.8, color="black")
                    ax.add_feature(cfeature.BORDERS, linewidth=0.6, edgecolor="black")
                    ax.set_xticks([])
                    ax.set_yticks([])
                    for spine in ax.spines.values():
                        spine.set_visible(False)
                    if hasattr(ax, "outline_patch"):
                        ax.outline_patch.set_visible(False)

                    if title == "Assimilation increment":
                        mesh_diff = ax.pcolormesh(
                            obs_mean["longitude"],
                            obs_mean["latitude"],
                            data,
                            cmap=cmap_diff,
                            vmin=increment_vmin,
                            vmax=increment_vmax,
                            shading="auto",
                        )
                    else:
                        mesh_no2 = ax.pcolormesh(
                            obs_mean["longitude"],
                            obs_mean["latitude"],
                            data,
                            cmap=cmap_no2,
                            norm=no2_norm,
                            shading="auto",
                        )
                    ax.set_title(title, fontsize=12, pad=6)

                if mesh_no2 is not None:
                    fig.canvas.draw()
                    if add_difference_panel:
                        bbox_no2 = Bbox.union([axs[0, 0].get_position(), axs[0, 1].get_position(), axs[1, 0].get_position()])
                        bbox_diff = axs[1, 1].get_position()
                        cbar_width = min(bbox_no2.width, bbox_diff.width)
                        cbar_height = 0.035
                        cbar_pad = 0.020
                        no2_left = (bbox_no2.x0 + bbox_no2.x1) / 2.0 - cbar_width / 2.0
                        diff_left = (bbox_diff.x0 + bbox_diff.x1) / 2.0 - cbar_width / 2.0
                        cbar_bottom = min(bbox_no2.y0, bbox_diff.y0) - cbar_pad - cbar_height

                        no2_cax = fig.add_axes([no2_left, cbar_bottom, cbar_width, cbar_height])
                        no2_ticks = np.linspace(no2_vmin, no2_vmax, 6)
                        cbar = fig.colorbar(
                            mesh_no2,
                            cax=no2_cax,
                            orientation="horizontal",
                            extend="both",
                            boundaries=no2_bounds,
                            ticks=no2_ticks,
                            spacing="proportional",
                        )
                        cbar.outline.set_edgecolor("black")
                        cbar.outline.set_linewidth(0.6)
                        cbar.ax.tick_params(labelsize=9, length=2, width=0.6, colors="black")
                        cbar.ax.set_xticklabels([f"{v/1e15:.0f}" for v in no2_ticks])
                        cbar.set_label("NO$_2$ column [molecules/cm$^2$]  (×10$^{15}$)", fontsize=11, color="black")

                        if mesh_diff is not None:
                            diff_cax = fig.add_axes([diff_left, cbar_bottom, cbar_width, cbar_height])
                            cbar_diff = fig.colorbar(
                                mesh_diff,
                                cax=diff_cax,
                                orientation="horizontal",
                                extend="both",
                                ticks=[-1.5e15, -1.0e15, -5e14, 0.0, 5e14, 1.0e15, 1.5e15],
                            )
                            cbar_diff.outline.set_edgecolor("black")
                            cbar_diff.outline.set_linewidth(0.6)
                            cbar_diff.ax.tick_params(labelsize=9, length=2, width=0.6, colors="black")
                            cbar_diff.set_label("Analysis - Forecast [molecules/cm$^2$]", fontsize=11, color="black")
                    else:
                        bbox_no2 = Bbox.union([ax.get_position() for ax in flat_axes if ax.get_visible()])
                        cbar_width = bbox_no2.width * 0.55
                        cbar_height = 0.035
                        cbar_pad = 0.020
                        no2_left = (bbox_no2.x0 + bbox_no2.x1) / 2.0 - cbar_width / 2.0
                        cbar_bottom = bbox_no2.y0 - cbar_pad - cbar_height
                        no2_cax = fig.add_axes([no2_left, cbar_bottom, cbar_width, cbar_height])
                        no2_ticks = np.linspace(no2_vmin, no2_vmax, 6)
                        cbar = fig.colorbar(
                            mesh_no2,
                            cax=no2_cax,
                            orientation="horizontal",
                            extend="both",
                            boundaries=no2_bounds,
                            ticks=no2_ticks,
                            spacing="proportional",
                        )
                        cbar.outline.set_edgecolor("black")
                        cbar.outline.set_linewidth(0.6)
                        cbar.ax.tick_params(labelsize=9, length=2, width=0.6, colors="black")
                        cbar.ax.set_xticklabels([f"{v/1e15:.0f}" for v in no2_ticks])
                        cbar.set_label("NO$_2$ column [molecules/cm$^2$]  (×10$^{15}$)", fontsize=11, color="black")

                    fig.suptitle(
                        f"Monthly mean ({month_period}) | Variable: NO$_2$ Column. "
                        f"Area: Italy, Experiment {self.cfg.experiment_name}",
                        fontsize=12,
                        y=0.98,
                    )
                    fig.subplots_adjust(left=0.03, right=0.97, top=0.90, bottom=0.24)

                    out_monthly = Path(self.cfg.output_dir) / f"{self.cfg.experiment_name}_CAMS_style_{month}_monthly_mean.png"
                    fig.savefig(out_monthly, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
                    plt.close(fig)
                    print(f"CAMS-style monthly mean saved to: {out_monthly}")
            except Exception as e:
                print(f"Failed CAMS-style monthly mean for month {month}: {e}")

        print(f"CAMS-style maps saved (last): {last_orbit_path}")
        return last_orbit_path


    def plot_domain_timeseries_pwp(self, nc_paths: List[Path]) -> Optional[Path]:
        df = self._build_domain_timeseries(nc_paths)
        if df is None or df.empty:
            return None

        daily = df.sort_index().resample("D").mean(numeric_only=True)
        if daily.empty:
            return None

        units = self.cfg.csv_and_plots_units
        ylabel = f"NO$_2$ total column [{units}]"

        # --- COLORS (match Sentinel style) ---
        bg_color = "#0c1c44"
        grid_color = "white"
        obs_color = "#a6ff00"   # neon green
        prior_color = "#4cc9f0" # cyan
        post_color = "#ff6b6b"  # soft red

        fig = plt.figure(figsize=(11, 5), facecolor=bg_color)
        ax = fig.add_subplot(111)
        ax.set_facecolor(bg_color)

        # --- LINES ---
        ax.plot(daily.index, daily["obs_s5p_tropomi"],
                marker="o", markersize=3, linewidth=1.8,
                color=obs_color,
                label="Satellite observations (Sentinel-5P)")

        ax.plot(daily.index, daily["prior_ens_mean"],
                marker="o", markersize=3, linewidth=1.5,
                color=prior_color,
                label="Model forecast (background)")

        ax.plot(daily.index, daily["posterior_ens_mean"],
                marker="o", markersize=3, linewidth=1.8,
                color=post_color,
                label="Assimilated analysis")

        # --- AXES STYLE ---
        ax.set_title("Daily-average NO$_2$ total column over Italy",
                    fontsize=14, color="white", pad=10)

        ax.set_ylabel(ylabel, color="white")
        ax.set_xlabel("Time (UTC)", color="white")

        ax.tick_params(colors="white")
        ax.grid(True, linestyle="--", alpha=0.2, color=grid_color)

        # remove spines
        for spine in ax.spines.values():
            spine.set_visible(False)

        # --- LEGEND ---
        leg = ax.legend(frameon=False, fontsize=9)
        for text in leg.get_texts():
            text.set_color("white")

        fig.autofmt_xdate()

        # # =====================================
        # # 📍 INSET MAP (kAIROS domain)
        # # =====================================
        # if ccrs is not None:
        #     # Make room on the right
        #     fig.subplots_adjust(right=0.78)

        #     # Add inset OUTSIDE main plot
        #     axins = fig.add_axes([0.80, 0.55, 0.18, 0.30], projection=ccrs.PlateCarree())

        #     axins.coastlines(resolution="10m", color="white", linewidth=0.5)
        #     axins.add_feature(cfeature.BORDERS, linewidth=0.3, edgecolor="white")

        #     # wider Europe view
        #     axins.set_extent([-10, 30, 30, 55], crs=ccrs.PlateCarree())

        #     # draw domain box
        #     axins.plot(
        #         [self.cfg.lonmin, self.cfg.lonmax, self.cfg.lonmax, self.cfg.lonmin, self.cfg.lonmin],
        #         [self.cfg.latmin, self.cfg.latmin, self.cfg.latmax, self.cfg.latmax, self.cfg.latmin],
        #         color=obs_color,
        #         linewidth=1.5
        #     )

        #     axins.set_title("kAIROS domain", color="white", fontsize=8)

        #     axins.tick_params(labelsize=6, colors="white")

        # # =====================================

        fig.tight_layout()

        out_path = Path(self.cfg.output_dir) / f"{self.cfg.experiment_name}_PWP_timeseries.png"
        fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
        plt.close(fig)

        print(f"PWP plot saved to: {out_path}")
        return out_path

    def plot_domain_timeseries(self, nc_paths: List[Path]) -> Optional[Path]:
        """Create a daily-mean time series plot from domain-aggregated values."""
        df = self._build_domain_timeseries(nc_paths)
        if df is None or df.empty:
            return None

        daily = df.sort_index().resample("D").mean(numeric_only=True)
        if daily.empty:
            return None

        units = self.cfg.csv_and_plots_units
        ylabel = f"Vertical Column Density NO₂ [{units}]"

        fig, ax = plt.subplots(figsize=(10, 4.5))
        ax.plot(daily.index, daily["obs_s5p_tropomi"], "o-", color="black", linewidth=1.5, markersize=3, label="Observations")
        ax.plot(daily.index, daily["prior_ens_mean"], "o-", color="tab:blue", linewidth=1.5, markersize=3, label="Background")
        ax.plot(daily.index, daily["posterior_ens_mean"], "o-", color="tab:red", linewidth=1.5, markersize=3, label="Analysis")

        ax.set_title("Daily-average Total Column NO₂")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        ax.legend(frameon=True)
        fig.autofmt_xdate()
        fig.tight_layout()

        out_path = Path(self.cfg.output_dir) / f"{self.cfg.experiment_name}_daily_timeseries.png"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"Domain time series plot saved to: {out_path}")
        return out_path

    def plot_domain_distributions(self, nc_paths: List[Path]) -> Optional[Path]:
        """Create an aggregated log-x histogram for obs/prior/posterior over all grid cells and times."""
        arrays = self._collect_domain_distributions(nc_paths)
        if arrays is None:
            return None
        obs, prior, post, time_range_label = arrays
        if obs.size == 0 or prior.size == 0 or post.size == 0:
            return None

        units = self.cfg.csv_and_plots_units
        xlabel = f"NO₂ [{units}]"

        positive = np.concatenate([obs, prior, post])
        positive = positive[np.isfinite(positive) & (positive > 0)]
        if positive.size == 0:
            return None

        vmin = float(np.nanmin(positive))
        vmax = float(np.nanmax(positive))
        bins = np.logspace(np.log10(vmin), np.log10(vmax), 45)

        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.hist(obs, bins=bins, alpha=0.5, color="gray", label="Observation",
                weights=np.ones_like(obs) * 100.0 / max(obs.size, 1))
        ax.hist(prior, bins=bins, alpha=0.5, color="tab:blue", label="Background",
                weights=np.ones_like(prior) * 100.0 / max(prior.size, 1))
        ax.hist(post, bins=bins, alpha=0.5, color="tab:red", label="Analysis",
                weights=np.ones_like(post) * 100.0 / max(post.size, 1))

        ax.set_xscale("log")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Frequency (%)")
        title = f"Total Column NO₂ Distributions, {time_range_label}, {self.cfg.experiment_name}"
        ax.set_title(title)
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(frameon=True)
        fig.tight_layout()

        out_path = Path(self.cfg.output_dir) / f"{self.cfg.experiment_name}_distributions.png"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"Domain distributions plot saved to: {out_path}")
        return out_path

    def plot_domain_ensemble_members(self, nc_paths: List[Path]) -> Optional[Path]:
        """
        Plot domain-mean ensemble members (prior/post) as spaghetti lines, plus
        the prior/post ensemble means and observations.
        """
        df_members = self._build_domain_members_timeseries(nc_paths)
        if df_members is None or df_members.empty:
            print("No ensemble member variables found; skipping member plot.")
            return None

        daily = df_members.sort_index().resample("D").mean(numeric_only=True)
        if daily.empty:
            return None

        units = self.cfg.csv_and_plots_units
        ylabel = f"Total Column NO₂ [{units}]"

        prior_member_cols = [c for c in daily.columns if c.startswith("prior_ens_member_")]
        post_member_cols = [c for c in daily.columns if c.startswith("posterior_ens_member_")]

        fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(11, 7), sharex=True)

        # Prior panel
        for c in prior_member_cols:
            ax0.plot(daily.index, daily[c], color="tab:blue", alpha=0.25, linewidth=0.9)
        ax0.plot(daily.index, daily["prior_ens_mean"], color="tab:blue", linewidth=2.0, label="Background mean")
        ax0.plot(daily.index, daily["obs_s5p_tropomi"], "o-", color="black", linewidth=1.2, markersize=2.5, label="Observations")
        ax0.set_title("Daily-average Total Column NO₂ (Background members)")
        ax0.set_ylabel(ylabel)
        ax0.grid(True, alpha=0.3)
        ax0.legend(frameon=True)

        # Posterior panel
        for c in post_member_cols:
            ax1.plot(daily.index, daily[c], color="tab:red", alpha=0.25, linewidth=0.9)
        ax1.plot(daily.index, daily["posterior_ens_mean"], color="tab:red", linewidth=2.0, label="Analysis mean")
        ax1.plot(daily.index, daily["obs_s5p_tropomi"], "o-", color="black", linewidth=1.2, markersize=2.5, label="Observations")
        ax1.set_title("Daily-average Total Column NO₂ (Analysis members)")
        ax1.set_ylabel(ylabel)
        ax1.grid(True, alpha=0.3)
        ax1.legend(frameon=True)

        fig.autofmt_xdate()
        fig.tight_layout()

        out_path = Path(self.cfg.output_dir) / f"{self.cfg.experiment_name}_daily_members.png"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"Ensemble members plot saved to: {out_path}")
        return out_path

    def plot_maps(self, nc_paths: List[Path]) -> None:
        """Generates subplot spatial maps from unified NetCDF files."""
        if ccrs is None or cfeature is None:
            print("Cartopy not available; skipping spatial maps.")
            return
        if not nc_paths:
            return

        base_cmap = plt.get_cmap(self.cfg.cmap_name)
        diff_cmap = plt.get_cmap("RdBu_r")
        ratio_cmap = plt.get_cmap("viridis")

        log_scale_vars = [
            "obs_s5p_tropomi_errstd",
            "prior_ensemble_spread",
            "posterior_ensemble_spread",
            "prior_rmse",
            "posterior_rmse",
            "prior_totalspread",
            "posterior_totalspread",
        ]

        for path_nc in nc_paths:
            print(f"\nGenerating map for: {path_nc.name}")
            self._create_subplot_map(path_nc, base_cmap, diff_cmap, ratio_cmap, log_scale_vars)

    def _create_subplot_map(self, path_nc: Path, base_cmap, diff_cmap, ratio_cmap, log_scale_vars: List[str]) -> None:
        """Helper to create and save a grid figure containing subplots for all variables."""
        try:
            with xr.open_dataset(path_nc) as dataset:
                if 'time' in dataset.coords:
                    time_values = dataset.time.values
                    time_val = str(time_values[0] if np.size(time_values) else time_values)[:13].replace('T', ' ') + ':00'
                else:
                    time_val = path_nc.name

                vars_to_plot = [v for v in self.cfg.plot_variables if v in dataset.variables]

                num_vars = len(vars_to_plot)
                if num_vars == 0:
                    return

                cols = 4
                rows = (num_vars + cols - 1) // cols

                fig, axs = plt.subplots(
                    rows, cols,
                    figsize=(7 * cols, 6 * rows),
                    subplot_kw={"projection": ccrs.PlateCarree()}
                )

                fig.patch.set_facecolor("white")
                axs = np.atleast_1d(axs).flatten()

                fig.suptitle(f"{self.cfg.experiment_name} | Time: {time_val}", fontsize=20, y=1.02)

                for idx, var_name in enumerate(vars_to_plot):
                    ax = axs[idx]
                    vmin, vmax = self.cfg.range_conc.get(var_name, (0, 1))
                    da = dataset[var_name]
                    if "time" in da.dims:
                        da = da.isel(time=0)

                    ax.set_title(var_name, fontsize=14, pad=10)
                    ax.set_extent(
                        [self.cfg.lonmin, self.cfg.lonmax, self.cfg.latmin, self.cfg.latmax],
                        ccrs.PlateCarree()
                    )

                    ax.coastlines(resolution=self.cfg.map_resolution)
                    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
                    ax.add_feature(cfeature.LAND, facecolor="white")
                    ax.add_feature(cfeature.OCEAN, facecolor="white")

                    is_diff = "diff_pct" in var_name
                    is_ratio = "ratio" in var_name

                    if is_diff:
                        img = ax.pcolormesh(
                            dataset["longitude"], dataset["latitude"], da,
                            alpha=0.8, cmap=diff_cmap, vmin=vmin, vmax=vmax,
                            shading="auto", zorder=1
                        )
                        cbar_label_used = "Relative Diff [%]"

                    elif is_ratio:
                        img = ax.pcolormesh(
                            dataset["longitude"], dataset["latitude"], da,
                            alpha=0.8, cmap=ratio_cmap, vmin=vmin, vmax=vmax,
                            shading="auto", zorder=1
                        )
                        cbar_label_used = "RMSE / TOTALSPREAD [-]"

                    elif var_name in log_scale_vars:
                        if vmin <= 0:
                            vmin = 5e13
                        norm = mcolors.LogNorm(vmin=vmin, vmax=vmax)
                        img = ax.pcolormesh(
                            dataset["longitude"], dataset["latitude"], da,
                            alpha=0.8, cmap=base_cmap, norm=norm,
                            shading="auto", zorder=1
                        )
                        cbar_label_used = self.cfg.cbar_label

                    else:
                        img = ax.pcolormesh(
                            dataset["longitude"], dataset["latitude"], da,
                            alpha=0.8, cmap=base_cmap, vmin=vmin, vmax=vmax,
                            shading="auto", zorder=1
                        )
                        cbar_label_used = self.cfg.cbar_label

                    cbar = plt.colorbar(
                        img, ax=ax, orientation="horizontal",
                        pad=0.15, fraction=0.05, extend="both"
                    )
                    cbar.set_label(cbar_label_used, fontsize=12)

                    x_ticks = np.linspace(self.cfg.lonmin, self.cfg.lonmax, 7)
                    y_ticks = np.linspace(self.cfg.latmin, self.cfg.latmax, 7)
                    ax.set_xticks(np.floor(x_ticks))
                    ax.set_yticks(np.floor(y_ticks))
                    ax.set_xticklabels(np.floor(x_ticks), rotation=45, fontsize=10)
                    ax.set_yticklabels(np.floor(y_ticks), fontsize=10)
                    ax.set_xlabel("Longitude", fontsize=12)
                    ax.set_ylabel("Latitude", fontsize=12)

                for ax in axs[num_vars:]:
                    ax.axis("off")

                fig.tight_layout(pad=2.0)
                out_path = path_nc.with_suffix(".png")
                plt.savefig(out_path, dpi=200, bbox_inches="tight", facecolor=fig.get_facecolor())
                print(f"Map successfully saved to: {out_path}")

        except Exception as e:
            print(f"Failed to plot {path_nc}: {e}")
        finally:
            plt.close("all")

    def _build_domain_timeseries(self, nc_paths: List[Path]) -> Optional[pd.DataFrame]:
        if not nc_paths:
            return None

        records = []
        for path_nc in sorted(nc_paths):
            try:
                with xr.open_dataset(path_nc) as ds:
                    time_val = self._extract_time(ds, path_nc)
                    if time_val is None:
                        continue

                    obs = self._domain_weighted_mean(ds, "obs_s5p_tropomi")
                    prior_mean = self._domain_weighted_mean(ds, "prior_ensemble_mean")
                    post_mean = self._domain_weighted_mean(ds, "posterior_ensemble_mean")
                    prior_spread = self._domain_weighted_mean(ds, "prior_ensemble_spread")
                    post_spread = self._domain_weighted_mean(ds, "posterior_ensemble_spread")

                records.append({
                    "time": time_val,
                    "obs_s5p_tropomi": obs,
                    "prior_ens_mean": prior_mean,
                    "posterior_ens_mean": post_mean,
                    "prior_ens_spread": prior_spread,
                    "posterior_ens_spread": post_spread,
                })
            except Exception as e:
                print(f"Failed to read {path_nc}: {e}")

        if not records:
            return None

        df = pd.DataFrame.from_records(records).set_index("time").sort_index()

        if self.cfg.csv_and_plots_units == "mol/m2":
            for col in df.columns:
                df[col] = df[col] / self.cfg.conversion_factor

        return df

    def _build_domain_members_timeseries(self, nc_paths: List[Path]) -> Optional[pd.DataFrame]:
        if not nc_paths:
            return None

        records = []
        member_name_re = re.compile(r"^(prior|posterior)_ens_member_(\d{2})$")

        for path_nc in sorted(nc_paths):
            try:
                with xr.open_dataset(path_nc) as ds:
                    time_val = self._extract_time(ds, path_nc)
                    if time_val is None:
                        continue

                    base = {
                        "time": time_val,
                        "obs_s5p_tropomi": self._domain_weighted_mean(ds, "obs_s5p_tropomi"),
                        "prior_ens_mean": self._domain_weighted_mean(ds, "prior_ensemble_mean"),
                        "posterior_ens_mean": self._domain_weighted_mean(ds, "posterior_ensemble_mean"),
                    }

                    member_vars = [v for v in ds.variables if member_name_re.match(v)]
                    if not member_vars:
                        continue

                    for v in sorted(member_vars):
                        base[v] = self._domain_weighted_mean(ds, v)

                records.append(base)
            except Exception as e:
                print(f"Failed to read {path_nc}: {e}")

        if not records:
            return None

        df = pd.DataFrame.from_records(records).set_index("time").sort_index()

        if self.cfg.csv_and_plots_units == "mol/m2":
            for col in df.columns:
                df[col] = df[col] / self.cfg.conversion_factor

        return df

    def _collect_domain_distributions(
        self, nc_paths: List[Path]
    ) -> Optional[Tuple[np.ndarray, np.ndarray, np.ndarray, str]]:
        if not nc_paths:
            return None

        obs_all: List[np.ndarray] = []
        prior_all: List[np.ndarray] = []
        post_all: List[np.ndarray] = []
        times: List[pd.Timestamp] = []

        for path_nc in sorted(nc_paths):
            try:
                with xr.open_dataset(path_nc) as ds:
                    time_val = self._extract_time(ds, path_nc)
                    if time_val is not None:
                        times.append(time_val)

                    obs = ds.get("obs_s5p_tropomi")
                    prior = ds.get("prior_ensemble_mean")
                    post = ds.get("posterior_ensemble_mean")
                    if obs is None or prior is None or post is None:
                        continue

                    obs_all.append(np.asarray(obs.values).ravel())
                    prior_all.append(np.asarray(prior.values).ravel())
                    post_all.append(np.asarray(post.values).ravel())
            except Exception as e:
                print(f"Failed to read {path_nc}: {e}")

        if not obs_all or not prior_all or not post_all:
            return None

        obs = np.concatenate(obs_all)
        prior = np.concatenate(prior_all)
        post = np.concatenate(post_all)

        mask = np.isfinite(obs) & np.isfinite(prior) & np.isfinite(post)
        obs = obs[mask]
        prior = prior[mask]
        post = post[mask]

        if self.cfg.csv_and_plots_units == "mol/m2":
            obs = obs / self.cfg.conversion_factor
            prior = prior / self.cfg.conversion_factor
            post = post / self.cfg.conversion_factor

        label = "all-times"
        if times:
            tmin = min(times)
            tmax = max(times)
            if (tmin.year, tmin.month) == (tmax.year, tmax.month):
                label = f"{tmin.month:02d}/{tmin.year}"
            else:
                label = f"{tmin.date()}–{tmax.date()}"

        return obs, prior, post, label

    @staticmethod
    def _extract_time(ds: xr.Dataset, path_nc: Path) -> Optional[pd.Timestamp]:
        if "time" in ds.coords:
            tvals = ds["time"].values
            if np.size(tvals) == 0:
                return None
            t0 = tvals[0] if np.size(tvals) > 0 else tvals
            return pd.to_datetime(t0)

        # fallback: parse from filename "<exp>_<YYYYMMDDHH>.nc"
        stem = path_nc.stem
        parts = stem.split("_")
        for part in reversed(parts):
            if part.isdigit() and len(part) == 10:
                return pd.to_datetime(part, format="%Y%m%d%H")
        return None

    @staticmethod
    def _domain_weighted_mean(ds: xr.Dataset, varname: str) -> float:
        if varname not in ds.variables:
            return float("nan")
        da = ds[varname]
        if "time" in da.dims:
            da = da.isel(time=0)
        spatial_dims = list(da.dims)
        if "nobs" in ds.variables and spatial_dims:
            weights = ds["nobs"]
            if "time" in weights.dims:
                weights = weights.isel(time=0)
            try:
                w = weights.where(np.isfinite(da))
                val = da.weighted(w).mean(dim=spatial_dims, skipna=True).values
                return float(val)
            except Exception:
                pass

        # fallback: unweighted mean
        return float(da.mean(dim=spatial_dims, skipna=True).values)


# ==========================================
# 4. MAIN PIPELINE EXECUTION
# ==========================================
def main():
    parser = argparse.ArgumentParser(description="Post-process DART obs_seq files into gridded NetCDF, CSV, and plots.")
    parser.add_argument("--workdir", type=Path, default=Config.workdir, help="Directory containing timestamp subfolders (YYYYMMDDHH).")
    parser.add_argument("--sub-dir-name", type=str, default=Config.sub_dir_name, help="Name of per-timestamp output subdirectory.")
    parser.add_argument("--output-dir", type=Path, default=None, help="Directory for CSV/plot outputs (defaults to workdir).")
    parser.add_argument("--skip-maps", action="store_true", help="Skip spatial map plots (Cartopy required).")
    parser.add_argument("--skip-diagnostic-dashboard", action="store_true", help="Skip the large multi-panel diagnostic dashboard maps.")
    parser.add_argument("--skip-csv", action="store_true", help="Skip writing the domain time series CSV.")
    parser.add_argument("--skip-ts-plot", action="store_true", help="Skip the domain daily time series plot.")
    parser.add_argument("--skip-dist-plot", action="store_true", help="Skip the domain distributions plot.")
    parser.add_argument("--skip-member-plot", action="store_true", help="Skip the domain ensemble members time series plot.")
    parser.add_argument("--csv-plots-units", choices=["mol/m2", "molecules/cm2"], default=Config.csv_and_plots_units, help="Units for CSV/plots.")
    parser.add_argument("--only-timestamp", type=str, default=None, help="Process only one timestamp folder (YYYYMMDDHH).")
    parser.add_argument("--debug-timestamp", type=str, default=None, help="Enable debug for a specific timestamp folder (YYYYMMDDHH).")
    parser.add_argument("--debug-obs-index", type=int, default=0, help="Observation index within the .final file to debug (0-based).")
    parser.add_argument("--debug-dump-obs", action="store_true", help="Write a debug text dump of the selected OBS block.")
    parser.add_argument("--debug-break", action="store_true", help="Call breakpoint() at the selected OBS block.")
    args = parser.parse_args()

    cfg = Config(
        workdir=args.workdir,
        sub_dir_name=args.sub_dir_name,
        output_dir=args.output_dir,
        csv_and_plots_units=args.csv_plots_units,
        debug_timestamp=args.debug_timestamp,
        debug_obs_index=args.debug_obs_index,
        debug_dump_obs=args.debug_dump_obs,
        debug_break=args.debug_break,
        only_timestamp=args.only_timestamp,
    )

    print("--- Starting Post-Processing Pipeline ---")

    processor = DataProcessor(cfg)
    nc_paths = processor.process_all_directories()

    visualizer = Visualizer(cfg)
    if not args.skip_csv:
        visualizer.export_domain_timeseries_csv(nc_paths)
    if not args.skip_ts_plot:
        visualizer.plot_domain_timeseries_pwp(nc_paths)
    if not args.skip_dist_plot:
        visualizer.plot_domain_distributions(nc_paths)
    if not args.skip_member_plot:
        visualizer.plot_domain_ensemble_members(nc_paths)
    if not args.skip_maps and not args.skip_diagnostic_dashboard:
        visualizer.plot_maps(nc_paths)

    if not args.skip_maps:
        visualizer.plot_cams_style_maps(
            nc_paths,
            month="2025-12",
            add_difference_panel=True,
            add_monthly_mean=True,
            separate_figures=True,
        )

    print("\n--- Pipeline Completed Successfully! ---")


if __name__ == "__main__":
    main()
