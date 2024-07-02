# Author: Alessandro D'Ausilio
# email: alessandro.dausilio@suez.com


# Python version based from
# https://github.com/apmizzi/DART_Chem/blob/main/apm_run_scripts/RUN_PERT_CHEM/EMISS_PERT/perturb_chem_emiss_CORR_RT_MA.f90
# Evensen, G. (2003). "The Ensemble Kalman Filter: theoretical formulation and practical implementation." Ocean Dynamics, 53(4), 343-367.
# An ensemble assessment of regional ozone model uncertainty with an explicit error representation. Boynard 2021

# Organize imports, remove duplicates and unused ones
import os
import numpy as np
import xarray as xr
from pathlib import Path
from tqdm import tqdm
from typing import Tuple
from pydantic import BaseSettings  # Changed from BaseConfig to BaseSettings
from glob import glob

# Define constants
EARTH_RADIUS_KM = 6371.0


# Settings class using Pydantic BaseSettings
class Settings(BaseSettings):
    path_emissions: Path = Path("/home/dausilia/emission_check/")
    name_netcdfs: str = "HERMESv3_*.nc"
    dim_to_groud_ncs: str = "time"
    tau: int = 6  # Time decorrelation scale in hours
    alpha: float = np.exp(-1 / tau)  # Smoothing coefficient
    var: str = "veNO2"
    mems: int = 1
    corr_length_hz: float = 30000
    coee_lenght_vz: float = 100
    spread: float = 0.8
    corr_time: int = 6

    class Config:
        env_prefix = "EMISSION_"  # Allow overriding settings with environment variables


settings = Settings()


def get_vertical_correlation_matrix(
    nx: int, ny: int, nz: int, corr_length_vt: float
) -> np.ndarray:
    """
    Generate vertical correlation matrix.

    Args:
        nx (int): Number of x grid points
        ny (int): Number of y grid points
        nz (int): Number of z grid points
        corr_length_vt (float): Vertical correlation length

    Returns:
        np.ndarray: Vertical correlation matrix
    """
    A = np.zeros((nx, ny, nz, nz))
    for k in tqdm(range(nz), unit="first_z"):
        for l in range(nz):
            vcov = np.exp(-abs(k - l) / corr_length_vt)
            A[:, :, k, l] = calculate_A_element(k, l, vcov, A)
    return A


def calculate_A_element(k: int, l: int, vcov: float, A: np.ndarray) -> float:
    """Helper function to calculate an element of matrix A."""
    if k == 0 and l == 0:
        return 1.0
    elif k == 0 and l > 0:
        return 0.0
    elif k == 1:
        if l == 0:
            return vcov
        elif l == 1:
            return np.sqrt(1.0 - A[:, :, k, l - 1] ** 2)
        else:
            return 0.0
    elif k >= 2:
        if l == 0:
            return vcov
        elif l < k:
            sum_term = np.sum(A[:, :, l, :l] * A[:, :, k, :l], axis=-1)
            return (vcov - sum_term) / A[:, :, l, l] if A[:, :, l, l].any() else 0
        elif l == k:
            sum_term = np.sum(A[:, :, k, :l] ** 2, axis=-1)
            return np.sqrt(1.0 - sum_term)
    return 0.0


def get_dist(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Calculate distance between two points on Earth's surface."""
    lat1_rad, lon1_rad = np.radians(lat1), np.radians(lon1)
    lat2_rad, lon2_rad = np.radians(lat2), np.radians(lon2)

    dlat, dlon = lat2_rad - lat1_rad, lon2_rad - lon1_rad

    a = (
        np.sin(dlat / 2) ** 2
        + np.cos(lat1_rad) * np.cos(lat2_rad) * np.sin(dlon / 2) ** 2
    )
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

    return EARTH_RADIUS_KM * c


def apply_horizontal_correlations(
    field: np.ndarray,
    nx: int,
    ny: int,
    nz: int,
    corr_length_hz: float,
    grid_length: float,
    lat_grid: np.ndarray,
    lon_grid: np.ndarray,
) -> np.ndarray:
    """Apply horizontal correlations to the field."""
    ngrid_corr = int(np.ceil(corr_length_hz / grid_length)) + 1
    chem_fac = np.zeros_like(field)

    for imem in tqdm(range(settings.mems), unit="member"):
        for i in tqdm(range(nx), unit="xcell", leave=False):
            for j in range(ny):
                ii_str, ii_end = max(0, i - ngrid_corr), min(nx, i + ngrid_corr)
                jj_str, jj_end = max(0, j - ngrid_corr), min(ny, j + ngrid_corr)

                dist = np.sqrt(
                    (lat_grid[i, j] - lat_grid[ii_str:ii_end, jj_str:jj_end]) ** 2
                    + (lon_grid[i, j] - lon_grid[ii_str:ii_end, jj_str:jj_end]) ** 2
                )
                within_distance = dist <= corr_length_hz
                wgt = (
                    np.exp(-(dist**2) / (corr_length_hz**2))[:, :, np.newaxis]
                    * within_distance[:, :, np.newaxis]
                )

                for k in range(nz):
                    chem_fac[imem, i, j, k] = np.sum(
                        wgt
                        * field[imem, ii_str:ii_end, jj_str:jj_end, k][
                            :, :, np.newaxis
                        ],
                        axis=(0, 1),
                    )
                    weights_sum = np.sum(wgt[:, :])
                    if weights_sum != 0:
                        chem_fac[imem, i, j, k] /= weights_sum
    return chem_fac


def box_muller_random_field(nx: int, ny: int, nz: int, num_mem: int) -> np.ndarray:
    """Generate random field using Box-Muller transform."""
    u1 = np.random.rand(num_mem, nx, ny, nz)
    u2 = np.random.rand(num_mem, nx, ny, nz)
    r = np.sqrt(-2.0 * np.log(u1))
    theta = 2.0 * np.pi * u2
    return r * np.cos(theta)


def apply_vertical_correlation(pert_fields: np.ndarray, A: np.ndarray) -> np.ndarray:
    """Apply vertical correlation to perturbation fields."""
    return np.einsum("ijkl,mijl->mijk", A, pert_fields)


def recenter_and_rescale(field: np.ndarray, spread: float) -> np.ndarray:
    """Recenter and rescale the field."""
    mean = np.mean(field, axis=-1, keepdims=True)
    std = np.std(field, axis=-1, ddof=1, keepdims=True)
    return (field - mean) * (spread / std)


def perturb_emission():
    """Perturb emission data."""
    netcdfs = glob(str(settings.path_emissions / settings.name_netcdfs))
    for netcdf_emi in netcdfs:
        try:
            emission_dataset = xr.open_dataset(netcdf_emi)
        except FileNotFoundError:
            print(f"Error: Could not find file {netcdf_emi}")
            return

        nx, ny, nz = (
            len(emission_dataset[settings.var].lon),
            len(emission_dataset[settings.var].lat),
            len(emission_dataset[settings.var].z),
        )

        A = get_vertical_correlation_matrix(nx, ny, nz, settings.corr_length_hz)
        chem_fac_pr = None

        for i in range(1):
            random_field = box_muller_random_field(nx, ny, nz, settings.mems)

            grid_length = get_dist(
                emission_dataset[settings.var].lat[0],
                emission_dataset[settings.var].lat[1],
                emission_dataset[settings.var].lon[0],
                emission_dataset[settings.var].lat[1],
            )
            chem_fac = apply_horizontal_correlations(
                random_field,
                nx,
                ny,
                nz,
                settings.corr_length_hz,
                grid_length,
                *np.meshgrid(emission_dataset.lat, emission_dataset.lon),
            )
            chem_fac = apply_vertical_correlation(chem_fac, A)
            chem_fac = recenter_and_rescale(chem_fac, settings.spread)

            if chem_fac_pr is not None:
                alpha = np.exp(-1 / settings.corr_time)
                chem_fac = alpha * chem_fac_pr + np.sqrt(1 - alpha**2) * chem_fac

            for imem in range(settings.mems):
                chem_fac_t = np.transpose(chem_fac, axes=[imem, 3, 2, 1])
                emission_dataset[settings.var][i, :, :, :] *= np.exp(
                    chem_fac_t[imem, :, :, :]
                )
            chem_fac_pr = chem_fac

        try:
            for imem in range(settings.mems):
                path_filename = (
                    Path(netcdf_emi.strip(".nc")) / f"emi_{imem}" / f"_{imem}.nc"
                )
                if not os.path.exists(path_filename):
                    os.makedirs(path_filename)
                emission_dataset[settings.var].to_netcdf(path_filename)
        except Exception as e:
            print(f"Error saving netCDF file: {e}")


if __name__ == "__main__":
    perturb_emission()
