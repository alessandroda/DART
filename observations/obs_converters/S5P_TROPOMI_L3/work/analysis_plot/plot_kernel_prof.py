import xarray as xr
import matplotlib.pyplot as plt

ds = xr.open_dataset(
    "/mnt/mumbai_n4r5/dausilio/projects/DART/observations/obs_converters/S5P_TROPOMI_L3/work/S5p_NO2_16742.nc"
)
y = ds.pressure[1].values[:-1]
x = ds.kernel_trop[250].values

# Plot kernel values along pressure
plt.figure(figsize=(6, 8))  # Adjust figsize as needed
plt.plot(x, y, marker="o", linestyle="-", label="Kernel Values")  # Line-point plot
plt.ylabel("Pressure")
plt.xlabel("Kernel Values")
plt.gca().invert_yaxis()  # Invert y-axis for pressure
plt.xscale("log")  # Use logarithmic scale for x-axis
plt.legend(loc="upper right")  # Display legend based on pressure levels
plt.savefig("kernel.png", dpi=300)
plt.close()
