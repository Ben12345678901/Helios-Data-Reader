"""
Helios Data Plotting Utilities

Author: Xu Zhao
Affiliation: York Plasma Institute, School of Physics, Engineering and Technology,
             University of York, York YO10 5DD, United Kingdom
Email: xu.zhao@york.ac.uk
Last Updated: 2023-11-14

Description:
This module provides functions for loading and processing Helios data, as well as 
plotting functions for visualizing density, electron temperature, and radius evolution.
Each function is designed to produce publication-quality figures with consistent style.
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

def load_and_process_data(file_path):
    """
    Loads and processes data from the specified Helios file.

    Parameters:
    - file_path: str, path to the data file.

    Returns:
    - dict: A dictionary containing processed data arrays for time, radius, density, 
            electron density, temperature, pressure, velocity, and volume.
    """
    data = xr.open_dataset(file_path)

    # Extract variables and convert to numpy arrays
    time_whole = data['time_whole'].values
    zone_boundaries = data['zone_boundaries'].values
    mass_density = data['mass_density'].values
    elec_density = data['elec_density'].values
    ion_temperature = data['ion_temperature'].values
    elec_temperature = data['elec_temperature'].values
    zone_mass = data['zone_mass'].values
    pressure = data['ion_pressure'].values + data['elec_pressure'].values
    fluid_velocity = data['fluid_velocity'].values / 100000  # Convert to km/s
    volume = zone_mass / mass_density

    # Calculate time edges for pcolormesh
    time_diff = np.diff(time_whole) / 2
    time_edges = np.concatenate(([time_whole[0] - time_diff[0]], time_whole[:-1] + time_diff, [time_whole[-1] + time_diff[-1]]))

    # Calculate radius edges for pcolormesh
    radius_diff = np.diff(zone_boundaries, axis=0) / 2
    radius_edges = np.vstack((zone_boundaries[0, :] - radius_diff[0, :],
                              zone_boundaries[:-1, :] + radius_diff,
                              zone_boundaries[-1, :] + radius_diff[-1, :]))

    return {
        "time_whole": time_whole,
        "zone_boundaries": zone_boundaries,
        "mass_density": mass_density,
        "elec_density": elec_density,
        "ion_temperature": ion_temperature,
        "elec_temperature": elec_temperature,
        "zone_mass": zone_mass,
        "pressure": pressure,
        "fluid_velocity": fluid_velocity,
        "volume": volume,
        "time_edges": time_edges,
        "radius_edges": radius_edges
    }

def plot_radius_evolution(time, radius, figsize=(8.5 / 2.54, 8.5 / 1.618 / 2.54), line_width=0.5,
                          line_color='black', font_size=7, font_family='Arial',
                          dpi=200, border_width=0.5, tick_length=3, tick_width=0.5):
    """
    Plots the evolution of radius over time for each zone.

    Parameters:
    - time: 1D array-like, time points.
    - radius: 2D array-like, radius values for each time and zone.
    - figsize: tuple, figure size in inches, default (8.5 cm, 6.5 cm).
    - line_width: float, line width for the plot.
    - line_color: str, line color.
    - font_size: int, font size for labels and title.
    - font_family: str, font family for labels and title.
    - dpi: int, resolution of the figure in dots per inch.
    - border_width: float, width of the figure border.
    - tick_length: float, length of ticks on both axes.
    - tick_width: float, width of ticks on both axes.

    Returns:
    - ax: Matplotlib Axes object, for further customization.
    """
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    for i in range(radius.shape[1]):
        ax.plot(time, radius[:, i], color=line_color, lw=line_width)

    ax.set_xlabel("Time (ns)", fontsize=font_size, fontfamily=font_family)
    ax.set_ylabel(r"Radius ($\rm{\mu}$m)", fontsize=font_size, fontfamily=font_family)
    ax.set_title("Radius Evolution Over Time for Each Zone", fontsize=font_size, fontfamily=font_family)
    ax.tick_params(axis='both', which='major', labelsize=font_size, length=tick_length, width=tick_width)
    for spine in ax.spines.values():
        spine.set_linewidth(border_width)
    return ax

def plot_density_pcolormesh(time_edges, radius_edges, density, figsize=(8.5 / 2.54, 8.5 / 1.618 / 2.54),
                            cmap='jet', font_size=7, font_family='Arial',
                            dpi=200, border_width=0.5, tick_length=3, tick_width=0.5):
    """
    Plots a density pcolormesh with a color bar.

    Parameters:
    - time_edges: 1D array-like, edges for the time axis (scaled to ns).
    - radius_edges: 1D or 2D array-like, edges for the radius axis (scaled to microns).
    - density: 2D array-like, density values for each time and radius.
    - figsize, cmap, font_size, font_family, dpi, border_width, tick_length, tick_width.

    Returns:
    - ax: Matplotlib Axes object.
    """

    # Convert figsize from cm to inches
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    # Plot pcolormesh for density data
    cmesh = ax.pcolormesh(time_edges, radius_edges.T, density.T, shading='auto', cmap=cmap)

    # Add color bar
    cbar = fig.colorbar(cmesh, ax=ax)
    cbar.ax.tick_params(labelsize=font_size, length=tick_length, width=tick_width)
    cbar.outline.set_linewidth(border_width)  # Set colorbar border thickness

    # Set colorbar label at the top and horizontal
    cbar_label = r"$\rho$ (g/cc)"
    cbar.ax.text(0.5, 1.02, cbar_label, ha='center', va='bottom',
                 fontsize=font_size, fontfamily=font_family, transform=cbar.ax.transAxes)

    # Set labels and title
    ax.set_xlabel("Time (ns)", fontsize=font_size, fontfamily=font_family)
    ax.set_ylabel(r"Radius ($\rm{\mu}$m)", fontsize=font_size, fontfamily=font_family)
    ax.set_title("Density Plot", fontsize=font_size, fontfamily=font_family)

    # Set font size and family for tick labels
    ax.tick_params(axis='both', which='major', labelsize=font_size, length=tick_length, width=tick_width)
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontname(font_family)

    # Set border width for all spines
    for spine in ax.spines.values():
        spine.set_linewidth(border_width)

    # Return the Axes object for further customization
    return ax

def plot_electron_temp_pcolormesh(time_edges, radius_edges, elec_temp, figsize=(8.5 / 2.54, 8.5 / 1.618 / 2.54),
                                  cmap='jet', font_size=7, font_family='Arial',
                                  dpi=200, border_width=0.5, tick_length=3, tick_width=0.5):
    """
    Plots an electron temperature (Te) pcolormesh with a color bar.

    Parameters:
    - time_edges, radius_edges, elec_temp, figsize, cmap, font_size, font_family, dpi, border_width, tick_length, tick_width.

    Returns:
    - ax: Matplotlib Axes object.
    """
    # Convert figsize from cm to inches
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    # Plot pcolormesh for electron temperature data
    cmesh = ax.pcolormesh(time_edges, radius_edges.T, elec_temp.T, shading='auto', cmap=cmap, vmin=0.025, vmax=70)

    # Add color bar
    cbar = fig.colorbar(cmesh, ax=ax)
    cbar.ax.tick_params(labelsize=font_size, length=tick_length, width=tick_width)
    cbar.outline.set_linewidth(border_width)  # Set colorbar border thickness

    # Set colorbar label at the top and horizontal
    cbar_label = r"$T_e$ (eV)"
    cbar.ax.text(0.5, 1.02, cbar_label, ha='center', va='bottom',
                 fontsize=font_size, fontfamily=font_family, transform=cbar.ax.transAxes)

    # Set labels and title
    ax.set_xlabel("Time (ns)", fontsize=font_size, fontfamily=font_family)
    ax.set_ylabel(r"Radius ($\rm{\mu}$m)", fontsize=font_size, fontfamily=font_family)
    ax.set_title("Electron Temperature Contour", fontsize=font_size, fontfamily=font_family)

    # Set font size and family for tick labels
    ax.tick_params(axis='both', which='major', labelsize=font_size, length=tick_length, width=tick_width)
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontname(font_family)

    # Set border width for all spines
    for spine in ax.spines.values():
        spine.set_linewidth(border_width)

    # Return the Axes object for further customization
    return ax


#example test

file = r"C:\Users\benny\OneDrive\Desktop\Helios\SACLA_Imprint_2025\CH_GLUE_GOLD_Lif.exo"

data_dict = load_and_process_data(file)
print(data_dict.keys())

time = data_dict['time_whole']*1e9
radius = data_dict['zone_boundaries']*1e4

# Calculate time edges and radius edges for pcolormesh (required for accurate cell boundaries)
time_edges = data_dict['time_edges'] *1e9
radius_edges = data_dict['radius_edges']*1e4

density = data_dict['mass_density']
elec_temp = data_dict['elec_temperature']

ax = plot_density_pcolormesh(time_edges, radius_edges, np.log(density))
#ax.set_xlim(0,100)
ax.set_ylim(-10,810)
plt.show()

'''
ax = plot_electron_temp_pcolormesh(time_edges, radius_edges, elec_temp)
ax.set_xlim(0,5)
ax.set_ylim(-10,60)

ax = plot_radius_evolution(time, radius,line_width=0.5)
ax.set_xlim(0,5)
ax.set_ylim(-10,60)
plt.show()
'''
#####################################
'''
Modified by Ben Fisher
Date: 13/03/2024

Add in shock tracking using the highest density gradient

currently does not track the shock properly but does exactly what is has been told to do and tracks
the interface boundary where the density gradient is the highest

would work for single material shocks
'''

def track_shock_position(radius, density):
    """
    Tracks the shock position in each time slice by locating the maximum absolute gradient in density.
    
    Parameters:
    - radius: 2D numpy array containing either cell centers or zone boundaries (n_time x n_zone) 
              or (n_time x n_zone+1) respectively.
    - density: 2D numpy array (n_time x n_zone) containing density values.
    
    Returns:
    - shock_positions: 1D numpy array of shock positions (in microns) for each time slice.
    """
    shock_positions = []
    
    # Check if radius is provided as boundaries (n_zone+1) and compute centers if needed
    if radius.shape[1] == density.shape[1] + 1:
        # Compute the centers from boundaries
        radius_centers = 0.5 * (radius[:, :-1] + radius[:, 1:])
    else:
        radius_centers = radius
    
    # Loop over each time slice to compute the gradient and locate the shock front.
    for t in range(density.shape[0]):
        dens_slice = density[t, :]
        rad_slice = radius_centers[t, :]
        # Compute the gradient of density with respect to radius
        grad = np.gradient(dens_slice, rad_slice)
        # Identify the shock front as the position with maximum absolute gradient
        shock_index = np.argmax(np.abs(grad))
        shock_positions.append(rad_slice[shock_index])
        
    return np.array(shock_positions)


def compute_shock_velocity(time, shock_positions):
    """
    Computes the shock velocity as the derivative of shock position with respect to time.
    
    Parameters:
    - time: 1D array of time values (in ns).
    - shock_positions: 1D array of shock positions (in microns).
    
    Returns:
    - shock_velocity: 1D array of shock velocities.
    
    Note:
    - With time in nanoseconds and radius in microns, the derivative (micron/ns) is equivalent to km/s,
      because 1 micron/ns = 1e-6 m/1e-9 s = 1e3 m/s = 1 km/s.
    """
    shock_velocity = np.gradient(shock_positions, time)
    return shock_velocity

def plot_shock_tracking(time, shock_positions, shock_velocity, font_size=7, font_family='Arial', dpi=200):
    """
    Plots the shock position and shock velocity versus time.
    
    Parameters:
    - time: 1D array of time values (in ns).
    - shock_positions: 1D array of shock positions (in microns).
    - shock_velocity: 1D array of shock velocities (in km/s).
    - font_size: int, font size for labels.
    - font_family: str, font family for labels.
    - dpi: int, resolution of the figure.
    
    Returns:
    - fig, axs: Matplotlib Figure and Axes objects.
    """
    fig, axs = plt.subplots(2, 1, figsize=(8, 6), dpi=dpi, sharex=True)
    
    # Shock position plot
    axs[0].plot(time, shock_positions, marker='o', linestyle='-', color='blue')
    axs[0].set_ylabel("Shock Position ($\mu$m)", fontsize=font_size, fontfamily=font_family)
    axs[0].set_title("Shock Position vs Time", fontsize=font_size, fontfamily=font_family)
    axs[0].tick_params(axis='both', labelsize=font_size)
    
    # Shock velocity plot
    axs[1].plot(time, shock_velocity, marker='o', linestyle='-', color='red')
    axs[1].set_xlabel("Time (ns)", fontsize=font_size, fontfamily=font_family)
    axs[1].set_ylabel("Shock Velocity (km/s)", fontsize=font_size, fontfamily=font_family)
    axs[1].set_title("Shock Velocity vs Time", fontsize=font_size, fontfamily=font_family)
    axs[1].tick_params(axis='both', labelsize=font_size)
    
    plt.tight_layout()
    return fig, axs


def plot_density_with_shock(time_edges, radius_edges, density, time, shock_positions,
                            figsize=(8.5 / 2.54, 8.5 / 1.618 / 2.54),
                            cmap='jet', font_size=7, font_family='Arial',
                            dpi=200, border_width=0.5, tick_length=3, tick_width=0.5):
    """
    Plots a density pcolormesh and overlays the shock position on top.

    Parameters:
    - time_edges: 1D array, time cell edges.
    - radius_edges: 2D array, radius cell edges.
    - density: 2D array, density values.
    - time: 1D array, time values (e.g., in ns) corresponding to shock_positions.
    - shock_positions: 1D array, shock positions (in microns) for each time slice.
    - (Other plotting parameters)

    Returns:
    - ax: Matplotlib Axes object.
    """
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    # Plot density pcolormesh
    cmesh = ax.pcolormesh(time_edges, radius_edges.T, density.T, shading='auto', cmap=cmap)
    cbar = fig.colorbar(cmesh, ax=ax)
    cbar.ax.tick_params(labelsize=font_size, length=tick_length, width=tick_width)
    cbar.outline.set_linewidth(border_width)
    cbar_label = r"$\rho$ (g/cc)"
    cbar.ax.text(0.5, 1.02, cbar_label, ha='center', va='bottom',
                  fontsize=font_size, fontfamily=font_family, transform=cbar.ax.transAxes)

    # Set labels and title
    ax.set_xlabel("Time (ns)", fontsize=font_size, fontfamily=font_family)
    ax.set_ylabel(r"Radius ($\rm{\mu}$m)", fontsize=font_size, fontfamily=font_family)
    ax.set_title("Density Plot with Shock Tracking", fontsize=font_size, fontfamily=font_family)
    ax.tick_params(axis='both', which='major', labelsize=font_size, length=tick_length, width=tick_width)
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontname(font_family)
    for spine in ax.spines.values():
        spine.set_linewidth(border_width)

    # Overlay the shock front. Adjust the color and markers as needed.
    ax.plot(time, shock_positions, color='white', linestyle='-', marker='o',
            markersize=3, label="Shock Front")
    ax.legend(loc='upper right', fontsize=font_size)
    return ax


# --- Usage Example ---

# (Assuming you already loaded your data and defined 'time', 'radius', and 'density')
# Here, time and radius are already scaled: time in ns and radius in microns.

# Load your data
data_dict = load_and_process_data(file)
time = data_dict['time_whole']*1e9         # time in ns
density = data_dict['mass_density']         # density values
zone_boundaries = data_dict['zone_boundaries']  # zone boundaries, shape (n_time, n_zone)

# Optionally, scale zone_boundaries if needed, e.g., converting to microns
radius = zone_boundaries * 1e4

# Define the time range for analysis (in ns)
time_min = 0
time_max = 20

# Create a boolean mask for the specified time range
mask = (time >= time_min) & (time <= time_max)

# Trim the time, density, and zone_boundaries arrays using the mask.
time_trimmed = time[mask]
density_trimmed = density[mask, :]
zone_boundaries_trimmed = zone_boundaries[mask, :]

# If needed, compute radius from the trimmed zone boundaries (with scaling)
radius_trimmed = zone_boundaries_trimmed * 1e4

# Recompute time_edges for the trimmed time array.
time_diff_trim = np.diff(time_trimmed) / 2
time_edges_trimmed = np.concatenate((
    [time_trimmed[0] - time_diff_trim[0]],
    time_trimmed[:-1] + time_diff_trim,
    [time_trimmed[-1] + time_diff_trim[-1]]
))

# Recompute radius_edges from the trimmed zone boundaries.
# (This code assumes that your original radius_edges were computed similarly)
radius_diff_trim = np.diff(zone_boundaries_trimmed, axis=0) / 2
radius_edges_trimmed = np.vstack((
    zone_boundaries_trimmed[0, :] - radius_diff_trim[0, :],
    zone_boundaries_trimmed[:-1, :] + radius_diff_trim,
    zone_boundaries_trimmed[-1, :] + radius_diff_trim[-1, :]
))
# Scale radius_edges as needed:
radius_edges_trimmed = radius_edges_trimmed * 1e4

# Compute shock positions and velocity using the trimmed data.
shock_positions = track_shock_position(radius_trimmed, density_trimmed)
shock_velocity = compute_shock_velocity(time_trimmed, shock_positions)

# Plot the density with the shock overlay using the trimmed data.
ax = plot_density_with_shock(
    time_edges_trimmed,         # trimmed time edges for pcolormesh
    radius_edges_trimmed,       # trimmed radius edges
    np.log(density_trimmed),    # using trimmed density (or log-density, as desired)
    time_trimmed,               # trimmed time for shock overlay
    shock_positions
)
plt.show()

# Optionally, plot the shock tracking (position and velocity) on separate subplots.
fig, axs = plot_shock_tracking(time_trimmed, shock_positions, shock_velocity)
ax.set_ylim(-10,810)
plt.show()
