from Helios_Data_Reader_V2 import *

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

ax = plot_density_pcolormesh(time_edges, radius_edges, density)
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