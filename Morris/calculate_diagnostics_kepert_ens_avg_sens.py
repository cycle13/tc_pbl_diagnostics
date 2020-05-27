#!/Users/kmn182/anaconda3/envs/ENV1/bin/python
# ------------------------------------------------------------
#
#	Program Name: calculate_diagnostics_kepert_ens_avg_sens.py
#
#	Purpose: To calculate model diagnostics for a given 
#            time step. Output comes from the "Kepert" model and
#            simulates a tropical cyclone on an f-plane. Output
#            is located at the radius of maximum wind.
#
#	WARNING: This code assumes Python 3, so there may be some 
#            syntax errors if running with Python 2.
#
# -------------------------------------------------------------

# Step 1 (import necessary modules)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from netCDF4 import Dataset
import os
import sys
import csv
from output_functions import calc_p_sfc_ens_avg

# define environmental path variable
GOOGLE_ROOT = os.environ['GOOGLE_ROOT']

# define the font family
mpl.rc('font',family='Arial')

# define the configurations
#configs = ['RCE.QPC6.ne0np4tcfplane.ne15x8.exp011.001']
configs = [sys.argv[1]]

# define the time step labels
time_steps = range(1, 17)
time_steps = np.asarray(time_steps, dtype='str')
for t_i, time in enumerate(time_steps):
	time_steps[t_i] = time.zfill(2)
	
# how many ensemble members?
num_ensembles = int(sys.argv[2])

ensembles = range(1, int(num_ensembles) + 1)

# how many IC ensemble members?
num_IC_ens = 3

# Step 2 (read in the input files)

for c_i, config_base in enumerate(configs):
	for e_i, ens in enumerate(ensembles):

		ens_str = str(ens).rjust(3, '0')
		
		config = config_base + '.' + ens_str

		# define the input path
		input_path = '/Users/kmn182/ICS_scratch/output/ens_avg/' + config + '/'
	
		# define the output path
		output_path = GOOGLE_ROOT + '/CPT_Project/ANALYSIS/diagnostics/' + config + '/'
		if not os.path.exists(output_path):
			os.mkdir(output_path)
	
		# calculate surface pressure minimum and surface wind maximum for the configuration (all times included)
		min_sfc_p, max_sfc_wind = calc_p_sfc_ens_avg(config, num_IC_ens)

		# Step 2 (create a CSV file to store the model diagnostics)

		f = open(output_path + 'diagnostics_kepert_' + config + '_ens_avg.csv', mode='w')
		csv_writer = csv.writer(f, delimiter=',')
		csv_writer.writerow(['Time Step (days)', '', 'Radius Max Wind (m)', '', 'Surface Inflow Angle (deg)', '', 'Inflow Depth (m)', '', 'Surface Max Wind Speed (m/s)', '', 'Min Surface Pressure (hPa)', '', 'Level of Max Wind (m)'])

		for t_i, time_step in enumerate(time_steps):

			# Step 3 (read in the input file...assuming one at this point)
		
			# define file name
			file_name = config + '_ens_avg_' + time_step + '.nc'

			dataset = Dataset(input_path + file_name)

			rad_T = dataset.variables['rad_T'][:] # average temperature along the radial (in degrees K)
			rad_q = dataset.variables['rad_Q'][:] # average specific humidity along the radial (in kg/kg)
			rad_v_theta = dataset.variables['rad_v_theta'][:] # average tangential wind along the radial (in m/s)
			rad_v_rad = dataset.variables['rad_v_rad'][:] # average radial wind along the radial (in m/s)
			rad_w = dataset.variables['rad_W'][:] # average vertical velocity along the radial (in m/s)
			rad_kvm = dataset.variables['rad_KVM'][:] # average vertical diffusivity of heat/moisture (in m**2/s)
			rad_tke = dataset.variables['rad_TKE'][:] # average turbulent kinetic energy (in m**2/s)
			rad_lscale = dataset.variables['rad_LSCALE'][:] # average turbulent length scale (in m)
			rad_theta = dataset.variables['rad_THETA'][:] # average potential temperature (in K)
			levels = dataset.variables['lev'][:] # vertical z levels (in m)
			radius = dataset.variables['radius'][:] # horizontal distance from the storm center (in m)

			min_level = int(np.nanmin(levels))
			min_radius = int(np.nanmin(radius))
			max_level = int(np.nanmax(levels))
			max_radius = int(np.nanmax(radius))

			# Step 4 (calculate the model diagnostics)

			# find the radius of maximum wind (tangential)
			max_wind = np.nanmax(rad_v_theta) # maximum tangential wind
			[max_levels, max_radius] = np.where(rad_v_theta == max_wind)
			max_radius_index = max_radius[0]
			radius_of_max_wind = radius[max_radius_index]

			# find 2x radius of maximum wind
			if 2*radius_of_max_wind > np.nanmax(radius):
				max_radius_index_2x = np.where(radius == np.nanmax(radius))[0][0]
			else:
				max_radius_index_2x = np.where(radius == 2*radius_of_max_wind)[0][0]

			# find the surface inflow angle at the radius of maximum wind
			# first, find the lowest model level with non-NaN data
			for l_i, level in enumerate(levels):
				level_v_theta = rad_v_theta[l_i, max_radius_index]
				if np.ma.is_masked(level_v_theta):
					continue
				else:
					lowest_level_index = l_i
					lowest_level = level
					break
			sfc_v_thetas = rad_v_theta[lowest_level_index, :]
			sfc_v_rads = rad_v_rad[lowest_level_index, :]
			sfc_v_theta = rad_v_theta[lowest_level_index, max_radius_index]
			sfc_v_rad = rad_v_rad[lowest_level_index, max_radius_index] # these assume that the surface is the lowest model layer with data
			sfc_inflow_angle = np.arctan(sfc_v_rad / sfc_v_theta) # in radians
			sfc_inflow_angle = np.absolute(sfc_inflow_angle * (180. / np.pi))
		
			# calculate the level at which the maximum tangential wind speed occurs
			max_v_theta = np.nanmax(rad_v_theta[:, max_radius_index])
			max_v_theta_index = np.where(rad_v_theta[:, max_radius_index] == max_v_theta)[0][0]
			max_v_theta_level = levels[max_v_theta_index] # level in meters

			# find the inflow depth based on the surface radial flow at the 2x radius of maximum wind
			sfc_radial_flow = np.absolute(rad_v_rad[lowest_level_index, max_radius_index_2x])
			# loop through each level at the 2x radius of maximum wind to find the level at which inflow gets below 10% of surface radial flow
			for l_i, level in enumerate(levels):
				level_radial_flow = np.absolute(rad_v_rad[l_i, max_radius_index_2x])
				if not np.ma.is_masked(level_radial_flow):
					if level_radial_flow > (0.1 * sfc_radial_flow):
						inflow_depth = level
						continue
					else:
						break
					
			dataset.close()
		
			# Step 5 (output the results to a netCDF file)

			time = time_step
			csv_writer.writerow([time, '', str(radius_of_max_wind), '', str(sfc_inflow_angle), '', str(inflow_depth), '', str(max_sfc_wind[t_i]), '', str(min_sfc_p[t_i]), '', max_v_theta_level])
	
		f.close() # close and save the file
		print('The file has been written for ' + config + '!')
	




	

			
		





