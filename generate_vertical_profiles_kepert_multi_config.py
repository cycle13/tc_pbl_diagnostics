#!/Users/kmn182/anaconda3/envs/ENV1/bin/python
# ------------------------------------------------------------
#
#	Program Name: generate_vertical_profiles_kepert_multi_config.py
#
#	Purpose: To plot vertical profiles of model output for a given 
#            time step. Output comes from the "Kepert" model and
#            simulates a tropical cyclone on an f-plane. Output
#            is located at the radius of maximum wind. This is done
#            for multiple configurations.
#
#	WARNING: This code assumes Python 3, so there may be some 
#            syntax errors if running with Python 2.
#
# -------------------------------------------------------------

# Step 1 (import necessary modules)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.animation as anim
from netCDF4 import Dataset
import os
import sys
import csv
from palettable.cmocean.diverging import Curl_20_r
from palettable.cartocolors.diverging import Earth_7
from palettable.cmocean.sequential import Ice_20_r, Tempo_20, Dense_20
from palettable.cubehelix import perceptual_rainbow_16_r, classic_16_r, cubehelix1_16_r, cubehelix2_16_r, jim_special_16_r, cubehelix3_16_r, perceptual_rainbow_16_r
from palettable.cubehelix import Cubehelix
from obs_functions import vickery, get_NOAA_profiles, composite_NOAA_profiles, composite_NASA_profiles, composite_G_IV_profiles, composite_P_3_profiles, composite_NOAA_HRD_profiles, std_dev
from output_functions import calc_d_dz

# define environmental path variable (will need to be changed based on user's system)
GOOGLE_ROOT = os.environ['GOOGLE_ROOT']

# define the font family
mpl.rc('font',family='Arial')

# plot linear vertical axes?
plot_linear = True # if False, logarithmic

# the following Boolean variables are used later for plotting reference profiles

# plot the observed wind profiles?
plot_vickery = False

# plot a NOAA observed wind profile?
plot_NOAA = False

# plot Zhang et al. (2012) Km profile?
plot_Zhang = True

# plot composites of observed profiles?
plot_composite_NOAA = False
plot_composite_NASA = False
plot_composite_P_3 = False
plot_composite_G_IV = False
plot_composite_NOAA_HRD = False

# plot betacast profile?
plot_betacast = False

if plot_Zhang:
	zhang_levs = np.array([110., 190., 280., 380., 490., 750.]) # these are estimates from Figure 8a
	zhang_kvm = np.array([25., 50., 40., 35., 20., 25.])

if plot_betacast:
	date = '2016-10-07-21600'
	config = 'TC-betacast.hindcast.Matthew.expb101'
	file = open(GOOGLE_ROOT + '/CPT_Project/ANALYSIS/diagnostics/' + config + '/diagnostics_betacasts_' + config + '.csv', mode='r')
	csv_reader = csv.reader(file, delimiter=',')
	counter = 0
	dates = []
	rads_max_wind = []
	for row in csv_reader: # indexing may need to change based on where in CSV file radius of maximum wind is calculated
		if counter == 0:
			counter += 1
			continue
		else:
			counter += 1
			dates.append(str(row[0]))
			rads_max_wind.append(str(row[2]))

	dates = np.asarray(dates, dtype='str')
	rads_max_wind = np.asarray(rads_max_wind, dtype='float')
	
	date_i = np.where(dates == date)[0][0]
	rad_max_wind = rads_max_wind[date_i]
	
	# read in the radial average file
	dataset = Dataset('/Users/kmn182/ICS_scratch/betacasts/output_regridded/' + config + '/' + config + '_' + date + '_radial_avg_regridded.nc', 'r')
	radius = dataset.variables['radius'][:]
	rad_max_wind_i = np.where(radius == rad_max_wind)[0][0]
	beta_v_theta = dataset.variables['rad_v_theta'][:, rad_max_wind_i]
	beta_v_rad = dataset.variables['rad_v_rad'][:, rad_max_wind_i]
	beta_khv = dataset.variables['rad_KVM'][:, rad_max_wind_i]
	beta_kvm = c_k10 * beta_khv
	beta_lscale = dataset.variables['rad_LSCALE'][:, rad_max_wind_i]
	beta_levs = dataset.variables['lev'][:]
	beta_dv_dz = calc_d_dz(beta_v_theta, beta_levs)

# plot CM1 profile? 
plot_cm1 = True
if plot_cm1:
	date = '000120'
	config = 'cm1out'
	file = open(GOOGLE_ROOT + '/CPT_Project/ANALYSIS/diagnostics/' + config + '/diagnostics_cm1_' + config + '.csv', 'r')
	csv_reader = csv.reader(file, delimiter=',')
	counter = 0
	dates = []
	rads_max_wind = []
	for row in csv_reader:
		if counter == 0:
			counter += 1
			continue
		else:
			counter += 1
			dates.append(str(row[0]))
			rads_max_wind.append(str(row[2]))

	dates = np.asarray(dates, dtype='str')
	rads_max_wind = np.asarray(rads_max_wind, dtype='float')
	
	date_i = np.where(dates == date)[0][0]
	rad_max_wind = rads_max_wind[date_i] 
	
	# read in the CM1 radial average file
	dataset = Dataset('/Users/kmn182/ICS_scratch/CM1_data/' + config + '_' + date + '_radial_avg.nc', 'r')
	radius = dataset.variables['radius'][:]
	rad_max_wind_i = np.where(radius == rad_max_wind)[0][0]
	cm1_w = dataset.variables['rad_W'][:, rad_max_wind_i]
	cm1_q = dataset.variables['rad_Q'][:, rad_max_wind_i]
	cm1_theta = dataset.variables['rad_THETA'][:, rad_max_wind_i]
	cm1_v_theta = dataset.variables['rad_v_theta'][:, rad_max_wind_i]
	cm1_v_rad = dataset.variables['rad_v_rad'][:, rad_max_wind_i]
	cm1_kvm = dataset.variables['rad_KMV'][:, rad_max_wind_i]
	cm1_levs = dataset.variables['lev'][:] * 1000.
	cm1_rh = dataset.variables['rad_RH'][:, rad_max_wind_i]
	cm1_dv_dz = calc_d_dz(cm1_v_theta, cm1_levs)
	
	# now, variables from an additional azimuthal average file
	dataset = Dataset('/Users/kmn182/ICS_scratch/CM1_data/azim_avg.nc', 'r')
	azim_radius = dataset.variables['xh'][:]
	azim_cm1_levs = dataset.variables['zf'][:]  # in meters
	
	# find the radius of maximum wind (tangential)
	azim_cm1_spd = dataset.variables['wsp'][0, :, 0, :]
	max_wind = np.nanmax(azim_cm1_spd) # maximum tangential wind
	[max_levels, max_radius] = np.where(azim_cm1_spd == max_wind)
	max_radius_index = max_radius[0]
	azim_rad_max_wind = azim_radius[max_radius_index]
	
	diff_rad_max_wind = np.absolute(azim_radius - azim_rad_max_wind)
	rad_max_wind_i = np.where(diff_rad_max_wind == np.nanmin(diff_rad_max_wind))[0][0]
	cm1_stke = dataset.variables['stke'][0, :, 0, rad_max_wind_i]
	cm1_leffs = dataset.variables['leffs'][0, :, 0, rad_max_wind_i]
	cm1_leffv = dataset.variables['leffv'][0, :, 0, rad_max_wind_i]
	cm1_leffu = dataset.variables['leffu'][0, :, 0, rad_max_wind_i]
	cm1_wp = dataset.variables['wpf'][0, :, 0, rad_max_wind_i]
	cm1_up = dataset.variables['upf'][0, :, 0, rad_max_wind_i]
	cm1_vp = dataset.variables['vpf'][0, :, 0, rad_max_wind_i]
	cm1_leff = cm1_leffu
	cm1_kvm = dataset.variables['kmv'][0, :, 0, rad_max_wind_i]
	cm1_tke = cm1_stke + 0.5 * (cm1_wp**2 + cm1_vp**2 + cm1_up**2)

# if NOAA composites are desired, go read-in the files
if plot_composite_NOAA:
	u_wind_profile, v_wind_profile, w_wind_profile, tan_wind_profile, rad_wind_profile, spd_profile, z_profile, p_profile, T_profile, q_profile, pot_T_profile, max_heights, inflow_angles, inflow_depths = composite_NOAA_profiles(0, 1000, 40., 80., storm_exceptions=False, chosen_storms=['Katrina', 'Isabel', 'Ivan', 'Emily', 'Rita', 'Wilma', 'Dean', 'Felix'])

num_panels = 9 # how many panels should be created?

# define the configurations
num_configs = int(sys.argv[1])
configs = [sys.argv[i] for i in range(2, 2 + num_configs)]
c_k10 = 0.5 # assume this is constant, but it may not be depending on the configuration

# define a list of colors
custom_colormap = Cubehelix.make(start_hue=240, end_hue=-300, min_sat=1, max_sat=2.5, gamma=1, min_light=0.3, max_light=0.8, n=num_configs, reverse=False, name='custom_colormap')
colors = custom_colormap.hex_colors

# define the time step labels
time_steps = ['08_10_12']
	
num_cols = 3 # how many columns to use when plotting the panels?

# the following lists will store all variables for all time steps to be plotted
# consider these master lists
vars_to_plot_all_times = []
vars_to_plot_short_names_all_times = []
vars_to_plot_names_all_times = []
vars_to_plot_units_all_times = []
vars_to_plot_ticks_all_times = []
max_level_all_configs = []
levels_all_configs = []

# Step 2 (read in the input files)

for c_i, config in enumerate(configs):

	#define the input and output paths
	input_path = '/Users/kmn182/ICS_scratch/output/ens_avg/' + config + '/'
	output_path = GOOGLE_ROOT + '/CPT_Project/FIGURES/multi_config/'
	
	# check if the directory already exists 
	# if it does not, make the directory
	if not os.path.exists(output_path):
		os.mkdir(output_path)

	# for each time step, read in the input data from radially/azimuthally-averaged output files
	for t_i, time_step in enumerate(time_steps):
	
		# these lists store data for the particular time step (will be added to master list)
		vars_to_plot = []
		vars_to_plot_short_names = []
		vars_to_plot_names = []
		vars_to_plot_units = []
		vars_to_plot_bounds = []
		vars_to_plot_coarser_bounds = []
		vars_to_plot_ticks = []
	
		# define the input file name
		# assumed to be an ensemble average
		file_name = config + '_ens_avg_' + time_step + '.nc'
		
		if time_step == '08_10_12': # this file needs to be made separately
			file_name = config + '.8_10_12_ens_avg.nc'

		# read in the data
		dataset = Dataset(input_path + file_name)

		rad_T = dataset.variables['rad_T'][:] # average temperature along the radial (in degrees K)
		rad_q = dataset.variables['rad_Q'][:] # average specific humidity along the radial (in kg/kg)
		rad_v_theta = dataset.variables['rad_v_theta'][:] # average tangential wind along the radial (in m/s)
		rad_v_rad = dataset.variables['rad_v_rad'][:] # average radial wind along the radial (in m/s)
		rad_w = dataset.variables['rad_W'][:] # average vertical velocity along the radial (in m/s)
		rad_khv = dataset.variables['rad_KVM'][:] # eddy diffusivity for heat (in m**2/s) (disregard the odd netCDF var name)
		rad_kvm = c_k10 * rad_khv # how CLUBB calculates eddy viscosity
		rad_tke = dataset.variables['rad_TKE'][:] # average turbulence kinetic energy (in m**2/s)
		rad_lscale = dataset.variables['rad_LSCALE'][:] # average turbulent length scale (in m)
		rad_theta = dataset.variables['rad_THETA'][:] # average potential temperature (in K)
		rad_rh = dataset.variables['rad_RH'][:] # average relative humidity
		levels = dataset.variables['lev'][:] # vertical z levels (in m)
		radius = dataset.variables['radius'][:] # horizontal distance from the storm center (in m)

		min_level = int(np.nanmin(levels))
		min_radius = int(np.nanmin(radius))
		max_level = int(np.nanmax(levels))
		max_radius = int(np.nanmax(radius))

		# Step 3 (find the radius of maximum wind)

		max_wind = np.nanmax(rad_v_theta) # maximum tangential wind

		# now, find the radius at which this wind occurs
		# the vertical profiles are plotted at this radius
		[max_levels, max_radius] = np.where(rad_v_theta == max_wind)
		max_radius_index = max_radius[0] # this is the column index
		
		rad_dv_dz = calc_d_dz(rad_v_theta[:, max_radius_index], levels)
		
		# Step 4 (set up the panel plots with the appropriate variables)
		
		num_cols = 3 # how many columns to use when plotting the panels?
		
		# the arrays to be plotted
		vars_to_plot.extend([rad_dv_dz, rad_v_theta[:, max_radius_index], rad_v_rad[:, max_radius_index], rad_w[:, max_radius_index],\
		rad_theta[:, max_radius_index], rad_kvm[:, max_radius_index], rad_tke[:, max_radius_index], rad_lscale[:, max_radius_index], rad_q[:, max_radius_index]])
		vars_to_plot_all_times.append(vars_to_plot)
		
		# the shortened identifiers for each variable (used in figure file name)
		vars_to_plot_short_names.extend(['dv_dz', 'v_theta', 'v_rad', 'w', 'theta', 'kvm', 'tke', 'lscale', 'q'])
		vars_to_plot_short_names_all_times.append(vars_to_plot_short_names)
		
		# the long names for each variable (used in subplot titles)
		vars_to_plot_names.extend(['dv/dz', 'Tangential Velocity', 'Radial Velocity', 'Vertical Velocity', \
		'Potential Temperature', 'Vertical Diffusivity', 'Turbulent Kinetic Energy', 'Turbulent Length Scale', 'Specific Humidity'])
		vars_to_plot_names_all_times.append(vars_to_plot_names)
	
		# the units for each variable (used in horizontal axis label)
		vars_to_plot_units.extend(['1/s', 'm/s', 'm/s', 'm/s', 'degrees K', 'm**2/s', 'm**2/s', 'm', 'kg/kg'])
		vars_to_plot_units_all_times.append(vars_to_plot_units)
		
		# the ticks for each variable (also used to define the horizontal axis limits)
		vars_to_plot_ticks.extend([np.around(np.linspace(-0.02, 0.1, 7), decimals=2), np.arange(30, 110, 10), np.arange(-45, 10, 5), np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]), \
		np.arange(300, 312, 2), np.arange(0, 120, 20), np.arange(0, 50, 10), np.arange(0, 750, 150), np.array([0.015, 0.016, 0.017, 0.018, 0.019, 0.020])])
		vars_to_plot_ticks_all_times.append(vars_to_plot_ticks)
		
		# check to see if the number of panels equals the number of variables to plot
		if len(vars_to_plot) % num_panels != 0:
			print('The number of variables listed does not equal the number of panels plotted!\nTerminating program...')
			exit()
	
		# check to see if the number of panels is even
		if num_panels % num_cols != 0:
			print('The number of panels is not a multiple of the number of columns!\nTerminating program...')
			exit()

		# if continuing to this point, there are no errors as described above
		
	# record information for the maximum vertical level and all levels for the data
	max_level_all_configs.append(max_level)
	levels_all_configs.append(levels)

# Step 5 (plot the variables of interest for each desired time step)

for t_i, time_step in enumerate(time_steps): 

	f0 = plt.figure(figsize=(28,28))
			
	ax0 = f0.add_subplot(111)

	ax0.axis('off')
	pos = ax0.get_position() # obtain the position of the axis
	ax0.set_position([pos.x0, pos.y0, pos.width, pos.height*1.05]) # this will likely need to be adjusted accordingly
	ax0.set_title('Ensemble Average Model Output Variables\nat Time Step ' + time_step, fontsize=60)

	# the output file name (will be appended to)
	output_name = 'vertical_profiles_kepert_' + config + '_' + time_step + '_' + str(num_panels) + '_panels'

	for v_i, var in enumerate(vars_to_plot):

		# add to the output file name, which will eventually list all variables plotted
		output_name = output_name + '_' + vars_to_plot_short_names[v_i] 

		# for each variable, a panel needs to be made
		ax = f0.add_subplot(num_panels / num_cols, num_cols, v_i + 1)

		ax.set_title(vars_to_plot_names[v_i], fontsize=36)

		ax.tick_params(labelsize=24)
			
		ax.set_xlim(left=np.nanmin(vars_to_plot_ticks[v_i]), right=np.nanmax(vars_to_plot_ticks[v_i])) # the distance from the center of the cyclone
		ax.set_ylim(bottom=10., top=5000.) # the vertical distance above the surface 

		ax.set_xlabel(vars_to_plot_units[v_i], fontsize=28)
		
		if v_i == 0:
			ax.set_ylabel('Height (m)', fontsize=28)

		ax.set_xticks(vars_to_plot_ticks[v_i])
		ax.set_xticklabels(vars_to_plot_ticks[v_i])

		# the following depends on whether or not a log vertical axis is desired
		if plot_linear == False:
			ax.set_yscale('log')

			ax.set_yticks(np.array([10, 100, 1000, 5000], dtype='int'))
			ax.set_yticklabels(np.array([10, 100, 1000, 5000], dtype='int'))
		else:
			ax.set_ylim(bottom=0., top=1200.)
			
			ax.set_yticks(np.array([0, 200, 400, 600, 800, 1000, 1200], dtype='int'))
			ax.set_yticks(np.array([0, 200, 400, 600, 800, 1000, 1200], dtype='int'))
		
		# for each variable, plot the vertical profile for each configuration to be compared
		for c_i, config in enumerate(configs):
		
			# this is done to find the proper location of the configuration's data in the large array of all data
			config_index = c_i * len(time_steps) + t_i

			# plot the vertical profile for the variable
			plt.plot(vars_to_plot_all_times[config_index][v_i], levels_all_configs[c_i], linewidth=4, color=colors[c_i], label=config[-6:])
			
		# the following sequence of conditional statements plots various reference profiles (see Boolean definitions above)
		if vars_to_plot_short_names[v_i] == 'v_theta' and plot_vickery: # plots one set of observations from 2009
			vickery_input_path = '/Users/kmn182/ICS_work/CPT_Project/ANALYSIS/obs/'
			num_obs = 6
			obs_ids = ['2009_obs001', '2009_obs002', '2009_obs003', '2009_obs004', '2009_obs005', '2009_obs006']
			if len(obs_ids) != num_obs:
				print('The list of observation IDs is of length ' + str(len(obs_ids)) + '. It must be of length ' + str(num_obs) + '!')
			winds, heights = vickery(num_obs, obs_ids, vickery_input_path)
			for n in range(0, num_obs):
				plt.plot(winds[n], heights[n], linewidth=4, color='black')
		
		if vars_to_plot_short_names[v_i] == 'v_theta' and plot_NOAA: # plot the observed vertical tangential wind profiles
			profile_path = GOOGLE_ROOT + '/CPT_Project/ANALYSIS/dropsondes/data/2005/Katrina/ublox.qc.eol/P-3.43/D20050829_135021.044535124_PQC.eol.radazm.Wwind.npz'
			u_wind_profile, v_wind_profile, tan_wind_profile, rad_wind_profile, spd_profile, z_profile = get_NOAA_profiles(profile_path)
			plt.plot(tan_wind_profile, z_profile, linewidth=4, color='black')
		
		if vars_to_plot_short_names[v_i] == 'v_rad' and plot_NOAA: # plot the observed vertical tangential wind profiles
			profile_path = GOOGLE_ROOT + '/CPT_Project/ANALYSIS/dropsondes/data/2005/Katrina/ublox.qc.eol/P-3.43/D20050829_135021.044535124_PQC.eol.radazm.Wwind.npz'
			u_wind_profile, v_wind_profile, tan_wind_profile, rad_wind_profile, spd_profile, z_profile = get_NOAA_profiles(profile_path)
			plt.plot(rad_wind_profile, z_profile, linewidth=4, color='black')
			
		elif vars_to_plot_short_names[v_i] == 'v_rad' and plot_composite_NOAA: # plot the observed vertical tangential wind profiles
			plt.plot(rad_wind_profile, z_profile, linewidth=4, color='black')
			
		if vars_to_plot_short_names[v_i] == 'v_theta' and plot_composite_NOAA: # plot the observed vertical tangential wind profiles
			plt.plot(tan_wind_profile, z_profile, linewidth=4, color='black')
			
		elif vars_to_plot_short_names[v_i] == 'v_theta' and plot_composite_NASA: # plot the observed vertical tangential wind profiles
			u_wind_profile, v_wind_profile, spd_profile, z_profile, p_profile = composite_NASA_profiles(40., 80., storm_exceptions=True, chosen_storms=['Erika', 'Joaquin'])
			plt.plot(spd_profile, z_profile, linewidth=4, color='black')
			
		elif vars_to_plot_short_names[v_i] == 'v_theta' and plot_composite_P_3: # plot the observed vertical tangential wind profiles
			u_wind_profile, v_wind_profile, spd_profile, z_profile, p_profile = composite_P_3_profiles(0, 80.)
			plt.plot(spd_profile, z_profile, linewidth=4, color='black')
			
		elif vars_to_plot_short_names[v_i] == 'v_theta' and plot_composite_G_IV: # plot the observed vertical tangential wind profiles
			u_wind_profile, v_wind_profile, spd_profile, z_profile, p_profile = composite_G_IV_profiles(40., 80.)
			plt.plot(spd_profile, z_profile, linewidth=4, color='black')
			
		elif vars_to_plot_short_names[v_i] == 'v_theta' and plot_composite_NOAA_HRD: # plot the observed vertical tangential wind profiles
			u_wind_profile, v_wind_profile, spd_profile, z_profile, p_profile = composite_NOAA_HRD_profiles(40., 100.)
			plt.plot(spd_profile, z_profile, linewidth=4, color='black')
			
		if vars_to_plot_short_names[v_i] == 'theta' and plot_composite_NOAA: # plot the observed vertical tangential wind profiles
			plt.plot(pot_T_profile, z_profile, linewidth=4, color='black')
			
		if vars_to_plot_short_names[v_i] == 'w' and plot_composite_NOAA: # plot the observed vertical tangential wind profiles
			plt.plot(w_wind_profile, z_profile, linewidth=4, color='black')
		
		if vars_to_plot_short_names[v_i] == 'q' and plot_composite_NOAA: # plot the observed vertical tangential wind profiles
			plt.plot(q_profile, z_profile, linewidth=4, color='black')
		
		if vars_to_plot_short_names[v_i] == 'T' and plot_composite_NOAA: # plot the observed vertical tangential wind profiles
			plt.plot(T_profile, z_profile, linewidth=4, color='black')
			
		if vars_to_plot_short_names[v_i] == 'w' and plot_cm1: # plot CM1 vertical winds
			plt.plot(cm1_w, cm1_levs, linewidth=8, color='black', label='CM1')
			
		if vars_to_plot_short_names[v_i] == 'q' and plot_cm1: # plot CM1 specific humidity
			plt.plot(cm1_q, cm1_levs, linewidth=8, color='black', label='CM1')
		
		if vars_to_plot_short_names[v_i] == 'theta' and plot_cm1: # plot CM1 potential temperature
			plt.plot(cm1_theta, cm1_levs, linewidth=8, color='black', label='CM1')
			
		if vars_to_plot_short_names[v_i] == 'v_rad' and plot_cm1: # plot CM1 radial wind
			plt.plot(cm1_v_rad, cm1_levs, linewidth=8, color='black', label='CM1')
			
		if vars_to_plot_short_names[v_i] == 'v_theta' and plot_cm1: # plot CM1 tangential wind
			plt.plot(cm1_v_theta, cm1_levs, linewidth=8, color='black', label='CM1')
			
		if vars_to_plot_short_names[v_i] == 'rh' and plot_cm1: # plot CM1 relative humidity
			plt.plot(cm1_rh, cm1_levs, linewidth=8, color='black', label='CM1')
			
		if vars_to_plot_short_names[v_i] == 'dv_dz' and plot_cm1: # plot CM1 dv_dz
			plt.plot(cm1_dv_dz, cm1_levs, linewidth=8, color='black', label='CM1')
			
		# now, plot some values from the CM1 azimuthal average file
		if vars_to_plot_short_names[v_i] == 'kvm' and plot_cm1: # plot CM1 vertical diffusivity
			plt.plot(cm1_kvm, azim_cm1_levs, linewidth=8, color='black', label='CM1')
			
		if vars_to_plot_short_names[v_i] == 'tke' and plot_cm1: # plot CM1 turbulent kinetic energy
			plt.plot(cm1_tke, azim_cm1_levs, linewidth=8, color='black', label='CM1')
			
		if vars_to_plot_short_names[v_i] == 'lscale' and plot_cm1: # plot CM1 turbulence length scale
			plt.plot(cm1_leff, azim_cm1_levs, linewidth=8, color='black', label='CM1')
			
		# now, plot some betacast output
		if vars_to_plot_short_names[v_i] == 'v_rad' and plot_betacast: # plot betacast radial velocity
			plt.plot(beta_v_rad, beta_levs, linewidth=8, color='grey', label='betacast')
			
		if vars_to_plot_short_names[v_i] == 'v_theta' and plot_betacast: # plot betacast tangential velocity
			plt.plot(beta_v_theta, beta_levs, linewidth=8, color='grey', label='betacast')
			
		if vars_to_plot_short_names[v_i] == 'kvm' and plot_betacast: # plot betacast vertical diffusivity
			plt.plot(beta_kvm, beta_levs, linewidth=8, color='grey', label='betacast')
			
		if vars_to_plot_short_names[v_i] == 'lscale' and plot_betacast: # plot betacast turbulence length scale
			plt.plot(beta_lscale, beta_levs, linewidth=8, color='grey', label='betacast')
			
		if vars_to_plot_short_names[v_i] == 'dv_dz' and plot_betacast: # plot betacast dv_dz
			plt.plot(beta_dv_dz, beta_levs, linewidth=8, color='black', label='CM1')
			
		if vars_to_plot_short_names[v_i] == 'kvm' and plot_Zhang: # plot Zhang et al. (2012) vertical diffusivity
			plt.plot(zhang_kvm, zhang_levs, linewidth=8, color='green', label='Zhang')
		
		# plot the legend in the top right subplot
		if v_i == 2:	
			plt.legend(loc='lower right', fontsize=24)

	# Step 6 (save the figure)
	plt.savefig(output_path + output_name + '_multi_config.png', dpi=400, format='png')

	print('The figure file has been generated for time step ' + time_step + '!')








