#!/Users/kmn182/anaconda3/envs/ENV1/bin/python
# ------------------------------------------------------------
#
#	Program Name: generate_diagnostics_time_series_kepert_multi_ens.py
#
#	Purpose: To plot the time evolution of diagnostic model
#            output for multiple configurations.
#
#	Updated: 11:00 PM Monday October 14, 2019
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
from palettable.cubehelix import Cubehelix
from obs_functions import read_best_track, composite_NOAA_profiles
from output_functions import running_avg

# define environmental path variable
GOOGLE_ROOT = os.environ['GOOGLE_ROOT']

# define the font family
mpl.rc('font',family='Arial')

# define the configurations
num_ensembles = int(sys.argv[2])
config_base = sys.argv[1]

ensembles = range(1, num_ensembles+1)

# define a list of colors
custom_colormap = Cubehelix.make(start_hue=240, end_hue=-300, min_sat=1, max_sat=2.5, gamma=1, min_light=0.3, max_light=0.8, n=num_ensembles, reverse=False, name='custom_colormap')
colors = custom_colormap.hex_colors

# define the variables to be plotted
vars_to_plot_all_configs = []
vars_to_plot_short_names_all_configs = []
vars_to_plot_long_names_all_configs = []
vars_to_plot_units_all_configs = []
vars_to_plot_ticks_all_configs = []

# plot betacast profile?
plot_betacast = True
if plot_betacast:
	date = '2016-10-07-21600'
	config = 'TC-betacast.hindcast.Matthew.expb101'
	file = open(GOOGLE_ROOT + '/CPT_Project/ANALYSIS/diagnostics/' + config + '/diagnostics_betacasts_' + config + '.csv', mode='r')
	csv_reader = csv.reader(file, delimiter=',')
	beta_dates = [] # the time step (in days), first column
	beta_rads_max_wind = [] # radius of maximum wind (in m), second column
	beta_sfc_inflow_angs = [] # surface inflow angle (in degrees), third column
	beta_inflow_depths = [] # inflow depth (in m), fourth column
	beta_max_sfc_winds = [] # maximum surface wind speed (extrapolated to 10m, in m/s), fifth column
	beta_min_sfc_ps = [] # minimum surface pressure (in hPa), sixth column
	beta_height_max_winds = [] # height of maximum tangential wind (in m), seventh column
	
	skipped_cols = True # are there column skips?
	
	if skipped_cols:
		skip_factor = 2
	else:
		skip_factor = 1
	
	counter = 0 # this is a counter to determine whether or not we are analyzing the first row (which has headers)
	for row in csv_reader:
		if counter == 0:
			counter += 1
			continue
		else:
			beta_dates.append(str(row[skip_factor * 0]))
			beta_rads_max_wind.append(float(row[skip_factor * 1]))
			beta_sfc_inflow_angs.append(float(row[skip_factor * 2]))
			beta_inflow_depths.append(float(row[skip_factor * 3]))
			beta_max_sfc_winds.append(float(row[skip_factor * 4]))
			beta_min_sfc_ps.append(float(row[skip_factor * 5]))
			beta_height_max_winds.append(float(row[skip_factor * 6]))
			counter += 1
			
	beta_dates = np.asarray(beta_dates, dtype='str')
	
	date_i = np.where(beta_dates == date)[0][0]
	beta_rad_max_wind = beta_rads_max_wind[date_i] 
	beta_sfc_inflow_ang = beta_sfc_inflow_angs[date_i]
	beta_inflow_depth = beta_inflow_depths[date_i]
	beta_max_sfc_wind = beta_max_sfc_winds[date_i]
	beta_min_sfc_p = beta_min_sfc_ps[date_i]
	beta_height_max_wind = beta_height_max_winds[date_i]

# plot CM1 profile? 
plot_cm1 = True
if plot_cm1:
	date = '000120'
	config = 'cm1out'
	file = open(GOOGLE_ROOT + '/CPT_Project/ANALYSIS/diagnostics/' + config + '/diagnostics_cm1_' + config + '.csv', 'r')
	csv_reader = csv.reader(file, delimiter=',')
	cm1_dates = [] # the time step (in days), first column
	cm1_rads_max_wind = [] # radius of maximum wind (in m), second column
	cm1_sfc_inflow_angs = [] # surface inflow angle (in degrees), third column
	cm1_inflow_depths = [] # inflow depth (in m), fourth column
	cm1_max_sfc_winds = [] # maximum surface wind speed (extrapolated to 10m, in m/s), fifth column
	cm1_min_sfc_ps = [] # minimum surface pressure (in hPa), sixth column
	cm1_height_max_winds = [] # height of maximum tangential wind (in m), seventh column
	
	skipped_cols = True # are there column skips?
	
	if skipped_cols:
		skip_factor = 2
	else:
		skip_factor = 1
	
	counter = 0 # this is a counter to determine whether or not we are analyzing the first row (which has headers)
	for row in csv_reader:
		if counter == 0:
			counter += 1
			continue
		else:
			cm1_dates.append(str(row[skip_factor * 0]))
			cm1_rads_max_wind.append(float(row[skip_factor * 1]))
			cm1_sfc_inflow_angs.append(float(row[skip_factor * 2]))
			cm1_inflow_depths.append(float(row[skip_factor * 3]))
			cm1_max_sfc_winds.append(float(row[skip_factor * 4]))
			cm1_min_sfc_ps.append(float(row[skip_factor * 5]))
			cm1_height_max_winds.append(float(row[skip_factor * 6]))
			counter += 1
	
	cm1_dates = np.asarray(cm1_dates, dtype='str')
	
	date_i = np.where(cm1_dates == date)[0][0]
	cm1_rad_max_wind = cm1_rads_max_wind[date_i] 
	cm1_sfc_inflow_ang = cm1_sfc_inflow_angs[date_i]
	cm1_inflow_depth = cm1_inflow_depths[date_i]
	cm1_max_sfc_wind = cm1_max_sfc_winds[date_i]
	cm1_min_sfc_p = cm1_min_sfc_ps[date_i]
	cm1_height_max_wind = cm1_height_max_winds[date_i]

# Step 2 (read in the input files)

for e_i, ens in enumerate(ensembles):

	config = config_base + '.' + str(ens).rjust(3, '0')

	# define the output path
	output_path = GOOGLE_ROOT + '/CPT_Project/FIGURES/multi_config/'
	
	# define the input path
	input_path = GOOGLE_ROOT + '/CPT_Project/ANALYSIS/diagnostics/' + config + '/'

	f = open(input_path + 'diagnostics_kepert_' + config + '_ens_avg.csv', mode='r')
	csv_reader = csv.reader(f, delimiter=',')
	
	time_step = [] # the time step (in days), first column
	rad_max_wind = [] # radius of maximum wind (in m), second column
	sfc_inflow_ang = [] # surface inflow angle (in degrees), third column
	inflow_depth = [] # inflow depth (in m), fourth column
	max_sfc_wind = [] # maximum surface wind speed (extrapolated to 10m, in m/s), fifth column
	min_sfc_p = [] # minimum surface pressure (in hPa), sixth column
	height_max_wind = [] # height of maximum tangential wind (in m), seventh column
	
	skipped_cols = True # are there column skips?
	
	if skipped_cols:
		skip_factor = 2
	else:
		skip_factor = 1
	
	counter = 0 # this is a counter to determine whether or not we are analyzing the first row (which has headers)
	for row in csv_reader:
		if counter == 0:
			counter += 1
			continue
		else:
			time_step.append(int(row[skip_factor * 0]))
			rad_max_wind.append(float(row[skip_factor * 1]))
			sfc_inflow_ang.append(float(row[skip_factor * 2]))
			inflow_depth.append(float(row[skip_factor * 3]))
			max_sfc_wind.append(float(row[skip_factor * 4]))
			min_sfc_p.append(float(row[skip_factor * 5]))
			height_max_wind.append(float(row[skip_factor * 6]))
			counter += 1
			
	# convert the lists into numpy arrays of the appropriate type
	time_step = np.asarray(time_step, dtype='int')
	rad_max_wind = np.asarray(rad_max_wind, dtype='float')
	sfc_inflow_ang = np.asarray(sfc_inflow_ang, dtype='float')
	inflow_depth = np.asarray(inflow_depth, dtype='float')
	max_sfc_wind = np.asarray(max_sfc_wind, dtype='float')
	min_sfc_p = np.asarray(min_sfc_p, dtype='float')
	height_max_wind = np.asarray(height_max_wind, dtype='float')
	
	max_sfc_wind = max_sfc_wind * (0.94 / 0.85)
	
	# convert these arrays to weighted averages
	num_days = 3
	rad_max_wind = running_avg(rad_max_wind, num_days)
	sfc_inflow_ang = running_avg(sfc_inflow_ang, num_days)
	inflow_depth = running_avg(inflow_depth, num_days)
	max_sfc_wind = running_avg(max_sfc_wind, num_days)
	min_sfc_p = running_avg(min_sfc_p, num_days)
	height_max_wind = running_avg(height_max_wind, num_days)
	
	vars_to_plot = [rad_max_wind, sfc_inflow_ang, inflow_depth, max_sfc_wind, min_sfc_p, height_max_wind]
	vars_to_plot_all_configs.extend(vars_to_plot)
	vars_to_plot_short_names = ['rad_max_wind', 'sfc_inflow_ang', 'inflow_depth', 'max_sfc_wind', 'p_sfc_min', 'height_max_wind']
	vars_to_plot_short_names_all_configs.extend(vars_to_plot_short_names)
	vars_to_plot_long_names = ['Radius of Maximum Wind', 'Surface Inflow Angle', 'Inflow Depth', 'Maximum Surface Wind Speed', 'Minimum Surface Pressure', 'Height of Maximum Wind']
	vars_to_plot_long_names_all_configs.extend(vars_to_plot_long_names)
	vars_to_plot_units = ['m/s', 'degrees', 'm', 'm/s', 'hPa', 'm']
	vars_to_plot_units_all_configs.extend(vars_to_plot_units)
	vars_to_plot_ticks = [np.arange(0, 350, 50), np.arange(0, 50, 10), np.arange(0, 6000, 1000), np.arange(0, 110, 10), np.arange(860, 1020, 20), np.arange(0, 3500, 500)]
	vars_to_plot_ticks_all_configs.extend(vars_to_plot_ticks)
	
# Step 3 (plot the time series of the diagnostic quantity)

# first, calculate the vertical profile data
u_wind_profile, v_wind_profile, w_wind_profile, tan_wind_profile, rad_wind_profile, spd_profile, z_profile, p_profile, T_profile, q_profile, pot_T_profile, max_heights, inflow_angles, inflow_depths = composite_NOAA_profiles(0., 1000., 40., 80.)

for v_i, var in enumerate(vars_to_plot): # loop through each variable to plot

	f0 = plt.figure(figsize=(20,20))
	
	plt.rc('xtick',labelsize=32)
	plt.rc('ytick',labelsize=32)
	
	# v_i index can be used for ticks and labels because these are all repeated for each configuration
	
	file_name = 'diagnostics_time_series_' + config_base + '_' + vars_to_plot_short_names[v_i] + '_multi_ens.png'
	
	plt.title('Ensemble Average Time Evolution\nof ' + vars_to_plot_long_names[v_i], fontsize=42)
	plt.xlim(np.nanmin(time_step), np.nanmax(time_step))
	plt.ylim(np.nanmin(vars_to_plot_ticks[v_i]), np.nanmax(vars_to_plot_ticks[v_i]))
	plt.xlabel('Time (days)', fontsize=32)
	plt.ylabel(vars_to_plot_long_names[v_i] + ' (' + vars_to_plot_units[v_i] + ')', fontsize=32)
	plt.xticks(np.arange(np.nanmin(time_step), np.nanmax(time_step) + 1, 1), np.arange(np.nanmin(time_step), np.nanmax(time_step) + 1, 1))
	plt.yticks(vars_to_plot_ticks[v_i], vars_to_plot_ticks[v_i])
	
	for c_i, config in enumerate(ensembles):
	
		config_index = int(c_i * int(len(vars_to_plot_all_configs)) / int(num_ensembles) + v_i)
	
		plt.plot(time_step, vars_to_plot_all_configs[config_index], linewidth=6, color=colors[c_i], label=str(config).rjust(3,'0'))
		
	if vars_to_plot_long_names[v_i] == 'Height of Maximum Wind':
		pct50 = np.nanmedian(max_heights)
		pct25 = np.nanpercentile(max_heights, 25)
		pct75 = np.nanpercentile(max_heights, 75)
		# min = np.nanmin(max_heights)
# 		max = np.nanmax(max_heights)
		plt.plot(time_step, np.ones([len(time_step)]) * pct50, linewidth=6, color='black')
		plt.plot(time_step, np.ones([len(time_step)]) * pct25, linewidth=6, color='black')
		plt.plot(time_step, np.ones([len(time_step)]) * pct75, linewidth=6, color='black')
		plt.plot(time_step, np.ones([len(time_step)]) * beta_height_max_wind, '--', linewidth=6, color='red', label='betacast')
		plt.plot(time_step, np.ones([len(time_step)]) * cm1_height_max_wind, linewidth=6, color='red', label='CM1')
		# plt.plot(time_step, np.ones([len(time_step)]) * min, linewidth=6, color='red')
# 		plt.plot(time_step, np.ones([len(time_step)]) * max, linewidth=6, color='red')
		plt.fill_between(time_step, pct25, pct50, facecolor='grey')
		plt.fill_between(time_step, pct50, pct75, facecolor='grey')
	elif vars_to_plot_long_names[v_i] == 'Surface Inflow Angle':
		pct50 = np.nanmedian(inflow_angles)
		pct25 = np.nanpercentile(inflow_angles, 25)
		pct75 = np.nanpercentile(inflow_angles, 75)
		# min = np.nanmin(inflow_angles)
# 		max = np.nanmax(inflow_angles)
# 		print(min, max)
		plt.plot(time_step, np.ones([len(time_step)]) * pct50, linewidth=6, color='black')
		plt.plot(time_step, np.ones([len(time_step)]) * pct25, linewidth=6, color='black')
		plt.plot(time_step, np.ones([len(time_step)]) * pct75, linewidth=6, color='black')
		plt.plot(time_step, np.ones([len(time_step)]) * beta_sfc_inflow_ang, '--', linewidth=6, color='red', label='betacast')
		plt.plot(time_step, np.ones([len(time_step)]) * cm1_sfc_inflow_ang, linewidth=6, color='red', label='CM1')
		# plt.plot(time_step, np.ones([len(time_step)]) * min, linewidth=6, color='red')
# 		plt.plot(time_step, np.ones([len(time_step)]) * max, linewidth=6, color='red')
		plt.fill_between(time_step, pct25, pct50, facecolor='grey')
		plt.fill_between(time_step, pct50, pct75, facecolor='grey')
	elif vars_to_plot_long_names[v_i] == 'Inflow Depth':
		pct50 = np.nanmedian(inflow_depths)
		pct25 = np.nanpercentile(inflow_depths, 25)
		pct75 = np.nanpercentile(inflow_depths, 75)
		# min = np.nanmin(inflow_depths)
# 		max = np.nanmax(inflow_depths)
# 		print(min, max)
		plt.plot(time_step, np.ones([len(time_step)]) * pct50, linewidth=6, color='black')
		plt.plot(time_step, np.ones([len(time_step)]) * pct25, linewidth=6, color='black')
		plt.plot(time_step, np.ones([len(time_step)]) * pct75, linewidth=6, color='black')
		plt.plot(time_step, np.ones([len(time_step)]) * beta_inflow_depth, '--', linewidth=6, color='red', label='betacast')
		plt.plot(time_step, np.ones([len(time_step)]) * cm1_inflow_depth, linewidth=6, color='red', label='CM1')
	# 	plt.plot(time_step, np.ones([len(time_step)]) * min, linewidth=6, color='red')
# 		plt.plot(time_step, np.ones([len(time_step)]) * max, linewidth=6, color='red')
		plt.fill_between(time_step, pct25, pct50, facecolor='grey')
		plt.fill_between(time_step, pct50, pct75, facecolor='grey')
		
	plt.legend(loc='upper right', fontsize=24, ncol=4)
	
	plt.savefig(output_path + file_name, dpi=400, format='png')
	
	plt.close(f0)
	
	print('The figure has been generated for ' + vars_to_plot_long_names[v_i] + '!')
	
	# Step 4 (plot the wind vs. pressure)
	
	# if the variable to plot is surface pressure, also plot wind speed as a function of pressure
	if vars_to_plot_short_names[v_i] == 'p_sfc_min':
		
		f1 = plt.figure(figsize=(20,20))
	
		plt.rc('xtick',labelsize=32)
		plt.rc('ytick',labelsize=32)
	
		# v_i index can be used for ticks and labels because these are all repeated for each configuration
	
		file_name = 'diagnostics_time_series_' + config_base + '_pressure_wind_multi_config.png'
	
		plt.title('Ensemble Average Wind vs. Pressure', fontsize=42)
		plt.xlim(np.nanmin(vars_to_plot_ticks[v_i]), np.nanmax(vars_to_plot_ticks[v_i])) # be careful that the prior index corresponds to wind!
		plt.ylim(np.nanmin(vars_to_plot_ticks[v_i-1]), np.nanmax(vars_to_plot_ticks[v_i-1]))
		plt.xlabel(vars_to_plot_long_names[v_i] + ' (' + vars_to_plot_units[v_i] + ')', fontsize=32)
		plt.ylabel(vars_to_plot_long_names[v_i-1] + ' (' + vars_to_plot_units[v_i-1] + ')', fontsize=32)
		plt.xticks(vars_to_plot_ticks[v_i], vars_to_plot_ticks[v_i])
		plt.yticks(vars_to_plot_ticks[v_i-1], vars_to_plot_ticks[v_i-1])
		
		for c_i, config in enumerate(ensembles):
	
			config_index = int(c_i * int(len(vars_to_plot_all_configs)) / int(num_ensembles) + v_i)
			
			# read-in the best-track pressure and wind observations
			txt_file_name = '/Users/kmn182/Documents/CPT_Project/ANALYSIS/obs/NHC_best_track/hurdat2-1851-2018-120319.txt'
			max_spd, min_ps = read_best_track(txt_file_name)
			
			if c_i == 0:
				plt.plot(min_ps, max_spd, '.', color='green', markersize=10, alpha=0.5, label='obs')
				plt.plot(cm1_min_sfc_ps, cm1_max_sfc_winds, '.', markersize=30, color='black', label='CM1')
				plt.plot(beta_min_sfc_ps, beta_max_sfc_winds, '.', markersize=30, color='grey', label='betacast')
		
			plt.plot(vars_to_plot_all_configs[config_index], vars_to_plot_all_configs[config_index-1], '.', markersize=30, color=colors[c_i], label=str(config).rjust(3,'0'))
			
		plt.legend(loc='upper right', fontsize=24, ncol=4)
	
		plt.savefig(output_path + file_name, dpi=400, format='png')
	
		plt.close(f1)
	
		print('The figure has been generated for ' + vars_to_plot_long_names[v_i] + '!')
	
	
	
	
	
	
	
	
	
	 