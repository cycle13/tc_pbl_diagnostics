#!/Users/kmn182/anaconda3/envs/ENV1/bin/python
# ------------------------------------------------------------
#
#	Program Name: output_functions.py
#
#	Purpose: Functions that aide in plotting model output
#
#	Updated: 2:00 PM Thursday December 5, 2019
#
#	WARNING: This code assumes Python 3, so there may be some 
#            syntax errors if running with Python 2.
#
# -------------------------------------------------------------

import numpy as np
from netCDF4 import Dataset
import os
import sys
import types
from scipy import ndimage

# "center_cyclone" function
# inputs:
# initial_cyclone (floats): the 2D array of the cyclone (any variable)
# lats (floats): the 1D array of latitudes
# lons (floats): the 1D array of longitudes
# outputs:
# centered_cyclone (floats): the 2D array representing the cyclone
# this assumes that the origin of the plane is (0,0)

def center_cyclone( ps, initial_cyclone, lats, lons ):
	
	# import the necessary modules
	from netCDF4 import Dataset
	import csv
	import numpy as np
	from scipy import ndimage
	
	# find the location of the variable field's centroid
	i,j = ndimage.measurements.center_of_mass(initial_cyclone)
	min_ps = np.nanmin(ps)
	# min_lat_index = int(np.around(i)) # a bit of a misnomer
# 	min_lon_index = int(np.around(j))
	min_lat_index = np.where(ps == min_ps)[0][0]
	min_lon_index = np.where(ps == min_ps)[1][0]
	zero_lat_index = np.where(lats==0.)[0][0]
	zero_lon_index = np.where(lons==0.)[0][0]
	lat_index_offset = -1 * (min_lat_index - zero_lat_index)
	lon_index_offset = -1 * (min_lon_index - zero_lon_index)
	
	# shift the initial cyclone's data to the new centered array
	centered_cyclone = np.ones([len(lats), len(lons)]) * np.nan
	for lat_i, lat in enumerate(lats):
		for lon_i, lon in enumerate(lons):
			if lat_i - lat_index_offset >= len(lats) or lon_i - lon_index_offset >= len(lons) or lat_i - lat_index_offset < 0 or lon_i - lon_index_offset < 0:
				continue
			centered_cyclone[lat_i, lon_i] = initial_cyclone[lat_i - lat_index_offset, lon_i - lon_index_offset] # note that boundaries may have NaN values
	
	# return the centered cyclone
	return centered_cyclone
	
# calc_p_sfc_ind_ens function
# --------------------
# inputs:
# config = string defining the configuration
# num_ensembles = integer number of ensembles
# --------------------
# outputs:
# p_sfc_min = time series of minimum surface pressure (in hPa), with an added dimension of ensemble members
# sfc_wind_max = time series of maximum surface wind (in m/s), with an added dimension of ensemble members
# --------------------
def calc_p_sfc_ind_ens(config, num_ensembles):

	# read in the necessary modules
	import numpy as np
	from netCDF4 import Dataset
	
	times = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16']
	
	ensembles = range(1, num_ensembles + 1)
	
	p_sfc_min = np.ones([len(times), num_ensembles]) * np.nan
	sfc_wind_max = np.ones([len(times), num_ensembles]) * np.nan
	
	for t_i, time in enumerate(times):
		for e_i, ensemble in enumerate(ensembles):
		
			# define the file path
			input_path = '/Users/kmn182/ICS_scratch/output/' + config + '.' + str(ensemble).rjust(3, '0') + '/run'
			file_name = config + '.' + str(ensemble).rjust(3, '0') + '.cam.h1.0001-01-' + time + '-00000.nc_regrid.nc'
			full_path = input_path + '/' + file_name
		
			# read in the surface pressure data from the netCDF file
			dataset = Dataset(full_path)
			p_sfc = dataset['PS'][0, :, :] # time x lat x lon # take the first time step 
			u = dataset['U'][0, -1, :, :] # time x level x lat x lon # take the first time step and lowest model level
			v = dataset['V'][0, -1, :, :] 
			spd = np.sqrt(u**2 + v**2)
			spd = spd * 0.85
		
			p_sfc_min[t_i, e_i] = np.nanmin(p_sfc) / 100. # the minimum surface pressure (converted from Pa to hPa)
			sfc_wind_max[t_i, e_i] = np.nanmax(spd) # the maximum surface wind speed
			
			dataset.close()
	
	# return the minimum surface pressure
	return p_sfc_min, sfc_wind_max
	
# calc_p_sfc_ens_avg function
# --------------------
# inputs:
# config = string defining the configuration
# num_ensembles = integer number of ensembles
# --------------------
# outputs:
# p_sfc_min = time series of minimum surface pressure (in hPa)...ensemble averaged
# sfc_wind_max = time series of maximum surface wind (in m/s)...ensemble averaged
# --------------------
def calc_p_sfc_ens_avg(config, num_ensembles):

	# read in the necessary modules
	import numpy as np
	from netCDF4 import Dataset
	
	times = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16']
	
	letters = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z']
	ensembles = letters[0:int(num_ensembles)]
	
	p_sfc_min = np.ones([len(times), num_ensembles]) * np.nan
	sfc_wind_max = np.ones([len(times), num_ensembles]) * np.nan
	
	for t_i, time in enumerate(times):
		for e_i, ensemble in enumerate(ensembles):
		
			# define the file path
			input_path = '/Users/kmn182/ICS_scratch/output/' + config + str(ensemble) + '/run'
			file_name = config + str(ensemble) + '.cam.h1.0001-01-' + time + '-00000.nc_regrid.nc'
			full_path = input_path + '/' + file_name
		
			# read in the surface pressure data from the netCDF file
			if config == 'RCE.QPC6.ne0np4tcfplane.ne15x8.exp712' and str(ensemble) == '16':
				continue
			print(config, str(ensemble))
			dataset = Dataset(full_path)
			dims = dataset['PS'][:].shape
			if dims[0] == 0:
				p_sfc = np.nan
				u = np.nan
				v = np.nan
			else:
				p_sfc = dataset['PS'][0, :, :] # time x lat x lon # take the first time step 
				u = dataset['U'][0, -1, :, :] # time x level x lat x lon # take the first time step and lowest model level
				v = dataset['V'][0, -1, :, :] 
			spd = np.sqrt(u**2 + v**2)
			spd = spd * 0.85 # this is an adjustment to account for moving from 60m to 10m
			
			print(t_i, e_i)
		
			p_sfc_min[t_i, e_i] = np.nanmin(p_sfc) / 100. # the minimum surface pressure (converted from Pa to hPa)
			sfc_wind_max[t_i, e_i] = np.nanmax(spd) # the maximum surface wind speed
			
			dataset.close()
			
	# calculate the ensemble average
	p_sfc_min = np.nanmean(p_sfc_min, axis=1)
	sfc_wind_max = np.nanmean(sfc_wind_max, axis=1)
	
	# return the minimum surface pressure
	return p_sfc_min, sfc_wind_max
	
# calc_p_sfc_cm1 function
# --------------------
# inputs:
# config = string defining the configuration
# times = string list of time labels
# --------------------
# outputs:
# p_sfc_min = time series of minimum surface pressure (in hPa)...ensemble averaged
# sfc_wind_max = time series of maximum surface wind (in m/s)...ensemble averaged
# --------------------
def calc_p_sfc_cm1(config, times):

	# read in the necessary modules
	import numpy as np
	from netCDF4 import Dataset
	
	p_sfc_min = np.ones([len(times)]) * np.nan
	sfc_wind_max = np.ones([len(times)]) * np.nan
	
	for t_i, time in enumerate(times):
		
		# define the file path
		input_path = '/Users/kmn182/ICS_scratch/CM1_data'
		file_name = config + '_' + time + '.nc'
		full_path = input_path + '/' + file_name
	
		# read in the surface pressure data from the netCDF file
		dataset = Dataset(full_path)
		p_sfc = dataset['psfc'][:, :] # lat x lon # take the first time step 
		u10 = dataset['u10'][:, :] # lat x lon # take the first time step and lowest model level
		v10 = dataset['v10'][:, :] 
		spd = np.sqrt(u10**2 + v10**2)
	
		p_sfc_min[t_i] = np.nanmin(p_sfc) / 100. # the minimum surface pressure (converted from Pa to hPa)
		sfc_wind_max[t_i] = np.nanmax(spd) # the maximum surface wind speed
		
		dataset.close()
	
	# return the minimum surface pressure
	return p_sfc_min, sfc_wind_max
	
# calc_p_sfc_betacast function
# --------------------
# inputs:
# config = string defining the configuration
# times = string list of time labels
# --------------------
# outputs:
# p_sfc_min = time series of minimum surface pressure (in hPa)...ensemble averaged
# sfc_wind_max = time series of maximum surface wind (in m/s)...ensemble averaged
# --------------------
def calc_p_sfc_betacast(config, times):

	# read in the necessary modules
	import numpy as np
	from netCDF4 import Dataset
	
	p_sfc_min = np.ones([len(times)]) * np.nan
	sfc_wind_max = np.ones([len(times)]) * np.nan
	
	for t_i, time in enumerate(times):
		
		# define the file path
		input_path = '/Users/kmn182/ICS_scratch/betacasts/output_regridded/' + config
		file_name = config + '.cam.h0.' + time + '.nc_regrid_centered.nc'
		full_path = input_path + '/' + file_name
	
		# read in the surface pressure data from the netCDF file
		dataset = Dataset(full_path)
		p_sfc = dataset['PSL'][:, :] # lat x lon # take the first time step 
		spd = dataset['U10'][:, :]
	
		p_sfc_min[t_i] = np.nanmin(p_sfc) / 100. # the minimum surface pressure (converted from Pa to hPa)
		sfc_wind_max[t_i] = np.nanmax(spd) # the maximum surface wind speed
		
		dataset.close()
	
	# return the minimum surface pressure
	return p_sfc_min, sfc_wind_max
	
# running_avg function
# --------------------
# inputs:
# time_series = array of strings representing a time series of the given variable
# num_days = integere number of days to perform the centered running average
# --------------------
# outputs:
# time_series_running_avg = the time series with a centered running average applied
# --------------------
def running_avg(time_series, num_days):

	# verify that the number of days chosen for the running average is not too large
	if num_days > len(time_series)/3.:
		print('The chosen number of days is too large!')
		exit()
		
	if num_days%2 == 0:
		print('The chosen number of days is even! Centered averaging not possible!')
		exit()
		
	# create the time series running average array
	time_series_running_avg = np.ones([len(time_series)]) * np.nan
		
	# loop through each element and calculate an x-day running average
	for v_i, val in enumerate(time_series):
		if v_i >= (num_days-1)/2 and v_i <= len(time_series) - 1 - (num_days-1)/2: # centered average
			time_series_running_avg[v_i] = np.nanmean(time_series[int(v_i - (num_days-1)/2):int(v_i + (num_days-1)/2 + 1)])
		elif v_i < (num_days-1)/2: # forward-skewed average
			#time_series_running_avg[v_i] = np.nanmean(time_series[0:num_days])
			time_series_running_avg[v_i] = time_series[v_i]
		elif v_i > len(time_series) - 1 - (num_days-1)/2: # backward-skewed average
			#time_series_running_avg[v_i] = np.nanmean(time_series[-1*num_days:])
			time_series_running_avg[v_i] = time_series[v_i]

	# output the running average of the time series
	return time_series_running_avg