#!/Users/kmn182/anaconda3/envs/ENV1/bin/python
# ------------------------------------------------------------
#
#	Program Name: obs_functions
#
#	Purpose: Functions that aide in plotting observed vertical wind profiles
#
#	Updated: 2:00 PM Tuesday December 3, 2019
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

# "vickery" function
# inputs:
# num_obs (integer): the number of observed profiles
# obs_ids (list of strings): the IDs of the observation (typically a time stamp or a time stamp with an observation number
# input_path (string): the input path for the profiles
# outputs:
# wind_profiles (list of arrays of floats): a list of wind profiles for each observation
# height_profiles (list of arrays of floats): a list of height profiles for each observation
# notes: 
# assumes input is CSV file of the following type "vickery" + obs_id + ".csv"
# also assumes two columns (row 0 = tangential wind, row 1 = height) with headers
def vickery( num_obs, obs_ids, input_path ):
	# import the necessary modules
	import csv
	import numpy as np
	
	# define the lists
	wind_profiles = []
	height_profiles = []
	time_stamps = []
	
	# read-in the vertical profiles
	for o_i, obs in enumerate(obs_ids):
		file_path = input_path + 'vickery' + obs + '.csv'
		csv_file = open(file_path, mode='r')
		csv_reader = csv.reader(csv_file, delimiter=',')
		wind_profile = []
		height_profile = []
		counter = 0
		for row in csv_reader:
			if len(row) != 2:
				break
			if counter == 0:
				counter += 1 
				continue
			else:
				counter += 1
				wind_profile.append(float(row[0]))
				height_profile.append(float(row[1]))
		wind_profile = np.asarray(wind_profile, dtype='float')
		height_profile = np.asarray(height_profile, dtype='float')
		
		# add the vertical profiles to the lists of vertical profiles
		wind_profiles.append(wind_profile)
		height_profiles.append(height_profile)
	
	# return the complete lists of profiles to be plotted
	return wind_profiles, height_profiles
	
# "get_NOAA_profiles" function
# inputs:
# profile_path = a string containing the path to the NOAA vertical wind profile
# outputs:
# u_wind_profile = a numpy array of floats containing the vertical profile of zonal wind
# v_wind_profile = a numpy array of floats containing the vertical profile of meridional wind
# tan_wind_profile = a numpy array of floats containing the vertical profile of tangential wind
# rad_wind_profile = a numpy array of floats containing the vertical profile of radial wind
# spd_profile = a numpy array of floats containing the vertical profile of total wind speed
# z_profile = a numpy array of floats containing the vertical height profile
# notes:
# the profile path must be to a numpy binary file
# otherwise, an error will occur
def get_NOAA_profiles(profile_path):
	
	# read in the numpy binary file
	file = np.load(profile_path)
	
	# extract the relevant profiles
	u_wind_profile = file['u_wind']
	v_wind_profile = file['v_wind']
	tan_wind_profile = file['tan_wind']
	rad_wind_profile = file['rad_wind']
	spd_profile = file['wind_speed']
	z_profile = file['z']
	
	# return the complete lists of profiles to be plotted
	return u_wind_profile, v_wind_profile, tan_wind_profile, rad_wind_profile, spd_profile, z_profile
	
# "get_NASA_profiles" function
# inputs:
# profile_path = a string containing the path to the NOAA vertical wind profile
# outputs:
# u_wind_profile = a numpy array of floats containing the vertical profile of zonal wind
# v_wind_profile = a numpy array of floats containing the vertical profile of meridional wind
# spd_profile = a numpy array of floats containing the vertical profile of total wind speed
# z_profile = a numpy array of floats containing the vertical height profile
# notes:
# the profile path must be to a numpy binary file
# otherwise, an error will occur
def get_NASA_profiles(profile_path):
	
	# read in the numpy binary file
	file = np.load(profile_path)
	
	# extract the relevant profiles
	u_wind_profile = file['u_wind']
	v_wind_profile = file['v_wind']
	spd_profile = file['wind_speed']
	z_profile = file['z']
	
	# return the complete lists of profiles to be plotted
	return u_wind_profile, v_wind_profile, spd_profile, z_profile
	
# "composite_NOAA_profiles" function
# inputs:
# min_radius = a float representing the minimum radius for the dropsonde profile (km))
# max_radius = a float representing the maximum radius for the dropsonde profile (km)
# min_max_speed = a float representing the minimum maximum wind speed for the dropsonde profile (m/s)
# max_max_speed = a float representing the maximum maximum wind speed for the dropsonde profile (m/x)
# year_exceptions = a boolean representing whether or not subsetting by year will occur (optional, default is False)
# storm_exceptions = a boolean representing whether or not subsetting by storm will occur (optional, default is False)
# chosen_years = if year_exceptions is true, an array of strings containing the years to include in the composite (default is None)
# chosen_storms = if storm_exceptions is true, an array of strings containing the named storms to include in the composite (default is None)
# outputs:
# u_wind_profile = a numpy array of floats containing the composite vertical profile of zonal wind (m/s)
# v_wind_profile = a numpy array of floats containing the composite vertical profile of meridional wind (m/s)
# tan_wind_profile = a numpy array of floats containing the composite vertical profile of tangential wind
# rad_wind_profile = a numpy array of floats containing the composite vertical profile of radial wind
# spd_profile = a numpy array of floats containing the composite vertical profile of total wind speed
# z_profile = a numpy array of floats containing the height (m)
# p_profile = a numpy array of floats containing the composite vertical pressure profile (hPa)
# T_profile = a numpy array of floats containing the composite vertical temperature profile (K)
# q_profile = a numpy array of floats containing the specific humidity profile (kg/kg)
# pot_T_profile = a numpy array of floats containing the potential temperature profile (K)
# notes:
# the profile path must be to a numpy binary file
# otherwise, an error will occur
def composite_NOAA_profiles(min_radius, max_radius, min_max_speed, max_max_speed, year_exceptions=False, storm_exceptions=False, chosen_years=None, chosen_storms=None):

	# load any necessary additional modules
	import os
	import sys
	import types
	from output_functions import running_avg
	import metpy.calc as calc
	
	# define environmental path variable
	GOOGLE_ROOT = os.environ['GOOGLE_ROOT']
	
	# define the file path
	path_first_part = GOOGLE_ROOT + '/CPT_Project/ANALYSIS/obs/dropsondes_1996-2012/'

	years = range(1996, 2013) # the years on record
	
	# define the uniform pressures over which interpolation will occur
	p_interp = np.linspace(1000., 500., 501)
	
	# define the summation array and the counter (for averaging purposes)
	u_sum = np.zeros([len(p_interp)]) 
	v_sum = np.zeros([len(p_interp)]) 
	v_tan_sum = np.zeros([len(p_interp)]) 
	v_rad_sum = np.zeros([len(p_interp)])
	spd_sum = np.zeros([len(p_interp)])
	z_sum = np.zeros([len(p_interp)])
	w_sum = np.zeros([len(p_interp)])
	T_sum = np.zeros([len(p_interp)])
	q_sum = np.zeros([len(p_interp)])
	pot_T_sum = np.zeros([len(p_interp)])
	counts = np.zeros([len(p_interp)])
	
	# create a list of height maxima, inflow angle, and inflow depth
	max_heights = []
	inflow_angles = []
	inflow_depths = []

	# loop through the subdirectories to obtain all relevant dropsonde profiles

	for y_i, year in enumerate(years):
	
		if year_exceptions:
			if year not in chosen_years:
				continue
				
		#print('Moving forward with ' + str(year))
	
		updated_path = path_first_part + str(year)
	
		# obtain the names of the storms (this is the next file level)
		storms = os.listdir(updated_path)
	
		for s_i, storm in enumerate(storms):
	
			if storm == '.DS_Store':
				continue
				
			if storm_exceptions:
				if storm not in chosen_storms:
					continue
	
			updated_path = path_first_part + str(year) + '/' + storm
		
			# obtain the sensor types 
			sensor_types = os.listdir(updated_path)
		
			for s_t_i, sensor_type in enumerate(sensor_types):
		
				if sensor_type == '.DS_Store':
					continue
		
				updated_path = path_first_part + str(year) + '/' + storm + '/' + sensor_type
			
				# obtain the aircraft type
				aircraft_types = os.listdir(updated_path)
			
				for a_t_i, aircraft_type in enumerate(aircraft_types):
			
					if aircraft_type == '.DS_Store':
						continue
					
					updated_path = path_first_part + str(year) + '/' + storm + '/' + sensor_type + '/' + aircraft_type
				
					# obtain the individual dropsonde files
					files = os.listdir(updated_path)
				
					for f_i, file in enumerate(files):
					
						# Step 4 (read in the file)
					
						if file == '.DS_Store':
							continue
							
						if file.endswith('.npz'): # numpy files only
	
							# extract the relevant profiles
							file_path = path_first_part + str(year) + '/' + storm + '/' + sensor_type + '/' + aircraft_type + '/' + file
							file = np.load(file_path)
							ps = file['p']
							u = file['u_wind']
							v = file['v_wind']
							w = file['w_wind']
							v_tan = file['tan_wind']
							v_rad = file['rad_wind']
							spd = file['wind_speed']
							z = file['z']
							r = file['r']
							T = file['T']
							q = file['q']
							pot_T = file['pot_T']
							
							# calculate the median radius of the observations
							r_median = np.nanmedian(r)
							
							if np.nanmax(spd) < min_max_speed or np.nanmax(spd) > max_max_speed: # filter out weaker cyclones
								continue
							
							if np.isnan(r_median):
								continue
								
							if r_median > max_radius or r_median < min_radius: # does the radius fit within the acceptable bounds
								continue
								
							# calculate the observed diagnostics
							# height of maximum wind
							max_v_tan = np.nanmax(v_tan)
							if np.isnan(max_v_tan):
								continue
							max_i = np.where(v_tan == max_v_tan)[0][0]
							max_height = z[max_i]
							max_heights.append(max_height)
							
							# "surface" inflow angle (actually, approximate to 60 m...the lowest model level)
							diffs_from_60 = np.absolute(z - 60.)
							min_diff = np.nanmin(np.absolute(diffs_from_60))
							if np.isnan(min_diff):
								continue
							min_diff_i = np.where(diffs_from_60 == min_diff)[0][0]
							v_tan_60 = v_tan[min_diff_i]
							v_rad_60 = v_rad[min_diff_i]
							inflow_angle = np.arctan(v_rad_60 / v_tan_60)
							inflow_angle = np.absolute(inflow_angle * 180. / np.pi)
							inflow_angles.append(inflow_angle)
							
							# "inflow" depth (height at which radial inflow drops below 10% of "surface" radial inflow)
							rad_thresh = np.absolute(v_rad_60 * 0.1)
							below_thresh_i = np.where(np.absolute(v_rad) < rad_thresh)[0]
							if len(below_thresh_i) == 0:
								continue
							inflow_depth_i = np.nanmin(below_thresh_i)
							inflow_depth = z[inflow_depth_i]
							inflow_depths.append(inflow_depth)
								
							# if getting to this point, then the dropsonde profile fits within the specifications
							
							# now, bin the vertical profiles based on the pressure level of each observation
							for p_i, p in enumerate(ps):
								diffs = np.absolute(p_interp - p)
								if np.isnan(p) or np.isnan(z[p_i]) or np.isnan(u[p_i]) or np.isnan(v[p_i]) or np.isnan(v_tan[p_i]) or np.isnan(v_rad[p_i]) or np.isnan(spd[p_i]) or np.isnan(w[p_i]) or \
								np.isnan(q[p_i]) or np.isnan(pot_T[p_i]) or np.isnan(T[p_i]):
									continue
								# if getting to this point, no invalid values
								min_diff = np.nanmin(diffs)
								m_d_i = np.where(diffs == min_diff)[0][0] # min_diff_index
								u_sum[m_d_i] += u[p_i]
								v_sum[m_d_i] += v[p_i]
								v_tan_sum[m_d_i] += v_tan[p_i]
								v_rad_sum[m_d_i] += v_rad[p_i]
								spd_sum[m_d_i] += spd[p_i]
								w_sum[m_d_i] += w[p_i]
								z_sum[m_d_i] += z[p_i]
								T_sum[m_d_i] += T[p_i]
								q_sum[m_d_i] += q[p_i]
								pot_T_sum[m_d_i] += pot_T[p_i]
								counts[m_d_i] += 1.
						
	# calculate the composite profiles
	u_wind_profile = u_sum / counts
	v_wind_profile = v_sum / counts
	w_wind_profile = w_sum / counts
	tan_wind_profile = v_tan_sum / counts
	rad_wind_profile = v_rad_sum / counts
	spd_profile = spd_sum / counts
	z_profile = z_sum / counts
	T_profile = T_sum / counts
	q_profile = q_sum / counts
	pot_T_profile = pot_T_sum / counts
	p_profile = p_interp
	
	# make the diagnostics lists into numpy arrays
	max_heights = np.asarray(max_heights, dtype='float')
	inflow_angles = np.asarray(inflow_angles, dtype='float')
	inflow_depths = np.asarray(inflow_depths, dtype='float')
	
	# calculate a running average for the vertical profiles
	num = 11
	u_wind_profile = running_avg(u_wind_profile, num)
	v_wind_profile = running_avg(v_wind_profile, num)
	w_wind_profile = running_avg(w_wind_profile, 101)
	tan_wind_profile = running_avg(tan_wind_profile, num)
	rad_wind_profile = running_avg(rad_wind_profile, num)
	spd_profile = running_avg(spd_profile, num)
	z_profile = running_avg(z_profile, num)
	p_profile = running_avg(p_profile, num)
	T_profile = running_avg(T_profile, num)
	q_profile = running_avg(q_profile, num)
	pot_T_profile = running_avg(pot_T_profile, num)
	
	# return the complete lists of profiles to be plotted
	return u_wind_profile, v_wind_profile, w_wind_profile, tan_wind_profile, rad_wind_profile, spd_profile, z_profile, p_profile, T_profile, q_profile, pot_T_profile, max_heights, inflow_angles, inflow_depths
	
# "composite_NASA_profiles" function
# inputs:
# min_max_speed = a float representing the minimum maximum wind speed for the dropsonde profile (m/s)
# max_max_speed = a float representing the maximum maximum wind speed for the dropsonde profile (m/s)
# year_exceptions = a boolean representing whether or not subsetting by year will occur (optional, default is False)
# storm_exceptions = a boolean representing whether or not subsetting by storm will occur (optional, default is False)
# chosen_years = if year_exceptions is true, an array of strings containing the years to include in the composite (default is None)
# chosen_storms = if storm_exceptions is true, an array of strings containing the named storms to include in the composite (default is None)
# outputs:
# u_wind_profile = a numpy array of floats containing the composite vertical profile of zonal wind (m/s)
# v_wind_profile = a numpy array of floats containing the composite vertical profile of meridional wind (m/s)
# spd_profile = a numpy array of floats containing the composite vertical profile of total wind speed
# z_profile = a numpy array of floats containing the height (m)
# p_profile = a numpy array of floats containing the composite vertical pressure profile (hPa)
# notes:
# the profile path must be to a numpy binary file
# otherwise, an error will occur
def composite_NASA_profiles(min_max_speed, max_max_speed, year_exceptions=False, storm_exceptions=False, chosen_years=None, chosen_storms=None):

	# load any necessary additional modules
	import os
	import sys
	import types
	from output_functions import running_avg
	
	# define environmental path variable
	GOOGLE_ROOT = os.environ['GOOGLE_ROOT']
	
	# define the file path
	path_first_part = GOOGLE_ROOT + '/CPT_Project/ANALYSIS/obs/NASA-2015/'
	
	# define the uniform pressures over which interpolation will occur
	p_interp = np.linspace(1000., 500., 501)
	
	# define the summation array and the counter (for averaging purposes)
	u_sum = np.zeros([len(p_interp)]) 
	v_sum = np.zeros([len(p_interp)]) 
	spd_sum = np.zeros([len(p_interp)])
	z_sum = np.zeros([len(p_interp)])
	counts = np.zeros([len(p_interp)])

	# loop through the subdirectories to obtain all relevant dropsonde profiles
	
	updated_path = path_first_part

	# obtain the names of the storms (this is the next file level)
	storms = os.listdir(updated_path)

	for s_i, storm in enumerate(storms):

		if storm == '.DS_Store':
			continue
			
		if storm_exceptions:
			if storm[7:-9] not in chosen_storms:
				continue

		updated_path = path_first_part + storm
			
		# obtain the individual dropsonde files
		files = os.listdir(updated_path)
	
		for f_i, file in enumerate(files):
		
			# Step 4 (read in the file)
		
			if file == '.DS_Store':
				continue
				
			if file.endswith('.npz'): # numpy files only

				# extract the relevant profiles
				file_path = path_first_part + storm + '/' + file
				file = np.load(file_path)
				ps = file['p']
				u = file['u_wind']
				v = file['v_wind']
				spd = file['wind_speed']
				z = file['z']
				
				# at this point, the only values remaining are those that have valid data for both pressure and the variable of interest
				
				if np.nanmax(spd) < min_max_speed or np.nanmax(spd) > max_max_speed: # filter out weaker cyclones
					continue
				
				# if getting to this point, then the dropsonde profile fits within the specifications
				
				# now, bin the vertical profiles based on the pressure level of each observation
				for p_i, p in enumerate(ps):
					diffs = np.absolute(p_interp - p)
					if np.isnan(p) or np.isnan(z[p_i]) or np.isnan(u[p_i]) or np.isnan(v[p_i]) or np.isnan(spd[p_i]):
						continue
					# if getting to this point, no invalid values
					min_diff = np.nanmin(diffs)
					m_d_i = np.where(diffs == min_diff)[0][0] # min_diff_index
					u_sum[m_d_i] += u[p_i]
					v_sum[m_d_i] += v[p_i]
					spd_sum[m_d_i] += spd[p_i]
					z_sum[m_d_i] += z[p_i]
					counts[m_d_i] += 1.
						
	# calculate the composite profiles
	counts = np.where(counts==0, np.nan, counts)
	u_wind_profile = u_sum / counts
	v_wind_profile = v_sum / counts
	spd_profile = spd_sum / counts
	z_profile = z_sum / counts
	p_profile = p_interp
	
	# calculate a running average for the vertical profiles
	num = 11
	u_wind_profile = running_avg(u_wind_profile, num)
	v_wind_profile = running_avg(v_wind_profile, num)
	spd_profile = running_avg(spd_profile, num)
	z_profile = running_avg(z_profile, num)
	p_profile = running_avg(p_profile, num)
	
	# return the complete lists of profiles to be plotted
	return u_wind_profile, v_wind_profile, spd_profile, z_profile, p_profile
	
# "composite_radiosonde_profiles" function
# inputs:
# min_max_speed = a float representing the minimum maximum wind speed for the dropsonde profile (m/s)
# max_max_speed = a float representing the maximum maximum wind speed for the dropsonde profile (m/s)
# year_exceptions = a boolean representing whether or not subsetting by year will occur (optional, default is False)
# site_exceptions = a boolean representing whether or not subsetting by observation site will occur (optional, default is False)
# chosen_years = if year_exceptions is true, an array of strings containing the years to include in the composite (default is None)
# chosen_sites = if storm_exceptions is true, an array of strings containing the 4-character observation sites to include in the composite (default is None)
# outputs:
# u_wind_profile = a numpy array of floats containing the composite vertical profile of zonal wind (m/s)
# v_wind_profile = a numpy array of floats containing the composite vertical profile of meridional wind (m/s)
# spd_profile = a numpy array of floats containing the composite vertical profile of total wind speed
# z_profile = a numpy array of floats containing the height (m)
# p_profile = a numpy array of floats containing the composite vertical pressure profile (hPa)
# notes:
# the profile path must be to a numpy binary file
# otherwise, an error will occur
def composite_radiosonde_profiles(min_max_speed, max_max_speed, year_exceptions=False, site_exceptions=False, chosen_years=None, chosen_sites=None):

	# load any necessary additional modules
	import os
	import sys
	import types
	from output_functions import running_avg
	
	# define environmental path variable
	GOOGLE_ROOT = os.environ['GOOGLE_ROOT']
	
	# define the file path
	path_first_part = GOOGLE_ROOT + '/CPT_Project/ANALYSIS/obs/radiosondes_NWS/'
	
	# define the uniform pressures over which interpolation will occur
	p_interp = np.linspace(1000., 500., 501)
	
	# define the summation array and the counter (for averaging purposes)
	u_sum = np.zeros([len(p_interp)]) 
	v_sum = np.zeros([len(p_interp)]) 
	spd_sum = np.zeros([len(p_interp)])
	z_sum = np.zeros([len(p_interp)])
	counts = np.zeros([len(p_interp)])

	# loop through the subdirectories to obtain all relevant dropsonde profiles
	
	updated_path = path_first_part
			
	# obtain the individual dropsonde files
	files = os.listdir(updated_path)

	for f_i, file in enumerate(files):
	
		# determine the site
		site = file[0:4]
		if site_exceptions:
			if site not in chosen_sites:
				continue
	
		# Step 4 (read in the file)
	
		if file == '.DS_Store':
			continue
			
		if file.endswith('.npz'): # numpy files only

			# extract the relevant profiles
			file_path = path_first_part + file
			file = np.load(file_path)
			ps = file['p']
			u = file['u_wind']
			v = file['v_wind']
			spd = file['wind_speed']
			z = file['z']
			
			# at this point, the only values remaining are those that have valid data for both pressure and the variable of interest
			
			if np.nanmax(spd) < min_max_speed or np.nanmax(spd) > max_max_speed: # filter out weaker cyclones
				continue
			
			# if getting to this point, then the dropsonde profile fits within the specifications
			
			# now, bin the vertical profiles based on the pressure level of each observation
			for p_i, p in enumerate(ps):
				diffs = np.absolute(p_interp - p)
				if np.isnan(p) or np.isnan(z[p_i]) or np.isnan(u[p_i]) or np.isnan(v[p_i]) or np.isnan(spd[p_i]):
					continue
				# if getting to this point, no invalid values
				min_diff = np.nanmin(diffs)
				m_d_i = np.where(diffs == min_diff)[0][0] # min_diff_index
				u_sum[m_d_i] += u[p_i]
				v_sum[m_d_i] += v[p_i]
				spd_sum[m_d_i] += spd[p_i]
				z_sum[m_d_i] += z[p_i]
				counts[m_d_i] += 1.
						
	# calculate the composite profiles
	counts = np.where(counts==0, np.nan, counts)
	u_wind_profile = u_sum / counts
	v_wind_profile = v_sum / counts
	spd_profile = spd_sum / counts
	z_profile = z_sum / counts
	p_profile = p_interp
	
	# calculate a running average for the vertical profiles
	num = 11
	u_wind_profile = running_avg(u_wind_profile, num)
	v_wind_profile = running_avg(v_wind_profile, num)
	spd_profile = running_avg(spd_profile, num)
	z_profile = running_avg(z_profile, num)
	p_profile = running_avg(p_profile, num)
	
	# return the complete lists of profiles to be plotted
	return u_wind_profile, v_wind_profile, spd_profile, z_profile, p_profile
	
# "composite_G_IV_profiles" function
# inputs:
# min_max_speed = a float representing the minimum maximum wind speed for the dropsonde profile (m/s)
# max_max_speed = a float representing the maximum maximum wind speed for the dropsonde profile (m/s)
# date_exceptions = a boolean representing whether or not subsetting by date (yyyymmdd) will occur (optional, default is False)
# storm_exceptions = a boolean representing whether or not subsetting by storm will occur (optional, default is False)
# chosen_dates = if year_exceptions is true, an array of strings containing the dates to include in the composite (default is None)
# chosen_storms = if storm_exceptions is true, an array of strings containing the storm to include in the composite (default is None)
# outputs:
# u_wind_profile = a numpy array of floats containing the composite vertical profile of zonal wind (m/s)
# v_wind_profile = a numpy array of floats containing the composite vertical profile of meridional wind (m/s)
# spd_profile = a numpy array of floats containing the composite vertical profile of total wind speed
# z_profile = a numpy array of floats containing the height (m)
# p_profile = a numpy array of floats containing the composite vertical pressure profile (hPa)
# notes:
# the profile path must be to a numpy binary file
# otherwise, an error will occur
def composite_G_IV_profiles(min_max_speed, max_max_speed, date_exceptions=False, storm_exceptions=False, chosen_dates=None, chosen_storms=None):

	# load any necessary additional modules
	import os
	import sys
	import types
	from output_functions import running_avg
	
	# define environmental path variable
	GOOGLE_ROOT = os.environ['GOOGLE_ROOT']
	
	# define the file path
	path_first_part = GOOGLE_ROOT + '/CPT_Project/ANALYSIS/obs/G-IV_2015-2016/'
	
	# define the uniform pressures over which interpolation will occur
	p_interp = np.linspace(1000., 500., 501)
	
	# define the summation array and the counter (for averaging purposes)
	u_sum = np.zeros([len(p_interp)]) 
	v_sum = np.zeros([len(p_interp)]) 
	spd_sum = np.zeros([len(p_interp)])
	z_sum = np.zeros([len(p_interp)])
	counts = np.zeros([len(p_interp)])

	# loop through the subdirectories to obtain all relevant dropsonde profiles
	
	updated_path = path_first_part
	
	# obtain the subdirectories
	subdirs = os.listdir(updated_path)
	
	for s_i, subdir in enumerate(subdirs):
	
		if not subdir.endswith('.frd'):
			continue
	
		if date_exceptions:
			if subdir[0:8] not in chosen_dates:
				continue
		
		updated_path = path_first_part + subdir
			
		# obtain the individual dropsonde files
		files = os.listdir(updated_path)

		for f_i, file in enumerate(files):
	
			# Step 4 (read in the file)
	
			if file == '.DS_Store':
				continue
			
			if file.endswith('.npz'): # numpy files only

				# extract the relevant profiles
				file_path = path_first_part + subdir + '/' + file
				file = np.load(file_path)
				ps = file['p']
				u = file['u_wind']
				v = file['v_wind']
				spd = file['wind_speed']
				z = file['z']
			
				# at this point, the only values remaining are those that have valid data for both pressure and the variable of interest
			
				if np.nanmax(spd) < min_max_speed or np.nanmax(spd) > max_max_speed: # filter out weaker cyclones
					continue
			
				# if getting to this point, then the dropsonde profile fits within the specifications
			
				# now, bin the vertical profiles based on the pressure level of each observation
				for p_i, p in enumerate(ps):
					diffs = np.absolute(p_interp - p)
					if np.isnan(p) or np.isnan(z[p_i]) or np.isnan(u[p_i]) or np.isnan(v[p_i]) or np.isnan(spd[p_i]):
						continue
					# if getting to this point, no invalid values
					min_diff = np.nanmin(diffs)
					m_d_i = np.where(diffs == min_diff)[0][0] # min_diff_index
					u_sum[m_d_i] += u[p_i]
					v_sum[m_d_i] += v[p_i]
					spd_sum[m_d_i] += spd[p_i]
					z_sum[m_d_i] += z[p_i]
					counts[m_d_i] += 1.
						
	# calculate the composite profiles
	counts = np.where(counts==0, np.nan, counts)
	u_wind_profile = u_sum / counts
	v_wind_profile = v_sum / counts
	spd_profile = spd_sum / counts
	z_profile = z_sum / counts
	p_profile = p_interp
	
	# calculate a running average for the vertical profiles
	num = 11
	u_wind_profile = running_avg(u_wind_profile, num)
	v_wind_profile = running_avg(v_wind_profile, num)
	spd_profile = running_avg(spd_profile, num)
	z_profile = running_avg(z_profile, num)
	p_profile = running_avg(p_profile, num)
	
	# return the complete lists of profiles to be plotted
	return u_wind_profile, v_wind_profile, spd_profile, z_profile, p_profile
	
# "composite_P_3_profiles" function
# inputs:
# min_max_speed = a float representing the minimum maximum wind speed for the dropsonde profile (m/s)
# max_max_speed = a float representing the maximum maximum wind speed for the dropsonde profile (m/s)
# date_exceptions = a boolean representing whether or not subsetting by date (yyyymmdd) will occur (optional, default is False)
# storm_exceptions = a boolean representing whether or not subsetting by storm will occur (optional, default is False)
# chosen_dates = if year_exceptions is true, an array of strings containing the dates to include in the composite (default is None)
# chosen_storms = if storm_exceptions is true, an array of strings containing the storm to include in the composite (default is None)
# outputs:
# u_wind_profile = a numpy array of floats containing the composite vertical profile of zonal wind (m/s)
# v_wind_profile = a numpy array of floats containing the composite vertical profile of meridional wind (m/s)
# spd_profile = a numpy array of floats containing the composite vertical profile of total wind speed
# z_profile = a numpy array of floats containing the height (m)
# p_profile = a numpy array of floats containing the composite vertical pressure profile (hPa)
# notes:
# the profile path must be to a numpy binary file
# otherwise, an error will occur
def composite_P_3_profiles(min_max_speed, max_max_speed, date_exceptions=False, storm_exceptions=False, chosen_dates=None, chosen_storms=None):

	# load any necessary additional modules
	import os
	import sys
	import types
	from output_functions import running_avg
	
	# define environmental path variable
	GOOGLE_ROOT = os.environ['GOOGLE_ROOT']
	
	# define the file path
	path_first_part = GOOGLE_ROOT + '/CPT_Project/ANALYSIS/obs/P-3_2015-2016/'
	
	# define the uniform pressures over which interpolation will occur
	p_interp = np.linspace(1000., 500., 501)
	
	# define the summation array and the counter (for averaging purposes)
	u_sum = np.zeros([len(p_interp)]) 
	v_sum = np.zeros([len(p_interp)]) 
	spd_sum = np.zeros([len(p_interp)])
	z_sum = np.zeros([len(p_interp)])
	counts = np.zeros([len(p_interp)])

	# loop through the subdirectories to obtain all relevant dropsonde profiles
	
	updated_path = path_first_part
	
	# obtain the subdirectories
	subdirs = os.listdir(updated_path)
	
	for s_i, subdir in enumerate(subdirs):
	
		if not subdir.endswith('.frd'):
			continue
	
		if date_exceptions:
			if subdir[0:8] not in chosen_dates:
				continue
		
		updated_path = path_first_part + subdir
			
		# obtain the individual dropsonde files
		files = os.listdir(updated_path)

		for f_i, file in enumerate(files):
	
			# Step 4 (read in the file)
	
			if file == '.DS_Store':
				continue
			
			if file.endswith('.npz'): # numpy files only

				# extract the relevant profiles
				file_path = path_first_part + subdir + '/' + file
				file = np.load(file_path)
				ps = file['p']
				u = file['u_wind']
				v = file['v_wind']
				spd = file['wind_speed']
				z = file['z']
			
				# at this point, the only values remaining are those that have valid data for both pressure and the variable of interest
			
				if np.nanmax(spd) < min_max_speed or np.nanmax(spd) > max_max_speed: # filter out weaker cyclones
					continue
			
				# if getting to this point, then the dropsonde profile fits within the specifications
			
				# now, bin the vertical profiles based on the pressure level of each observation
				for p_i, p in enumerate(ps):
					diffs = np.absolute(p_interp - p)
					if np.isnan(p) or np.isnan(z[p_i]) or np.isnan(u[p_i]) or np.isnan(v[p_i]) or np.isnan(spd[p_i]):
						continue
					# if getting to this point, no invalid values
					min_diff = np.nanmin(diffs)
					m_d_i = np.where(diffs == min_diff)[0][0] # min_diff_index
					u_sum[m_d_i] += u[p_i]
					v_sum[m_d_i] += v[p_i]
					spd_sum[m_d_i] += spd[p_i]
					z_sum[m_d_i] += z[p_i]
					counts[m_d_i] += 1.
						
	# calculate the composite profiles
	counts = np.where(counts==0, np.nan, counts)
	u_wind_profile = u_sum / counts
	v_wind_profile = v_sum / counts
	spd_profile = spd_sum / counts
	z_profile = z_sum / counts
	p_profile = p_interp
	
	# calculate a running average for the vertical profiles
	num = 11
	u_wind_profile = running_avg(u_wind_profile, num)
	v_wind_profile = running_avg(v_wind_profile, num)
	spd_profile = running_avg(spd_profile, num)
	z_profile = running_avg(z_profile, num)
	p_profile = running_avg(p_profile, num)
	
	# return the complete lists of profiles to be plotted
	return u_wind_profile, v_wind_profile, spd_profile, z_profile, p_profile
	
# "composite_NOAA_HRD_profiles" function
# inputs:
# min_max_speed = a float representing the minimum maximum wind speed for the dropsonde profile (m/s)
# max_max_speed = a float representing the maximum maximum wind speed for the dropsonde profile (m/s)
# date_exceptions = a boolean representing whether or not subsetting by date (yyyymmdd) will occur (optional, default is False)
# storm_exceptions = a boolean representing whether or not subsetting by storm will occur (optional, default is False)
# chosen_dates = if year_exceptions is true, an array of strings containing the dates to include in the composite (default is None)
# chosen_storms = if storm_exceptions is true, an array of strings containing the storm to include in the composite (default is None)
# outputs:
# u_wind_profile = a numpy array of floats containing the composite vertical profile of zonal wind (m/s)
# v_wind_profile = a numpy array of floats containing the composite vertical profile of meridional wind (m/s)
# spd_profile = a numpy array of floats containing the composite vertical profile of total wind speed
# z_profile = a numpy array of floats containing the height (m)
# p_profile = a numpy array of floats containing the composite vertical pressure profile (hPa)
# notes:
# the profile path must be to a numpy binary file
# otherwise, an error will occur
def composite_NOAA_HRD_profiles(min_max_speed, max_max_speed, date_exceptions=False, storm_exceptions=False, chosen_dates=None, chosen_storms=None):

	# load any necessary additional modules
	import os
	import sys
	import types
	from output_functions import running_avg
	
	# define environmental path variable
	GOOGLE_ROOT = os.environ['GOOGLE_ROOT']
	
	# define the file path
	path_first_part = GOOGLE_ROOT + '/CPT_Project/ANALYSIS/obs/NOAA_HRD_2017-2018/'
	
	# define the uniform pressures over which interpolation will occur
	p_interp = np.linspace(1000., 500., 501)
	
	# define the summation array and the counter (for averaging purposes)
	u_sum = np.zeros([len(p_interp)]) 
	v_sum = np.zeros([len(p_interp)]) 
	spd_sum = np.zeros([len(p_interp)])
	z_sum = np.zeros([len(p_interp)])
	counts = np.zeros([len(p_interp)])

	# loop through the subdirectories to obtain all relevant dropsonde profiles
	
	updated_path = path_first_part
	
	# obtain the subdirectories
	storms = os.listdir(updated_path)
	
	for s_i, storm in enumerate(storms):
	
		if storm == '.DS_Store':
			continue
			
		if storm_exceptions:
			if storm not in chosen_storms:
				continue
			
		updated_path = path_first_part + storm
		
		flights = os.listdir(updated_path)
		
		for fl_i, flight in enumerate(flights):
		
			if flight == '.DS_Store':
				continue
	
			if not flight.endswith('_FRD'):
				continue
	
			if date_exceptions:
				if flight[0:8] not in chosen_dates:
					continue
		
			updated_path = path_first_part + storm + '/' + flight
			
			# obtain the individual dropsonde files
			files = os.listdir(updated_path)

			for f_i, file in enumerate(files):
	
				# Step 4 (read in the file)
	
				if file == '.DS_Store':
					continue
			
				if file.endswith('.npz'): # numpy files only

					# extract the relevant profiles
					file_path = path_first_part + storm + '/' + flight + '/' + file
					file = np.load(file_path)
					ps = file['p']
					u = file['u_wind']
					v = file['v_wind']
					spd = file['wind_speed']
					z = file['z']
			
					# at this point, the only values remaining are those that have valid data for both pressure and the variable of interest
			
					if np.nanmax(spd) < min_max_speed or np.nanmax(spd) > max_max_speed: # filter out weaker cyclones
						continue
			
					# if getting to this point, then the dropsonde profile fits within the specifications
			
					# now, bin the vertical profiles based on the pressure level of each observation
					for p_i, p in enumerate(ps):
						diffs = np.absolute(p_interp - p)
						if np.isnan(p) or np.isnan(z[p_i]) or np.isnan(u[p_i]) or np.isnan(v[p_i]) or np.isnan(spd[p_i]):
							continue
						# if getting to this point, no invalid values
						min_diff = np.nanmin(diffs)
						m_d_i = np.where(diffs == min_diff)[0][0] # min_diff_index
						u_sum[m_d_i] += u[p_i]
						v_sum[m_d_i] += v[p_i]
						spd_sum[m_d_i] += spd[p_i]
						z_sum[m_d_i] += z[p_i]
						counts[m_d_i] += 1.
						
	# calculate the composite profiles
	counts = np.where(counts==0, np.nan, counts)
	u_wind_profile = u_sum / counts
	v_wind_profile = v_sum / counts
	spd_profile = spd_sum / counts
	z_profile = z_sum / counts
	p_profile = p_interp
	
	# calculate a running average for the vertical profiles
	num = 11
	u_wind_profile = running_avg(u_wind_profile, num)
	v_wind_profile = running_avg(v_wind_profile, num)
	spd_profile = running_avg(spd_profile, num)
	z_profile = running_avg(z_profile, num)
	p_profile = running_avg(p_profile, num)
	
	# return the complete lists of profiles to be plotted
	return u_wind_profile, v_wind_profile, spd_profile, z_profile, p_profile
	
# "std_dev" function
# inputs:
# config: string that defines the experiment's configuration
# num_ensembles: the integer number of ensemble members
# time_step: a string representing the time step in days
# var_name: a string representing the variable's name
# outputs:
# std: a 2D array of standard deviations for that time step
def std_dev( config, num_ensembles, time_step, var_name ):
	# import the necessary modules
	import numpy as np
	from netCDF4 import Dataset
	
	# read in the first ensemble member file for this time step
	file_path = '/Users/kmn182/ICS_scratch/output/' + config + '.001/' + config + '.001_' + time_step + '.nc'
	dataset = Dataset(file_path)
	var = dataset.variables[var_name][:]
	dims = var.shape # the dimensions of the variable file (should be levels x radii)
	num_levels = dims[0]
	num_rad = dims[1]
	
	# make an array to store all ensemble data (should be size levels x radii x ensemble members)
	array = np.ones([num_levels, num_rad, num_ensembles]) * np.nan
	array[:, :, 0] = var # the first time step
	
	# read in the ensemble files for the remaining ensembles
	for n in range(2, num_ensembles + 1):
		file_path = '/Users/kmn182/ICS_scratch/output/' + config + '.' + str(n).rjust(3, '0') + '/' + config + '.' + str(n).rjust(3, '0') + '_' + time_step + '.nc'
		dataset = Dataset(file_path)
		var = dataset.variables[var_name][:]
		array[:, :, n-1] = var
		
	# at this point, a full 3D array of the model output variable is complete
	# now, calculate standard deviation over axis 2
	std = np.nanstd(array, axis=2) 
	
	dataset.close() # close the dataset
	
	# return the standard deviation for the configuration and time step
	return std
	
# "read_best_track function
# inputs:
# file_path: string that represents the full file path of the text or text-like file
# outputs:
# min_ps: a numpy array of floats containing each minimum surface pressure observation (in hPa)
# max_spd: a numpy array of floats containing each maximum sustained wind observation (in knots)
def read_best_track( file_path ):
	# import the necessary modules
	import numpy as np
	
	# read-in the file
	file = open(file_path, 'r')
	lines = file.readlines()
	
	# create the speed and pressure lists
	max_spd = []
	min_ps = []
	
	# read in the values line by line
	for line in lines:
		split_line = line.split(', ')
		
		if len(split_line) < 4: # those header lines that do not have numerical data have only three elements
			continue
		elif split_line[6] == '-999' or split_line[7] == '-999': # is the speed (element 6) or pressure (element 7) invalid?
			continue
		
		# if getting to this point, the line has valid data
		
		# add to the lists
		max_spd.append(float(split_line[6]))
		min_ps.append(float(split_line[7]))
		
	# make the lists into numpy arrays
	max_spd = np.asarray(max_spd, dtype='float')
	min_ps = np.asarray(min_ps, dtype='float')
	
	# convert the speed from knots into m/s
	max_spd = max_spd * 0.514
	
	# output the results
	return max_spd, min_ps
	

	
