#!/Users/kmn182/anaconda3/envs/ENV1/bin/python
# ------------------------------------------------------------
#
#	Program Name: generate_betacast_modeled_centers.py
#
#	Purpose: To get the modeled cyclone center
#
#	WARNING: This code assumes Python 3, so there may be some 
#            syntax errors if running with Python 2.
#
# -------------------------------------------------------------

# Step 1 (import necessary modules)

import numpy as np
from netCDF4 import Dataset
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
from palettable.cmocean.diverging import Curl_20_r
from palettable.cartocolors.diverging import Earth_7
from palettable.cmocean.sequential import Ice_20_r, Tempo_20, Dense_20
from palettable.cubehelix import perceptual_rainbow_16_r, classic_16_r, cubehelix1_16_r, cubehelix2_16_r, jim_special_16_r, cubehelix3_16_r, perceptual_rainbow_16_r
from palettable.cubehelix import Cubehelix

input_path = '/Users/kmn182/ICS_scratch/output/'
ouput_path = '/Users/kmn182/Documents/CPT_Project/FIGURES/'

# define the font family
mpl.rc('font',family='Arial')

configs = ['RCE.QPC6.ne0np4tcfplane.ne15x8.exp709', 'RCE.QPC6.ne0np4tcfplane.ne15x8.exp710'] 
num_configs = len(configs)

# define a list of colors
custom_colormap = Cubehelix.make(start_hue=240, end_hue=-300, min_sat=1, max_sat=2.5, gamma=1, min_light=0.3, max_light=0.8, n=num_configs, reverse=False, name='custom_colormap')
colors = custom_colormap.hex_colors

time_steps = ['01']
time_steps = time_steps + list(range(2,17))
time_steps = np.asarray(time_steps, dtype='str')
for t_i, time_step in enumerate(time_steps):
	if t_i == 0:
		continue
	time_steps[t_i] = time_step.rjust(2,'0')
	
ensembles = ['001']
ensembles = ensembles + list(range(2,21))
ensembles = np.asarray(ensembles, dtype='str')
for e_i, ens in enumerate(ensembles):
	if e_i == 0:
		continue
	ensembles[e_i] = ens.rjust(3,'0')
	
num_ensembles = len(ensembles)
num_time_steps = len(time_steps)

total_points = num_ensembles * num_time_steps

CDs = np.ones([num_configs, total_points]) * np.nan
CKs = np.ones([num_configs, total_points]) * np.nan

# Step 2 (create the figure...not currently used)
# see below

# Step 3 (calculate Ck and Cd for each configuration)

for c_i, config in enumerate(configs):
	
	ouput_path = '/Users/kmn182/Documents/CPT_Project/FIGURES/multi_config/'
	counter = 0

	for t_i, time_step in enumerate(time_steps):
	
		for e_i, ens in enumerate(ensembles):
		
			input_path = '/Users/kmn182/ICS_scratch/output/' + config + '.' + ens + '/run/'

			# read in the file
			file_name = config + '.' + ens + '.cam.h3.0001-01-' + time_step + '-00000.nc_regrid.nc'
			dataset = Dataset(input_path + file_name, 'r+')
		
			# define constants
			Rd = 287. # gas constant
			cp = 1004. # heat capacity
			Lv = 2.5104e6 # latent heat of vaporization
			p0 = 100000. # standard pressure in Pa

			# read in required variables
			shf = dataset.variables['SHFLX'][:] # sensible heat flux
			lhf = dataset.variables['LHFLX'][:] # latent heat flux
			tbot = dataset.variables['TBOT'][:] # temperature at bottom
			ubot = dataset.variables['UBOT'][:] # U wind at bottom
			ustar = dataset.variables['USTAR'][:] # friction velocity
			vbot = dataset.variables['VBOT'][:] # V wind at bottom
			ps = dataset.variables['PS'][:] # time x lat x lon
			wbot = dataset.variables['QBOT'][:] # mixing ratio at bottom
			sst = dataset.variables['SST'][:] # sea-surface temperature
			taux = dataset.variables['TAUX'][:] # tau in X direction
			tauy = dataset.variables['TAUY'][:] # tau in Y direction
			a = dataset.variables['hyam'][:] # hybrid A levels
			b = dataset.variables['hybm'][:] # hybrid B levels
			
			dataset.close()

			# print some SST output
			print("SST -- max: " + str(np.nanmax(sst)) + "   min: " + str(np.nanmin(sst)) + "   mean: " + str(np.nanmean(sst)))

			# calculate derived variables
			tau = np.sqrt(taux**2 + tauy**2)
			dims = ps.shape # time x lats x lons
			fullpres = np.ones([dims[0], len(a), dims[1], dims[2]]) * np.nan
			for k_i, k in enumerate(a):
				fullpres[:, k_i, :, :] = a[k_i] * p0 + b[k_i] * ps # like NCL function pres_hybrid_ccm
			rho = fullpres[:, 29, :, :] / (Rd * tbot)
			windbot = np.sqrt(vbot**2 + ubot**2)
			theta = tbot * (p0 / fullpres[:, 29, :, :])**0.286

			# calculate DELQ
			qbot = wbot / (1 + wbot) # same as mixhum_convert in NCL

			# from Hobbs 1977
			esat = 33.8639 * ((0.00738 * (sst - 273.15) + 0.8072)**8 - 0.000019 * np.absolute(1.8 * (sst - 273.15) + 48.) + 0.0013)
			wsat = 0.622 * esat / (ps / 100.)
			qsat = wsat / (1 + wsat)
		
			# method from CESM
			#qsat2 = 640380. / np.exp(5107.4 / sst)
			#qsat2 := 0.98 * qsat2 / rho

			delq = qsat - qbot

			# calculate DELT
		
			delt = sst - tbot

			# calculate exchange coefficients

			CH = shf / (rho * cp * windbot * delt)
			CD = (ustar / windbot)**2
			CQ = lhf / (rho * Lv * windbot * delq)

			#CK = (shf+lhf)/(rho*Lv*windbot*delq+rho*cp*windbot*delt)
			CK = (CH * shf + CQ * lhf) / (shf + lhf)
			
			ratio = CK / CD
			
			CKs[c_i, counter] = np.nanmedian(CK)
			CDs[c_i, counter] = np.nanmedian(CD)
			
			# Step 4 (add to the existing netCDF file)

			# Create coordinate variables
			# C_D = dataset.createVariable('CD',np.float64,('time', 'lat', 'lon'))
# 			C_H = dataset.createVariable('CH',np.float64,('time', 'lat', 'lon'))
# 			C_Q = dataset.createVariable('CQ',np.float64,('time', 'lat', 'lon'))
# 			C_K = dataset.createVariable('CK',np.float64,('time', 'lat', 'lon'))
# 			RATIO = dataset.createVariable('RATIO',np.float64,('time', 'lat', 'lon'))
# 			WINDBOT = dataset.createVariable('WINDBOT',np.float64,('time', 'lat', 'lon'))
# 
# 			# Fill the variables with data
# 			C_D[:, :, :] = CD
# 			C_H[:, :, :] = CH
# 			C_Q[:, :, :] = CQ
# 			C_K[:, :, :] = CK
# 			RATIO[:, :, :] = ratio
# 			WINDBOT[:, :, :] = windbot
# 			
# 			dataset.close()

f0 = plt.figure(figsize=(10,10))

plt.title('Ck/Cd Ratio', fontsize=32)
plt.xlabel('Drag Coefficient', fontsize=24)
plt.ylabel('Entropy Coefficient', fontsize=24)
plt.xlim(0, 0.002)
plt.ylim(0, 0.002)
plt.xticks(np.linspace(0., 0.002, 5), np.linspace(0., 0.002, 5))
plt.yticks(np.linspace(0., 0.002, 5), np.linspace(0., 0.002, 5))

# plot the baseline
plt.plot(np.linspace(0., 0.002, 5), np.linspace(0., 0.002, 5), linewidth=4, color='black')

plt.legend(loc='upper right')
plt.savefig(output_path + config + '_ck_cd.png', dpi=400, format='png')

for c_i, config in enumerate(configs):
	plt.plot(CDs[c_i,:], CKs[c_i,:], '.', markersize=30, color=colors[c_i], label=config)

