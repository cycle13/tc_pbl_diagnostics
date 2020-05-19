#!/Users/kmn182/anaconda3/envs/ENV1/bin/python
# ------------------------------------------------------------
#
#	Program Name: generate_sfc_exch_coeff_plots.py
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

# define the font family
mpl.rc('font',family='Arial')

# define the configurations to be plotted
configs = ['RCE.QPC6.ne0np4tcfplane.ne15x8.exp715', 'RCE.QPC6.ne0np4tcfplane.ne15x8.exp718', 'RCE.QPC6.ne0np4tcfplane.ne15x8.exp721', 'RCE.QPC6.ne0np4tcfplane.ne15x8.exp722', 'RCE.QPC6.ne0np4tcfplane.ne15x8.exp723'] 
num_configs = len(configs)

# define a list of colors
custom_colormap = Cubehelix.make(start_hue=240, end_hue=-300, min_sat=1, max_sat=2.5, gamma=1, min_light=0.3, max_light=0.8, n=num_configs, reverse=False, name='custom_colormap')
colors = custom_colormap.hex_colors

# define the input and output paths
input_path = '/Users/kmn182/ICS_work/CPT_Project/ANALYSIS/cd_ck/'
output_path = '/Users/kmn182/Documents/CPT_Project/FIGURES/multi_config/'

# Step 2 (plot the points for each configuration and time step for Ck vs. Cd)

# set up the figure
f0 = plt.figure(figsize=(12,12))

plt.title('Ck/Cd Ratio', fontsize=60)
plt.xlabel('Drag Coefficient', fontsize=32)
plt.ylabel('Enthalpy Coefficient', fontsize=32)
plt.xlim(0, 0.020)
plt.ylim(0, 0.020)
plt.xticks(np.linspace(0., 0.020, 6), np.linspace(0., 0.020, 6), fontsize=24)
plt.yticks(np.linspace(0., 0.020, 6), np.linspace(0., 0.020, 6), fontsize=24)

# plot the baseline 1-to-1 ratio
plt.plot(np.linspace(0., 0.020, 6), np.linspace(0., 0.020, 6), linewidth=4, color='black')

for c_i, config in enumerate(configs): # loop through all configurations
	f = np.load(input_path + config + '_cd_ck.npz')
	CD = f['CDs'] # drag coefficient
	CK = f['CKs'] # enthalpy coefficient
	plt.plot(CD, CK, '.', markersize=20, color=colors[c_i], label=config)

# finalize the figure
plt.legend(loc='upper left', fontsize=18)
plt.savefig(output_path + config + '_ck_cd.png', dpi=400, format='png')

plt.close()

# Step 3 (plot Ck vs. U10)

f0 = plt.figure(figsize=(12,12))

plt.title('Ck vs. U10', fontsize=60)
plt.xlabel('10-meter Wind Speed (m/s)', fontsize=32)
plt.ylabel('Enthalpy Coefficient', fontsize=32)
plt.xlim(0, 80)
plt.ylim(0, 0.005)
plt.xticks(np.linspace(0., 80., 5), np.asarray(np.linspace(0., 80., 5), dtype='int'), fontsize=24)
plt.yticks(np.linspace(0., 0.005, 6), np.linspace(0., 0.005, 6), fontsize=24)

for c_i, config in enumerate(configs):
	f = np.load(input_path + config + '_cd_ck.npz')
	CK_4D = np.ravel(f['CK_4D'])
	U10_4D = np.ravel(f['U10_4D'])
	plt.plot(U10_4D, CK_4D, '.', markersize=10, color=colors[c_i], alpha=0.5, label=config)

plt.legend(loc='lower right', fontsize=18)
plt.savefig(output_path + config + '_ck_u10.png', dpi=400, format='png')

plt.close()

# Step 4 (plot the points Cd vs. U10)

f0 = plt.figure(figsize=(12,12))

plt.title('Cd vs. U10', fontsize=60)
plt.xlabel('10-meter Wind Speed (m/s)', fontsize=32)
plt.ylabel('Drag Coefficient', fontsize=32)
plt.xlim(0, 80)
plt.ylim(0, 0.005)
plt.xticks(np.linspace(0., 80., 5), np.asarray(np.linspace(0., 80., 5), dtype='int'), fontsize=24)
plt.yticks(np.linspace(0., 0.005, 6), np.linspace(0., 0.005, 6), fontsize=24)

for c_i, config in enumerate(configs):
	f = np.load(input_path + config + '_cd_ck.npz')
	CD_4D = np.ravel(f['CD_4D'])
	U10_4D = np.ravel(f['U10_4D'])
	Umps = np.where(U10_4D >= 33., 33., U10_4D)
	CD_LY = (0.0027 / Umps) + 0.000142 + (0.0000764 * Umps) - (3.14807e-13 * (Umps**6))
	plt.plot(U10_4D, CD_4D, '.', markersize=10, color=colors[c_i], alpha=0.5, label=config)
	if c_i == len(configs)-1:
		plt.plot(U10_4D, CD_LY, '.', markersize=10, color='black', alpha=0.5, label='LY09')
	else:
		plt.plot(U10_4D, CD_LY, '.', markersize=10, color='black', alpha=0.5)

plt.legend(loc='lower right', fontsize=18)
plt.savefig(output_path + config + '_cd_u10.png', dpi=400, format='png')

plt.close()

print('The figures have been generated!')



