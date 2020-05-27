#!/bin/bash
# ------------------------------------------------------------
#
#	Program Name: run_kepert_analysis_sens.sh
#
#	Purpose: Call python scripts to calculate model diagnostics,
#            plot cross sections, and plot vertical profiles. 
#            Output comes from CESM tropical cyclone model
#            runs on an f-plane. 
#
#	WARNING: This code assumes Python 3, so there may be some 
#            syntax errors if running with Python 2.
#
# -------------------------------------------------------------

# Step 1 (set the root path)

export GOOGLE_ROOT=/Users/kmn182/Documents

# Step 2 (move to the proper directory)

cd /Users/kmn182/Documents/CPT_Project/CODE/Morris

# Step 3 (define configurations)

configs=( "RCE.QPC6.ne0np4tcfplane.ne15x8.expS02" )
num_ensembles=40 # assumed to be labelled (001, 002, ...)

# Step 4 (run the analysis code in succession)

for c in "${configs[@]}"
do	
	# calculate the diagnostics
	# echo "Now starting analysis for configuration $c."
# 	echo "------------------------"
# 	python calculate_diagnostics_kepert_ens_avg_sens.py ${c} ${num_ensembles}
# 	echo "---------------------------------"
# 	echo "Diagnostics have been calculated for the ensemble average!"
	
	# Now run the plotting scripts
	python generate_vertical_profiles_kepert_multi_ens.py ${c} ${num_ensembles} 
	echo "--------------------------------------"
	echo "Vertical profiles have been generated for the ensemble averages!"
	echo "--------------------------------------"
	python generate_diagnostics_time_series_kepert_multi_ens.py ${c} ${num_ensembles}
	echo "--------------------------------------"
	echo "Diagnostics time series have been generated for the ensemble averages!"
	echo "--------------------------------------"
	#python calc_sfc_exch_coeff.py ${num_configs} ${all_configs[*]}
	#echo " -------------------------------------"
	#echo "Surface exchange coefficients have been calculated and plotted!"
	#echo "--------------------------------------"
done

echo "END OF ANALYSIS"
echo "--------------------------------------"

