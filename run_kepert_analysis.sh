#!/bin/bash
# ------------------------------------------------------------
#
#	Program Name: run_kepert_analysis.sh
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

cd /Users/kmn182/Documents/CPT_Project/CODE

# Step 3 (define configurations)

calc_diagnostics="False" # should diagnostics time series be created?
diagnostics_configs=( "RCE.QPC6.ne0np4tcfplane.ne15x8.exp722" "RCE.QPC6.ne0np4tcfplane.ne15x8.exp723" ) # these are the configurations that need diagnostics calculated

plotted_configs=( "RCE.QPC6.ne0np4tcfplane.ne15x8.exp715" "RCE.QPC6.ne0np4tcfplane.ne15x8.exp718" "RCE.QPC6.ne0np4tcfplane.ne15x8.exp721" "RCE.QPC6.ne0np4tcfplane.ne15x8.exp722" "RCE.QPC6.ne0np4tcfplane.ne15x8.exp723" ) # the configurations to be plotted
num_configs=${#plotted_configs[@]}
num_ensembles=20 # assumed to be labelled (001, 002, ...)

# Step 4 (run the analysis)

# first, calculate the diagnostics if necessary
if [ ${calc_diagnostics} == "True" ]
then
	for c in "${diagnostics_configs[@]}"
	do	
		echo "Now starting analysis for configuration $c."
		echo "------------------------"
		python calculate_diagnostics_kepert_ens_avg.py $c ${num_ensembles}
		echo "---------------------------------"
		echo "Diagnostics have been calculated for the ensemble average!"
		echo "---------------------------------"
	done
fi

# now, run scripts to plot data for multiple configurations
python generate_vertical_profiles_kepert_multi_config.py ${num_configs} ${plotted_configs[*]} 
echo "--------------------------------------"
echo "Vertical profiles have been generated for the ensemble averages!"
echo "--------------------------------------"
python generate_diagnostics_time_series_kepert_multi_config.py ${num_configs} ${plotted_configs[*]}
echo "--------------------------------------"
echo "Diagnostics time series have been generated for the ensemble averages!"
echo "--------------------------------------"


