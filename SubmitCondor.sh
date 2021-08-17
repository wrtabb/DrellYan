#!/bin/bash

lepType=ELE
sampleType=SAMPLE_LL
runType=TEST
root_file=(
	DYLL_M10to50_EE
	DYLL_M50to100_EE
	DYLL_M100to200_EE
	DYLL_M200to400_EE
	DYLL_M400to500_EE
	DYLL_M500to700_EE
	DYLL_M700to800_EE
	DYLL_M800to1000_EE
	DYLL_M1000to1500_EE
	DYLL_M1500to2000_EE
	DYLL_M2000to3000_EE
)

for index in ${!root_file[*]}; do
	echo "Beginning to process ${root_file[$index]}"
	condor_submit \
		arg1=${runType} \
		arg2=${sampleType} \
		arg3=${lepType} \
		arg4=${root_file[$index]} \
		arg5=${disk_needed[$index]} \
		condor_control.condor
done
