#!/bin/bash

lepType=(
	ELE
	#MUON
	)
sampleType=SAMPLE_LL
runType=TEST

for index in ${!lepType[*]}; do
	echo "Beginning to process ${lepType[$index]}"
	condor_submit \
		arg1=${runType} \
		arg2=${sampleType} \
		arg3=${lepType[$index]} \
		condor_control.condor
done
