#!/bin/bash

lepType=(
	ELE
#	MUON
)
sampleType=SAMPLE_LL
runType=V2P6
root_file=(
#	DYLL_M10to50
#	DYLL_M50to100
#	DYLL_M100to200
#	DYLL_M200to400
#	DYLL_M400to500
#	DYLL_M500to700
#	DYLL_M700to800
#	DYLL_M800to1000
#	DYLL_M1000to1500
#	DYLL_M1500to2000
#	DYLL_M2000to3000
	DATA_RunB
	DATA_RunC
	DATA_RunD
	DATA_RunE
	DATA_RunF
	DATA_RunG
	DATA_RunHver2
	DATA_RunHver3
)

for index1 in ${!lepType[*]}; do
	for index2 in ${!root_file[*]}; do
		echo "Beginning to process ${root_file[$index2]} for ${lepType[$index1]}"
		condor_submit \
			arg1=${runType} \
			arg2=${sampleType} \
			arg3=${lepType[$index1]} \
			arg4=${root_file[$index2]} \
			condor_control.condor
	done
done
