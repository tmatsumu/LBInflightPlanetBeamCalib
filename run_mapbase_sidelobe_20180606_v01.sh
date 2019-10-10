#!/bin/sh

dir=/Users/tomotake_matsumura/Qsync/Projects/Entries/20180305_BeamCalibration/20130801_PlanetCal_v2/
dir_out=${dir}/output/ 
dir_src=${dir}/src/

D_mm[0]=400
D_mm[1]=400
D_mm[2]=400
D_mm[3]=400
D_mm[4]=400
D_mm[5]=400
D_mm[6]=400
D_mm[7]=400
D_mm[8]=400

nu_obs[0]=40
nu_obs[1]=50
nu_obs[2]=60
nu_obs[3]=68
nu_obs[4]=78
nu_obs[5]=89
nu_obs[6]=100
nu_obs[7]=119
nu_obs[8]=140

src_planet[0]='Jupiter'
src_planet[1]='Mars'
src_planet[2]='Saturn'

#option=ellipticalGaussian
option=airyfunction

for j in {0..8}; do
	for k in {0..2}; do
	    python ${dir_src}/main_mapbase_sidelobe_v8.py ${D_mm[j]} ${nu_obs[$j]} ${src_planet[$k]} ${dir_out} ${option}
	done
done

