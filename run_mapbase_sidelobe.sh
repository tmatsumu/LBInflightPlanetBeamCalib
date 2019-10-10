#!/bin/sh

D_mm[0]=250
D_mm[1]=300
D_mm[2]=350
D_mm[3]=400
D_mm[4]=450
D_mm[5]=500
D_mm[6]=550
D_mm[7]=600

nu_obs[0]=60
nu_obs[1]=78
nu_obs[2]=100
nu_obs[3]=140
nu_obs[4]=195
nu_obs[5]=280

src_planet[0]='Mars'
src_planet[1]='Jupiter'
src_planet[2]='Saturn'

for i in {0..7}; do
    for j in {0..5}; do
	for k in {0..2}; do
	    python main_mapbase_sidelobe_v4.py ${D_mm[$i]} ${nu_obs[$j]} ${src_planet[$k]}
	done
    done
done