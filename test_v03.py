import numpy as np
import pylab as py
import sys
import glob
import lib_planetcal as lib_p
import lib_lbv24 as lbv24

pi = np.pi
h = 6.67e-34
c = 2.99792458e8 
kb = 1.38e-23
radeg = (180./pi)

nu_obsGHz = 100.

x0_rad,y0_rad,sigma_x_rad,sigma_y_rad,phi,samplerate = lbv24.gen_cambus(nu_obsGHz)
res_rad = sigma_x_rad*0.1 # res is in unit of rad
map_width_rad = sigma_x_rad*10. # unit in rad

par_in = np.array([x0_rad,y0_rad,sigma_x_rad,sigma_y_rad,phi,1.])

elip = lib_p.ellipticalGaussian()
elip.resol_rad = res_rad
elip.map_width_rad = map_width_rad
elip.par = par_in

edgetaper = 0.5
radius = 200.e-3
nu_obs = nu_obsGHz * 1.e9

theta, out = elip.gen_flatTruncGauss_2D(nu_obs,edgetaper,radius)

print theta, out
print len(theta), len(out)

py.plot(theta, out)
py.semilogy()
py.show()