import numpy as np
import pylab as py
import sys
import glob
import lib_planetcal as lib_p
print ''

pi = np.pi
radeg = (180./pi)
res = .05/radeg
map_width_rad = 3./radeg
noise = 5./(res*radeg*60.)*np.sqrt(300.)

nu_obs = 140.e9
Dapt_mm = 300.
src_planet = 'Jupiter'
x0,y0,sigma_x,sigma_y,phi,samplerate = lib_p.gen_InputParams(nu_obs,Dapt_mm)
beam_solid_str = lib_p.beam_solidangle_AirySymmetric(nu_obs,Dapt_mm*1e-3)

Amp = lib_p.planet_info(src_planet,nu_obs,beam_solid_str)
par_in = np.array([x0,y0,sigma_x,sigma_y,phi,Amp])
elip = lib_p.ellipticalGaussian()
elip.resol_rad = res
elip.map_width_rad = map_width_rad
elip.par = par_in
elip.beam_sigma = sigma_x
X, Y, MAP_S = elip.gen_flatellip_map()

num = len(MAP_S[0])
MAP_N = noise*np.random.normal(0.,1.,[num,num])
MAP_SN = MAP_S+MAP_N

num_bin = 100
Bkk_S, kx, ky = lib_p.fft_2d(res,res,MAP_S)
Bk_S_kr, Bk_S_mean, Bk_S_std, Bk_S_med = lib_p.cal_Bk(num_bin,kx,ky,Bkk_S)

Bkk_N, kx, ky = lib_p.fft_2d(res,res,MAP_N)
Bk_N_kr, Bk_N_mean, Bk_N_std, Bk_N_med = lib_p.cal_Bk(num_bin,kx,ky,Bkk_N)

Bkk_SN, kx, ky = lib_p.fft_2d(res,res,MAP_SN)
Bk_SN_kr, Bk_SN_mean, Bk_SN_std, Bk_SN_med = lib_p.cal_Bk(num_bin,kx,ky,Bkk_SN)


py.figure(0, figsize=(18,10))

py.subplot(231)
X = X*radeg
Y = Y*radeg
Z = MAP_S
xmin, xmax, ymin, ymax = np.min(X), np.max(X), np.min(Y), np.max(Y)
extent = xmin, xmax, ymin, ymax
im1 = py.imshow(Z, extent=extent)
py.colorbar()

py.subplot(232)
Z = MAP_N
xmin, xmax, ymin, ymax = np.min(X), np.max(X), np.min(Y), np.max(Y)
extent = xmin, xmax, ymin, ymax
im1 = py.imshow(Z, extent=extent)
py.colorbar()

py.subplot(233)
Z = MAP_SN
xmin, xmax, ymin, ymax = np.min(X), np.max(X), np.min(Y), np.max(Y)
extent = xmin, xmax, ymin, ymax
im1 = py.imshow(Z, extent=extent)
py.colorbar()

py.subplot(234)
xmin, xmax, ymin, ymax = np.min(kx), np.max(kx), np.min(ky), np.max(ky)
extent = xmin, xmax, ymin, ymax
im1 = py.imshow(np.abs(Bkk_S), extent=extent)
py.colorbar()

py.subplot(235)
xmin, xmax, ymin, ymax = np.min(kx), np.max(kx), np.min(ky), np.max(ky)
extent = xmin, xmax, ymin, ymax
im1 = py.imshow(np.abs(Bkk_N), extent=extent)
py.colorbar()

py.subplot(236)
py.errorbar(np.array(Bk_S_kr), Bk_S_mean/max(Bk_S_mean), Bk_S_std/max(Bk_S_mean), fmt='o')
#py.semilogx()
#py.errorbar(Bk_N_kr, Bk_N_mean/max(Bk_N_mean), Bk_N_std/max(Bk_N_mean),fmt='o')
#py.semilogx()
py.errorbar(Bk_SN_kr, Bk_SN_mean/max(Bk_SN_mean), Bk_SN_std/max(Bk_SN_mean),fmt='.')
py.semilogy()
py.ylim([1e-12,1.1])

py.savefig( 'png/'+src_planet+'_'+str(nu_obs*1e-9)+'_'+str(Dapt_mm)+'.png' )
np.savez( 'png/'+src_planet+'_'+str(nu_obs*1e-9)+'_'+str(Dapt_mm) )
#py.show()
