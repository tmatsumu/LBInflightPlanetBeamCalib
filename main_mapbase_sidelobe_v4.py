import numpy as np
import pylab as py
import sys
import glob
import lib_planetcal as lib_p
print ''

def noise_info(nu_obs):
    if nu_obs == '60': NETarr = 10.3; det_num = 304.
    if nu_obs == '78': NETarr = 6.5; det_num = 304.
    if nu_obs == '100': NETarr = 4.7; det_num = 304.
    if nu_obs == '140': NETarr = 3.7; det_num = 370.
    if nu_obs == '195': NETarr = 3.1; det_num = 370.
    if nu_obs == '280': NETarr = 3.8; det_num = 370.
    return NETarr*np.sqrt(det_num)

print '---'
print sys.argv[1]
print sys.argv[2]
print sys.argv[3]
print '---'

pi = np.pi
radeg = (180./pi)
#res = .05/radeg
map_width_rad = 3./radeg
#noise = 5./(res*radeg*60.)*np.sqrt(300.)*0.5

nu_obs = float(sys.argv[2]) #140.e9
nu_obs = nu_obs*1.e9
Dapt_mm = float(sys.argv[1])  #300.
src_planet = sys.argv[3]  #'Jupiter'

x0,y0,sigma_x,sigma_y,phi,samplerate = lib_p.gen_InputParams(nu_obs,Dapt_mm)
beam_solid_str = lib_p.beam_solidangle_AirySymmetric(nu_obs,Dapt_mm*1e-3)

res = sigma_x*0.1
noise = np.sqrt(0.5)*noise_info(sys.argv[2])/(res*radeg*60.)
map_width_rad = sigma_x*20.

Amp = lib_p.planet_info(src_planet,nu_obs,beam_solid_str)
par_in = np.array([x0,y0,sigma_x,sigma_y,phi,Amp])


elip = lib_p.ellipticalGaussian()
elip.resol_rad = res
elip.map_width_rad = map_width_rad
elip.par = par_in
# elip.beam_sigma = sigma_x
# X, Y, MAP_S = elip.gen_flatellip_map()
X, Y, MAP_S = elip.gen_flatairy_map(nu_obs,Dapt_mm/2.*1e-3)
MAP_S = MAP_S*Amp

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
py.xlabel('$\\theta_x$ [degs]')
py.ylabel('$\\theta_y$ [degs]')
py.title('Map space beam, '+src_planet)

py.subplot(232)
Z = MAP_N
xmin, xmax, ymin, ymax = np.min(X), np.max(X), np.min(Y), np.max(Y)
extent = xmin, xmax, ymin, ymax
im1 = py.imshow(Z, extent=extent)
py.colorbar()
py.xlabel('$\\theta_x$ [degs]')
py.ylabel('$\\theta_y$ [degs]')
py.title('Map space noise')

py.subplot(233)
Z = MAP_SN
xmin, xmax, ymin, ymax = np.min(X), np.max(X), np.min(Y), np.max(Y)
extent = xmin, xmax, ymin, ymax
im1 = py.imshow(Z, extent=extent)
py.colorbar()
py.xlabel('$\\theta_x$ [degs]')
py.ylabel('$\\theta_y$ [degs]')
py.title('Map space beam+noise')

py.subplot(234)
xmin, xmax, ymin, ymax = np.min(kx), np.max(kx), np.min(ky), np.max(ky)
extent = xmin, xmax, ymin, ymax
im1 = py.imshow(np.abs(Bkk_S), extent=extent)
cbar=py.colorbar(im1)
#cbar.ax.ticklabel_format(style='sci', scilimits=(0,0)) 
py.xlabel('$k_x$')
py.ylabel('$k_y$')
py.title('Fourier space S+N')

#py.subplot(235)
#xmin, xmax, ymin, ymax = np.min(kx), np.max(kx), np.min(ky), np.max(ky)
#extent = xmin, xmax, ymin, ymax
#im1 = py.imshow(np.abs(Bkk_N), extent=extent)
#py.colorbar()
#py.xlabel('$\\theta_x$')
#py.ylabel('$\\theta_y$')
#py.title('Fourier space noise')

#++++++++
ind = np.isnan(Bk_SN_mean)
tmp = 0.
tmp_i = 0.
for i in range(len(Bk_SN_kr)):
    if ((ind[i] == False) & ((i>80) & (i<97))):
        tmp += Bk_SN_mean[i]/max(Bk_SN_mean)
        tmp_i += 1.
if tmp_i == 0: tmp_i = 1.e-30
noise_floor = tmp/tmp_i
#np.savez( 'npy/'+src_planet+'_'+str(nu_obs*1e-9)+'GHz_'+str(Dapt_mm)+'mm',
#          noise=noise, noise_floor=noise_floor, Dapt_mm=Dapt_mm, nu_obs=nu_obs, src_planet=src_planet,
#          help='noise, noise_floor, D_mm, nu_obs, src_planet')

py.subplot(235)
py.errorbar(np.array(Bk_S_kr), Bk_S_mean/max(Bk_S_mean), Bk_S_std/max(Bk_S_mean), fmt='o')
#py.semilogx()
#py.errorbar(Bk_N_kr, Bk_N_mean/max(Bk_N_mean), Bk_N_std/max(Bk_N_mean),fmt='o')
#py.semilogx()
py.errorbar(Bk_SN_kr, Bk_SN_mean/max(Bk_SN_mean), Bk_SN_std/max(Bk_SN_mean),fmt='.')
py.text(200,0.1, 'noise floor= %02.2e' % (noise_floor))
py.text(200,0.03, '$D$ = %02.2d mm' % (Dapt_mm))
py.text(200,0.01, '$\\nu$ = %02.2d GHz' % (nu_obs*1e-9))
py.semilogy()
py.ylim([1e-8,1.1])
py.xlim([0,1000])
py.xlabel('$|k|$')
py.ylabel('$B_k$')
py.title('Beam in Fourier space')

#++++++++
py.subplot(236)
#py.errorbar(np.array(Bk_S_kr), (Bk_SN_mean/max(Bk_SN_mean) - Bk_S_mean/max(Bk_S_mean))/(Bk_S_mean/max(Bk_S_mean)), Bk_S_std/max(Bk_S_mean), fmt='o')
py.plot(np.array(Bk_S_kr), (Bk_SN_mean/max(Bk_SN_mean) - Bk_S_mean/max(Bk_S_mean))/(Bk_S_mean/max(Bk_S_mean)), 'o')
#py.semilogx()
#py.errorbar(Bk_N_kr, Bk_N_mean/max(Bk_N_mean), Bk_N_std/max(Bk_N_mean),fmt='o')
#py.semilogx()
#py.errorbar(Bk_SN_kr, Bk_SN_mean/max(Bk_SN_mean), Bk_SN_std/max(Bk_SN_mean),fmt='.')
py.semilogy()
py.ylim([1e-6,1e4])
py.xlim([0,1000])
py.xlabel('$|k|$')
py.ylabel('$\\Delta B_k/B_k$')
py.title('Beam in Fourier space')

py.savefig( 'png/'+src_planet+'_'+str(nu_obs*1e-9)+'GHz_'+str(Dapt_mm)+'mm.png' )
#np.savez( 'png/'+src_planet+'_'+str(nu_obs*1e-9)+'_'+str(Dapt_mm) )
np.savez( 'npy/'+src_planet+'_'+str(nu_obs*1e-9)+'GHz_'+str(Dapt_mm)+'mm',
          noise=noise, noise_floor=noise_floor, Dapt_mm=Dapt_mm, nu_obs=nu_obs, src_planet=src_planet,
          help='noise, noise_floor, D_mm, nu_obs, src_planet')

#py.show()
