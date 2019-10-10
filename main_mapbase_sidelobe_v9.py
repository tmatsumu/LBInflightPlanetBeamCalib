import numpy as np
import pylab as py
import sys
import glob
import lib_planetcal as lib_p
from scipy.optimize import curve_fit
import lib_lbv24 as lbv24


'''
The values below are updated as of 2018-03-05 (PhaseA1)

2018-June-6: v8
The beam shape is now taken from the analysitcala expression written as eq. 3.22 in the H. Imada's thesis.
checked to see if the result is consistent with the past results

2018-June-18: v9
take argument in order to sweep the parameter spaces
'''

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

pi = np.pi
radeg = (180./pi)
#res = .05/radeg
#map_width_rad = 2./radeg
#noise = 5./(res*radeg*60.)*np.sqrt(300.)*0.5

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Telescope = sys.argv[1]
#FreqWafer = sys.argv[2]
#diameter_mm = float(sys.argv[3])
#src_planet = sys.argv[4]  #'Jupiter'
option = 'TruncGaussian'
dir_out = '/Users/tomotake_matsumura/Qsync/Projects/Entries/20180305_BeamCalibration/20130801_PlanetCal_v2/output/'

def run_mainloadSN(Telescope,FreqWafer,diameter_mm,src_planet,option_noise=''):

  #Telescope = 'LFT'
  #FreqWafer = 'p60'
  Fnum, BeamWaistFact, radius = lbv24.Telescope_info(Telescope)
  #radius = 400.e-3/2.
  radius = diameter_mm*1e-3/2.
  Dapt_mm = radius*2.*1e3
  if Telescope == 'LFT':
    nu_obsGHz = lbv24.nuobsGHz_FreqWafer_LFT(FreqWafer)
    Dpixmm = lbv24.Dpixmm_info_PhaseA1_LFTv24(FreqWafer)
    beam_arcmin = lbv24.beamarcmin_info_PhaseA1_LFT(nu_obsGHz)
  if Telescope == 'reflect_HFT_MF':
    nu_obsGHz = lbv24.nuobsGHz_FreqWafer_HFT(FreqWafer)
    Dpixmm = lbv24.Dpixmm_info_PhaseA1_HFTv24(FreqWafer)
    beam_arcmin = lbv24.beamarcmin_info_PhaseA1_HFT(nu_obsGHz)
  nu_obs = nu_obsGHz*1.e9
  edgetaper = lbv24.CalcSEdgeTaper_v24(nu_obs,Fnum,Dpixmm,BeamWaistFact)

  #src_planet = 'Jupiter'
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  print ''
  print '--beam properties--'
  #x0_rad,y0_rad,sigma_x_rad,sigma_y_rad,phi,samplerate = lib_p.gen_InputParams_PhaseA1_LFT(nu_obs,Dapt_mm)
  #x0_rad,y0_rad,sigma_x_rad,sigma_y_rad,phi,samplerate = lbv24.gen_cambus(nu_obsGHz)
  x0_rad,y0_rad,sigma_x_rad,sigma_y_rad,phi,samplerate = lbv24.gen_cambus(beam_arcmin)
  #print sigma_x_rad, sigma_x_rad*0.1
  res_rad = sigma_x_rad*0.1 # res is in unit of rad
  map_width_rad = sigma_x_rad*10. # unit in rad

  # x0, y0, sigma_x, sigma_y are in unit of rad
  par_in = np.array([x0_rad,y0_rad,sigma_x_rad,sigma_y_rad,phi,1.])

  elip = lib_p.ellipticalGaussian()
  elip.resol_rad = res_rad
  elip.map_width_rad = map_width_rad
  elip.par = par_in

  if option == 'TruncGaussian':
      X, Y, MAP_S = elip.gen_flatTruncGauss_map(nu_obs,edgetaper,radius)
      tmp_theta, tmp_out = elip.gen_flatTruncGauss_2D(nu_obs,edgetaper,radius)
      beam_solid_str = 2.*pi*np.sum(np.sin(tmp_theta)*tmp_out)*(tmp_theta[1]-tmp_theta[0])
      print 'beam solid angle', beam_solid_str
  #beam_solid_str = np.sum(MAP_S)
  #beam_solid_str = lib_p.beam_solidangle_AirySymmetric(nu_obs,Dapt_mm*1e-3)
  # x0, y0, sigma_x, sigma_y are in unit of rad

  print 'input:', src_planet,nu_obs,beam_solid_str
  Amp = lib_p.planet_info(src_planet,nu_obs,beam_solid_str)
  MAP_S = MAP_S*Amp

  print 'amp', Amp
  print '--beam properties--'
  print ''
  #sys.exit()
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  #noise = np.sqrt(0.5)*noise_info(sys.argv[2])/(res*radeg*60.)
  #NET_per_det = noise_info_PhaseA1_LFT(sys.argv[2])
  #NET_per_det = noise_info_PhaseA1_HFT(sys.argv[2])
  if Telescope == 'LFT':
    NET_arr, numDet = lbv24.noise_info_PhaseA1_LFT(nu_obsGHz)
  if Telescope == 'reflect_HFT_MF':
    NET_arr, numDet = lbv24.noise_info_PhaseA1_HFT(nu_obsGHz)

  if option_noise == 'det': NET_per_det = NET_arr * np.sqrt(numDet)
  if option_noise == 'band': NET_per_det = NET_arr 
  obs_eff = 0.85*0.85
  time_interval_sec = 3.*obs_eff*3600.*24.*365.
  time2spent_in_one_pixel = res_rad**2/(4.*pi)*time_interval_sec
  noise = NET_per_det/np.sqrt(time2spent_in_one_pixel)
  print ''
  print '--noise properties--'
  print 'NET_per_det', NET_per_det
  print 'res_rad, res_rad**2, time2spent_in_one_pixel', res_rad, res_rad*radeg, res_rad**2, time2spent_in_one_pixel
  print 'noise', noise
  print '--noise properties--'
  print ''
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  #MAP_S = MAP_S*Amp

  num = len(MAP_S[0])
  MAP_N = noise*np.random.normal(0.,1.,[num,num])
  MAP_SN = MAP_S+MAP_N

  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  num_bin = 100
  Bkk_S, kx, ky = lib_p.fft_2d(res_rad,res_rad,MAP_S)
  Bk_S_kr, Bk_S_mean, Bk_S_std, Bk_S_med = lib_p.cal_Bk(num_bin,kx,ky,Bkk_S)

  Bkk_N, kx, ky = lib_p.fft_2d(res_rad,res_rad,MAP_N)
  Bk_N_kr, Bk_N_mean, Bk_N_std, Bk_N_med = lib_p.cal_Bk(num_bin,kx,ky,Bkk_N)

  Bkk_SN, kx, ky = lib_p.fft_2d(res_rad,res_rad,MAP_SN)
  Bk_SN_kr, Bk_SN_mean, Bk_SN_std, Bk_SN_med = lib_p.cal_Bk(num_bin,kx,ky,Bkk_SN)

  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  py.figure(0, figsize=(18,10))

  py.subplot(231)
  X_deg = X*radeg
  Y_deg = Y*radeg
  Z = MAP_S
  xmin, xmax, ymin, ymax = np.min(X_deg), np.max(X_deg), np.min(Y_deg), np.max(Y_deg)
  extent = xmin, xmax, ymin, ymax
  im1 = py.imshow(Z, extent=extent)
  py.colorbar()
  py.xlabel('$\\theta_x$ [degs]')
  py.ylabel('$\\theta_y$ [degs]')
  py.title('Map space beam, '+src_planet)

  py.subplot(232)
  Z = MAP_N
  xmin, xmax, ymin, ymax = np.min(X_deg), np.max(X_deg), np.min(Y_deg), np.max(Y_deg)
  extent = xmin, xmax, ymin, ymax
  im1 = py.imshow(Z, extent=extent)
  py.colorbar()
  py.xlabel('$\\theta_x$ [degs]')
  py.ylabel('$\\theta_y$ [degs]')
  py.title('Map space noise')

  #
  data = lib_p.meshgrid2array_Z(MAP_SN)
  data_err = noise
  def fit_elipbeam_2D( (X_rad,Y_rad), x0_rad, y0_rad, sigma_x_rad, sigma_y_rad, theta_rad, amplitude):
      #par_init,x,y,data,noise_sigma,g):
      elip = lib_p.ellipticalGaussian()
      elip.resol_rad = res_rad
      elip.map_width_rad = map_width_rad
      elip.par = np.array([x0_rad, y0_rad, sigma_x_rad, sigma_y_rad, theta_rad, amplitude])
      Z = amplitude * lib_p.ellipGauss(X_rad,Y_rad,x0_rad, y0_rad, sigma_x_rad, sigma_y_rad, theta_rad)
      z = lib_p.meshgrid2array_Z(Z)
      return np.array(z)

  par_in = np.array([0.,0.,sigma_x_rad,sigma_y_rad,phi,np.max(MAP_SN)])
  par_out, par_err = curve_fit(fit_elipbeam_2D, (X,Y), data, sigma=data_err, p0=par_in)
  print '-----------------------------------------'
  print '(x0,y0) arcsec = ', par_out[0]*radeg*3600., par_out[1]*radeg*3600.
  print 'del (x0,y0) arcsec =', np.sqrt(np.diag(par_err))[0]*radeg*3600., np.sqrt(np.diag(par_err))[1]*radeg*3600.
  print 'FWHM_x, FWHM_y= (', par_out[2]*radeg*60*np.sqrt(8.*np.log(2.)), par_out[3]*radeg*60*np.sqrt(8.*np.log(2.)), ') arcmin'
  print 'theta', par_out[4]
  print 'amplitude=', par_out[5]
  print np.sqrt(np.diag(par_err))[0]*3600., np.sqrt(np.diag(par_err))[1]*3600., 'arcsec', np.sqrt(np.diag(par_err))[2:]
  print np.sqrt(np.diag(par_err))/par_out
  print '-----------------------------------------'

  py.subplot(233)
  Z = MAP_SN
  xmin, xmax, ymin, ymax = np.min(X_deg), np.max(X_deg), np.min(Y_deg), np.max(Y_deg)
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
  #sigma_tmp = 69./60./radeg/np.sqrt(8.*np.log(2.))
  #py.plot(np.array(Bk_S_kr), np.exp(-np.array(Bk_S_kr)**2*sigma_tmp**2))
  #py.semilogx()
  #py.errorbar(Bk_N_kr, Bk_N_mean/max(Bk_N_mean), Bk_N_std/max(Bk_N_mean),fmt='o')
  #py.semilogx()
  py.errorbar(np.array(Bk_SN_kr), Bk_SN_mean/max(Bk_SN_mean), Bk_SN_std/max(Bk_SN_mean),fmt='.')
  py.errorbar(np.array(Bk_N_kr), Bk_N_mean/max(Bk_SN_mean), Bk_N_std/max(Bk_SN_mean),fmt='.')
  py.text(200,0.1, 'noise floor= %02.2e' % (noise_floor))
  py.text(200,0.03, '$D$ = %02.2d mm' % (Dapt_mm))
  py.text(200,0.01, '$\\nu$ = %02.2d GHz' % (nu_obs*1e-9))
  py.semilogy()
  py.ylim([1e-6,1.1])
  py.xlim([0,1000])
  py.xlabel('$|k|$')
  py.ylabel('$B_k$')
  py.title('Beam in Fourier space')

  #++++++++
  py.subplot(236)
  py.errorbar(np.array(Bk_S_kr), Bk_S_mean/max(Bk_S_mean), Bk_S_std/max(Bk_S_mean), fmt='o')
  py.errorbar(np.array(Bk_SN_kr), Bk_SN_mean/max(Bk_SN_mean), Bk_SN_std/max(Bk_SN_mean),fmt='.')
  #py.errorbar(np.array(Bk_S_kr), (Bk_SN_mean/max(Bk_SN_mean) - Bk_S_mean/max(Bk_S_mean))/(Bk_S_mean/max(Bk_S_mean)), Bk_S_std/max(Bk_S_mean), fmt='o')
  #py.plot(np.array(Bk_S_kr), (Bk_SN_mean/max(Bk_SN_mean) - Bk_S_mean/max(Bk_S_mean))/(Bk_S_mean/max(Bk_S_mean)), 'o')
  #py.semilogx()
  #py.errorbar(Bk_N_kr, Bk_N_mean/max(Bk_N_mean), Bk_N_std/max(Bk_N_mean),fmt='o')
  #py.semilogx()
  #py.errorbar(Bk_SN_kr, Bk_SN_mean/max(Bk_SN_mean), Bk_SN_std/max(Bk_SN_mean),fmt='.')
  #py.semilogy()
  py.ylim([1e-6,1.1])
  py.xlim([1,10000])
  py.loglog()
  py.xlabel('$|k|$')
  py.ylabel('$\\Delta B_k/B_k$')
  py.title('Beam in Fourier space')

  if option=='TruncGaussian':
    py.savefig( dir_out+'/TruncGaussian/png/'+src_planet+'_'+str(nu_obs*1e-9)+'GHz_'+str(Dapt_mm)+'mm.png' )
    np.savez( dir_out+'/TruncGaussian/npy/'+src_planet+'_'+str(nu_obs*1e-9)+'GHz_'+str(Dapt_mm)+'mm', \
            noise=noise, noise_floor=noise_floor, Dapt_mm=Dapt_mm, nu_obs=nu_obs, src_planet=src_planet, \
            par_out=par_out, par_err=par_err, \
            help='noise, noise_floor, D_mm, nu_obs, src_planet')
  py.clf()

  if option=='airyfunction':
      py.savefig( dir_out+'/AiryFunction/png/'+src_planet+'_'+str(nu_obs*1e-9)+'GHz_'+str(Dapt_mm)+'mm.png' )
  #np.savez( 'png/'+src_planet+'_'+str(nu_obs*1e-9)+'_'+str(Dapt_mm) )
      np.savez( dir_out+'/AiryFunction/npy/'+src_planet+'_'+str(nu_obs*1e-9)+'GHz_'+str(Dapt_mm)+'mm', \
                noise=noise, noise_floor=noise_floor, Dapt_mm=Dapt_mm, nu_obs=nu_obs, src_planet=src_planet, \
                par_out=par_out, par_err=par_err, \
                help='noise, noise_floor, D_mm, nu_obs, src_planet')

  if option=='ellipticalGaussian':
      py.savefig( dir_out+'/ellipticalGaussian/png/'+src_planet+'_'+str(nu_obs*1e-9)+'GHz_'+str(Dapt_mm)+'mm.png' )
  #np.savez( 'png/'+src_planet+'_'+str(nu_obs*1e-9)+'_'+str(Dapt_mm) )
      np.savez( dir_out+'/ellipticalGaussian/npy/'+src_planet+'_'+str(nu_obs*1e-9)+'GHz_'+str(Dapt_mm)+'mm', \
                noise=noise, noise_floor=noise_floor, Dapt_mm=Dapt_mm, nu_obs=nu_obs, src_planet=src_planet, \
                par_out=par_out, par_err=par_err, \
                help='noise, noise_floor, D_mm, nu_obs, src_planet')

  return noise_floor
