import numpy as np
import pylab as py
import sys

pi = np.pi
radeg = (180./pi)

nu_obs = float(sys.argv[2]) #140.e9
nu_obs = nu_obs*1.e9
Dapt_mm = float(sys.argv[1])  #300.
src_planet = sys.argv[3]  #'Jupiter'
dir_out = sys.argv[4]
option = sys.argv[5]

if option=='airyfunction':
	out = np.load(dir_out+'/AiryFunction/npy/'+src_planet+'_'+str(nu_obs*1e-9)+'GHz_'+str(Dapt_mm)+'mm.npz')
if option=='ellipticalGaussian':
	out = np.load(dir_out+'/ellipticalGaussian/npy/'+src_planet+'_'+str(nu_obs*1e-9)+'GHz_'+str(Dapt_mm)+'mm.npz')
par_out = out['par_out']
par_err = np.sqrt(np.diag(out['par_err']))
nu_obs = out['nu_obs']
src_planet = out['src_planet']
Dapt_mm = out['Dapt_mm']

#print 'Band[GHz], Dapt[mm], src, (x0,y0)[arcmin], frac.FWHM_x[%], frac.FWHM_y[%] '
print '%s %1.0d   %1.2f   %s   %1.3f   %1.3f   %1.3e   %1.3e' % \
	(option, nu_obs*1e-9, Dapt_mm, src_planet, par_err[0]*radeg*3600., par_err[1]*radeg*3600., \
		par_err[2]/par_out[2], par_err[3]/par_out[3])

