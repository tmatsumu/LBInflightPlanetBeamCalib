import numpy as np
import pylab as py
import sys

dir = '/Users/tomotake_matsumura/Documents/Projects/litebird/20130801_PlanetCal_v2/npy/'

src_planet = ['Jupiter', 'Saturn', 'Mars']
nu_obs = ['60.0', '78.0', '100.0', '140.0', '195.0', '280.0']
Dapt_mm = ['250.0', '300.0', '350.0', '400.0', '450.0', '500.0', '550.0', '600.0']
fmt = ['or','ob','og','om','oc','ok']
num_nu = len(nu_obs)
num_Dapt = len(Dapt_mm)

out_noisefloor = []
out_Daptmm = []
out_nuobs = []
out_noisein = []
out_srcplanet = []

i_src_planet = 0
i_nu_obs = 1
for i in range(0,num_nu):
    for j in range(0,num_Dapt):
        for k in range(i_src_planet,i_src_planet+1):
            out = np.load(dir+src_planet[k]+'_'+nu_obs[i]+'GHz'+'_'+Dapt_mm[j]+'mm.npz')
            out_noisefloor.append(out['noise_floor'])
            out_Daptmm.append(out['Dapt_mm'])
            out_nuobs.append(out['nu_obs'])
            out_noisein.append(out['noise'])
            out_srcplanet.append(out['src_planet'])

#    print out_Daptmm,out_noisefloor
    py.plot(out_Daptmm,out_noisefloor,fmt[i])
py.title(src_planet[k])
py.xlabel('Aperture diameter [mm]')
py.xlim([200,650])
py.semilogy()

py.show()

