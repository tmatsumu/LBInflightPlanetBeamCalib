import numpy as np
import pylab as py
import sys
import glob
import lib_planetcal as lib_p
import scipy.special as sp

def flattruncGauss(alpha, beta):
	sumout = np.zeros(len(beta))
	m_max = 10
	n_max = 10
	for m in range(1,m_max):
		for n in range(1,n_max):
			sumout = sumout + (2.*alpha)**(m+n)*sp.jv(m,beta)/beta**m * sp.jv(n,beta)/beta**n
	print sumout
	out = np.exp(-2.*alpha)/(1.-np.exp(-alpha))**2 * sumout
	return out

c = 2.99792458e8 
#R = Dapt / 2.
pi = np.pi
radeg = (180./pi)

#freq = 
#wavelength = c/freq
#Te = 0.
wavelength = 300.e-6
R = 5000.e-3

theta = np.linspace(1,40.,1000) /3600.*pi/180.

#alpha = np.log(Te)/(-2.)
alpha = 1.e-10
beta = 2.*pi/wavelength * R * np.sin(theta)

out = flattruncGauss(alpha, beta)

print len(theta), len(beta), len(out)
py.plot(theta*radeg*3600., out)
print theta, theta*radeg*3600., out
py.semilogy()
py.ylim([1e-6,1])

py.show()


