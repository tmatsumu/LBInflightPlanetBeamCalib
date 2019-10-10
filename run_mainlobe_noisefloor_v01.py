import numpy as np
import pylab as py
import lib_lbv24 as lbv24
import main_mapbase_sidelobe_v9 as v9

src_planet_arr = ['Jupiter', 'Saturn', 'Mars']
psm_arr = ['r','g','b']
dir_out = '/Users/tomotake_matsumura/Qsync/Projects/Entries/20180305_BeamCalibration/20130801_PlanetCal_v2/'

option_det = 'band'

#++++++++++++++++++++++++++++++++++++++++++++++++++++++
py.figure(1)
num_planet = len(src_planet_arr)
for j in range(0,num_planet):
	src_planet = src_planet_arr[j]
	print ''
	print '#++++++++++++++++++++++++++++++++++++++++++++++++++'
	print src_planet
	#++++++++++++++++++++++++++++++++++++++++++++++++++++++
	Telescope = 'LFT'
	FreqWafer_arr = lbv24.Gen_ArrFreqWafer(Telescope)
	num = len(FreqWafer_arr)
	diameter_mm = 400.

	band_LFT = np.zeros(num)
	LFTnoisefloor_arr = np.zeros((2,num))
	for i in range(0,num):
		band_LFT[i] = lbv24.nuobsGHz_FreqWafer_LFT(FreqWafer_arr[i])
		if option_det=='det':
			LFTnoisefloor_arr[0,i] = v9.run_mainloadSN(Telescope,FreqWafer_arr[i],diameter_mm,src_planet,option_noise='det')
		if option_det=='band':
			LFTnoisefloor_arr[1,i] = v9.run_mainloadSN(Telescope,FreqWafer_arr[i],diameter_mm,src_planet,option_noise='band')

	#++++++++++++++++++++++++++++++++++++++++++++++++++++++
	Telescope = 'reflect_HFT_MF'
	FreqWafer_arr = lbv24.Gen_ArrFreqWafer(Telescope)
	num = len(FreqWafer_arr)
	diameter_mm = 300.

	band_HFT = np.zeros(num)
	HFTnoisefloor_arr = np.zeros((2,num))
	for i in range(0,num):
		band_HFT[i] = lbv24.nuobsGHz_FreqWafer_HFT(FreqWafer_arr[i])
		print band_HFT[i]
		if option_det=='det':
			HFTnoisefloor_arr[0,i] = v9.run_mainloadSN(Telescope,FreqWafer_arr[i],diameter_mm,src_planet,option_noise='det')
		if option_det=='band':
			HFTnoisefloor_arr[1,i] = v9.run_mainloadSN(Telescope,FreqWafer_arr[i],diameter_mm,src_planet,option_noise='band')

	#++++++++++++++++++++++++++++++++++++++++++++++++++++++
	py.figure(1)
	if option_det=='det':
		py.plot(band_LFT,np.log10(LFTnoisefloor_arr[0])*10,'.'+psm_arr[j],label='LFT detector each')
		py.plot(band_HFT,np.log10(HFTnoisefloor_arr[0])*10,'.'+psm_arr[j],label='HFT detector each')

	if option_det=='band':
		py.plot(band_LFT,np.log10(LFTnoisefloor_arr[1])*10,'^'+psm_arr[j],label='LFT band average, '+src_planet)
		py.plot(band_HFT,np.log10(HFTnoisefloor_arr[1])*10,'v'+psm_arr[j],label='HFT band average, '+src_planet)

	print '#++++++++++++++++++++++++++++++++++++++++++++++++++'
	print ''

py.xlim([0,500])
py.ylim([-70,-20])

py.xlabel('Frequency [GHz]')
py.ylabel('Noise floor [dB]')
py.legend(loc='best')
py.grid()
py.savefig(dir_out+'/output/TruncGaussian/ConfigVer20180422_noisefloor.png')
