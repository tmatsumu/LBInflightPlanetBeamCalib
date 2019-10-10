import numpy as np

pi = np.pi
radeg = (180./pi)
c = 2.99792458e8 

def Gen_ArrFreqWafer(Telescope):
    if Telescope == 'LFT':
        return ['p40','p50','p60','p68a','p68b','p78a','p78b','p89a','p89b','p100','p119','p140']
    if Telescope == 'reflect_MFT':
        return ['p100','p119','p140','p166','p195','p235','p280','p337','p402']
    if Telescope == 'reflect_HFT':
        return ['p100','p119','p140','p166','p195','p235','p280','p337','p402']

def nuobsGHz_FreqWafer(FreqWafer,Telescope):
    if Telescope == 'LFT':
        if FreqWafer == 'p40': return 40.
        if FreqWafer == 'p50': return 50.
        if FreqWafer == 'p60': return 60.
        if FreqWafer == 'p68a': return 68.
        if FreqWafer == 'p68b': return 68.
        if FreqWafer == 'p78a': return 78
        if FreqWafer == 'p78b': return 78.
        if FreqWafer == 'p89a': return 89.
        if FreqWafer == 'p89b': return 89.
        if FreqWafer == 'p100': return 100.
        if FreqWafer == 'p119': return 119.
        if FreqWafer == 'p140': return 140.
    if Telescope == 'MFT':
        if FreqWafer == 'p40': return 40.
        if FreqWafer == 'p50': return 50.
        if FreqWafer == 'p60': return 60.
        if FreqWafer == 'p68a': return 68.
        if FreqWafer == 'p68b': return 68.
        if FreqWafer == 'p78a': return 78
        if FreqWafer == 'p78b': return 78.
        if FreqWafer == 'p89a': return 89.
        if FreqWafer == 'p89b': return 89.
        if FreqWafer == 'p100': return 100.
        if FreqWafer == 'p119': return 119.
        if FreqWafer == 'p140': return 140.
    if Telescope == 'HFT':
        if FreqWafer == 'p40': return 40.
        if FreqWafer == 'p50': return 50.
        if FreqWafer == 'p60': return 60.
        if FreqWafer == 'p68a': return 68.
        if FreqWafer == 'p68b': return 68.
        if FreqWafer == 'p78a': return 78
        if FreqWafer == 'p78b': return 78.
        if FreqWafer == 'p89a': return 89.
        if FreqWafer == 'p89b': return 89.
        if FreqWafer == 'p100': return 100.
        if FreqWafer == 'p119': return 119.
        if FreqWafer == 'p140': return 140.


def nuobsGHz_FreqWafer_LFT(FreqWafer):
    if FreqWafer == 'p40': return 40.
    if FreqWafer == 'p50': return 50.
    if FreqWafer == 'p60': return 60.
    if FreqWafer == 'p68a': return 68.
    if FreqWafer == 'p68b': return 68.
    if FreqWafer == 'p78a': return 78
    if FreqWafer == 'p78b': return 78.
    if FreqWafer == 'p89a': return 89.
    if FreqWafer == 'p89b': return 89.
    if FreqWafer == 'p100': return 100.
    if FreqWafer == 'p119': return 119.
    if FreqWafer == 'p140': return 140.

def nuobsGHz_FreqWafer_MFT(FreqWafer):
    if FreqWafer == 'p100': return 100.
    if FreqWafer == 'p119': return 119.
    if FreqWafer == 'p140': return 140.
    if FreqWafer == 'p166': return 166.
    if FreqWafer == 'p195': return 195.
    if FreqWafer == 'p235': return 235
    if FreqWafer == 'p280': return 280.
    if FreqWafer == 'p337': return 337.
    if FreqWafer == 'p402': return 402.

def nuobsGHz_FreqWafer_HFT(FreqWafer):
    if FreqWafer == 'p100': return 100.
    if FreqWafer == 'p119': return 119.
    if FreqWafer == 'p140': return 140.
    if FreqWafer == 'p166': return 166.
    if FreqWafer == 'p195': return 195.
    if FreqWafer == 'p235': return 235
    if FreqWafer == 'p280': return 280.
    if FreqWafer == 'p337': return 337.
    if FreqWafer == 'p402': return 402.

def Dpixmm_info_PhaseA1_LFTv27(FreqWafer):
    if FreqWafer == 'p40': return 30.
    if FreqWafer == 'p50': return 30.
    if FreqWafer == 'p60': return 30.
    if FreqWafer == 'p68a': return 30.
    if FreqWafer == 'p68b': return 18.
    if FreqWafer == 'p78a': return 30
    if FreqWafer == 'p78b': return 18.
    if FreqWafer == 'p89a': return 30.
    if FreqWafer == 'p89b': return 18.
    if FreqWafer == 'p100': return 18.
    if FreqWafer == 'p119': return 18.
    if FreqWafer == 'p140': return 18.

def Dpixmm_info_PhaseA1_MFTv27(FreqWafer):
    if FreqWafer == 'p100': return 12.
    if FreqWafer == 'p119': return 12.
    if FreqWafer == 'p140': return 12.
    if FreqWafer == 'p166': return 12.
    if FreqWafer == 'p195': return 12.
    if FreqWafer == 'p235': return 12
    if FreqWafer == 'p280': return 5.2
    if FreqWafer == 'p337': return 5.2
    if FreqWafer == 'p402': return 5.2

def Dpixmm_info_PhaseA1_HFTv27(FreqWafer):
    if FreqWafer == 'p100': return 12.
    if FreqWafer == 'p119': return 12.
    if FreqWafer == 'p140': return 12.
    if FreqWafer == 'p166': return 12.
    if FreqWafer == 'p195': return 12.
    if FreqWafer == 'p235': return 12
    if FreqWafer == 'p280': return 5.2
    if FreqWafer == 'p337': return 5.2
    if FreqWafer == 'p402': return 5.2

def Telescope_info(LMHFT):
	if LMHFT == 'LFT': Fnum = 3.5; BeamWaistFact = 2.6; radius = 400.e-3/2.
    if LMHFT == 'MFT': Fnum = 3.5; BeamWaistFact = 2.6; radius = 300.e-3/2.
    if LMHFT == 'HFT': Fnum = 3.5; BeamWaistFact = 2.6; radius = 200.e-3/2.
	return Fnum, BeamWaistFact, radius

def noise_info_PhaseA1_LFT(nu_obsGHz):
    if nu_obsGHz == 40: NETarr = 17.3; det_num = 42.
    if nu_obsGHz == 50: NETarr = 9.4; det_num = 56.
    if nu_obsGHz == 60: NETarr = 9.7; det_num = 42.
    if nu_obsGHz == 68: NETarr = 5.4; det_num = 170.
    if nu_obsGHz == 78: NETarr = 4.9; det_num = 156.
    if nu_obsGHz == 89: NETarr = 4.0; det_num = 170.
    if nu_obsGHz == 100: NETarr = 4.8; det_num = 114.
    if nu_obsGHz == 119: NETarr = 3.8; det_num = 114.
    if nu_obsGHz == 140: NETarr = 3.6; det_num = 114.
    return NETarr, det_num

def noise_info_PhaseA1_HFT(nu_obsGHz):
    if nu_obsGHz == 100: NETarr = 4.6; det_num = 222.
    if nu_obsGHz == 119: NETarr = 4.1; det_num = 148.
    if nu_obsGHz == 140: NETarr = 2.9; det_num = 222.
    if nu_obsGHz == 166: NETarr = 3.4; det_num = 148.
    if nu_obsGHz == 195: NETarr = 2.8; det_num = 222.
    if nu_obsGHz == 235: NETarr = 3.8; det_num = 148.
    if nu_obsGHz == 280: NETarr = 4.4; det_num = 338.
    if nu_obsGHz == 337: NETarr = 5.5; det_num = 338.
    if nu_obsGHz == 402: NETarr = 9.4; det_num = 338.
    return NETarr, det_num

def beamarcmin_info_PhaseA1_LFT(nu_obsGHz):
    if nu_obsGHz == 40: beamarcmin = 69.2
    if nu_obsGHz == 50: beamarcmin = 56.9
    if nu_obsGHz == 60: beamarcmin = 49.0
    if nu_obsGHz == 68: beamarcmin = 40.8
    if nu_obsGHz == 78: beamarcmin = 36.1
    if nu_obsGHz == 89: beamarcmin = 32.3
    if nu_obsGHz == 100: beamarcmin = 27.7
    if nu_obsGHz == 119: beamarcmin = 23.7
    if nu_obsGHz == 140: beamarcmin = 20.7
    return beamarcmin

def beamarcmin_info_PhaseA1_HFT(nu_obsGHz):
    if nu_obsGHz == 100: beamarcmin = 37.0
    if nu_obsGHz == 119: beamarcmin = 31.6
    if nu_obsGHz == 140: beamarcmin = 27.6
    if nu_obsGHz == 166: beamarcmin = 24.2
    if nu_obsGHz == 195: beamarcmin = 21.7
    if nu_obsGHz == 235: beamarcmin = 19.6
    if nu_obsGHz == 280: beamarcmin = 13.2
    if nu_obsGHz == 337: beamarcmin = 11.2
    if nu_obsGHz == 402: beamarcmin = 9.7
    return beamarcmin

def gen_cambus(beam_arcmin):
    x0 = 0.
    y0 = 0.
    sigma_x = (beam_arcmin/np.sqrt(8.*np.log(2.))/60./radeg)
    sigma_y = (beam_arcmin/np.sqrt(8.*np.log(2.))/60./radeg)
    phi = (45./radeg)
    samplerate = 10. # in unit of Hz
    return x0, y0, sigma_x, sigma_y, phi, samplerate

def CalcSEdgeTaper_v24(Freq,Fnum,Dpixmm,BeamWaistFact):
	wavelength = c/Freq
	out = np.exp(-pi**2/2.*(Dpixmm*1e-3/BeamWaistFact/Fnum/wavelength)**2)
	return out

def CalcSpillOver_v24(Freq,Fnum,Dpixmm,BeamWaistFact):
	wavelength = c/Freq
	SpillOver = 1.-np.exp(-pi**2/2.*(Dpixmm*1e-3/BeamWaistFact/Fnum/wavelength)**2)
	return SpillOver

def CalcEdgeTaper_dB(Te):
	Te_dB = np.log(Te)*10.
	return Te_dB
