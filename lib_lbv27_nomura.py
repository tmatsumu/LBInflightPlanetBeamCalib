import numpy as np

pi = np.pi
radeg = (180./pi)
c = 2.99792458e8 

def Gen_ArrFreqWafer(Telescope):
    if Telescope == 'LFT':
        return ['p40','p50','p60','p68a','p68b','p78a','p78b','p89a','p89b','p100','p119','p140']
    if Telescope == 'MFT':
        return ['p100','p119','p140','p166','p195']
    if Telescope == 'HFT':
        return ['p195','p235','p280','p337','p402']

def nuobsGHz_FreqWafer(FreqWafer,Telescope):
    if Telescope == 'LFT':
        if FreqWafer == 'p40': return 40.
        if FreqWafer == 'p50': return 50.
        if FreqWafer == 'p60': return 60.
        if FreqWafer == 'p68a': return 68.
        if FreqWafer == 'p68b': return 68.
        if FreqWafer == 'p78a': return 78.
        if FreqWafer == 'p78b': return 78.
        if FreqWafer == 'p89a': return 89.
        if FreqWafer == 'p89b': return 89.
        if FreqWafer == 'p100': return 100.
        if FreqWafer == 'p119': return 119.
        if FreqWafer == 'p140': return 140.

    if Telescope == 'MFT':
        if FreqWafer == 'p100': return 100.
        if FreqWafer == 'p119': return 119.
        if FreqWafer == 'p140': return 140.
        if FreqWafer == 'p166': return 166.
        if FreqWafer == 'p195': return 195.
    
    if Telescope == 'HFT':
        if FreqWafer == 'p195': return 195.
        if FreqWafer == 'p235': return 235.
        if FreqWafer == 'p280': return 280.
        if FreqWafer == 'p337': return 337.
        if FreqWafer == 'p402': return 402.

def Dpixmm_info_PhaseA2_v27(FreqWafer,Telescope):
    if Telescope == 'LFT':
        if FreqWafer == 'p40': return 23.6
        if FreqWafer == 'p50': return 23.6
        if FreqWafer == 'p60': return 23.6
        if FreqWafer == 'p68a': return 23.6
        if FreqWafer == 'p68b': return 15.6
        if FreqWafer == 'p78a': return 23.6
        if FreqWafer == 'p78b': return 15.6
        if FreqWafer == 'p89a': return 23.6
        if FreqWafer == 'p89b': return 15.6
        if FreqWafer == 'p100': return 15.6
        if FreqWafer == 'p119': return 15.6
        if FreqWafer == 'p140': return 15.6

    if Telescope == 'MFT':
        if FreqWafer == 'p100': return 11.6
        if FreqWafer == 'p119': return 11.6
        if FreqWafer == 'p140': return 11.6
        if FreqWafer == 'p166': return 11.6
        if FreqWafer == 'p195': return 11.6
    
    if Telescope == 'HFT':
        if FreqWafer == 'p195': return 6.6
        if FreqWafer == 'p235': return 6.6
        if FreqWafer == 'p280': return 6.6
        if FreqWafer == 'p337': return 6.6
        if FreqWafer == 'p402': return 5.7
    
def Telescope_info(LMHFT):
    if LMHFT == 'LFT': Fnum = 3.0; BeamWaistFact = 2.75; radius = 400.e-3/2.
    if LMHFT == 'MFT': Fnum = 2.2; BeamWaistFact = 2.75; radius = 300.e-3/2.
    if LMHFT == 'HFT': Fnum = 2.2; BeamWaistFact = 3.1; radius = 200.e-3/2.
    return Fnum, BeamWaistFact, radius

def noise_info_PhaseA2(nu_obsGHz,Telescope):
    if Telescope == 'LFT':
        if nu_obsGHz == 40: NETarr = 21.32; det_num = 64.
        if nu_obsGHz == 50: NETarr = 13.81; det_num = 64.
        if nu_obsGHz == 60: NETarr = 11.09; det_num = 64.
        if nu_obsGHz == 68: NETarr = 6.82; det_num = 208.
        if nu_obsGHz == 78: NETarr = 5.57; det_num = 208.
        if nu_obsGHz == 89: NETarr = 4.8; det_num = 208.
        if nu_obsGHz == 100: NETarr = 5.59; det_num = 144.
        if nu_obsGHz == 119: NETarr = 4.17; det_num = 144.
        if nu_obsGHz == 140: NETarr = 3.92; det_num = 144.
        return NETarr, det_num

    if Telescope == 'MFT':
        if nu_obsGHz == 100: NETarr = 4.38; det_num = 366.
        if nu_obsGHz == 119: NETarr = 2.76; det_num = 488.
        if nu_obsGHz == 140: NETarr = 2.98; det_num = 366.
        if nu_obsGHz == 166: NETarr = 2.61; det_num = 488.
        if nu_obsGHz == 195: NETarr = 3.57; det_num = 366.
        return NETarr, det_num

    if Telescope == 'HFT':
        if nu_obsGHz == 195: NETarr = 5.05; det_num = 254.
        if nu_obsGHz == 235: NETarr = 5.21; det_num = 254.
        if nu_obsGHz == 280: NETarr = 6.92; det_num = 254.
        if nu_obsGHz == 337: NETarr = 10.22; det_num = 254.
        if nu_obsGHz == 402: NETarr = 23.34; det_num = 338.
        return NETarr, det_num

def beamarcmin_info_PhaseA2(nu_obsGHz,Telescope):
    if Telescope == 'LFT':
        if nu_obsGHz == 40: beamarcmin = 69.3
        if nu_obsGHz == 50: beamarcmin = 56.8
        if nu_obsGHz == 60: beamarcmin = 49.0
        if nu_obsGHz == 68: beamarcmin = 41.6
    #    if nu_obsGHz == 68: beamarcmin = 44.5
        if nu_obsGHz == 78: beamarcmin = 36.9
    #    if nu_obsGHz == 78: beamarcmin = 40.0
        if nu_obsGHz == 89: beamarcmin = 33.0
    #    if nu_obsGHz == 89: beamarcmin = 36.7
        if nu_obsGHz == 100: beamarcmin = 30.2
        if nu_obsGHz == 119: beamarcmin = 26.3
        if nu_obsGHz == 140: beamarcmin = 23.7
        return beamarcmin

    if Telescope == 'MFT':
        if nu_obsGHz == 100: beamarcmin = 37.8
        if nu_obsGHz == 119: beamarcmin = 33.6
        if nu_obsGHz == 140: beamarcmin = 30.8
        if nu_obsGHz == 166: beamarcmin = 28.9
        if nu_obsGHz == 195: beamarcmin = 28.0
        return beamarcmin

    if Telescope == 'HFT':
        if nu_obsGHz == 195: beamarcmin = 28.6
        if nu_obsGHz == 235: beamarcmin = 24.7
        if nu_obsGHz == 280: beamarcmin = 22.5
        if nu_obsGHz == 337: beamarcmin = 20.9
        if nu_obsGHz == 402: beamarcmin = 17.9
        return beamarcmin

def gen_cambus(beam_arcmin):
    x0 = 0.
    y0 = 0.
    sigma_x = (beam_arcmin/np.sqrt(8.*np.log(2.))/60./radeg)
    sigma_y = (beam_arcmin/np.sqrt(8.*np.log(2.))/60./radeg)
    phi = (45./radeg)
    samplerate = 10. # in unit of Hz
    return x0, y0, sigma_x, sigma_y, phi, samplerate

def CalcSEdgeTaper_v27(Freq,Fnum,Dpixmm,BeamWaistFact):
    wavelength = c/Freq
    out = np.exp(-pi**2/2.*(Dpixmm*1e-3/BeamWaistFact/Fnum/wavelength)**2)
    return out

def CalcSpillOver_v27(Freq,Fnum,Dpixmm,BeamWaistFact):
    wavelength = c/Freq
    SpillOver = 1.-np.exp(-pi**2/2.*(Dpixmm*1e-3/BeamWaistFact/Fnum/wavelength)**2)
    return SpillOver

def CalcEdgeTaper_dB(Te):
    Te_dB = np.log(Te)*10.
    return Te_dB
