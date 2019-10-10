import numpy as np

pi = np.pi
c = 2.99792458e8 

def nuobs_FreqWafer_LFT(FreqWafer):
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

def nuobs_FreqWafer_HFT(FreqWafer):
    if FreqWafer == 'p100': return 100.
    if FreqWafer == 'p119': return 119.
    if FreqWafer == 'p140': return 140.
    if FreqWafer == 'p166': return 166.
    if FreqWafer == 'p195': return 195.
    if FreqWafer == 'p235': return 235
    if FreqWafer == 'p280': return 280.
    if FreqWafer == 'p337': return 337.
    if FreqWafer == 'p402': return 402.

def Dpixmm_info_PhaseA1_LFTv24(FreqWafer):
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

def Telescope_info(LHFT):
	if LHFT == 'LFT': Fnum = 3.5; BeamWaistFact = 2.6
	if LHFT == 'reflect_HFT_MF': Fnum = 3.5; BeamWaistFact = 2.6
	if LHFT == 'refract_HFT_HF': Fnum = 2.2; BeamWaistFact = 3.1
	return Fnum, BeamWaistFact

def CalcSpillOver_v24(Freq,Fnum,Dpixmm,BeamWaistFact):
	wavelength = c/Freq
	out = 1.-np.exp(-pi**2/2.*(Dpixmm*1e-3/BeamWaistFact/Fnum/wavelength)**2)
	return out

def CalcSEdgeTaper_v24(Freq,Fnum,Dpixmm,BeamWaistFact):
	wavelength = c/Freq
	out = np.exp(-pi**2/2.*(Dpixmm*1e-3/BeamWaistFact/Fnum/wavelength)**2)
	return out

def CalcEdgeTaper_dB(Te):
	Te = CalcSpillOver_v24(nu_obs*1.e9,Fnum,Dpixmm,BeamWaistFact)
	Te_dB = np.log(Te)*10.
	return Te_dB

FreqWafer = 'p40'
nu_obs = nuobs_FreqWafer_LFT(FreqWafer)
Dpixmm = Dpixmm_info_PhaseA1_LFTv24(FreqWafer)
Fnum, BeamWaistFact = Telescope_info('LFT')

Te = CalcSpillOver_v24(nu_obs*1.e9,Fnum,Dpixmm,BeamWaistFact)
print np.log10(Te)*10.





