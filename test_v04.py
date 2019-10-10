import numpy as np
import lib_planetcal as lib_p

src_planet = 'Jupiter'
nu_obs = 40.e9
beam_solid_str = 2.481e-4

print lib_p.planet_info(src_planet,nu_obs,beam_solid_str)