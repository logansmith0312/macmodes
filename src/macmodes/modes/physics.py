import numpy as np
from macmodes.constants import SECONDS_PER_YEAR, ETA, R_CMB

def get_delta(period_years: int) -> float:
    '''given mode period in years returns skin layer depth'''
    freq = 2 * np.pi/(period_years*SECONDS_PER_YEAR)
    delta = np.sqrt(2*ETA/freq)
    return delta/R_CMB

# Linearly varying Buoyancy frequency
def get_N_linear(r, N_max, H):
    return N_max * ( (r-(1-H)) / H )
