import numpy as np

# CMB radius
#R_CMB_DIM: float = 3.485e6
R_CMB_DIM: float = 3.48e6
R_CMB: float = 1.0 # non dimenaional CMB radius
# Time constants
SECONDS_PER_DAY: float = 86400.0
DAYS_PER_YEAR: float = 365.25
SECONDS_PER_YEAR: float = SECONDS_PER_DAY*DAYS_PER_YEAR
OMEGA: float = 2*np.pi / SECONDS_PER_DAY
# E&M/fluid core stuff
MU0: float = 4*np.pi*1e-7
#RHO_CORE: float = 1.1e4
RHO_CORE: float = 1e4
COND_CORE: float = 1e6
ETA: float = 1/(COND_CORE*MU0)
BR_SCALE: float = R_CMB_DIM*OMEGA*np.sqrt(RHO_CORE*MU0)
#BR_SCALE: float = np.sqrt(OMEGA*RHO_CORE*MU0*ETA)

EK_PM: float = ETA/(OMEGA*R_CMB_DIM**2)

COND_MANTLE: float = 1e8
CP: float = COND_MANTLE*MU0*R_CMB_DIM*OMEGA

# path constants
import os

SRC_DIR  = os.path.dirname(os.path.abspath(__file__)) #macmodes/src
MACMODES_DIR = os.path.dirname(os.path.dirname(SRC_DIR)) #macmodes
