import numpy as np
from macmodes.constants import R_CMB

def uniform_grid(Ntot: int, r0: float) -> np.ndarray[float]:
    '''Uniform grid'''
    return np.linspace(r0, R_CMB, Ntot)

def cheb_grid(Ntot: int, r0: float) -> np.ndarray[float]:
    '''Chebyshev nodes of the second kind over [r0, R_CMB]'''
    
    Nmax: int = Ntot - 1
    r: np.ndarray = np.zeros(Ntot)
    
    for k in range(1,Nmax):
        r[k] = 0.5*(+r0)-0.5*(R_CMB-r0)*np.cos(k*np.pi/Nmax)
    r[0] = r0
    r[Nmax] = R_CMB
    return r

def get_grid(Ntot: int, r0: float, gridflg: int) -> np.ndarray[float]:
    '''Defines either uniform or chebyshev grid'''
    if gridflg == 0:
        r = uniform_grid(Ntot,r0,R_CMB)
    elif gridflg == 1:
        r = cheb_grid(Ntot,r0,R_CMB)
    else:
        raise ValueError(f"Invalid value for gridflg '{gridflg}'")
    return r

     
def get_grid_params(Ntot: int, r0: float,r: np.ndarray) -> list[np.ndarray[float]]:
    '''Defines necessary derivative coefficients'''
    Nmax: int = Ntot - 1
        
    dr = np.zeros(Ntot) # Radial stepsizes
    r_1 = np.zeros(Ntot) 
    r_2 = np.zeros(Ntot)

    for k in range(1, Ntot):
        dr[k] = r[k] - r[k-1]
        r_1[k] = 1/r[k]
        r_2[k] = 1/r[k]**2
    if r0 == 0.:
        r_1[0] = 0.
        r_2[0] = 0.
    else:
        r_1[0] = 1./r0
        r_2[0] = 1./r0**2
    dr[0] = dr[1]

    # Radial derivative functions
    # d/dr (first derivative)
    dr1a = np.zeros(Ntot) 
    dr1b = np.zeros(Ntot)
    dr1c = np.zeros(Ntot)
    
    # d2/dr2 (second derivative)
    dr2a = np.zeros(Ntot)
    dr2b = np.zeros(Ntot)
    dr2c = np.zeros(Ntot) 
    
    # Interior points
    for k in range(1, Nmax):
        dr1a[k] = -dr[k+1] / (dr[k] * (dr[k]+dr[k+1]))
        dr1b[k] = (dr[k+1] - dr[k]) / (dr[k]*dr[k+1])
        dr1c[k] = dr[k] / (dr[k+1] * (dr[k] + dr[k+1]))

        dr2a[k] = 2. / (dr[k] * (dr[k] + dr[k+1]))
        dr2b[k] = -2. / (dr[k] * dr[k+1])
        dr2c[k] = 2./(dr[k+1] * (dr[k] + dr[k+1]))

    # Inner boundary point
    dr1a[0] = - (2.0 * dr[1] + dr[2]) / (dr[1] * (dr[1] + dr[2]))
    dr1b[0] =   (dr[1] + dr[2]) / (dr[1] * dr[2])
    dr1c[0] = - dr[1] / (dr[2] * (dr[1] + dr[2]))

    # These are not used
    dr2a[0] =  2.0 / (dr[1] * (dr[1] + dr[2]))
    dr2b[0] = -2.0 / (dr[1] * dr[2])
    dr2c[0] =  2.0 / (dr[2] * (dr[1] + dr[2]))


    # Radial derivative functions: outer boundary point
    dr1a[Nmax]= (2.*dr[Nmax] + dr[Nmax-1])/( dr[Nmax]*(dr[Nmax]+dr[Nmax-1]))
    dr1b[Nmax]=-(dr[Nmax]+dr[Nmax-1])/(dr[Nmax]*dr[Nmax-1])
    dr1c[Nmax]= dr[Nmax]/( dr[Nmax-1]*(dr[Nmax]+dr[Nmax-1]))
    
    # These are not used
    dr2a[Nmax] = 2./(dr[Nmax] * (dr[Nmax]+dr[Nmax-1]))
    dr2b[Nmax] =-2./(dr[Nmax]*dr[Nmax-1])
    dr2c[Nmax] = 2./(dr[Nmax-1] * (dr[Nmax]+dr[Nmax-1]))
    
    return dr, r_1, r_2, dr1a, dr1b, dr1c, dr2a, dr2b, dr2c
    

