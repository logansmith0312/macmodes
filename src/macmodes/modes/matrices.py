import numpy as np
from macmodes.constants import R_CMB, EK_PM, CP
from macmodes.modes.physics import get_delta, get_N_linear

def c(ell: int) -> float:
    '''Coriolis coefficient'''
    return  np.sqrt( ell**2 / ((2*ell-1) * (2*ell+1)) )


# Matrices
def setup_mats(r, grid_params, config):

    n_rad: int = config.n_rad
    ell_max: int = config.ell_max
    H: float = config.H
    N_max: float = config.N_max
    Br_cmb: float = config.Br_cmb
    Ek: float = config.Ek
    delta_per_est: int = config.delta_per_est

    r0: float = R_CMB - H
    delta = get_delta(delta_per_est)

    dr, r_1, r_2, dr1a, dr1b, dr1c, dr2a, dr2b, dr2c = grid_params
    Nmax: int = n_rad - 1
    ell_max_h: int = ell_max//2

    N: np.ndarray[float] = get_N_linear(r, N_max, H)

    A = np.zeros((5*ell_max_h*n_rad, 5*ell_max_h*n_rad), dtype=complex)
    B = np.zeros((5*ell_max_h*n_rad, 5*ell_max_h*n_rad), dtype=complex)
    
    for nl in range(ell_max_h):
        ello = 2*nl+1 # odd degree spherical harmonic
        ellpo = ello*(ello+1) # odd degree factor ell*(ell+1)
        elle = 2*nl+2 # even degree spherical harmonic
        ellpe = elle*(elle+1) # even degree factor ell*(ell+1)
        
        for k in range(1, Nmax):
            kz = k + nl*n_rad # zeta index
            kw = kz + ell_max_h*n_rad # W index
            kwl = kw - n_rad # W left spherical harmonic index
            kwr = kw  # W right spherical harmonic index
            kd = kw + ell_max_h*n_rad # D index
            kp = kd + ell_max_h*n_rad # Z index
            kpl = kp # Z left spherical harmonic index
            kpr = kp + n_rad # Z right spherical harmonic index
            kt = kp + ell_max_h*n_rad
            
            # zeta equations
            B[kz, kz] = -1j
            
            # Diffusion
            A[kz, kz-1] = Ek*dr2a[k]
            A[kz, kz] = Ek*(dr2b[k] - ellpe*r_2[k])
            A[kz, kz+1] = Ek*dr2c[k]
            
            # Buoyancy
            A[kz, kd] = N[k]**2*ellpe*r_2[k]
            
            # Coriolis
            factCora = 2*(elle - 1)/elle
            clea = c(elle)
            A[kz, kpl-1] = - factCora * clea * dr1a[k]
            A[kz, kpl] = - factCora * clea * (dr1b[k] - elle*r_1[k])
            A[kz, kpl+1] = - factCora * clea * dr1c[k]
            
            factCorb = 2*(elle+2)/(elle+1)
            if elle != ell_max: # Truncate: do not include ell_max+1 spherical harmonic part when ell=ell_max
                cleb = c(elle+1)
                A[kz, kpr-1] = - factCorb * cleb * dr1a[k]
                A[kz, kpr] = - factCorb * cleb * (dr1b[k] + (elle + 1)*r_1[k])
                A[kz, kpr+1] = - factCorb * cleb * dr1c[k]
                
            # W equations (no time derivative part)
            
            # Laplacian part
            A[kw, kw-1] = dr2a[k]
            A[kw, kw] = dr2b[k] - ellpe*r_2[k]
            A[kw, kw+1] = dr2c[k]
            
            # zeta part
            A[kw, kz] = -1.
            
            # D equations 
            B[kd, kd] = -1j
            A[kd, kw] = 1.
            
            # Z equations
            B[kp, kp] = -1j
            
            # Diffusion
            A[kp, kp-1] = Ek*dr2a[k]
            A[kp, kp] = Ek*(dr2b[k] - ellpo*r_2[k])
            A[kp, kp+1] = Ek*dr2c[k]
            
            # Coriolis
            factCora = 2*(ello - 1)/ello
            if ello != 1: # Truncate: do not include ell=0 spherical harmonic part when ell=1
                cloa = c(ello)
                A[kp, kwl-1] = factCora * cloa * dr1a[k]
                A[kp, kwl] = factCora * cloa * (dr1b[k] - ello*r_1[k])
                A[kp, kwl+1] = factCora * cloa * dr1c[k]
                
            factCorb = 2*(ello+2)/(ello+1)
            clob = c(ello+1)
            A[kp, kwr-1] = factCorb * clob * dr1a[k]
            A[kp, kwr] = factCorb * clob * (dr1b[k] + (ello+1)*r_1[k])
            A[kp, kwr+1] = factCorb * clob * dr1c[k]
            
            # Lorentz
            A[kp, kt-1] = Br_cmb * dr1a[k]
            A[kp, kt] = Br_cmb * dr1b[k]
            A[kp, kt+1] = Br_cmb * dr1c[k]
            
            # tau equations
            B[kt, kt] = -1j
            
            # Diffusion
            A[kt, kt-1] = EK_PM*dr2a[k]
            A[kt, kt] = EK_PM*(dr2b[k]-ellpo*r_2[k]) 
            A[kt, kt+1] = EK_PM*dr2c[k]
            
            # Advection
            A[kt, kp-1] = Br_cmb*dr1a[k]
            A[kt, kp] = Br_cmb*dr1b[k]
            A[kt, kp+1] = Br_cmb*dr1c[k]
            
        
        # Boundary conditions at r=r0
        k=0
        kz = nl*n_rad
        kw = kz + ell_max_h*n_rad
        kd = kw + ell_max_h*n_rad
        kp = kd + ell_max_h*n_rad
        kt = kp + ell_max_h*n_rad
        
        # free slip vtheta
        A[kz, kz] = 1.*dr[k]**2 
        A[kz, kw+1] = -1. + (1.-dr[k]*r_1[k]) / (1.+dr[k]*r_1[k])
        # no normal flow vr
        A[kw, kw] = 1. 
        A[kd, kd] = 1. 
        # free slip vphi
        A[kp, kp] = dr1a[k] - 2*r_1[k] 
        A[kp, kp+1] = dr1b[k]
        A[kp, kp+2] = dr1c[k]
        # pseudo vacuum bphi
        #A[kt, kt] = 1. 
        # match potential field of inner core for purely diffusive magnetic field
        A[kt, kt] = r_1[k] - dr1a[k] + (1-1j)/delta
        A[kt, kt+1] = -dr1b[k]
        A[kt, kt+2] = -dr1c[k]
        A[kt, kp] = -Br_cmb/EK_PM
        
        
        # Boundary conditions at r=R_CMB
        k=Nmax
        kz = (nl+1)*n_rad - 1
        kw = kz + ell_max_h*n_rad
        kd = kw + ell_max_h*n_rad
        kp = kd + ell_max_h*n_rad
        kt = kp + ell_max_h*n_rad

        # free slip vtheta
        A[kz, kz] = 1.*dr[k]**2
        A[kz, kw-1] = -1. + (1.+dr[k]*r_1[k]) / (1.-dr[k]*r_1[k])
        # no normal flow vr
        A[kw, kw] = 1.
        A[kd, kd] = 1.
        # free slip vphi
        A[kp, kp] = dr1a[k] - 2*r_1[k]
        A[kp, kp-1] = dr1b[k]
        A[kp, kp-2] = dr1c[k]
        #pseudo vacuum bphi
        A[kt, kt] = 1.
        # For EM coupling at mantle:
        #A[kt, kp] = CP*Br_cmb 
        
    return A, B

