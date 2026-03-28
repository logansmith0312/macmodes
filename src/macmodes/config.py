from dataclasses import dataclass
from macmodes.constants import EK_PM

@dataclass
class ModeConfig:

    n_rad: int = 101 # number of radial grid points
    n_theta: int = 100 # number of colatitudinal grid points
    ell_max: int = 4 # max spherical harmonic degree

    H: float = 140e3/R_CMB # dimensionless stratified layer thickness
    N_max: float = 0.84 # max value of dimensionles buoyancy frequency
    Br_cmb: float = 0.62e-3/BR_SCALE # dimesionless strength of radial magnetic
                                     # field at cmb
    Ek: float = EK_PM # ekman number (require Ek <= EK_PM for magnetic diffusion
                      # to dominate

    delta_per_est: int # estimated period to compute guess of skin layer depth 
    gridflg: int = 0 # 0 for uniform, 1 for chebyshev

    solver_target_per: int = 30 # mode period in years as target for solver
    nmodes: int = 30 # number of modes to compute and return










    A, B = setup_mats(Nrad, r0, rmax, grid_params, Ek, Ntildemax, ellmax, Br, Ek_Pm, Cp, delta_tilde)
