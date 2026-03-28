from dataclasses import dataclass
from macmodes.constants import R_CMB, BR_SCALE, EK_PM

@dataclass
class ModeConfig:

    n_rad: int = 3 # number of radial grid points
    n_theta: int = 3 # number of colatitudinal grid points
    ell_max: int = 2 # max spherical harmonic degree

    H: float = 140e3/R_CMB # dimensionless stratified layer thickness
    N_max: float = 0.84 # max value of dimensionles buoyancy frequency
    Br_cmb: float = 0.62e-3/BR_SCALE # dimesionless strength of radial magnetic
                                     # field at cmb
    Ek: float = EK_PM # ekman number (require Ek <= EK_PM for magnetic diffusion
                      # to dominate

    delta_per_est: int = 50 # estimated period to compute guess of skin layer depth 
    gridflg: int = 0 # 0 for uniform, 1 for chebyshev
    iterate_delta_flg: int = 0 # number of iterations of solving eigenvalue
                               # problem and recomputing delta

    solver_target_per: int = 30 # mode period in years as target for solver
    n_modes: int = 30 # number of modes to compute and return

    def __post_init__(self):
        if self.n_rad < 3:
            raise ValueError(f"n_rad must be >= 3, got {self.n_rad}")
        if self.n_theta < 3:
            raise ValueError(f"n_theta must be >= 3, got {self.n_theta}")
        if not (self.ell_max % 2) == 0:
            raise ValueError(f"ell_max must be even, got {self.ell_max}")


configtest = ModeConfig() #tests defaults

