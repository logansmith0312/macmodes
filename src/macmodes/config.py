from dataclasses import dataclass
from macmodes.constants import R_CMB_DIM, BR_SCALE, EK_PM

@dataclass
class ModeParams:

    n_rad: int = 3 # number of radial grid points
    n_theta: int = 3 # number of colatitudinal grid points
    ell_max: int = 2 # max spherical harmonic degree

    H: float = 140e3/R_CMB_DIM # dimensionless stratified layer thickness
    N_max: float = 0.84 # max value of dimensionles buoyancy frequency
    Br_cmb: float = 0.62e-3/BR_SCALE # dimesionless strength of radial magnetic
                                     # field at cmb
    Ek: float = EK_PM # ekman number (require Ek <= EK_PM for magnetic diffusion
                      # to dominate

    delta_per_est: int = 50 # estimated period to compute guess of skin layer depth 

    n_modes: int = 30 # number of modes to compute and return
    st_target_per_yrs: float = 100 # frequency that shift invert will target

    def __post_init__(self):
        if self.n_rad < 3:
            raise ValueError(f"n_rad must be >= 3, got {self.n_rad}")
        if self.n_theta < 3:
            raise ValueError(f"n_theta must be >= 3, got {self.n_theta}")
        if not (self.ell_max % 2) == 0:
            raise ValueError(f"ell_max must be even, got {self.ell_max}")

@dataclass
class ModeFlags:
    gridflg: int = 0 # 0 for uniform, 1 for chebyshev
    inner_bdry_flg: int = 0 # 0 for pseudo vacuum, 1 to match IC potential    
    outer_bdry_flg: int = 0 # 0 for pseudo vacuum, 1 for EM coupling with mantle

    def __post_init__(self):
        if (self.gridflg < 0) or (self.gridflg > 1):
            raise ValueError(f"gridflg must be 0 or 1, got {gridflg}")
        if (self.inner_bdry_flg < 0) or (self.inner_bdry_flg > 1):
            raise ValueError(f"inner_bdry_flg must be 0 or 1, got {inner_bdry_flg}")
        if (self.outer_bdry_flg < 0) or (self.outer_bdry_flg > 1):
            raise ValueError(f"outer_bdry_flg must be 0 or 1, got {outer_bdry_flg}")

@dataclass
class PlotConfig:
    mode_start: int
    mode_end: int
    save_flg: int

paramstest = ModeParams() #tests defaults
flagstest = ModeFlags()

