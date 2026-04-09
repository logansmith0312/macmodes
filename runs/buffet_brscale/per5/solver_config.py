# Template config file

from macmodes.config import ModeParams, ModeFlags
from macmodes.constants import R_CMB_DIM, BR_SCALE, EK_PM

params = ModeParams(
        n_rad = 150, # number of radial grid points
        n_theta = 100, # number of colatitudinal grid points
        ell_max = 30, # max spherical harmonic degree

        H = 140e3/R_CMB_DIM, # dimensionless stratified layer thickness
        N_max = 0.84, # max value of dimensionles buoyancy frequency
        Br_cmb = 0.62e-3/BR_SCALE, # dimesionless strength of radial magnetic
        # field at cmb
        Ek = EK_PM, # ekman number (require Ek <= EK_PM for magnetic diffusion
        # to dominate

        delta_per_est = 100, # estimated period to compute guess of skin layer depth

        n_modes = 5, # number of modes to compute and return
        st_target_per_yrs = 5, # frequency that shift invert will target
        )

flags = ModeFlags(
        gridflg = 1, # 0 for uniform, 1 for chebyshev
        inner_bdry_flg = 1, # 0 for pseudo-vacuum, 1 to match IC potential
        outer_bdry_flg = 0, # 0 for pseudo-vacuum, 1 for mantle EM coupling
        )
