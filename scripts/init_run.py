import os
from macmodes.modes.grids import get_grid, get_grid_params
from macmodes.io.save import save_npz
from macmodes.config import ModeConfig
from macmodes.constants import R_CMB, BR_SCALE, EK_PM

def main():
    ''' set parameters by adjusting the ModeConfig instance, and saves the config object
    and gridparams list to specified directory'''

    run_dir_name = 'test_dir'

    config = ModeConfig(
        n_rad = 100, # number of radial grid points
        n_theta = 100, # number of colatitudinal grid points
        ell_max = 2, # max spherical harmonic degree

        H = 140e3/R_CMB, # dimensionless stratified layer thickness
        N_max = 0.84, # max value of dimensionles buoyancy frequency
        Br_cmb = 0.62e-3/BR_SCALE, # dimesionless strength of radial magnetic
        # field at cmb
        Ek = EK_PM, # ekman number (require Ek <= EK_PM for magnetic diffusion
        # to dominate

        delta_per_est = 50, # estimated period to compute guess of skin layer depth
        gridflg = 0, # 0 for uniform, 1 for chebyshev
        iterate_delta_flg = 0, # number of iterations of solving eigenvalue
                               # problem and recomputing delta

        n_modes = 10, # number of modes to compute and return
        )

    n_rad: int = config.n_rad
    H: float = config.H
    gridflg: int = config.gridflg
    iterate_delta_flg: int = config.iterate_delta_flg
    n_modes: int = config.n_modes

    r0: float = R_CMB - H

    #create grid
    r: np.ndarray[float] = get_grid(n_rad, r0, gridflg)

    print("compute_modes: entering get_grid_params")
    grid_params: tuple[np.ndarray[float]]  = get_grid_params(n_rad, r0, r)

    # save r, grid_params and config into specified run_dir_name

    script_dir = os.path.dirname(os.path.abspath(__file__))
    root_dir = os.path.dirname(script_dir)

    npz_dir = os.path.join(root_dir, "runs", run_dir_name, "npzfiles")
    save_npz(os.path.join(npz_dir, "run_params"), r=r, grid_params=grid_params, config=config)

if __name__=="__main__":
    main()
