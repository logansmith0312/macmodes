import numpy as np
from macmodes.constants import R_CMB

def compute_modes(config):
    '''computes mac modes'''

    n_rad: int = config.n_rad
    H: float = config.H
    gridflg: int = config.gridflg
    solver_target_per: int = config.solver_target_per
    n_modes: int = config.n_modes

    r0: float = R_CMB - H
    solver_target_freq = 2*np.pi/(solver_target_per*yr*rot) # solver target in 
                                                            # nondimensional frequency
    #create grid
    r: np.ndarray = get_grid(n_rad, r0, gridflg)
    grid_params: list[np.ndarray]  = get_grid_params(n_rad, r0, r)

    # initialize matrices
    A: np.ndarray
    B: np.ndarray
    A, B = setup_mats(r, grid_params, config)

    # solve eigenvalue problem (include delta iteration in solver routine)
    # perhaps set flag for if I want delta_iteration
     

    # sort the modes 

    # save raw data to npz in macmodes/data

