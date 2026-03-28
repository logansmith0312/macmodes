import numpy as np
from macmodes.constants import R_CMB, SECONDS_PER_YEAR, OMEGA
from macmodes.modes.grids import get_grid, get_grid_params
from macmodes.modes.matrices import setup_mats

def solve_eig(A: np.ndarray[complex], B: np.ndarray[complex], ):
    ''' solves eigenvalue problem Ax +omega*Bx = 0'''


def compute_modes(config):
    '''computes mac modes'''

    n_rad: int = config.n_rad
    H: float = config.H
    gridflg: int = config.gridflg
    iterate_delta_flg: int = config.iterate_delta_flg
    solver_target_per: int = config.solver_target_per
    n_modes: int = config.n_modes

    r0: float = R_CMB - H
    solver_target_freq: float = 2*np.pi/(solver_target_per*SECONDS_PER_YEAR*OMEGA) # solver target in 
                                                            # nondimensional frequency
    #create grid
    r: np.ndarray[float] = get_grid(n_rad, r0, gridflg)
    grid_params: list[np.ndarray[float]]  = get_grid_params(n_rad, r0, r)

    # initialize matrices
    A: np.ndarray[complex]
    B: np.ndarray[complex]
    A, B = setup_mats(r, grid_params, config)

    # solve eigenvalue problem (include delta iteration in solver routine)
    # perhaps set flag for if I want delta_iteration

    eigvals, eigvecs = solve_eig(A, B, )

     

    # sort the modes 

    # save raw data to npz in macmodes/data

