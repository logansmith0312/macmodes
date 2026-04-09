import os
import sys
from macmodes.modes.grids import get_grid, get_grid_params
from macmodes.io.save import save_npz
from macmodes.io.load import load_npz, load_solver_config, load_plot_config
from macmodes.config import ModeParams, ModeFlags
from macmodes.constants import R_CMB, BR_SCALE, EK_PM, MACMODES_DIR

def main():
    ''' set parameters by adjusting the ModeParams and ModeFlags instance, and 
    saves the config objects, grid, gridparams list'''


    if os.path.isfile("solver_config.py"):
        params,flags = load_solver_config()
    else:
        raise RuntimeError(
                "solver_config.py not found in current run directory"
                )

    r0: float = R_CMB - params.H
    
    #create grid
    r: np.ndarray[float] = get_grid(params.n_rad, r0, flags.gridflg)
    
    #print('r: ')
    #print(r)

    print("build_params: entering get_grid_params")
    grid_params: tuple[np.ndarray[float]]  = get_grid_params(params.n_rad, r0, r)

    #print('grid params: ')    
    #print(grid_params)
    # save r, grid_params and config into npzfiles/run_params.npz 

    save_npz(os.path.join("npzfiles", "run_config"), r=r, grid_params=grid_params, params=params, flags=flags)

if __name__=="__main__":
    main()
