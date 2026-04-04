import os
import petsc4py
import sys
from macmodes.modes.solver import compute_modes
from macmodes.io.save import save_npz, save_mode_table
from macmodes.io.load import load_npz, load_object
from macmodes.config import ModeConfig
from macmodes.constants import R_CMB, BR_SCALE, EK_PM

petsc4py.init(sys.argv)
from petsc4py import PETSc

def main():

    # read in r, grid_params, config
    current_dir = os.getcwd()
    npz_dir = os.path.join(current_dir, "npzfiles")
    run_params_path = os.path.join(npz_dir, "run_params.npz")
    if os.path.isfile(run_params_path):
        with load_npz(run_params_path) as run_params:
            r = load_object(run_params, 'r')
            grid_params = load_object(run_params, 'grid_params') 
            config = load_object(run_params, 'config')
    else:
        raise RuntimeError(f"run_params not found in {current_dir}/npzfiles, ensure you are in runs/your_run when you call compute_modes.py, and execute init_run.py first")
    print(type(r))
    print(type(grid_params))
    print(type(config))
    print("")
    print("===================================================================")
    print("COMPUTE_MODES (MAIN): Routine to compute MAC Modes of oscillation")  
    print("in Boussinesq, perfectly conducting, rotating fluid confined to")
    print("spherical shell.")
    print("")
    print("parameters/flags:")
    for key, val in vars(config).items():
        print(f"\t\t{key}: {val}")
    print(f"\t\tmatrix size: {5*config.ell_max*config.n_rad//2}")
    print("===================================================================")
    print("")

    modes = compute_modes(r, grid_params, config)

    # save modes instance to current directory
    save_npz(os.path.join(npz_dir, "modes"), modes=modes)

    # save readable table of mode data
    table_path = os.path.join(current_dir, "mode_info.dat")
    save_mode_table(table_path, modes)
    
    # print to terminal
    with open(table_path, 'r') as file:
        for line in file:
            print(line, end='')

if __name__ == "__main__":
    main()
