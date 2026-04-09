import os
import petsc4py
import sys
import numpy as np
from macmodes.modes.solver import compute_modes
from macmodes.io.save import save_npz, save_mode_table
from macmodes.io.load import load_npz, load_object
from macmodes.config import ModeParams, ModeFlags
from macmodes.constants import R_CMB, BR_SCALE, EK_PM, SECONDS_PER_YEAR, OMEGA

petsc4py.init(sys.argv)
from petsc4py import PETSc
from mpi4py import MPI

def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # read in r, grid_params, config objects
    current_dir = os.getcwd()
    npz_dir = os.path.join(current_dir, "npzfiles")
    run_config_path = os.path.join(npz_dir, "run_config.npz")

    # load only on rank 0
    if rank == 0:
        if os.path.isfile(run_config_path):
            with load_npz(run_config_path) as run_params:
                r = load_object(run_params, 'r')
                grid_params = load_object(run_params, 'grid_params') 
                params = load_object(run_params, 'params')
                flags = load_object(run_params, 'flags')
        else:
            raise RuntimeError(
                    f"run_config not found in {npz_dir}," 
                    "ensure you are in runs/your_run and ran init_run.py first"
                    )
    else:
        r = grid_params = params = flags = None

    # broad cast to all ranks
    r = comm.bcast(r, root=0)
    grid_params = comm.bcast(grid_params, root=0)
    params = comm.bcast(params, root=0)
    flags = comm.bcast(flags, root=0)

    if rank == 0:
        print("")
        print("===================================================================")
        print("COMPUTE_MODES (MAIN): Routine to compute MAC Modes of oscillation")  
        print("in Boussinesq, perfectly conducting, rotating fluid confined to")
        print("spherical shell.")
        print("")
        print("parameters:")
        for key, val in vars(params).items():
            if isinstance(val, int):
                print(f"\t\t{key}: {val}")
            else:
                print(f"\t\t{key}: {val:.5e}")
        print(f"\t\tst_target: {2*np.pi/(params.st_target_per_yrs*SECONDS_PER_YEAR*OMEGA):.5e}")
        print(f"\t\tmatrix size: {5*params.ell_max*params.n_rad//2}")
        print("flags:")
        for key, val in vars(flags).items():
            print(f"\t\t{key}: {val}")
        print("===================================================================")
        print("")

    modes = compute_modes(r, grid_params, params, flags)

    if rank == 0:
        if modes.n_modes > 0:
            # save modes instance to current directory
            save_npz(os.path.join(npz_dir, "modes"), modes=modes)

            # save readable table of mode data
            table_path = os.path.join(current_dir, "mode_info.dat")
            save_mode_table(table_path, modes)
            
            # print to terminal
            with open(table_path, 'r') as file:
                for line in file:
                    print(line, end='')
        else:
            print("compute_modes: no modes to save")

if __name__ == "__main__":
    main()
