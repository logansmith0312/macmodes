#===============================================================================
mode_start = 1
mode_end = 2
#===============================================================================

import os
import numpy as np
from macmodes.analysis.plotting import contour_plot_mac
from macmodes.io.load import load_npz, load_object
from macmodes.io.save import save_plot

def main():

    # load in grid, run_config and modes
    current_dir = os.getcwd()
    npz_dir = os.path.join(current_dir, "npzfiles")
    run_config_path = os.path.join(npz_dir, "run_config.npz")
    modes_path = os.path.join(npz_dir, "modes.npz")

    if os.path.isfile(run_config_path):
        with load_npz(run_config_path) as run_params:
            r = load_object(run_params, 'r')
            grid_params = load_object(run_params, 'grid_params')
            params = load_object(run_params, 'params')
            flags = load_object(run_params, 'flags')
    else:
        raise RuntimeError(
                f"run_config not found in {npz_dir},"
                "ensure you are in runs/your_run and ran build_config first"
                )

    if os.path.isfile(modes_path):
        with load_npz(modes_path) as mode_data:
            modes = load_object(mode_data, 'modes')
    else:
        raise RuntimeError(
                f"modes.npz not found in {npz_dir}," 
                "ensure you are in runs/your_run and ran compute_modes first"
                )
    mode_inds = np.arange(mode_start, mode_end + 1)

    for ind in mode_inds:
        figs = contour_plot_mac(ind, r, grid_params, params.n_rad, params.n_theta, params.ell_max, modes.evecs)
        
        names = [f'v_re', f'v_im', f'b_re', f'b_im']
        plt_dir = os.path.join(current_dir, "plots", f"mode_{ind}")

        print(f"save_plots: saving mode={ind}, per={modes.pers_yrs[ind-1]:.5f} ")
        for k in range(len(figs)):
            save_plot(os.path.join(plt_dir, names[k]),figs[k])


if __name__ == "__main__":
    main()
