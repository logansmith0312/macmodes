import os
import sys
import matplotlib.pyplot as plt
from macmodes.io.load import load_npz, load_object
from macmodes.analysis.plotting import contour_plot_mac
def main():

    if len(sys.argv) > 1:
        mode_ind=int(sys.argv[1])
    else:
        raise RuntimeError(
                "No mode specified, call should be e.g. mode_look 4"
                )

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

    figs = contour_plot_mac(mode_ind, r, grid_params, params.n_rad, params.n_theta, params.ell_max, modes.evecs)
    plt.show()


if __name__ == "__main__":
    main()
