import os
import numpy as np

def save_npz(path, **arrays) -> None:
    '''Save numpy arrays and metadata to a compressed .npz file.'''
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.savez(path, **arrays, allow_pickle=True)

def save_mode_table(path, modes) -> None:
    '''creates and saves table of mode info at "path" ''' 
    os.makedirs(os.path.dirname(path), exist_ok=True)

    with open(path, "w") as f:
        f.write('mode' + '      ' + 'real' + '      ' + 'imag' + '      ' + 'per(yr)' + '      ' + 'Q\n')
        for ii in range(0, modes.nconv):
            f.write('{0}      {1:.7f}      {2:.7f}      {3:.7f}      {4:.3f}\n'\
                    .format(ii+1, modes.freqs[ii], modes.damps[ii], modes.pers_yrs[ii], modes.qs[ii]))


        
