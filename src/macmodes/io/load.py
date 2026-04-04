import os
import numpy as np

def load_npz(path):
    '''Load a .npz file and return a dict-like object.'''
    return np.load(path, allow_pickle=True)

def load_object(data, key):
    '''np.load stores classes as a 0-dim np.ndarray, so this checks if the data[key]
    is a 0-dim array, and if so unpacks with the item() method, else just returns the
    array itself'''

    obj = data[key]

    if isinstance(obj, np.ndarray) and obj.shape == ():
        return obj.item()   # stored python object
    else:
        return obj          # real numpy array
