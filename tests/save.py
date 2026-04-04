import os
from macmodes.io import save, load
import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
test_path = os.path.join(current_dir, "test.npz")
test_arr = np.ones(5)
print("test_arr:")
print(test_arr)
save.save_npz(test_path, test_arr=test_arr)

test_arr = load.load_npz(test_path)
print("test_arr loaded:")
print(test_arr['test_arr'])
print(type(test_arr['test_arr']))

