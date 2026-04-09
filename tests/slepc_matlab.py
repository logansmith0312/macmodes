import numpy as np

path = '/mnt/c/Users/logansmith/research/macmodes/Code/Python/MAC Code/code/matrices.npz'
with np.load(path) as data:
    A_matlab = data['A']
    B_matlab = data['B']

path = '/home/logansmith/code/macmodes/runs/mat_test/matrices_slepc.npz'
with np.load(path) as data:
    A_slepc = data['A']
    B_slepc = data['B']

print(f"A_matlab size: {A_matlab.shape}")
print(f"B_matlab size: {B_matlab.shape}")
print(f"A_slepc size: {A_slepc.shape}")
print(f"B_slepc size: {B_slepc.shape}")


mask = ~np.isclose(A_matlab, A_slepc)
indices = np.argwhere(mask)
print(indices)

for inds in indices:
    print(f"matlab: {A_matlab[inds[0], inds[1]]}")
    print(f"slepc: {A_slepc[inds[0], inds[1]]}")



