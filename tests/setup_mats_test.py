from macmodes.config import ModeConfig
from macmodes.modes.solver import compute_modes

config = ModeConfig()

A, B = compute_modes(config)

print("size A:")
print(A.shape)
print("size B:")
print(B.shape)
