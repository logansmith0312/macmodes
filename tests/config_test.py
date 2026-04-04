import os
from macmodes.io.save import save_npz
from macmodes.io.load import load_npz, load_object
from macmodes.config import ModeConfig



config = ModeConfig()
current_dir = os.path.dirname(os.path.abspath(__file__))
config_path = os.path.join(current_dir, "config.npz")
save_npz(config_path, config=config)

with load_npz(config_path) as data:
    config = load_object(data, "config")

print(config)
print("Testing config object:")
for key, val in vars(config).items():
    print(f"{key}: {val}")
