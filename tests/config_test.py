from macmodes.config import ModeConfig


config = ModeConfig()

print("Testing config object:")
for key, val in vars(config).items():
    print(f"{key}: {val}")
