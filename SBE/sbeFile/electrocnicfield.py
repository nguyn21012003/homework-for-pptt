import numpy as np
from constant import t_max, t0, delta_t, E0, dt

# Time span
tSpan = np.arange(t0, t_max, dt)

# Calculate the Gaussian pulse
E_t = E0 * np.exp(-(tSpan**2) / (delta_t**2))

# Save the data to a text file
with open("gaussian_pulse_data.txt", "w") as f:
    for t, E in zip(tSpan, E_t):
        f.write(f"{t} {E}\n")  # Save time in fs and E(t)
