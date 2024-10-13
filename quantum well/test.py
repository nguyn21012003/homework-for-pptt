import numpy as np


def generate_intervals(z0=45, step=np.pi / 2):
    intervals = []
    x_vals = np.arange(0, z0, step)  # Now x ranges up to z0

    for x in x_vals:
        # Avoid points where tan(x) diverges (odd multiples of pi/2)
        if np.isclose(np.mod(x, np.pi / 2), 0, atol=1e-5) and not np.isclose(np.mod(x, np.pi), 0, atol=1e-5):
            continue

        # Ensure the term inside the sqrt is valid
        if x == 0 or (z0 / x) < 1:  # Avoid division by zero and invalid sqrt
            continue

        # Evaluate the function tan(x) - sqrt((z0 / x)^2 - 1)
        sqrt_term = np.sqrt((z0 / x) ** 2 - 1)
        tan_x = np.tan(x)
        result = tan_x - sqrt_term

        intervals.append((x, result))

    return intervals


z0 = 45  # Updated z0 to 45

intervals = generate_intervals(z0)
print(intervals)
