import numpy as np
from math import sqrt

L = 1
t = 3
rho = 0.01
Tension = 40
c = sqrt(Tension / rho)

dt = 0.0001
dx = 0.01
beta = (c * dt) / dx
print(beta)
X = np.arange(0, L + dx, dx)
T = np.arange(0, t + dt, dt)
U = []
for _ in range(t + 1):
    row = []
    for _ in range(t + 1):
        row.append(0.0)
    U.append(row)

U = np.zeros((len(X) - 1, len(T) - 1))
print(U)
