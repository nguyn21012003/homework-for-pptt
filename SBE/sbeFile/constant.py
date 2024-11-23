from math import sqrt, pi


"""
Physical constants
"""


N = 100  # Số phương trình
hbar = 658.5  # meV.fs
chi_0 = 0.5
E_R = 4.2  # meV
delta_t = 50  # femtosecond ### Bề rộng xung laser
Delta_0 = 50  # meV
t_max = 1200  # femtosecond
t0 = -3 * delta_t
dt = 2  # femtosecond
e_max = 300  # meV
delta_e = e_max / N
T2 = 50  # femtosecond
a0 = 125  # ban kinh borh
E0 = 100
domega = 200 / N
Egap = E0  # eV
omega0 = Egap / hbar

constantE = delta_e * sqrt(delta_e)
C0 = (delta_e * sqrt(delta_e)) / (2 * (pi**2) * E_R ** (3 / 2) * a0**3)
