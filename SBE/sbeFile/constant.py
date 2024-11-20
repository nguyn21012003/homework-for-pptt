from math import sqrt, pi


"""
Physical constants
"""


N = 100
hbar = 658.5  # meV.fs
chi_0 = 0.0001
E_R = 4.2  # meV
delta_t = 50  # femtosecond ### Bề rộng xung laser
Delta_0 = 50  # meV
t_max = 1000  # femtosecond
t0 = -3 * delta_t
dt = 2  # femtosecond
e_max = 300  # meV
delta_e = e_max / N
T2 = 200  # femtosecond
a0 = 125  # ban kinh borh
E0 = 100
domega = 200 / N
Egap = E0  # eV
omega0 = Egap / hbar

constantE = delta_e * sqrt(delta_e)
C0 = (delta_e * sqrt(delta_e) * 8 * sqrt(8)) / (2 * pi**2 * E_R ** (3 / 2))


omega_min = -250
omega_max = -omega_min
