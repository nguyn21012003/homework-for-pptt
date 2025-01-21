from math import sqrt, pi


"""
Physical constants
"""


N = 100  # Số phương trình
hbar = 658.5  # meV.fs
chi_0 = 0.005  ### 0.01 lên 0.1 lên 0.5, 1
E_R = 4.2  # meV
delta_t = 150  # femtosecond ### Bề rộng xung laser
Delta_0 = 5.0  # meV ∆0 được gọi là năng lượng trội của photon
### delta_t = Delta_0 = 100 de ve Nt Pt
t_max = 1000  # femtosecond
t0 = -3 * delta_t
dt = 2  # femtosecond
e_max = 300  # meV
delta_e = e_max / N
T2 = 210  # femtosecond
a0 = 125  # ban kinh borh
E0 = chi_0  # meV/cm
domega = 200 / N
Egap = E0  # eV
omega0 = Egap / hbar
gamma = 6.5 * 1e-20


C0 = delta_e * sqrt(delta_e) * (1e24) / (2 * pi**2 * E_R ** (3 / 2) * a0**3)
