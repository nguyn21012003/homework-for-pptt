import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import hbar

# Input variables
hbar = 658.5  # meV fs
N = 100
E_max = 300  # meV
dE = E_max / N

d0 = 30  # meV
delta_t = 25  # fs
dt = 2  # fs
tmax = 500  # fs
t0 = -3 * delta_t
T2 = 200  # fs
E_R = 4.2
chi_0 = 0.5


def RK4(F, t, fp, dt, N, omega_R, E_n, dE, d0, T2):
    k1 = dt * F(N, omega_R, E_n, fp, dE, d0, T2)
    k2 = dt * F(N, omega_R, E_n, fp + k1 * 0.5, dE, d0, T2)
    k3 = dt * F(N, omega_R, E_n, fp + k2 * 0.5, dE, d0, T2)
    k4 = dt * F(N, omega_R, E_n, fp + k3, dE, d0, T2)
    fp = fp + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return fp


# Def fj and fn equation
def fp_init(N, f_e, f_h, p_n):
    fp = np.zeros((2, N), dtype="complex")
    for n in range(N):
        fp[0][n] = f_e + 1j * f_h
        fp[1][n] = p_n
    return fp


# Define g function
def g(n, n1, dE):
    return (1 / np.sqrt(n * dE)) * np.log(abs((np.sqrt(n) + np.sqrt(n1)) / (np.sqrt(n) - np.sqrt(n1))))


# Define energy function
def E_n(N, E_R, dE, g, fp):
    E_n_values = np.zeros(N)
    for n in range(N):
        for n1 in range(N):
            E_n_values[n] += (np.sqrt(E_R) / np.pi) * dE * g(n + 1, n1 + 1, dE) * (np.real(fp[0, n1]) + np.imag(fp[0, n1]))
    return E_n_values


# Define omega_R function
def omega_R(N, delta_t, chi_0, t, E_R, dE, g, fp):
    omega_R_values = np.zeros(N, dtype="complex")
    for n in range(N):
        for n1 in range(N):
            omega_R_values[n] += (1 / hbar) * (
                0.5 * (hbar * np.sqrt(np.pi) / delta_t) * chi_0 * np.exp(-(t**2) / delta_t**2)
                + (np.sqrt(E_R) / np.pi) * dE * g(n + 1, n1 + 1, dE) * fp[1, n1]
            )
    return omega_R_values


# Define sbe function
def sbe(N, omega_R, E_n, fp, dE, d0, T2):
    F = np.zeros((2, N), dtype="complex")
    for n in range(N):
        F[0][n] = -2 * np.imag(omega_R[n] * np.conjugate(fp[1, n]))
        F[1][n] = -(1j / hbar) * (n * dE - d0 - E_n[n]) * fp[1, n] + 1j * (1 - np.real(fp[0, n]) - np.imag(fp[0, n])) * omega_R[n] - (fp[1, n] / T2)
    return F


def main():
    # Initial conditions
    f0_e, f0_h, p0 = 0, 0, 0
    fp = fp_init(N, f0_e, f0_h, p0)

    results = []
    nt = int(tmax - t0 / dt)

    for i in range(nt):
        t = t0 + i * dt

        # Calculate E_n and omega_R at each time step
        E_n_values = E_n(N, E_R, dE, g, fp)
        omega_R_values = omega_R(N, delta_t, chi_0, t, E_R, dE, g, fp)

        # Update fp using RK4
        fp = RK4(sbe, t, fp, dt, N, omega_R_values, E_n_values, dE, d0, T2)
        results.append(fp)

    with open("sbe.txt", "w") as sbe_w:
        s0 = "{0:^10} {delim} {1:^10} {delim} {2:^10} {delim} {3:^10} {delim} {4:^10} \n"
        sbe_w.write(s0.format("E", "t", "f_e", "f_h", "p_n", delim="|"))
        for n in range(nt):
            s1 = "{0:^10} {delim} {1:^10} {delim} {2:^10} {delim} {3:^10} {delim} {4:^10} \n"
            sbe_w.write(s1.format(n * dE, t0 + n * dt, np.real(fp[0][n]), np.imag(fp[0][n]), fp[1][n], delim="|"))


main()
