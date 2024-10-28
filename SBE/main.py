import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import hbar
from numpy import real as RE, imag as IM, conjugate, sqrt, log, pi, exp


def RK4(F, t, fp, dt):
    k1 = dt * F(t, fp)
    k2 = dt * F(t + dt / 2, fp + k1 * dt / 2)
    k3 = dt * F(t + dt / 2, fp + k2 * dt / 2)
    k4 = dt * F(t + dt, fp + dt * k3)
    fp = fp + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return fp


# Def fj and fn equations
def fp(N, icInput):
    f_e, f_h, p_n = icInput
    fArr = np.zeros([2, N], dtype="complex")
    for n in range(1, N + 1):
        fArr[0][n] = f_e + 1j * f_h
        fArr[1][n] = p_n * n
    return fArr


# Time derivative of fj and fn (SBE equations)
def F(N, omega_R, E_n, fp, dE, d0, T2):
    F = np.zeros((2, N))
    for n in range(1, N + 1):
        F[0][N] = -2 * IM(omega_R * conjugate(fp[1][n]))
        F[1][N] = -(1j / hbar) * (n * dE - d0 - E_n) * fp[1][n] + 1j * (1 - RE(fp[0][n]) - IM(fp[0][n])) * omega_R - (fp[1][n] / T2)
    return F


# Some input equations
def g(n, n1, dE):
    return (1 / sqrt(n * dE)) * log(abs((sqrt(n) + sqrt(n1)) / (sqrt(n) - sqrt(n1))))


def E_n(N, E_R, dE, g, fp):
    for n in range(1, N + 1):
        for n1 in range(1, N + 1):
            E_n = (sqrt(E_R) / pi) * dE * g(n, n1, dE) * (RE(fp[0][n1]) + IM(fp[0][n1]))
    return E_n


def omega_R(N, delta_t, chi_0, t, E_R, dE):
    for n in range(1, N + 1):
        for n1 in range(1, N + 1):
            omega_R = (1j / hbar) * (
                1 / 2 * (hbar * sqrt(pi) / delta_t) * chi_0 * exp(-(t**2) / delta_t**2) + (sqrt(E_R) / pi) * dE * g(n, n1, dE) * fp[1][n1]
            )
    return omega_R


def input():
    # Input variables
    N = 100
    E_max = 300  # meV
    dE = E_max / N
    delta_t = 25  # fs
    d0 = 30  # meV
    dt = 2  # fs
    t_max = 300  # fs
    t0 = -3 * delta_t

    # Initial conditions
    icInput = np.zeros(3)
    print(icInput[0])
    print(icInput[1])
    print(icInput[2])


def main():
    input()


if __name__ == "__main__":
    main()
