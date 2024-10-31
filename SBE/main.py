import numpy as np
import matplotlib.pyplot as plt

# from scipy.constants import hbar
from numpy import real as RE, imag as IM, conjugate, sqrt, log, pi, exp
from numpy import typing as npt


# Input variables
a0 = 125
N = 100
E_max = 300  # meV
dE = E_max / N
E_R = 4.2
delta_t = 25  # fs
del0 = 30  # meV
dt = 2  # fs
t_max = 500  # fs
t0 = -3 * delta_t
hbar = 658.5  # MeV fs
chi_0 = 0.5
T2 = 200  # fs

# evẻi won í khong
# Y = np.zeros([2, N + 1], dtype="complex")
# for n in range(1, N + 1):
#    Y[0][n] = 0 + 1j * 0  # f_e +i f_h
#    Y[1][n] = 0 + 1j * 0  # p_n


def rk4(dF: npt.NDArray, tn: npt.NDArray, yn: npt.NDArray, h: float) -> npt.NDArray:
    # tn í value in tSpan  array
    k1 = dF(tn, yn)
    k2 = dF(tn + 0.5 * h, yn + 0.5 * h * k1)
    k3 = dF(tn + 0.5 * h, yn + 0.5 * h * k2)
    k4 = dF(tn + h, yn + h * k3)
    return yn + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def g(n, n1):
    return (1 / sqrt(n * dE)) * log(abs((sqrt(n) + sqrt(n1)) / (sqrt(n) - sqrt(n1))))


def E_n(g, f_e, f_h):
    E = np.zeros(N + 1)
    for n in range(1, N + 1):
        for n1 in range(1, N + 1):
            if n1 != n:
                E[n] += (sqrt(E_R) / pi) * dE * g(n, n1) * (2 * f_e[n1])
    # f_e[n1] + f_h[n1]
    return E


def omega_R(t, g, p_n):

    OMEGA = np.zeros(N + 1, dtype="complex")
    for n in range(1, N + 1):
        for n1 in range(1, N + 1):
            if n1 != n:
                OMEGA[n] += (1 / hbar) * (
                    0.5 * (hbar * sqrt(pi) / delta_t) * chi_0 * exp(-(t**2) / delta_t**2) + (sqrt(E_R) / pi) * dE * g(n, n1) * p_n[n1]
                )

    return OMEGA


def dF(t, Y):

    F = np.zeros([2, N + 1], dtype="complex")

    f_e = RE(Y[0])
    f_h = IM(Y[0])
    p_n = Y[1]
    E = E_n(g, f_e, f_h)
    OMEGA = omega_R(t, g, p_n)
    for n in range(1, N + 1):
        F[0][n] = -2 * IM(OMEGA[n] * conjugate(p_n[n]))  # phuc = thuc -> f_h =0
        F[1][n] = -(1j / hbar) * (n * dE - del0 - E[n]) * p_n[n] + 1j * (1 - 2 * f_e[n]) * OMEGA[n] - (p_n[n] / T2)

    return F


def solve_sys_ode(dF: npt.NDArray, h: float, rk4: npt.NDArray, N: int) -> npt.NDArray:
    y = []  # để lưu giá trị của Y
    Y = np.zeros([2, N + 1], dtype=complex)
    tSpan = np.linspace(t0, t_max, N)
    print(tSpan)
    for tn in tSpan:
        for j in range(N):
            Y = rk4(dF, tn, Y, h)  # lặp trong yn
            y.append(Y)
    return tSpan, y


def main():

    t, solution = solve_sys_ode(dF, dt, rk4, N)
    fe = []
    pn = []
    for i in range(1, N + 1):
        fe.append(RE(solution[i][0][N]))
        pn.append(abs(solution[i][1]) ** 2)
    plt.plot(t, solution)
    plt.show()


if __name__ == "__main__":
    main()
