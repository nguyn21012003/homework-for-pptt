import numpy as np
import matplotlib.pyplot as plt

# from scipy.constants import hbar
from numpy import real as RE, imag as IM, conjugate, sqrt, log, pi, exp
from numpy import typing as npt
from time import time
import cProfile

# Input variables
N = 100
hbar = 658.5  # MeV fs
chi_0 = 0.5
E_R = 4.2
delta_t = 200  # fs
delta_0 = 100  # meV
dt = 2  # fs
t_max = 500  # fs
t0 = -3 * delta_t
e_max = 300  # meV
delta_e = e_max / N
T2 = 200  # fs
a0 = 125


def rk4(dF: npt.NDArray, tn: npt.NDArray, yn: npt.NDArray, h: float) -> npt.NDArray:
    # tn í value in tSpan  array
    k1 = dF(tn, yn)
    k2 = dF(tn + 0.5 * h, yn + 0.5 * h * k1)
    k3 = dF(tn + 0.5 * h, yn + 0.5 * h * k2)
    k4 = dF(tn + h, yn + h * k3)
    return yn + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def g(n, n1):
    return (1 / sqrt(n * delta_e)) * log(abs((sqrt(n) + sqrt(n1)) / (sqrt(n) - sqrt(n1))))


def E_n(n, g, f_e, f_h):
    E = 0
    for n1 in range(1, N + 1):
        if n1 == n:
            continue
        E += (sqrt(E_R) / pi) * delta_e * g(n, n1) * (f_e[n1] + f_h[n1])
    return E


def omega_R(n, t, g, p_n):

    sum1 = 0
    for n1 in range(1, N + 1):
        if n1 == n:
            continue
        sum1 += g(n, n1) * p_n[n1]
    OMEGA = (0.5 * (hbar * sqrt(pi) / delta_t) * chi_0 * exp(-(t**2) / delta_t**2) + (sqrt(E_R) / pi) * delta_e * sum1) / hbar

    return OMEGA


def dF(t, Y):

    F = np.zeros((2, N + 1), dtype="complex")

    print(RE(Y[0]))
    for n in range(1, N + 1):
        E = E_n(n, g, RE(Y[0]), IM(Y[0]))
        OMEGA = omega_R(n, t, g, Y[1])

        F[0][n] = -2 * IM(OMEGA * conjugate(Y[1][n])) + 1j * -2 * IM(OMEGA * conjugate(Y[1][n]))
        F[1][n] = -(1j / hbar) * (n * delta_e - delta_0 - E) * Y[1][n] + 1j * (1 - RE(Y[0])[n] - IM(Y[0])[n]) * OMEGA - (Y[1][n] / T2)

    return F


def solve_sys_ode(dF: npt.NDArray, h: float, rk4: npt.NDArray, N: int) -> npt.NDArray:

    Y = np.zeros((2, N + 1), dtype="complex")

    tSpan = np.linspace(t0, t_max, N)

    Polarization = np.zeros(N)
    Energy = np.zeros(N)
    fe = []
    p = []
    for ti in tSpan:
        Y = rk4(dF, ti, Y, h)
        fe_t = []
        p_t = []
        for n in range(0, N):
            fe_t.append(RE(Y[0][n]))
            p_t.append((abs(Y[1][n])))
            Polarization[n] = delta_e * sqrt(delta_e) * sqrt(n) * RE(Y[0][n])
            Energy[n] = n * delta_e
        fe.append(fe_t)
        p.append(p_t)

    return tSpan, Energy, fe, p


def main():
    start = time()

    t, Energy, fe, p = solve_sys_ode(dF, dt, rk4, N)
    E, T = np.meshgrid(Energy, t)
    fe = np.array(fe)
    p = np.array(p)

    end = time()
    print(end - start)

    fig = plt.figure(figsize=(16, 8))

    ax1 = fig.add_subplot(121, projection="3d")
    ax1.plot_wireframe(T, E, fe, cmap="viridis")
    ax1.set_ylabel(r"Năng lượng MeV")
    ax1.set_title(r"Hàm phân bố $f_e$ theo thời gian và năng lượng")

    ax2 = fig.add_subplot(122, projection="3d")
    ax2.plot_wireframe(T, E, p, cmap="plasma")
    ax2.set_ylabel(r"Năng lượng MeV")
    ax2.set_title(r"Hàm phân bố $f_p$ theo thời gian và năng lượng")

    # ax3 = fig.add_subplot(133, projection="3d")
    # ax3.plot_surface(T, E, P, cmap="copper")

    plt.savefig("delta20fs")

    plt.show()


if __name__ == "__main__":
    main()
