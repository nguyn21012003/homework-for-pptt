import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, log, pi, exp

from numpy import real as RE, imag as IM, conjugate
from numpy import typing as npt
from time import time
from matplotlib.backends.backend_pdf import PdfPages
from tqdm import tqdm

# Input variables
N = 100
hbar = 658.5  # meV.fs
chi_0 = 0.5
E_R = 4.2  # meV
delta_t = 20  # femtosecond
delta_0 = -10  # meV
t_max = 500  # femtosecond
t0 = -3 * delta_t
dt = 2  # femtosecond
e_max = 300  # meV
delta_e = e_max / N
T2 = 200  # femtosecond
a0 = 125  # ban kinh borh

# Egap = 1.77e3  # eV
# omega0 = Egap / hbar

constant = delta_e * sqrt(delta_e)


def rk4(dF, tn, yn, h):
    # tn is value in tSpan  array
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

    summ = 0
    for n1 in range(1, N + 1):
        if n1 == n:
            continue
        summ += g(n, n1) * p_n[n1]
    OMEGA = (0.5 * (hbar * sqrt(pi) / delta_t) * chi_0 * exp(-(t**2) / delta_t**2) + (sqrt(E_R) / pi) * delta_e * summ) / hbar

    return OMEGA


def dF(t, Y):

    F = np.zeros((2, N + 1), dtype=complex)

    # print(RE(Y[0]))
    for n in range(1, N + 1):
        E = E_n(n, g, RE(Y[0]), IM(Y[0]))
        OMEGA = omega_R(n, t, g, Y[1])

        F[0][n] = -2 * IM(OMEGA * conjugate(Y[1][n]))
        F[1][n] = -(1j / hbar) * (n * delta_e - delta_0 - E) * Y[1][n] + 1j * (1 - RE(Y[0])[n] - IM(Y[0])[n]) * OMEGA - (Y[1][n] / T2) * 0

    return F


def solve_sys_ode(dF: npt.NDArray, h: float, rk4: npt.NDArray, N: int) -> npt.NDArray:

    Y = np.zeros((2, N + 1), dtype=complex)

    tSpan = np.linspace(t0, t_max, N)
    # tSpan = np.arange(t0, t_max, 100)
    print(tSpan)
    Polarization = np.zeros(N)
    NumberDensity = np.zeros(N)
    EnergyEps = np.zeros(N)
    fe = []
    p = []
    for ti in tqdm(range(len(tSpan)), desc="Processing", unit="step "):

        Y = rk4(dF, tSpan[ti], Y, h)
        fe_t = []
        p_t = []
        NumberDensity_sum = 0
        Polarization_sum = 0
        for n in range(N):
            fe_t.append(RE(Y[0][n]))
            p_t.append((abs(Y[1][n])))

            EnergyEps[n] = (n + 1) * delta_e

            NumberDensity_sum += sqrt(n + 1) * RE(Y[0][n])
            Polarization_sum += abs(constant * sqrt(n + 1) * RE(Y[1][n]))

        NumberDensity[ti] = NumberDensity_sum
        Polarization[ti] = Polarization_sum

        fe.append(fe_t)
        p.append(p_t)

    return tSpan, EnergyEps, fe, p, Polarization, NumberDensity


def Et(arg):
    pass


def FourierTrans(arg):
    pass


def multipage(filename, figs=None, dpi=200):

    pp = PdfPages(filename)
    if figs is None:
        figs = [plt.figure(n) for n in plt.get_fignums()]
    for fig in figs:
        fig.savefig(pp, format="pdf")
    pp.close()


def main():

    start = time()

    t, EnergyEps, fe, p, Polarization, NumberDensity = solve_sys_ode(dF, dt, rk4, N)

    end = time()

    print(end - start)
    E, T = np.meshgrid(EnergyEps, t)

    fe = np.array(fe)
    p = np.array(p)

    color = "inferno"
    ombre = "copper"
    transparent = 0.8

    fig1 = plt.figure(figsize=(16, 10))
    fig2 = plt.figure(figsize=(16, 10))
    # fig3 = plt.figure(figsize=(16, 10))

    ax1 = fig1.add_subplot(121, projection="3d")
    ax1.plot_wireframe(T, E, fe)
    ax1.set_ylabel(r"Năng lượng MeV")
    ax1.set_xlabel(r"Thời gian t(fs)")
    ax1.set_title(r"Hàm phân bố $f_e$ theo thời gian và năng lượng")

    ax2 = fig1.add_subplot(122, projection="3d")
    ax2.plot_wireframe(T, E, p)
    ax2.set_ylabel(r"Năng lượng MeV")
    ax2.set_xlabel(r"Thời gian t(fs)")
    ax2.set_title(r"Hàm phân bố $f_p$ theo thời gian và năng lượ9ng")

    ax3 = fig2.add_subplot(121)
    ax3.plot(t, Polarization, color="purple")
    ax3.set_ylabel(r"Phân cực toàn phần")
    ax3.grid()
    ax3.set_xlabel(r"Thời gian t(fs)")
    ax3.set_title(r"Hàm phân cực $|P_n(t)|$ theo thời gian và năng lượng")

    ax4 = fig2.add_subplot(122)
    ax4.plot(t, NumberDensity, color="purple")
    ax4.grid()
    ax4.set_ylabel(r"Mật độ toàn phần")
    ax4.set_xlabel(r"Thời gian t(fs)")
    ax4.set_title(r"Hàm phân bố toàn phần $N(t)$ theo thời gian và năng lượng")

    # ax5 = fig3.add_subplot(122)
    # ax5.plot(t, Energy)

    # ax3 = fig.add_subplot(133, projection="3d")
    # ax3.plot_surface(T, E, P, cmap="copper")

    multipage(f"{delta_t}&{delta_0}fs with N={N}&{chi_0}.pdf")

    plt.show()


if __name__ == "__main__":
    main()
