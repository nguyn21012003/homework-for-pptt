import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, log, pi

from numpy import real as RE, imag as IM, conjugate
from numpy import typing as npt
from time import time
from matplotlib.backends.backend_pdf import PdfPages
from tqdm import tqdm
import csv, os

from sbeFile.constant import N, hbar, chi_0, E_R, delta_t, Delta_0, t_max, t0, dt, delta_e, E0, constantE, C0


# Input variables

# print(C0)
#### Nếu mà ∆_0 là hệ số detunning(hệ số vượt rào cao thì sẽ dẫn đến xung(trường điện từ) đi sâu vô hơn, nhưng chi_0 càng bé là cường độ xung càng bé dẫn đến bị tắt dần càng nhanh )


def rk4(dF, tn, yn, h, T2):
    # tn is value in tSpan  array
    k1 = dF(tn, yn, T2)
    k2 = dF(tn + 0.5 * h, yn + 0.5 * h * k1, T2)
    k3 = dF(tn + 0.5 * h, yn + 0.5 * h * k2, T2)
    k4 = dF(tn + h, yn + h * k3, T2)

    return yn + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def g(n, n1):
    return (1 / sqrt(n * delta_e)) * log(abs((sqrt(n) + sqrt(n1)) / (sqrt(n) - sqrt(n1))))


def E_n(n, g, f_e, f_h):

    E = 0
    for n1 in range(N + 1):
        if n1 == n:
            continue
        E += (sqrt(E_R) / pi) * delta_e * g(n, n1) * (f_e[n1] + f_h[n1])

    return E


def omega_R(n, t, g, p_n):

    summ = 0
    for n1 in range(N + 1):
        if n1 == n:
            continue
        summ += g(n, n1) * p_n[n1]
    OMEGA = (0.5 * (hbar * sqrt(pi) / delta_t) * chi_0 * np.exp(-(t**2) / delta_t**2) + (sqrt(E_R) / pi) * delta_e * summ) / hbar

    return OMEGA


def dF(t, Y, T2):

    F = np.zeros([2, N + 1], dtype="complex")

    # print(RE(Y[0]))
    for n in range(1, N + 1):
        E = E_n(n, g, RE(Y[0]), IM(Y[0]))
        OMEGA = omega_R(n, t, g, Y[1])

        F[0][n] = -2 * IM(OMEGA * conjugate(Y[1][n]))
        F[1][n] = -(1j / hbar) * (n * delta_e - Delta_0 - E) * Y[1][n] + 1j * (1 - RE(Y[0])[n] - IM(Y[0])[n]) * OMEGA - (Y[1][n] / T2)

    return F


def solve_sys_ode(dF: npt.NDArray, dt: float, rk4: npt.NDArray, N: int) -> npt.NDArray:

    # tSpan = np.linspace(t0, t_max, N)
    tSpan = np.arange(t0, t_max, dt)
    # print(tSpan)
    # print((t_max - t0) / len(tSpan))

    Y_T2Fix = np.zeros([2, N + 1], dtype="complex")
    Polarization_T2Fix = np.zeros(len(tSpan))
    NumberDensity_T2Fix = np.zeros(len(tSpan))
    EnergyEps_T2Fix = np.zeros(N)
    fe_T2Fix = []
    p_T2Fix = []
    # print(len(Y))
    T2Fix = 210
    for ti in tqdm(range(len(tSpan)), desc="Processing", unit="step", ascii=" #"):
        gamma = 6.5e-20
        Y_T2Fix = rk4(dF, tSpan[ti], Y_T2Fix, dt, T2Fix)
        # print(RE(Y[0]))
        fe_t_T2Fix = []
        p_t_T2Fix = []
        NumberDensity_sum_T2Fix = 0
        Polarization_sum_T2Fix = 0
        for n in range(N + 1):
            if RE(Y_T2Fix[0][n]) != 0.0:
                fe_t_T2Fix.append(RE(Y_T2Fix[0][n]))
                p_t_T2Fix.append((abs(Y_T2Fix[1][n])))
            # print(fe_t)

            EnergyEps_T2Fix[n - 1] = n * delta_e

            NumberDensity_sum_T2Fix += C0 * sqrt(n + 1) * RE(Y_T2Fix[0][n])
            Polarization_sum_T2Fix += constantE * sqrt(n + 1) * RE(Y_T2Fix[1][n])

        NumberDensity_T2Fix[ti] = NumberDensity_sum_T2Fix
        Polarization_T2Fix[ti] = Polarization_sum_T2Fix
        # print(Polarization)

        fe_T2Fix.append(fe_t_T2Fix)
        p_T2Fix.append(p_t_T2Fix)

    # listEom = []
    listAlpha_T2Fix = []
    omega_T2Fix = np.linspace(-200, 200, len(tSpan))
    for omega_i_T2Fix in range(len(omega_T2Fix)):
        # omega[omega_i] += omega0 + omega_i * domega
        Piomega_T2Fix = 0
        Eiomega_T2Fix = 0
        for i in range(0, len(tSpan)):
            Piomega_T2Fix += tSpan[i] * Polarization_T2Fix[i] * np.exp(1j * omega_T2Fix[omega_i_T2Fix] * tSpan[i] / hbar)
            Eiomega_T2Fix += E0 * tSpan[i] * np.exp(1j * omega_T2Fix[omega_i_T2Fix] * tSpan[i] / hbar) * np.exp(-(tSpan[i] ** 2) / (delta_t**2))

        alpha_T2Fix = IM(Piomega_T2Fix / Eiomega_T2Fix)
        listAlpha_T2Fix.append(alpha_T2Fix)

    np.array(listAlpha_T2Fix)

    ####################################################################### thay đổi T2 theo thời gian và N(t)
    tSpan = np.arange(t0, t_max, dt)
    # print(tSpan)
    # print((t_max - t0) / len(tSpan))

    Y = np.zeros([2, N + 1], dtype="complex")
    Polarization = np.zeros(len(tSpan))
    NumberDensity = np.zeros(len(tSpan))
    EnergyEps = np.zeros(N)
    fe = []
    p = []
    T2 = 200

    for ti in tqdm(range(len(tSpan)), desc="Processing", unit="step", ascii=" #"):
        gamma = 6.5e-20
        Y = rk4(dF, tSpan[ti], Y, dt, T2)
        # print(RE(Y[0]))
        fe_t = []
        p_t = []
        NumberDensity_sum = 0
        Polarization_sum = 0
        for n in range(N + 1):
            if RE(Y[0][n]) != 0.0:
                fe_t.append(RE(Y[0][n]))
                p_t.append((abs(Y[1][n])))
            # print(fe_t)

            EnergyEps[n - 1] = n * delta_e

            NumberDensity_sum += C0 * sqrt(n + 1) * RE(Y[0][n])
            Polarization_sum += constantE * sqrt(n + 1) * RE(Y[1][n])

        NumberDensity[ti] = NumberDensity_sum
        Polarization[ti] = Polarization_sum
        T2 = 1 / (1 / T2 + gamma * NumberDensity_sum)
        # print(Polarization)

        fe.append(fe_t)
        p.append(p_t)

    listEom = []
    listAlpha = []
    omega = np.linspace(-250, 250, len(tSpan))
    for omega_i in range(len(omega)):
        # omega[omega_i] += omega0 + omega_i * domega
        Piomega = 0
        Eiomega = 0
        for i in range(0, len(tSpan)):
            Piomega += tSpan[i] * Polarization[i] * np.exp(1j * omega[omega_i] * tSpan[i] / hbar)
            Eiomega += E0 * tSpan[i] * np.exp(1j * omega[omega_i] * tSpan[i] / hbar) * np.exp(-(tSpan[i] ** 2) / (delta_t**2))

        alpha = IM(Piomega / Eiomega)
        listAlpha.append(alpha)

    np.array(listAlpha)

    return tSpan, EnergyEps, fe, p, Polarization, NumberDensity, listAlpha, listAlpha_T2Fix  # , listEom


def multipage(filename, figs=None):

    pp = PdfPages(filename)
    if figs is None:
        figs = [plt.figure(n) for n in plt.get_fignums()]
    for fig in figs:
        fig.savefig(pp, format="pdf")
    pp.close()


def createFolder(folderName):
    if not os.path.exists(folderName):
        os.makedirs(folderName)
    else:
        pass


def saveLog(T, E, fe, p, Polarization, NumberDensity, t):
    file = "SBE.data.txt"
    with open(file, "w") as writefile:
        header = ["t"]
        writer = csv.DictWriter(writefile, fieldnames=header, delimiter="\t")
        writer.writeheader()
        print(len(t))


def plot(T, E, fe, p, Polarization, NumberDensity, t, Alpha, Alpha_T2Fix):
    fig1 = plt.figure(figsize=(16, 10))
    fig2 = plt.figure(figsize=(16, 10))
    fig3 = plt.figure(figsize=(16, 10))

    color = "inferno"
    ombre = "copper"
    transparent = 0.8

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

    omega = np.linspace(-250, 250, len(t))
    ax5 = fig3.add_subplot(121)
    ax5.plot(omega, Alpha)

    ax6 = fig3.add_subplot(122)
    ax6.plot(omega, Alpha_T2Fix)

    # ax3 = fig.add_subplot(133, projection="3d")
    # ax3.plot_surface(T, E, P, cmap="copper")

    folderName = "data"
    createFolder(folderName)
    pdf_path = os.path.join(folderName, f"{delta_t}&{Delta_0}fs with N={N}&{chi_0}.pdf")

    multipage(pdf_path)

    plt.show()


def main():

    start = time()

    t, EnergyEps, fe, p, Polarization, NumberDensity, listAlpha, listAlpha_T2Fix = solve_sys_ode(dF, dt, rk4, N)

    end = time()

    print(end - start)

    E, T = np.meshgrid(EnergyEps, t)
    fe = np.array(fe)
    p = np.array(p)

    # print(fe.shape)
    # print(T.shape)
    # print(E.shape)
    # print(T)
    # print(len(fe))

    with open("data.txt", "w") as writefile:
        header = ["t", "Eps", "fe"]
        writer = csv.DictWriter(writefile, fieldnames=header)
        writer.writeheader()
        for i in range(len(T)):
            writer.writerow(
                {
                    "t": f"{i}",
                    "Eps": f"{E[i]}",
                    "fe": f"{fe[i]}",
                }
            )
    plot(T, E, fe, p, Polarization, NumberDensity, t, listAlpha, listAlpha_T2Fix)
    saveLog(T, E, fe, p, Polarization, NumberDensity, t)


if __name__ == "__main__":
    main()
