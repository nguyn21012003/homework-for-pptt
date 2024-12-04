import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, log, pi
from numpy import real as RE, imag as IM, conjugate

########################################################
from matplotlib.backends.backend_pdf import PdfPages
from tqdm import tqdm
import csv, os, subprocess
from time import time
from numpy import typing as npt

#### Local Library
from sbeFile.constant import N, hbar, chi_0, E_R, delta_t, Delta_0, t_max, t0, dt, delta_e, T2, E0, C0


# print(C0)
#### Nếu mà ∆_0 là hệ số detunning(hệ số vượt rào cao thì sẽ dẫn đến xung(trường điện từ) đi sâu vô hơn, nhưng chi_0 càng bé là cường độ xung càng bé dẫn đến bị tắt dần càng nhanh )


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
    OMEGA = (0.5 * (hbar * sqrt(pi) / delta_t) * chi_0 * np.exp(-(t**2) / delta_t**2) + (sqrt(E_R) / pi) * delta_e * summ) / hbar

    return OMEGA


def dF(t, Y):

    F = np.zeros([2, N + 1], dtype="complex")

    # print(RE(Y[0]))
    for n in range(1, N + 1):
        E = E_n(n, g, RE(Y[0]), IM(Y[0]))
        OMEGA = omega_R(n, t, g, Y[1])

        F[0][n] = -2 * IM(OMEGA * conjugate(Y[1][n]))
        F[1][n] = -(1j / hbar) * (n * delta_e - Delta_0 - E) * Y[1][n] + 1j * (1 - RE(Y[0])[n] - IM(Y[0])[n]) * OMEGA - (Y[1][n] / T2)

    return F


def solve_sys_ode(dF: npt.NDArray, dt: float, rk4: npt.NDArray, N: int) -> npt.NDArray:
    fileWriteSBE = f"{delta_t}&{Delta_0}fs with N={N}&{chi_0}.data.sbe.txt"
    fileWriteAbsortion = f"{delta_t}&{Delta_0}fs with N={N}&{chi_0}.data.absortion.txt"

    print(fileWriteSBE)
    print(fileWriteAbsortion)
    # tSpan = np.linspace(t0, t_max, N)
    tSpan = np.arange(t0, t_max, dt)
    omega = np.linspace(-100, 100, len(tSpan))
    # print(tSpan)
    # print((t_max - t0) / len(tSpan))

    Y = np.zeros([2, N + 1], dtype="complex")
    EnergyEps = np.zeros(N + 1)

    # Polarization = np.zeros(len(tSpan), dtype="complex")
    Polarization = np.zeros(len(tSpan))
    NumberDensity = np.zeros(len(tSpan))
    # listEom = []
    listAlpha = []
    fe = []
    p = []
    # print(len(Y))
    with open(fileWriteSBE, "w", newline="") as writefile:
        header = [
            f"{'Thoi gian':^9}",
            f"{'EnergyEpsilon':^13}",
            f"{'fe_t':^40}",
            f"{'p_t':^40}",
            f"{'Mat do toan phan':^40}",
            f"{'Phan cuc toan phan':^40}",
        ]
        writer = csv.DictWriter(writefile, fieldnames=header, delimiter="\t")
        writer.writeheader()
        prev_ti = None

        for ti in tqdm(range(len(tSpan)), desc="Processing", unit="step ", ascii=" #"):

            Y = rk4(dF, tSpan[ti], Y, dt)
            # print(RE(Y[0]))
            fe_t = []
            p_t = []
            NumberDensity_sum = 0
            Polarization_sum = 0
            for n in range(N + 1):
                fe_t.append(RE(Y[0][n]))
                p_t.append(abs((Y[1][n])))
                # print(fe_t)

                EnergyEps[n] = n * delta_e

                NumberDensity_sum += C0 * sqrt(n + 1) * RE(Y[0][n])
                Polarization_sum += delta_e * sqrt(delta_e) * sqrt(n + 1) * abs(Y[1][n])

            NumberDensity[ti] = NumberDensity_sum
            Polarization[ti] = Polarization_sum
            # print(Polarization)

            fe.append(fe_t)
            p.append(p_t)
            if prev_ti is not None and ti != prev_ti:
                writer.writerow({})  ### Để ghi file cách dòng

            for i in range(len(fe_t)):

                writer.writerow(
                    {
                        f"{'Thoi gian':^9}": f"{i:^9}",
                        f"{'EnergyEpsilon':^13}": f"{EnergyEps[i]:^13}",
                        f"{'fe_t':^40}": f"{fe_t[i]:^40}",
                        f"{'p_t':^40}": f"{p_t[i]:^40}",
                        f"{'Mat do toan phan':^40}": f"{NumberDensity[i]:^40}",
                        f"{'Phan cuc toan phan':^40}": f"{Polarization[i]:^40}",
                    }
                )
            prev_ti = ti
    with open(fileWriteAbsortion, "w", newline="") as writefileAbsortion:
        header = [
            f"{'Thoi gian':^9}",
            f"{'tan so omega':^25}",
            f"{'Polarization':^60}",
            f"{'NumberDensity':^45}",
            f"{'Alpha':^21}",
            f"{'RE_POmega':^45}",
            f"{'IM_POmega':^45}",
            f"{'RE_EOmega':^45}",
            f"{'IM_EOmega':^45}",
        ]
        writerfileAbsortion = csv.DictWriter(writefileAbsortion, fieldnames=header, delimiter="\t")
        writerfileAbsortion.writeheader()

        for omega_i in range(len(omega)):

            Piomega = 0
            Eiomega = 0
            for i in range(0, len(tSpan)):
                Piomega += dt * Polarization[i] * np.exp(1j * omega[omega_i] * tSpan[i] / hbar)
                Eiomega += E0 * dt * np.exp(1j * omega[omega_i] * tSpan[i] / hbar) * np.exp(-(tSpan[i] ** 2) / (delta_t**2))

            alpha = IM(Piomega / Eiomega)
            listAlpha.append(alpha)
            writerfileAbsortion.writerow(
                {
                    f"{'Thoi gian':^9}": f"{tSpan[omega_i]:^9}",
                    f"{'tan so omega':^25}": f"{omega[omega_i]:^25}",
                    f"{'Polarization':^60}": f"{Polarization[omega_i]:^60}",
                    f"{'NumberDensity':^45}": f"{NumberDensity[omega_i]:^45}",
                    f"{'Alpha':^21}": f"{listAlpha[omega_i]:^21}",
                    f"{'RE_POmega':^45}": f"{RE(Piomega):^45}",
                    f"{'IM_POmega':^45}": f"{IM(Piomega):^45}",
                    f"{'RE_EOmega':^45}": f"{RE(Eiomega):^45}",
                    f"{'IM_EOmega':^45}": f"{IM(Eiomega):^45}",
                }
            )

        np.array(listAlpha)

    return tSpan, EnergyEps, fe, p, Polarization, NumberDensity, listAlpha  # , listEom


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


def plot(T, E, fe, p, Polarization, NumberDensity, t, Pom):
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

    ax5 = fig3.subplots()
    omega = np.linspace(-250, 250, len(t))
    ax5.plot(omega, Pom)

    folderName = "data"
    createFolder(folderName)
    pdf_path = os.path.join(folderName, f"{delta_t}&{Delta_0}fs with N={N}&{chi_0}.pdf")

    multipage(pdf_path)

    plt.show()


def gnuPlot():
    """
    Không xài dòng nào thì # vào dòng đó
    """
    fileWriteSBE = f"{delta_t}&{Delta_0}fs with N={N}&{chi_0}.data.sbe.txt"
    fileWriteAbsortion = f"{delta_t}&{Delta_0}fs with N={N}&{chi_0}.data.absortion.txt"
    print(fileWriteAbsortion)
    with open("gnuPlot.gp", "w") as gnuplot:
        gnuplot.write(
            f"""

    #set multiplot layout 1,2


    set ylabel "y"
    set xlabel "x"
    set zlabel "z"


    set grid
    set key horiz

    #splot "{fileWriteSBE}" u 1:2:3 with lines tit "fe_t" 
    
    
    #set ylabel "y"
    #set xlabel "x"
    #set zlabel "z"
    #splot "{fileWriteSBE}" u 1:2:4 with lines tit "p_t"


    set xlabel "Omega"
    set ylabel "alpha(omega)"
    #plot "{fileWriteAbsortion}" u 2:3 with lines tit "Phân cực toàn phần"
    #plot "{fileWriteAbsortion}" u 2:4 with lines tit "Mật độ toàn phần"
    plot "{fileWriteAbsortion}" u 2:5 with lines tit "Phổ hấp thụ"


    #unset multiplot

    pause -1

"""
        )
    subprocess.run(["gnuplot", "gnuPlot.gp"])


def main():

    start = time()

    t, EnergyEps, fe, p, Polarization, NumberDensity, listPom = solve_sys_ode(dF, dt, rk4, N)

    end = time()

    print(end - start)

    E, T = np.meshgrid(EnergyEps, t)
    fe = np.array(fe)
    p = np.array(p)

    plot(T, E, fe, p, Polarization, NumberDensity, t, listPom)
    gnuPlot()


if __name__ == "__main__":
    main()
