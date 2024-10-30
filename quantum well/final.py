import numpy as np
from math import pi
import csv
from tabulate import tabulate
from numpy import sin, exp, cos, tan, sqrt
import matplotlib.pyplot as plt


def para(_):
    """Description

    Args:
        _ (int): Nhập số trường hợp để tính z0 từ bé tới lớn tương đương với thế nông, chiều dài giếng hẹp cho tới thế sâu, chiều dài giếng rộng,

    Returns:
        tuple: Chứa các parameter cần dùng để giải bài toán ví dụ như a,m,V0,z0,....
    """
    a0 = 1e-9  # m
    m = 0.067 * 9.11e-31  # kg
    Ce = 1.6e-19  # eV to Joule(J)
    hbar_ev = 6.582119569e-16  # eV*s

    hbar_Js = 1.054571817e-34  # J*s
    match _:
        case 1:  # z0 ~ 6.707666865904062
            V0 = 1
            a = 10 * a0
            z0 = sqrt(m * V0 * a * a / (2 * (hbar_ev) ** 2 * Ce))
            return a, m, Ce, V0, z0, hbar_ev, hbar_Js
        case 2:  # z0 ~ 12.779821793191344
            V0 = 3
            a = 11 * a0
            z0 = sqrt(m * V0 * a * a / (2 * (hbar_ev) ** 2 * Ce))
            return a, m, Ce, V0, z0, hbar_ev, hbar_Js
        case 3:  # z0 ~ 26.830667463616248
            V0 = 4
            a = 20 * a0
            z0 = sqrt(m * V0 * a * a / (2 * (hbar_ev) ** 2 * Ce))
            return a, m, Ce, V0, z0, hbar_ev, hbar_Js
        case 4:  # z0 ~ 42.42301016380012
            V0 = 10
            a = 20 * a0
            z0 = sqrt(m * V0 * a * a / (2 * (hbar_ev) ** 2 * Ce))
            return a, m, Ce, V0, z0, hbar_ev, hbar_Js
        case 5:  # z0 ~ 237.15183634105392
            V0 = 50
            a = 50 * a0
            z0 = sqrt(m * V0 * a * a / (2 * (hbar_ev) ** 2 * Ce))
            return a, m, Ce, V0, z0, hbar_ev, hbar_Js
        case 6:  # z0 ~ 1341.5333731808123
            V0 = 100
            a = 200 * a0
            z0 = sqrt(m * V0 * a * a / (2 * (hbar_ev) ** 2 * Ce))
            return a, m, Ce, V0, z0, hbar_ev, hbar_Js
        case 7:  # z0 ~ 42423.01016380013
            V0 = 0.1
            a = 0.2 * a0
            z0 = sqrt(m * V0 * a * a / (2 * (hbar_ev) ** 2 * Ce))
            return a, m, Ce, V0, z0, hbar_ev, hbar_Js
        case 8:  # z0 ~ 42423.01016380013
            V0 = 83e6
            a = 2.0e-15
            z0 = 0.094
            return a, m, Ce, V0, z0, hbar_ev, hbar_Js


a, m, Ce, V0, z0, hbar_ev, hbar_Js = para(int(input("Nhập trường hợp: ")))


def fzEven(z, z0):  # even_sol
    y1 = np.tan(z)

    y = np.sqrt((z0 / z) ** 2 - 1)
    return y1 - y


def readLog(file, solbisection, solnewton, solsecant, N):
    with open(file, "w", newline="") as writefile:
        header = ["n", "Bisection", "Newton Raphson", "Secant"]
        writer = csv.DictWriter(writefile, fieldnames=header)
        writer.writeheader()
        for i in range(N):
            writer.writerow(
                {
                    "n": i,
                    "Bisection": solbisection[i],
                    "Newton Raphson": solnewton[i],
                    "Secant": solsecant[i],
                }
            )
    with open(file, "r") as readfile:
        reader = csv.DictReader(readfile)
        rows = []
        for row in reader:
            rows.append(row)
        return tabulate(rows, headers="key", tablefmt="fancy_grid")


def tableLog(file, solbisection, solnewton, solsecant, N):
    with open(file, "w", newline="", encoding="utf-8") as writefile:
        header = [f"{'Khoảng':^6}", f"{'Bisection':^18}", f"{'Newton Raphson':^18}", f"{'Secant':^18}"]
        writer = csv.DictWriter(writefile, fieldnames=header, delimiter="|")
        writer.writeheader()
        for i in range(N):
            writer.writerow(
                {
                    f"{'Khoảng':^6}": f"{i:^6}",
                    f"{'Bisection':^18}": f"{solbisection[i]:^18}",
                    f"{'Newton Raphson':^18}": f"{solnewton[i]:^18}",
                    f"{'Secant':^18}": f"{solsecant[i]:^18}",
                }
            )


def plotWaveFunction(z, z0, a, V0, N, hbar_Js, m, Ce):

    x = np.linspace(-2, 2, N)

    kappa = []
    l = []
    energy = []

    for i in range(len(z)):
        E = (z[i] * hbar_Js) ** 2 / (2 * m * a**2) / Ce - V0
        l.append(z[i])  # công thức 2.155 slide # thứ nguyên là m^-1
        kappa.append(z[i] * tan(z[i]))  # công thức 2.156 slide # thứ nguyên là m^-1
        energy.append(E)  # công thức 2.173 Griffiths

    print(energy)

    def V(x, V0, a):
        Vp = np.zeros(len(x))
        for i in range(len(x)):
            if -1 <= x[i] <= 1:
                Vp[i] = -V0  # = 1 ev * 1.6*10^-19

        return Vp

    def Ψ(x, V0, l, kappa, a):

        # D = V0 / 60
        D = sqrt(1 / (1 + 1 / kappa))

        F = D * (cos(l)) / exp(-kappa)

        wave = np.zeros(len(x))

        for i in range(len(x)):
            if x[i] <= -1:
                wave[i] = F * exp(kappa * x[i])
            elif x[i] >= 1:
                wave[i] = F * exp(-kappa * x[i])
            else:
                wave[i] = D * cos(l * x[i])
        return wave

    fig, ax = plt.subplots(1, 3, figsize=(15, 7))
    U = V(x, V0, a)
    ax[2].plot(x, U)
    ax[1].grid(True)

    for i in range(len(l)):
        ax[1].plot(x, energy[i] + Ψ(x, V0, l[i], kappa[i], a))
        ax[2].axhline(y=energy[i], ls="--", xmin=0.1, xmax=0.9)
        ax[2].annotate(f"E= {energy[i]:10f} eV", xy=(x[-1], energy[i]), xycoords="data")

    z1 = np.arange(1e-15, 4 * np.pi, 1 / N)
    rhs = tan(z1)
    lhs = sqrt((z0 / z1) ** 2 - 1)

    ax[0].plot(z1, rhs)
    ax[0].plot(z1, lhs)
    ax[0].set_ylim(-0.2, 10 * z0)
    ax[0].grid(True)

    ax[0].legend([r"tan(z)", r"$\sqrt{\dfrac{z_0^2}{z^2} - 1}$"], loc="upper right")
    ax[1].legend([r"$\psi_1(x)$", r"$\psi_2(x)$", r"$\psi_3(x)$"], loc="upper right")

    ax[0].set_title("a)", y=0, pad=-25, verticalalignment="top")
    ax[1].set_title("b)", y=0, pad=-25, verticalalignment="top")
    ax[2].set_title("c)", y=0, pad=-25, verticalalignment="top")

    plt.savefig("Final.pdf")

    plt.show()


def bisection(f: float, a: float, b: float, eps: float, N: int, z0: float) -> float:
    """Description

    Sử dụng thuật toán Bisection để tìm nghiêm xấp xỉ của phương trình f(x) = 0 nằm trong [a,b].

    Parameters:
        f (function): function
                    Hàm cần tính toán.
        a,b (float): Giá trị đầu nằm trong khoảng [a,b].
        eps (float): Sai số
        N (int): Số vòng lặp
        z0 (float): Tham số.

    Returns:
        c: Nghiệm hay điểm cắt nhau hay điểm làm cho f(x) = 0. Nếu như hàm f(x) tại các điểm p0,p1,p_n cùng dấu nhau thì chúng ta được nghiệm không hội tụ với n vòng lặp
    """
    for i in range(N):
        c = (a + b) / 2
        if abs(f(c, z0)) < eps:
            break
        if f(a, z0) * f(c, z0) < 0:
            b = c
        elif f(c, z0) * f(b, z0) < 0:
            a = c
        else:
            break
    return c


def newton(f: float, p0: float, eps: float, N: int, z0: float) -> float:
    df = lambda z, z0: 1 / (cos(z) ** 2) + z0**2 / (z**3 * sqrt(z0**2 / z**2 - 1))
    """Description

    Sử dụng thuật toán Newton Raphson cho f(x) = 0 cho ra nghiệm xấp xỉ quanh p0

    Parameters:
        f (function): function
                    Hàm cần tính toán.
        p0 (float): Giá trị đầu nằm trong khoảng [a,b].
        eps (float): Sai số
        N (int): Số vòng lặp
        z0 (float): Tham số.

    Returns:
        p_n: Nghiệm hay điểm cắt nhau hay điểm làm cho f(x) = 0.
    """
    p0_n = p0
    for n in range(0, N + 1):

        p0_n = p0_n - f(p0_n, z0) / df(p0_n, z0)
        if df(p0_n, z0) == 0:
            break
        if abs(f(p0_n, z0)) < eps:
            break

    return p0_n


def secant(f: float, p0: float, p1: float, eps: float, N: int, z0: float) -> float:
    """Description

    Sử dụng thuật toán Secant cho f(x) = 0 cho ra nghiệm xấp xỉ trên đoạn [p0,p1].

    Parameters:
        f (function): function
                    Hàm cần tính toán.
        p0,p1 (float): Giá trị đầu nằm trong khoảng [a,b].
        eps (float): Sai số
        N (int): Số vòng lặp
        z0 (float): Tham số.

    Returns:
        p_n: Nghiệm hay điểm cắt nhau hay điểm làm cho f(x) = 0. Nếu như hàm f(x) tại các điểm p0,p1,p_n cùng dấu nhau thì chúng ta được nghiệm không hội tụ với n vòng lặp
    """

    p0_n = p0
    p1_n = p1
    for n in range(N + 1):
        p_n = p0_n - f(p0_n, z0) * (p1_n - p0_n) / (f(p1_n, z0) - f(p0_n, z0))

        if f(p0_n, z0) * f(p_n, z0) < 0:
            p0_n = p0_n
            p1_n = p_n
        elif f(p1_n, z0) * f(p_n, z0) < 0:
            p0_n = p_n
            p1_n = p1_n

        if abs(f(p_n, z0)) < eps:
            break

        if abs(p1_n - p0_n) < eps:
            break

        if f(p_n, z0) == 0:
            break
    return p_n


def generateIntervals(delta_inte):
    interval = {}

    key = 0

    while True:

        inter_i = []
        p0 = 0 + key * pi + delta_inte
        p1 = pi / 2 + key * pi - delta_inte

        if p1 > z0:
            p1 = z0 - delta_inte
            inter_i.append(float(p0))
            inter_i.append(float(p1))
            interval[key] = inter_i
            break

        if p1 == pi / 2 + key * pi:
            break

        inter_i.append(float(p0))
        inter_i.append(float(p1))
        interval[key] = inter_i
        key += 1

    check_inter = {}
    for key, value in interval.items():
        if value[0] <= value[1]:
            check_inter[key] = value
    print(f"Có {len(check_inter)} đoạn cho ra nghiệm z \n")
    return check_inter


def main():
    file_log = "datafinitewell.txt"
    file_table_log = "datafinitewell.table.txt"

    N = 10000
    eps = 1e-15
    delta_inte = 1e-4

    listSolSecant = []
    listSolNewton = []
    listSolBisection = []

    check_inter = generateIntervals(delta_inte)

    for i in check_inter:

        p0, p1 = check_inter[i]
        solBisection = bisection(fzEven, p0, p1, eps, N, z0)
        solNewton = newton(fzEven, p1, eps, N, z0)
        solSecant = secant(fzEven, p0, p1, eps, N, z0)

        listSolBisection.append(float(solBisection))
        listSolNewton.append(float(solNewton))
        listSolSecant.append(float(solSecant))

    print(check_inter, "\n")
    print(listSolBisection, "Bisection \n")
    print(listSolNewton, "Newton \n")
    print(listSolSecant, "Secant \n")
    print(z0)
    # print(delta_inte)

    tableLog(file_table_log, listSolBisection, listSolNewton, listSolSecant, len(listSolSecant))

    plotWaveFunction(listSolBisection, z0, a, V0, N, hbar_Js, m, Ce)


if __name__ == "__main__":
    main()

"""def fzOdd(z, z0):  # odd_sol
    y2 = -1 / np.tan(z)
    y = np.sqrt((z0 / z) ** 2 - 1)
    return y2 - y"""
