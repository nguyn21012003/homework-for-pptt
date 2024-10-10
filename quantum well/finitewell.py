import numpy as np
from numpy import sqrt, sin, cos, tan, exp
from math import pi
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider
import csv
from tabulate import tabulate
from scipy import optimize


a0 = 0.529e-9  ## Hang so Bohr
m = 9.31e-31
Ce = 1.6e-19  # MeV to Joule(J)
hbar_ev = 6.582119569e-16  # eV*s

hbar_Js = 1.054571817e-34

a = 1 * a0
V0 = 10
z0 = sqrt(0.067 * m * V0 * a * a / (2 * (hbar_ev) ** 2 * Ce))


"""def para(id_):
    global a, V0, z0
    match id_:
        case 1:  ########### Nông,hẹp => z0 bé ~ 2.4222738389475014
            a = 2 * a0  # meter
            V0 = 20 * Ce
            z0 = a / hbar * sqrt(2 * m * V0)
            print(z0)

        case 2:  ########### Sâu,rộng => z0 lớn ~ 7.266821516842505
            a = 6 * a0  # meter
            V0 = 20 * Ce
            z0 = a / hbar * sqrt(2 * m * V0)
            print(z0)

        case 3:  ########### Nông,hẹp => z0 bé ~ 0.012111369194737508
            a = 0.1 * a0  # meter
            V0 = 0.2 * Ce
            z0 = a / hbar * sqrt(2 * m * V0)
            print(z0)

        case 4:  ########### Sâu,rộng => z0 lớN ~ 135
            a = 50 * a0  # meter
            V0 = 100 * Ce
            z0 = a / hbar * sqrt(2 * m * V0)
            print(z0)

        case 5:
            a = 50 * a0  # meter
            V0 = 100 * Ce
            z0 = a / hbar * sqrt(2 * m * V0)
            print(z0)"""


def fz(z, z0):  # even_sol
    return tan(z) - sqrt((z0 / z) ** 2 - 1)


def bisection(f, a, b, N, eps):
    a = float(a)
    b = float(b)
    if a > b:
        a = b
        b = a

    na = np.zeros(N)
    nb = np.zeros(N)
    nc = np.zeros(N)
    na[0] = a
    nb[0] = b

    count = 0

    for i in range(N):
        if i + 1 < N:
            nc[i] = (na[i] + nb[i]) / 2

            if f(nc[i], z0) == 0:
                break

            if abs(f(nc[i], z0)) < eps:
                break

            if f(na[i], z0) * f(nc[i], z0) < 0:
                nb[i + 1] = nc[i]
                na[i + 1] = na[i]
                count += 1

            elif f(nc[i], z0) * f(nb[i], z0) < 0:
                nb[i + 1] = nb[i]
                na[i + 1] = nc[i]
                count += 1

    return (na, nb, nc), count


def df_newton(z, z0):
    return 1 / (cos(z) ** 2) + z0**2 / (z**2 * sqrt(z0**2 - z**2))


def newtonr(f, df, p0, N, eps):
    p = np.zeros(N)
    p[0] = p0
    count = 0
    listp = []
    for i in range(1, N + 1):
        if i + 1 < N:
            p[i] = p[i - 1] - f(p[i - 1], z0) / df(p[i - 1], z0)
            count += 1
            if abs(p[i] - p[i - 1]) <= eps:
                break
        else:
            print(f"Nghiem khong hoi tu voi {count} vong lap !")
            break

    # for i in range(len(p)):
    #    if p[i] != 0:
    #        listp.append(float(p[i]))

    return p, count


def df_secant(f, pn1, pn2):
    return (f(pn1, z0) - f(pn2, z0)) / (pn1 - pn2)


def secant(f, df, p0, p1, N, eps):
    p = np.zeros(N)
    p[0] = p0
    p[1] = p1
    count = 0
    listp = []
    for i in range(2, N + 2):
        if i + 1 < N:
            p[i] = p[i - 1] - f(p[i - 1], z0) / df(f, p[i - 1], p[i - 2])
            count += 1
            if abs(p[i] - p[i - 1]) <= eps:
                break
        else:
            print(f"Nghiem khong hoi tu voi {count} vong lap !")
            break

    # for i in range(len(p)):
    #    if p[i] != 0:
    #        listp.append(float(p[i]))

    return p, count


def read_log(file, solbisection, solnewton, solsecant, N):
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
        return tabulate(rows, headers="keys", tablefmt="fancy_grid")


def table_log(file, solbisection, solnewton, solsecant, N):
    with open(file, "w", newline="") as writefile:
        header = [f"{'n':^3}", f"{'Bisection':^18}", f"{'Newton Raphson':^18}", f"{'Secant':^18}"]
        writer = csv.DictWriter(writefile, fieldnames=header)
        writer.writeheader()
        for i in range(N):
            writer.writerow(
                {
                    f"{'n':^3}": f"{i:^3}",
                    f"{'Bisection':^18}": f"{solbisection[i]:^18}",
                    f"{'Newton Raphson':^18}": f"{solnewton[i]:^18}",
                    f"{'Secant':^18}": f"{solsecant[i]:^18}",
                }
            )


def V(x, V0, a):
    Vp = np.zeros(len(x))
    width = a
    height = V0
    for i in range(len(x)):
        if -width / 2 <= abs(x[i]) <= width / 2:
            Vp[i] = -V0

    return Vp


def sqWellsolution(z0: float, h_step: float) -> tuple:
    even_sol = np.array([], dtype=float)
    odd_sol = np.array([], dtype=float)
    income_wave = np.arange(h_step, z0, h_step)

    even = True
    for i in range(len(income_wave) - 1):
        if even == True:
            if fz(income_wave[i], z0) * fz(income_wave[i + 1], z0) < 0:
                even_sol = np.append(even_sol, optimize.brentq(fz, income_wave[i], income_wave[i + 1], args=z0))
            if even == False:
                if odd_sol(income_wave[i], z0) * odd_sol(income_wave[i + 1], z0) < 0:
                    odd_sol = np.append(odd_sol, optimize.brentq(odd_sol, income_wave[i], income_wave[i + 1], args=z0))

                even = True
    return even_sol, odd_sol


esol, osol = sqWellsolution(z0, 0.1)

el = esol * 2 / a
ol = osol * 2 / a

ekappa = 2 * esol / a * tan(esol)
okappa = abs(2 * osol) / a * tan(osol)

eneer_even = (2 * hbar_Js**2 * esol**2 / (m * a**2)) / Ce
eneer_odd = (2 * hbar_Js**2 * osol**2 / (m * a**2)) / Ce

print(z0, "nang luong")


def wavefunc(x, V0, el, ekappa, a):
    A = V0 / 20.0
    B = A * (cos(el * a / 2)) / exp(-ekappa * a / 2)

    wave = np.zeros(len(x))
    for i in range(len(x)):
        if x[i] <= -a / 2:
            wave[i] = B * exp(ekappa * x[i])
        elif x[i] >= a / 2:
            wave[i] = B * exp(-ekappa * x[i])
        else:
            wave[i] = A * cos(ekappa * x[i])
    return wave


def plot_data(x0, x1, z, z0, N):

    fig, axs = plt.subplots(2, 2, figsize=(15, 7))

    z1 = np.arange(0, 4 * pi, 1 / N)
    rhs = tan(z1)
    lhs = sqrt((z0 / z1) ** 2 - 1)

    axs[0, 0].plot(z1, rhs, "magenta")
    axs[0, 0].plot(z1, lhs, "navy")
    axs[0, 0].set_ylim(-0.2, z0)
    axs[0, 0].grid(True)

    axs[0, 0].axhline(0, color="red", linewidth=0.5)
    axs[0, 0].axvline(0, color="red", linewidth=0.5)

    x = np.linspace(-a, a, N)
    Vp = V(x, V0, a)
    axs[0, 1].plot(x, Vp)
    axs[0, 1].axhline(0, color="red", ls="--")
    axs[0, 1].set_xlabel("x", fontsize=14)
    axs[0, 1].set_ylabel("V(x)", fontsize=14)

    for i in range(len(el)):
        axs[0, 1].plot(x, eneer_even[i] + wavefunc(x, V0, el[i], ekappa[i], a))
        axs[0, 1].axhline(y=eneer_even[i], xmin=0.1666, xmax=0.8333)
        axs[0, 1].annotate("E= %3.2f eV" % eneer_even[i], xy=(x[-1], eneer_even[i]), xycoords="data")

    plt.show()


def main():

    # id_ = int(input("Nhập trường hợp muốn tính toán: "))

    # para(id_)

    file_log = "datafinitewell.txt"
    file_table_log = "datafinitewell.table.txt"

    N = 100
    eps = 1e-8

    interval = {
        "intv1": [0 + eps, pi / 2 - eps],
        "intv2": [pi / 2 + eps, pi - eps],
        "intv3": [pi + eps, 3 / 2 * pi - eps],
        "intv4": [3 / 2 * pi + eps, 4 / 2 * pi - eps],
        "intv5": [4 / 2 * pi + eps, 5 / 2 * pi - eps],
    }

    x0 = 0.1
    x1 = 1.5

    p0 = 1
    p1 = x1

    sol_bisection, count1 = bisection(fz, x0, x1, N, eps)
    na, nb, z_bisection = sol_bisection

    z_newton, count2 = newtonr(fz, df_newton, p0, N, eps)

    z_secant, count3 = secant(fz, df_secant, p0, p1, N, eps)

    table_log(file_table_log, z_bisection, z_newton, z_secant, N)
    read = read_log(file_log, z_bisection, z_newton, z_secant, N)
    # print(read)

    plot_data(x0, x1, z_newton, z0, N)
    sqWellsolution(z0, 0.001)


if __name__ == "__main__":
    main()
