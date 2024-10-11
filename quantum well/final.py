import numpy as np
from numpy import sqrt, sin, cos, tan, exp
from math import pi
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider
import csv
from tabulate import tabulate
from scipy import optimize

####################### Local

from para import case_
from bisection import bisection
from Secant import secant

a, m, Ce, V0, z0, hbar_ev, hbar_Js = case_(int(input("Nhập trường hợp: ")))
print(z0)


def fz_even(z, z0):  # even_sol
    return tan(z) - sqrt((z0 / z) ** 2 - 1)


def fz_odd(z, z0):  # odd_sol
    return -1 / tan(z) - sqrt((z0 / z) ** 2 - 1)


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
    even_solution = np.array([], dtype=float)
    odd_solution = np.array([], dtype=float)
    income_wave = np.arange(h_step, z0, h_step)

    even = True
    for i in range(len(income_wave) - 1):
        if even == True:
            if fz_even(income_wave[i], z0) * fz_even(income_wave[i + 1], z0) < 0:
                even_solution = np.append(even_solution, optimize.brentq(fz_even, income_wave[i], income_wave[i + 1], args=z0))
                even = False
        if even == False:
            if fz_odd(income_wave[i], z0) * fz_odd(income_wave[i + 1], z0) < 0:
                odd_solution = np.append(odd_solution, optimize.brentq(fz_odd, income_wave[i], income_wave[i + 1], args=z0))
                even = True
    return even_solution, odd_solution


esol, osol = sqWellsolution(z0, 0.1)

el = esol * 2 / a
ol = osol * 2 / a

ekappa = 2 * esol / a * tan(esol)
okappa = abs(2 * osol) / a * tan(osol)

eneer_even = (2 * hbar_Js**2 * esol**2 / (m * a**2)) / Ce
eneer_odd = (2 * hbar_Js**2 * osol**2 / (m * a**2)) / Ce

print(esol)
print(osol)


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
    lhs = sqrt((z0 / z) ** 2 - 1)

    axs[0, 0].plot(z1, rhs, "magenta")
    axs[0, 0].plot(z, lhs, "navy")
    axs[0, 0].set_ylim(-0.2, 10 * z0)
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
        axs[0, 1].plot(x, -(eneer_even[i] + wavefunc(x, V0, el[i], ekappa[i], a)))
        axs[0, 1].axhline(y=-eneer_even[i], xmin=0.1666, xmax=0.8333)
        axs[0, 1].annotate(f"E= {-eneer_even[i]:3.2f} eV", xy=(x[-1], eneer_even[i]))

    plt.show()


def main():
    file_log = "datafinitewell.txt"
    file_table_log = "datafinitewell.table.txt"

    N = 100
    eps = 1e-15

    x0 = 1
    x1 = 2

    interval = {
        0: [0 + eps, pi / 2 - eps],
    }

    k = 0
    while k < int(z0):
        inter_i = []
        p0 = pi / 2 + k * pi + 2 * eps
        p1 = pi / 2 + (k + 1) * pi - 2 * eps
        interval[k + 1] = inter_i
        if p1 >= z0:
            p1 = z0 - eps
            inter_i.append(p0)
            inter_i.append(p1)

            break
        inter_i.append(p0)
        inter_i.append(p1)
        k += 1

    # sol_bisection = bisection(fz_even, x0, x1, N, eps)
    # na, nb, z_bisection = sol_bisection
    sol_Secant = []
    for i in interval:
        p0, p1 = interval[i]
        ket_qua = secant(fz_even, p0, p1, eps, N, z0)
        sol_Secant.append(ket_qua)
    plot_data(x0, x1, sol_Secant, z0, N)


if __name__ == "__main__":
    main()
