import numpy as np
from numpy import sqrt, sin, cos, tan, exp
from math import pi
import csv
from tabulate import tabulate
from scipy import optimize

import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider

####################### Local

from para import case_
from bisection import bisection
from Secant import secant
from Newton import newton

a, m, Ce, V0, z0, hbar_ev, hbar_Js = case_(int(input("Nhập trường hợp: ")))


def fz_even(z, z0):  # even_sol
    y1 = np.tan(z)

    y = np.sqrt((z0 / z) ** 2 - 1)
    return y1 - y


def fz_odd(z, z0):  # odd_sol
    y2 = -1 / np.tan(z)
    y = np.sqrt((z0 / z) ** 2 - 1)
    return y2 - y


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


def main():
    file_log = "datafinitewell.txt"
    file_table_log = "datafinitewell.table.txt"

    N = 10000
    eps = 1e-15
    delta_inte = 0.1

    interval = {}
    print(z0)

    k = 0
    while k < int(z0):
        inter_i = []
        p0 = 0 + k * pi + delta_inte
        p1 = pi / 2 + k * pi - delta_inte
        if p1 > z0:
            p1 = z0 - delta_inte
            inter_i.append(float(p0))
            inter_i.append(float(p1))
            interval[k] = inter_i
            print(k)
            break

        if p1 == pi / 2 + k * pi:
            print("ngu")
            break
        inter_i.append(float(p0))
        inter_i.append(float(p1))
        interval[k] = inter_i
        k += 1
    check_inter = {}
    for k, v in interval.items():
        if v[0] <= v[1]:
            check_inter[k] = v
    print(check_inter)

    sol_Secant = []
    sol_Newton = []
    for i in check_inter:
        p0, p1 = check_inter[i]
        solSecant = secant(fz_even, p0, p1, eps, N, z0)
        solNewton = newton(fz_even, p1, eps, N, z0)
        sol_Secant.append(float(solSecant))
        sol_Newton.append(float(solNewton))
    print(sol_Secant)
    print(sol_Newton)


if __name__ == "__main__":
    main()
