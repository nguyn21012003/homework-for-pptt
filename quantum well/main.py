import numpy as np
from math import pi
import csv
from tabulate import tabulate

####################### Local

from para import para
from Secant import secant
from Newton import newton
from Bisection import bisection
from wavefunction import plotWaveFunction

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

    plotWaveFunction(listSolBisection, a, V0, N, hbar_Js, m, Ce)


if __name__ == "__main__":
    main()

"""def fzOdd(z, z0):  # odd_sol
    y2 = -1 / np.tan(z)
    y = np.sqrt((z0 / z) ** 2 - 1)
    return y2 - y"""
