import numpy as np
from numpy import sqrt, sin, cos, tan, exp
from math import pi
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider
import csv
from tabulate import tabulate

a0 = 0.529e-10  ## Hang so Bohr
m = 9.109e-31  # kg
Ce = 1.6e-19  # MeV to Joule(J)
hbar = 1.0546e-34  # J*s

a = None
V0 = None
z0 = None


def para(id):
    global a, V0, z0
    match id:
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

        case 3:  ########### Sâu,rộng => z0 lớN ~ 2.4222738389475014
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
            print(z0)


def fz(z):
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

            if f(nc[i]) == 0:
                break

            if abs(f(nc[i])) < eps:
                break

            if f(na[i]) * f(nc[i]) < 0:
                nb[i + 1] = nc[i]
                na[i + 1] = na[i]
                count += 1

            elif f(nc[i]) * f(nb[i]) < 0:
                nb[i + 1] = nb[i]
                na[i + 1] = nc[i]
                count += 1

    return (na, nb, nc), count


def df_newton(z):
    return 1 / (cos(z) ** 2) + (z0) ** 2 / (z**2 * sqrt(z0**2 - z**2))


def newtonr(f, df, p0, N, eps):
    p = np.zeros(N)
    p[0] = p0
    count = 0
    listp = []
    for i in range(1, N + 1):
        if i + 1 < N:
            p[i] = p[i - 1] - f(p[i - 1]) / df(p[i - 1])
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
    return (f(pn1) - f(pn2)) / (pn1 - pn2)


def secant(f, df, p0, p1, N, eps):
    p = np.zeros(N)
    p[0] = p0
    p[1] = p1
    count = 0
    listp = []
    for i in range(2, N + 2):
        if i + 1 < N:
            p[i] = p[i - 1] - f(p[i - 1]) / df(f, p[i - 1], p[i - 2])
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


def calc_psi(z, z0, a):

    k = sqrt(z0**2 - z**2) / a
    l = z / a

    Ek = -((k * hbar) ** 2) / (2 * m)
    El = ((l * hbar) ** 2) / (2 * m) - V0

    print(Ek)
    print(El)
    D = 1 / sqrt(a + 1 / k)
    F = exp(k * a) * cos(l * a) / sqrt(a + 1 / k)

    print(D, "\n")
    print(F)

    return D, F, k, l


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


def plot_data(x0, x1, z, z0, N):
    fig, axs = plt.subplots(2, 2, figsize=(15, 7))
    z = np.arange(0, 5 * pi, 1 / N)
    rhs = tan(z)
    lhs = sqrt((z0 / z) ** 2 - 1)
    axs[0, 0].plot(z, rhs)
    axs[0, 0].plot(z, lhs)
    axs[0, 0].set_ylim(-0.5, z0)
    axs[0, 0].grid(True)

    axs[0, 0].axhline(0, color="red", linewidth=0.5)
    axs[0, 0].axvline(0, color="red", linewidth=0.5)

    plt.show()


def fname(arg):
    pass


def main():

    id = int(input("Nhập trường hợp muốn tính toán: "))

    para(id)

    file_log = "datafinitewell.txt"
    file_table_log = "datafinitewell.table.txt"

    N = 100
    eps = 1e-14

    x0 = 0.1
    x1 = 1.5

    p0 = x0
    p1 = x1

    sol_bisection, count1 = bisection(fz, x0, x1, N, eps)
    na, nb, z_bisection = sol_bisection

    z_newton, count2 = newtonr(fz, df_newton, p0, N, eps)

    z_secant, count3 = secant(fz, df_secant, p0, p1, N, eps)

    table_log(file_table_log, z_bisection, z_newton, z_secant, N)
    read = read_log(file_log, z_bisection, z_newton, z_secant, N)
    print(read)

    calc_psi(z_newton, z0, a)

    plot_data(x0, x1, z_secant, z0, N)


if __name__ == "__main__":
    main()
