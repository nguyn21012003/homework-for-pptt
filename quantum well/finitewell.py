import numpy as np
from numpy import sqrt, sin, cos, tan
from math import pi
import matplotlib.pyplot as plt
from math import exp


def z0_cal():
    hbar = 1.05457182e-34
    m = 9.11e-31
    a = 0.5e-9
    V0 = 25 * 1.61e-19
    z0 = a / hbar * sqrt(2 * m * V0)
    return z0


def f(z):
    z0 = z0_cal()
    f = tan(z) - sqrt((z0 / z) ** 2 - 1)
    return f


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


def df_newton(x):
    return exp(x) - 4


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


def save_log(arg):
    pass


def plot_data(arg):
    pass


def main():
    N = 100
    x0 = 4
    x1 = 5
    eps = 10e-15

    p0 = x0
    p1 = x1

    sol_bisection, count1 = bisection(f, x0, x1, N, eps)
    na, nb, nc = sol_bisection
    print(nc, "\n")

    sol_newton, count2 = newtonr(f, df_newton, p0, N, eps)
    print(sol_newton, "\n")

    sol_secant, count3 = secant(f, df_secant, p0, p1, N, eps)
    print(sol_secant, "\n")


if __name__ == "__main__":
    main()
