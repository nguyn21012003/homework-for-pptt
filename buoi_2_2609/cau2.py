import numpy as np
import csv
from numpy import exp, pi
from math import log as ln, cos, sin
from tabulate import tabulate


def f(x):
    return cos(x) - x


def df_newton(x):
    return -sin(x) - 1


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


def main():
    a, b = 0, pi / 2

    p0 = b
    if not a <= p0 <= b:
        p0 = float(input("nhap lai p0: "))
    N = 100
    eps = 10e-8
    p, count = newtonr(f, df_newton, p0, N, eps)
    print(a, b, p)


if __name__ == "__main__":
    main()
