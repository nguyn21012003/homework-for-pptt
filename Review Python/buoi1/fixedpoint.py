import math
from math import exp as e, cos, sin, sqrt
import numpy as np
import csv


def fx1(x):
    return sqrt(10 - x**3) / 2


def fx2(x):
    return x - (x**3 + 4 * x**2 - 10)


def fx3(x):
    return sqrt(10 / x - 4 * x)


def fx4(x):
    return sqrt(10 / (x + 4))


def fx5(x):
    return x - (x**3 + 4 * x**2 - 10) / (3 * x**2 + 8 * x)


def fixedpoint(fx, p0, a, b, N, eps):
    if a <= p0 <= b:
        x = np.zeros(N)
        x[0] = p0
        count = 0
        for n in range(N):
            if n + 1 < N:
                try:
                    x[n + 1] = fx(x[n])
                    print(f"x[{n}] = {x[n]}")
                    if x[n] >= b * 10:
                        print(f"nghiem khong hoi tu voi {count} vong lap.")
                        return n, x
                    elif abs(x[n] - x[n - 1]) <= eps and x[n] <= b * 10:
                        print(f"nghiem hoi tu {x[count]} voi {count} vong lap.")
                        return n, x
                except (ZeroDivisionError, ValueError) as Z:
                    print(f"{Z}")
                    return n, x

                count += 1
            if n + 1 >= N:
                print(f"nghiem khong hoi tu voi {count} vong lap.")
                return n, x
    else:
        p0 = float(input("Nhập lại p0: "))


def main():
    p0 = float(input("Nhập p0: "))
    N = 20  ### so lan lap
    eps = 1e-5
    a, b = 1, 2
    sol1 = []
    sol2 = []
    sol3 = []
    sol4 = []
    sol5 = []

    n1, x1 = fixedpoint(fx1, p0, a, b, N, eps)
    print(x1)


if __name__ == "__main__":
    main()
