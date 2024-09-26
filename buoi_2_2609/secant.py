from numpy import exp
import numpy as np


def f(x):
    return exp(x) - 4 * x - 5


def df(f, pn2, pn1):
    df = (f(pn2) - f(pn1)) / (pn2 - pn1)
    return df


def secant(f, p0, p1, eps, N, df):

    p = np.zeros(N)
    p[0] = p0
    p[1] = p1
    count = 0
    for i in range(2, N):
        p[i] = p[i - 1] - f(p[i - 1]) / df(f, p[i - 2], p[i - 1])
        count += 1
        if abs(p[i] - p[i - 1]) <= eps:
            return p, count
        if i + 1 >= N:
            print(f"nghiem khong hoi tu voi {count} lap.")
            return p, i


def main():
    p0 = 5
    p1 = 6
    N = 100
    eps = 1e-5
    print(secant(f, p0, p1, eps, N, df))


if __name__ == "__main__":
    main()
