###########
# Newton Raphson
#################### f(p) = f(p0) + (p-p0)f'(p0)
#################### f(p) = 0(vi p la nghiem) nen 0 = f(p0) + (p-p0)f'(p0)
#################### tao chuoi pn = pn+1 - f(pn-1)


import numpy as np

# from sympy import diff, symbols
from numpy import exp


def f(x):
    return exp(x) - 4 * x - 5


def dfN(x):
    return exp(x) - 4


def newton_raph(f, df, N, p0, eps):
    p = np.zeros(N)
    p[0] = p0
    count = 0
    for i in range(1, N):
        p[i] = p[i - 1] - f(p[i - 1]) / df(p[i - 1])
        count += 1

        if abs(p[i] - p[i - 1]) <= eps:
            return p, count
        if i + 1 >= N:
            print(f"nghiem khong hoi tu voi {count} lap.")
            return p, i


##############################################################
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
    N = 100
    eps = 1e-15
    p0 = 0
    p1 = 4
    # p, i = newton_raph(f, df, N, p0, eps)
    # p, i = newton_raph(f1, df1, N, p0, eps)
    print(secant(f, p0, p1, eps, N, df))


if __name__ == "__main__":
    main()
