import numpy as np
from numpy import tan, sqrt
from math import pi


def function(z, z0, N):
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


def main():
    N = 100
    hbar = 1.05457182e-34
    m = 9.31e-31
    V0 = 32 * hbar**2 / (m * a**2)
    a = 1
    z0 = a / hbar * sqrt(2 * m * V0)
    function(z0, hbar, N)


if __name__ == "__main__":
    main()
