import numpy as np
from numpy import sqrt, sin, cos, tan, exp
from para import case_


def bisection(f, a, b, N, eps, cs):

    L, z0, _, _ = case_(cs)

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

    for i in range(N):
        if i + 1 < N:
            nc[i] = (na[i] + nb[i]) / 2

            if f(nc[i], z0) == 0:
                break

            if abs(f(nc[i], z0)) < eps:
                break

            if f(na[i], z0) * f(nc[i], z0) < 0:
                nb[i + 1] = nc[i]
                na[i + 1] = na[i]

            elif f(nc[i], z0) * f(nb[i], z0) < 0:
                nb[i + 1] = nb[i]
                na[i + 1] = nc[i]

    return (na, nb, nc)
