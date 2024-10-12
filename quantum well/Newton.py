import math
from numpy import sqrt, cos


df = lambda z, z0: 1 / (cos(z) ** 2) + z0**2 / (z**3 * sqrt(z0**2 / z**2 - 1))


def newton(f, p0, eps, N, z0):
    n = 1
    while n < N:
        p = p0 - f(p0, z0) / df(p0, z0)
        if abs(f(p0, z0)) < eps:
            break
        if abs(p - p0) < eps:
            break
        p0 = p
        n += 1
    return p
