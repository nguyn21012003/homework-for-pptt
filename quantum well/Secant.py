import math
import numpy


def secant(f, p0, p1, eps, N, z0):
    i = 2
    while i < N:
        f_p0 = f(p0, z0)
        f_p1 = f(p1, z0)

        p2 = p1 - f_p1 * (p1 - p0) / float(f_p1 - f_p0)

        if abs(f_p1) < eps:
            break

        if abs(p2 - p1) < eps:
            break

        p0, p1 = p1, p2
        i += 1

    return p1
