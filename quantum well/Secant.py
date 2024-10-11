# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 18:47:55 2024

@author: hoann
"""
import math
import numpy
from numpy import pi


def secant(f, p0, p1, e, nmax, z0):
    i = 2
    p = []
    while i < nmax:
        f_p0 = f(p0, z0)
        f_p1 = f(p1, z0)

        # Công thức Secant
        p2 = p1 - f_p1 * (p1 - p0) / float(f_p1 - f_p0)

        if abs(f_p1) < e:  # Nếu giá trị gần nghiệm
            break

        # Kiểm tra điều kiện dừng dựa trên sai số
        if abs(p2 - p1) < e:
            break

        p0, p1 = p1, p2
        i += 1

    return p1
