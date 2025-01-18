from numpy import sqrt
import numpy as np
import csv


a, b = 1, 2
N = 1000
epsilon = 10 ** (-6)
p0 = 1.2


def g1(x):
    return sqrt(10 - x**3) / 2


def g2(x):
    return x - (x**3 + 4 * x**2 - 10)


def g3(x):
    return sqrt(10 / (4 + x))


def g4(x):
    return x - (x**3 + 4 * x**2 - 10) / (3 * x**2 + 8 * x)


def fixpoint(g, a, b, N, epsilon, p0):
    x = np.zeros(N)
    x[0] = p0

    if a <= p0 <= b:
        for i in range(N):
            if i + 1 < N:
                x[i + 1] = g(x[i])
                if abs(x[i + 1] - x[i]) < epsilon:
                    break

    return x


x1 = fixpoint(g1, a, b, N, epsilon, p0)
x2 = fixpoint(g2, a, b, N, epsilon, p0)
x3 = fixpoint(g3, a, b, N, epsilon, p0)
x4 = fixpoint(g4, a, b, N, epsilon, p0)


with open("fileFixpoint.txt", "w", newline="") as writefile:
    header = ["n", "x1", "x2", "x3", "x4"]
    writer = csv.DictWriter(writefile, fieldnames=header)
    writer.writeheader()
    for i in range(N):
        writer.writerow(
            {
                "n": i,
                "x1": x1[i],
                "x2": x2[i],
                "x3": x3[i],
                "x4": x4[i],
            }
        )
