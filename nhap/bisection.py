import numpy as np
from numpy import exp,cos
import csv


def f(x):
    f = exp(x) - 2 - cos(exp(x)-2)

    return f


def bisection(f, a, b):
    c = (a + b) / 2
    N = 1000
    epsilon = 10**-6
    listc = np.zeros(N)
    lista = np.zeros(N)
    listb = np.zeros(N)
    lista[0] = a
    listb[0] = b

    for i in range(N):

        if i + 1 < N:
            listc[i] = (lista[i] + listb[i]) / 2
            fa = f(lista[i])
            fb = f(listb[i])
            fc = f(listc[i])
            if fa * fc < 0:
                lista[i + 1] = lista[i]
                listb[i + 1] = listc[i]
            elif fb * fc < 0:
                lista[i + 1] = listc[i]
                listb[i + 1] = listb[i]
            if abs(fc) < epsilon:
                break
    with open("data.txt", "w", newline="") as writefile:
        header = ["a", "b", "c"]
        writer = csv.DictWriter(writefile, fieldnames=header)
        writer.writeheader()
        for i in range(N):
            writer.writerow(
                {
                    "a": lista[i],
                    "b": listb[i],
                    "c": listc[i],
                }
            )


a = 1
b = 5
bisection(f, a, b)
