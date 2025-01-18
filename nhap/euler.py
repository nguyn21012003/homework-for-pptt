import numpy as np
import csv

y0 = 1 / 2
x0 = 0
xmax = 2
N = 20
h = (xmax - x0) / N


def f1(x, y):  #### tuong duong voi y'
    return y - x**2 + 1


### Thuật toán
def euler(f, x, y, h):
    y = y + h * f(x, y)
    return y


def rk4(f, x, y, h):
    k1 = f(x, y)
    k2 = f(x + 0.5 * h, y + 0.5 * h * k1)
    k3 = f(x + 0.5 * h, y + 0.5 * h * k2)
    k4 = f(x + h, y + h * k3)
    return y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def solve(f, x0, xmax, algorithm):
    x = np.arange(x0, xmax + h, h)
    y = np.zeros(len(x))
    y[0] = y0
    for i in range(len(x) - 1):
        y[i + 1] = algorithm(f, x[i], y[i], h)
        x0 += h

    return x, y


x1, y1 = solve(f1, x0, xmax, euler)
x2, y2 = solve(f1, x0, xmax, rk4)


with open("fileGhiODE.txt", "w", newline="") as file:
    header = ["x", "euler", "rk4"]
    writer = csv.DictWriter(file, fieldnames=header)
    writer.writeheader()
    for i in range(len(x1) - 1):
        writer.writerow(
            {
                "x": x1[i],
                "euler": y1[i],
                "rk4": y2[i],
            }
        )
