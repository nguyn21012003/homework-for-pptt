import numpy as np
from random import random, uniform
from math import cos, sqrt
import matplotlib.pyplot as plt
import time


def f(x):
    return x**2


def f1(x, y):
    return cos(x**4) + 3 * y**2


def monte(f, a, b, N, eps):
    count = 0
    S = np.zeros(N)
    S[0] = 0
    ans = 0
    list_x = np.zeros(N)
    list_y = np.zeros(N)
    for i in range(N):
        x = uniform(a, b)
        y = uniform(f(0), f(x))
        S[i] = f(x)
        count += 1
        ans += (b - a) * S[i] / N
        list_x[i] = x
        list_y[i] = y
    return ans, list_x, list_y


def riemann(f, x0, a, b, N, eps):
    count = 0
    h = (b - a) / N
    S = np.zeros(N)
    S[0] = 0
    ans = 0
    x = np.zeros(N)
    x[0] = x0
    for i in range(1, N):
        x[i] = a + i * h
        count += i
        S[i] += f(x[i])
        ans += (b - a) * S[i] / N
    return ans, x


def monte2D(f, a, b, c, d, N, eps):
    S = np.zeros(N)
    count = 0
    S[0] = 0
    ans = 0
    list_x = np.zeros(N)
    list_y = np.zeros(N)
    for i in range(N):
        x = uniform(a, b)
        y = uniform(c, d)
        S[i] += f(x, y)
        ans += (b - a) * (d - c) * S[i] / N
        list_x[i] += x
        list_y[i] += y
        count += i
    return ans, list_x, list_y


def riemann2D(f, x0, y0, a, b, c, d, N, eps):
    count = 0
    ans = 0
    hx = (b - a) / N
    hy = (d - c) / N
    S = np.zeros(N)
    S[0] = 0
    x = np.zeros(N)
    y = np.zeros(N)
    x[0] = x0
    y[0] = y0
    for i in range(1, N):
        x[i] = x0 + i * hx
        y[i] = y0 + i * hy
        count += i
        S[i] += f(x[i], y[i])
        ans += (b - a) * (d - c) * S[i] / N
    return ans, x, y


def plot(xmonte1D, ymonte1D, xrieman1D, N):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 7))

    x = np.linspace(0, 5, N)
    y = f(x)

    ax1.plot(x, y, color="r", label=r"y=$x^2$", lw=2)
    ax1.scatter(xmonte1D, ymonte1D, marker="o", s=10)
    ax2.bar(xrieman1D, ymonte1D)

    ax1.legend([r"y=$x^2$", "sỏi"])
    ax1.set_title(r"Monte-Carlo cho phương trình $y = x^2$")
    ax1.set_xlabel(r"x")
    ax1.set_ylabel(r"$f(x)$", rotation=0)

    # ax2

    fig.savefig("Monte-Carlo1D.pdf")

    plt.show()


def main():

    N = 10000
    eps = sqrt(N)

    #########################################
    a1, b1 = 0, 5
    # Câu 1: tính tích phân x = [0,5]  của s^2 và vẽ minh họa số sỏi ném trong miền tạo bởi f(x) và trục tung
    ans, list_x, list_y = monte(f, a1, b1, N, eps)
    #########################################
    # Câu 2:
    a, b, c, d = 4, 6, 0, 1
    x0 = a
    y0 = c
    ans1, list_x1, list_y1 = monte2D(f1, a, b, c, d, N, eps)
    ans2, list_x2 = riemann(f, x0, a, b, N, eps)
    ans3, list_x3, list_y3 = riemann2D(f1, x0, y0, a, b, c, d, N, eps)

    print(ans)
    plot(list_x, list_y, list_x2, N)


if __name__ == "__main__":
    main()
