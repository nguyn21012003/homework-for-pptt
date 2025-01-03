import numpy as np
from numpy import abs
import matplotlib.pyplot as plt


def Fx(x):
    return x**2 + 10 * np.sin(x)


def g(gamma, eta, vt, x):
    x = x - vt

    vt = gamma * vt + eta * (2 * x + 10 * np.cos(x))
    return vt


def Solution(N: int, eps, gamma, eta):

    x1 = [5]
    x2 = [-5]

    it = 0
    vt1 = 0
    vt2 = 0
    for i in range(0, N):
        vt1 = g(gamma, eta, vt1, x1[i])
        vt2 = g(gamma, eta, vt2, x2[i])
        x1.append(float(x1[i] - vt1))
        x2.append(float(x2[i] - vt2))

        if abs(x1[i] - x1[i - 1]) <= eps:
            it = i

            break
    return x1, x2


def main():
    N = 10000
    eps = 1e-5
    gamma = 0.9
    eta = 0.1
    x1, x2 = Solution(N, eps, gamma, eta)
    Farr1 = []
    Farr2 = []

    print((x2))
    for i in range(len(x1)):
        Farr1.append(Fx(x1[i]))
        Farr2.append(Fx(x2[i]))

    plt.plot(x1, Farr1)
    plt.plot(x2, Farr2)
    plt.show()


if __name__ == "__main__":
    main()
