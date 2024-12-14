import numpy as np
from numpy.random import rand
from numpy import abs
import matplotlib.pylab as plt


def lostFucntion(N, x, y):
    L = 1 / (2 * N) * abs(y - x)
    return L


def grad(N, x, y, w):
    gradL = 1 / N * x.T.dot(x.dot(w) - y)
    return gradL


def Solution(N, eps, gamma, x, y):

    w1 = np.zeros(1, N)
    vt1 = 0
    w = 0
    b = np.dot(x.T, y)
    A = np.dot(x.T, x)
    for i in range(0, N):
        vt1 = grad(N, x, y, w)
        w1.append(float(w1[i] - vt1))

        if abs(w1[i] - w1[i - 1]) <= eps:

            it = i
            break

    return w1


def generator(N):
    x = rand(N, 1)
    y = 3 + 4 * x + 0.1 * rand(N, 1)


def main():
    N = 1000
    eps = 1e-15
    gamma = 0.9
    x, y = generator(N)
    Solution(N, eps, gamma, x, y)
    plt.plot(x, y, "o", markersize="1")
    plt.show()


if __name__ == "__main__":
    main()
