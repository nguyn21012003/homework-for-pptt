from numpy import sin, exp, cos, tan, sqrt
import numpy as np
import matplotlib.pyplot as plt
import random


def plotWaveFunction(z, a, V0, N):

    x = np.linspace(-2 * a, 2 * a, N)
    x1 = np.linspace(-2 * a, 2 * a, N)

    x2 = np.linspace(a, 2 * a, N)
    x3 = np.linspace(-a, -2 * a, N)

    def V(x, V0, a):
        Vp = np.zeros(len(x))
        for i in range(len(x)):
            if -a <= x[i] <= a:
                Vp[i] = -V0
            else:
                Vp[i] = 0
        return Vp

    def F(z):
        return exp(z * tan(z)) * cos(z) / sqrt((z * tan(z) + 1) / (z * tan(z) / a))

    def D(z):
        return 1 / sqrt((z * tan(z) + 1) / (z * tan(z) / a))

    def Ψ2(x):
        wave = []
        for i in range(len(z)):
            f = F(z)[i]
            Z = z[i]
            wave.append(f * exp(-Z * tan(Z) * x / a))
        return np.array(wave)

    def Ψ1(x):
        wave2 = []
        for i in range(len(z)):
            d = D(z)[i]
            Z = z[i]
            wave2.append(d * cos(Z * x / a))
        return np.array(wave2)

    fig, ax = plt.subplots()
    U = V(x, V0, a)
    ax.plot(x, U)
    # ax.set_xlabel(r"$x$", fontsize=18)
    # ax.set_ylabel(r"$U(x)$", fontsize=18)
    #
    for i in range(len(z)):
        cc = (random.random(), random.random(), random.random())
        ax.plot(x1, Ψ1(x1)[i], color=cc, label=str(i))
        ax.plot(x2, Ψ2(x2)[i], color=cc, label=str(i))
    plt.grid()
    # plt.legend()
    plt.show()
