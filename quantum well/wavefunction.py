from numpy import sin, exp, cos, tan, sqrt
import numpy as np
import matplotlib.pyplot as plt


def plotWaveFunction(z, a, V0, N, hbar_Js, m, Ce):

    x = np.linspace(-2 * a, 2 * a, N)

    x1 = np.linspace(-2 * a, 2 * a, N)
    x2 = np.linspace(a, 2 * a, N)
    x3 = np.linspace(-a, -2 * a, N)

    kappa = []
    l = []
    energy = []

    for i in range(len(z)):
        l.append(z[i] / a)  # công thức 2.155 slide
        kappa.append(z[i] / a * tan(z[i]))  # công thức 2.156 slide
        energy.append(-(20 * (hbar_Js) ** 2 * z[i] ** 2) / (m * a**2) / Ce)  # công thức 2.173 Griffiths

    def V(x, V0, a):
        Vp = np.zeros(len(x))
        for i in range(len(x)):
            if -a <= x[i] <= a:
                Vp[i] = -V0
            else:
                Vp[i] = 0
        return Vp

    def Ψ(x, V0, l, kappa, a):

        D1 = V0 / 20
        D = sqrt(1 / (a + 1 / kappa))

        F = D * (cos(l * a)) / exp(-kappa * a)
        print(D)

        wave = np.zeros(len(x))

        for i in range(len(x)):
            if x[i] <= -a:
                wave[i] = F * exp(kappa * x[i])
            elif x[i] >= a:
                wave[i] = F * exp(-kappa * x[i])
            else:
                wave[i] = D * cos(l * x[i])
        return wave

    fig, ax = plt.subplots()
    U = V(x, V0, a)
    ax.plot(x, U)
    print(Ψ(x, V0, l[i], kappa[i], a))
    for i in range(len(l)):
        ax.plot(x, energy[i] + Ψ(x, V0, l[i], kappa[i], a))
        ax.axhline(y=0, ls="-.")
        ax.axhline(y=energy[i], ls="--", xmin=0.1666, xmax=0.8333)
        ax.annotate(f"E= {energy[i]:10f} eV", xy=(x[-1], energy[i]), xycoords="data")
    # ax.set_xlabel(r"$x$", fontsize=18)
    # ax.set_ylabel(r"$U(x)$", fontsize=18)
    #
    plt.grid()
    # plt.legend()
    plt.show()
