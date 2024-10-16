from numpy import sin, exp, cos, tan, sqrt
import numpy as np
import matplotlib.pyplot as plt


def plotWaveFunction(z, z0, a, V0, N, hbar_Js, m, Ce):

    x = np.linspace(-2 * a, 2 * a, N)

    kappa = []
    l = []
    energy = []

    for i in range(len(z)):
        l.append(z[i] / a)  # công thức 2.155 slide
        kappa.append(z[i] / a * tan(z[i]))  # công thức 2.156 slide
        energy.append(((hbar_Js * z[i]) ** 2) / (2 * m * a**2) / Ce)  # công thức 2.173 Griffiths

    def V(x, V0, a):
        Vp = np.zeros(len(x))
        for i in range(len(x)):
            if x[i] <= -a or x[i] >= a:
                Vp[i] = V0

        return Vp

    def Ψ(x, V0, l, kappa, a):

        D = V0 / 60
        # D1 = sqrt(1 / (a + 1 / kappa))

        F = D * (cos(l * a)) / exp(-kappa * a)

        wave = np.zeros(len(x))

        for i in range(len(x)):
            if x[i] <= -a:
                wave[i] = F * exp(kappa * x[i])
            elif x[i] >= a:
                wave[i] = F * exp(-kappa * x[i])
            else:
                wave[i] = D * cos(l * x[i])
        return wave

    fig, ax = plt.subplots(1, 2, figsize=(15, 7))
    U = V(x, V0, a)
    ax[1].plot(x, U)
    ax[1].set_title("b)")
    ax[1].grid(True)
    for i in range(len(l)):
        ax[1].plot(x, energy[i] + Ψ(x, V0, l[i], kappa[i], a))
        ax[1].axhline(y=energy[i], ls="--", xmin=0.2, xmax=0.8)
        ax[1].annotate(f"E= {energy[i]:10f} eV", xy=(x[-1], energy[i]), xycoords="data")

    z1 = np.arange(1e-15, 4 * np.pi, 1 / N)
    rhs = tan(z1)
    lhs = sqrt((z0 / z1) ** 2 - 1)
    ax[0].plot(z1, rhs)
    ax[0].plot(z1, lhs)
    ax[0].set_ylim(-0.2, 10 * z0)
    ax[0].grid(True)
    ax[0].set_title("a)")

    plt.savefig("Final.pdf")

    plt.show()
