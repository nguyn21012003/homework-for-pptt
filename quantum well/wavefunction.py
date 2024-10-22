from numpy import sin, exp, cos, tan, sqrt
import numpy as np
import matplotlib.pyplot as plt


def plotWaveFunction(z, z0, a, V0, N, hbar_Js, m, Ce):

    x = np.linspace(-2, 2, N)

    kappa = []
    l = []
    energy = []

    for i in range(len(z)):
        E = (z[i] * hbar_Js) ** 2 / (2 * m * a**2) / Ce - V0
        l.append(z[i])  # công thức 2.155 slide # thứ nguyên là m^-1
        kappa.append(z[i] * tan(z[i]))  # công thức 2.156 slide # thứ nguyên là m^-1
        energy.append(E)  # công thức 2.173 Griffiths

    print(energy)

    def V(x, V0, a):
        Vp = np.zeros(len(x))
        for i in range(len(x)):
            if -1 <= x[i] <= 1:
                Vp[i] = -V0  # = 1 ev * 1.6*10^-19

        return Vp

    def Ψ(x, V0, l, kappa, a):

        # D = V0 / 60
        D = sqrt(1 / (1 + 1 / kappa))

        F = D * (cos(l)) / exp(-kappa)

        wave = np.zeros(len(x))

        for i in range(len(x)):
            if x[i] <= -1:
                wave[i] = F * exp(kappa * x[i])
            elif x[i] >= 1:
                wave[i] = F * exp(-kappa * x[i])
            else:
                wave[i] = D * cos(l * x[i])
        return wave

    fig, ax = plt.subplots(1, 3, figsize=(15, 7))
    U = V(x, V0, a)
    ax[2].plot(x, U)
    ax[1].grid(True)

    for i in range(len(l)):
        ax[1].plot(x, energy[i] + Ψ(x, V0, l[i], kappa[i], a))
        ax[2].axhline(y=energy[i], ls="--", xmin=0.1, xmax=0.9)
        ax[2].annotate(f"E= {energy[i]:10f} eV", xy=(x[-1], energy[i]), xycoords="data")

    z1 = np.arange(1e-15, 4 * np.pi, 1 / N)
    rhs = tan(z1)
    lhs = sqrt((z0 / z1) ** 2 - 1)

    ax[0].plot(z1, rhs)
    ax[0].plot(z1, lhs)
    ax[0].set_ylim(-0.2, 10 * z0)
    ax[0].grid(True)

    ax[0].legend([r"tan(z)", r"$\sqrt{\dfrac{z_0^2}{z^2} - 1}$"], loc="upper right")
    ax[1].legend([r"$\psi_1(x)$", r"$\psi_2(x)$", r"$\psi_3(x)$"], loc="upper right")

    ax[0].set_title("a)", y=0, pad=-25, verticalalignment="top")
    ax[1].set_title("b)", y=0, pad=-25, verticalalignment="top")
    ax[2].set_title("c)", y=0, pad=-25, verticalalignment="top")

    plt.savefig("Final.pdf")

    plt.show()
