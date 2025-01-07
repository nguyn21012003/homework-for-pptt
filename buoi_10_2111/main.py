import numpy as np
import matplotlib.pyplot as plt
from numpy import sqrt, sin, pi, exp
from tqdm import tqdm


def VComputed(x, mu, sigma):

    V = -1e4 * exp(-((x - mu) ** 2) / (2 * sigma**2))

    return V


def plotSomething(x, t, args):
    Z = args
    X, Y = np.meshgrid(x, t)

    fig1 = plt.figure(figsize=(16, 10))
    ax1 = fig1.add_subplot(121, projection="3d")
    ax1.plot(Y, X, Z)
    ax2 = fig1.add_subplot(122)
    ax2.plot(x, args)

    plt.show()
    return None


def PsiArr(Nt, Nx, x):
    psi = np.zeros([Nt, Nx], dtype=complex)
    psi0 = sqrt(2) * sin(pi * x)
    psi[0] = psi0
    return psi


def solveEquation(Nt, Nx, dt, dx, psi, V):
    for t in tqdm(range(0, Nt - 1)):
        for i in range(1, Nx - 1):
            psi[t + 1][i] = psi[t][i] + 1j / 2 * dt / dx**2 * (psi[t][i + 1] - 2 * psi[t][i] + psi[t][i - 1]) - 1j * dt * V[i] * psi[t][i]

        norm = 0
        for i in range(1, Nx - 1):
            norm += (np.absolute(psi[t + 1][i]) ** 2) * dx

        for i in range(1, Nx - 1):
            psi[t + 1][i] = psi[t + 1][i] / norm

    return psi


def main():
    # Nx = 301
    # Nt = 5001
    dx = 0.01
    dt = 0.001
    L = 1  # met
    tmax = 3  # s
    # x = np.linspace(0, L, Nx)
    # t = np.linspace(0, tmax, Nt)
    x = np.arange(0, L + dx, dx)
    t = np.arange(0, tmax + dt, dt)

    mu = L / 2
    sigma = L / 20

    V_pot = VComputed(x, mu, sigma)
    # plotSomething(x, V_pot)
    psi = PsiArr(Nt, Nx, x)
    psiComputed = solveEquation(Nt, Nx, dt, dx, psi, V_pot)
    plotSomething(x, t, np.absolute(psiComputed[5000]) ** 2)


if __name__ == "__main__":
    main()
