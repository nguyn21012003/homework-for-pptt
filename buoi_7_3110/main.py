import numpy as np
from mpmath import *
from Jacobi import FArrMatrixJacobian
from solveLoop import solveLoop


from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np


def createMatrix(dim):

    n = int(dim / 40)

    AMatrix = np.zeros((dim, dim))
    BMatrix = np.zeros(dim)
    for i in range(n):
        BMatrix[n] = -100

    for i in range(dim):
        AMatrix[i, i] = -4

        if i - 1 >= 0:
            AMatrix[i, i - 1] = 1

        if i + 1 < dim:
            AMatrix[i, i + 1] = 1

        if i - n >= 0:
            AMatrix[i, i - n] = 1

        if i + n < dim:
            AMatrix[i, i + n] = 1

    return AMatrix, BMatrix


def plotSol(solution, x):
    pass


def main():
    N = 10000
    dim = 120
    x = np.full(dim, 1)
    AMatrix, BMatrix = createMatrix(dim)

    s, i = solveLoop(AMatrix, BMatrix, dim, FArrMatrixJacobian, N, x)

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    X, Y = np.linspace(0, 40, dim), np.linspace(0, 40, dim)

    ax.plot_wireframe(X, Y, s, rstride=10, cstride=10)
    plt.show()


if __name__ == "__main__":
    main()
