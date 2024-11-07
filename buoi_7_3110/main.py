import numpy as np
from Jacobi import FArrMatrixJacobian
from solveLoop import solveLoop


import matplotlib.pyplot as plt
import numpy as np

np.set_printoptions(precision=10, linewidth=120, edgeitems=12)


def createMatrix(subDiagonal, dim):

    AMatrix = np.zeros((dim, dim))

    BMatrix = np.full((dim, 1), -100)

    for i in range(dim):
        AMatrix[i, i] = -4

        if i - 1 >= 0:
            AMatrix[i, i - 1] = 1

        if i + 1 < dim:
            AMatrix[i, i + 1] = 1

        if i - subDiagonal >= 0:
            AMatrix[i, i - subDiagonal] = 1

        if i + subDiagonal < dim:
            AMatrix[i, i + subDiagonal] = 1

    return AMatrix, BMatrix


def plotSol(solution, x):
    pass


def main():
    N = 10000
    grid = 10
    m, n, subDiagonal = 2, 2, 3
    dim = m * subDiagonal + n + 1
    eigenFunction = np.full(dim, 1)
    AMatrix, BMatrix = createMatrix(subDiagonal, dim)

    print(AMatrix)
    print(BMatrix)
    u, step = solveLoop(AMatrix, BMatrix, dim, FArrMatrixJacobian, N, eigenFunction)
    print(u)

    v = np.zeros([m + 1, n + 1])
    v = np.reshape(u, (m + 1, n + 1))

    print(v)


if __name__ == "__main__":
    main()
