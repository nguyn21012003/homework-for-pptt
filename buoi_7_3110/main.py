import numpy as np
from mpmath import *
import matplotlib.pyplot as plt
from Jacobi import FArrMatrixJacobian
from solveLoop import solveLoop

dim = 35


def createMatrix():
    # or any other size you need

    # Initialize an n x n matrix with zeros
    AMatrix = np.zeros((dim, dim))
    BMatrix = np.zeros(dim)
    BMatrix[dim - 1] = -100
    # Use a for loop to set the values on the diagonals
    for i in range(dim):
        AMatrix[i, i] = -4  # Main diagonal
        if i > 0:
            AMatrix[i, i - 1] = 1  # Lower diagonal
        if i < dim - 1:
            AMatrix[i, i + 1] = 1  # Upper diagonal

    return AMatrix, BMatrix


def main():
    N = 100
    x = np.full(dim, 1)
    AMatrix, BMatrix = createMatrix()

    s = solveLoop(AMatrix, BMatrix, dim, FArrMatrixJacobian, N, x)

    a = np.linalg.solve(AMatrix, BMatrix)
    print(a)


if __name__ == "__main__":
    main()
