import numpy as np
from numpy import sqrt
from Jacobi import FArrMatrixGauss, FArrMatrixJacobian, solveLoop, saveLog


AMatrix = [
    [-1, 0, 0, 1 / sqrt(2), 1, 0, 0, 0],
    [0, -1, 0, 1 / sqrt(2), 0, 0, 0, 0],
    [0, 0, -1, 0, 0, 0, 1 / 2, 0],
    [0, 0, 0, -1 / sqrt(2), 0, sqrt(3) / 2, 0, 0],
    [0, 0, 0, 0, -1, 0, 0, 1],
    [0, 0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, -1 / sqrt(2), 0, -1, -1 / 2, 0],
    [0, 0, 0, 0, 0, 0, -sqrt(3) / 2, -1],
]
BMatrix = [0, 0, 0, 0, 0, 10000, 0, 0]


def main():
    dim = 8
    N = 100
    x = np.full(dim, 1)

    xByhand = []
    xMatrix, i = solveLoop(AMatrix, BMatrix, dim, FArrMatrixJacobian, N, x)
    xGauss, i = solveLoop(AMatrix, BMatrix, dim, FArrMatrixGauss, N, x)
    s = np.linalg.solve(AMatrix, BMatrix)
    print(s)


if __name__ == "__main__":
    main()
