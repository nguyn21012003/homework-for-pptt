import numpy as np
from mpmath import *
import matplotlib.pyplot as plt
from Jacobi import FArrMatrixJacobian
from solveLoop import solveLoop

AMatrix = [
    [0.25, 0, 0.25, 0],
    [0, 0.25, 0.25, 0],
    [0, 0.25, 0.25, 0],
    [0.25, 0, 0, 0.25],
]

BMatrix = [
    0,
    0,
    -0.25,
    -0.25,
]


def main():
    N = 1000
    dim = 4
    x = np.zeros(4)

    s = solveLoop(AMatrix, BMatrix, dim, FArrMatrixJacobian, N, x)
    print(s)


if __name__ == "__main__":
    main()
