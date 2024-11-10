import numpy as np
from numpy import typing as npt


def FArrMatrixGauss(
    AMatrix: npt.NDArray, BMatrix: npt.NDArray, dim: int, eigenFunction: npt.NDArray
) -> npt.NDArray:  ## Đưa mảng vào bằng ma trận sử dụng pp Gaussian-Seidel
    F = np.zeros(dim**2)

    for i in range(dim**2):
        sumAij_1 = 0
        sumAij_2 = 0
        for j in range(i):
            sumAij_1 += -AMatrix[i][j] * F[j]
        for j in range(i + 1, dim**2):
            sumAij_2 += -AMatrix[i][j] * eigenFunction[j]

        F[i] = (BMatrix[i] + sumAij_1 + sumAij_2) / AMatrix[i][i]

    return F
