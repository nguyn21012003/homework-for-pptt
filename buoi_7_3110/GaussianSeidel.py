import numpy as np
from numpy import typing as npt


def FArrMatrixGauss(
    AMatrix: npt.NDArray, BMatrix: npt.NDArray, dim: int, x: npt.NDArray
) -> npt.NDArray:  ## Đưa mảng vào bằng ma trận sử dụng pp Gaussian-Seidel
    F = np.zeros(dim)

    for i in range(dim):
        sumAij_1 = 0
        sumAij_2 = 0
        for j in range(i):
            sumAij_1 += -AMatrix[i][j] * F[j]
        for j in range(i + 1, dim):
            sumAij_2 += -AMatrix[i][j] * x[j]

        F[i] = (BMatrix[i] + sumAij_1 + sumAij_2) / AMatrix[i][i]

    return F
