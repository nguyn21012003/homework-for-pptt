import numpy as np
from numpy import abs
import csv
from numpy import typing as npt


def FArrMatrixJacobian(AMatrix, BMatrix, dim: int, x: npt.NDArray) -> npt.NDArray:  ## Đưa mảng vào bằng ma trận sử dụng pp Jacobian
    F = np.zeros(dim)

    for i in range(dim):
        sumAij = 0
        for j in range(dim):
            if j != i:
                sumAij += -AMatrix[i][j] * x[j]

        F[i] = 1 / AMatrix[i][i] * (sumAij + BMatrix[i])

    return F
