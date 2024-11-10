import numpy as np
from numpy import typing as npt


def FArrMatrixJacobian(
    AMatrix: npt.NDArray, BMatrix: npt.NDArray, dim: int, eigenFunction: npt.NDArray
) -> npt.NDArray:  ## Đưa mảng vào bằng ma trận sử dụng pp Jacobian
    F = np.zeros(dim**2)

    for i in range(dim**2):
        sumAij = 0
        for j in range(dim**2):
            if j != i:
                sumAij += -AMatrix[i][j] * (eigenFunction[j])
        F[i] = 1 / AMatrix[i][i] * (sumAij + BMatrix[i])

    return F
