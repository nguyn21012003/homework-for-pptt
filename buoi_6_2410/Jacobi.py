import numpy as np
from numpy import abs
import csv
from numpy import typing as npt

AConst = [
    [10, -1, 2, 0],
    [-1, 11, -2, 3],
    [2, -1, 10, -1],
    [0, 3, -1, 8],
]

BConst = [6, 25, -11, 15]


def FArr(dim: int, x: npt.NDArray) -> npt.NDArray:  ## Đưa mảng vào bằng tay
    F = np.zeros(dim)
    F[0] = 1 / 10 * x[1] - 1 / 5 * x[2] + 6 / 10
    F[1] = 1 / 11 * x[0] + 2 / 11 * x[2] - 3 / 11 * x[3] + 25 / 11
    F[2] = -2 / 10 * x[0] + 1 / 10 * x[1] + 1 / 10 * x[3] - 11 / 10
    F[3] = -3 / 8 * x[1] + 1 / 8 * x[2] + 15 / 8

    return F


def FArrMatrixJacobian(dim: int, x: npt.NDArray) -> npt.NDArray:  ## Đưa mảng vào bằng ma trận sử dụng pp Jacobian
    F = np.zeros(dim)

    for i in range(dim):
        sumAij = 0
        for j in range(dim):
            if j != i:
                sumAij += -AConst[i][j] * x[j]

        F[i] = 1 / AConst[i][i] * (sumAij + BConst[i])

    return F


def FArrMatrixGauss(dim: int, x: npt.NDArray) -> npt.NDArray:  ## Đưa mảng vào bằng ma trận sử dụng pp Gaussian-Seidel
    F = np.zeros(dim)

    for i in range(dim):
        sumAij_1 = 0
        sumAij_2 = 0
        for j in range(i):
            sumAij_1 += -AConst[i][j] * F[j]
        for j in range(i + 1, dim):
            sumAij_2 += -AConst[i][j] * x[j]

        F[i] = (BConst[i] + sumAij_1 + sumAij_2) / AConst[i][i]

    return F


def solveLoop(dim: int, F: npt.NDArray, N: int) -> list:
    x = np.zeros(4)
    listX = [x]  ### Lưu giá trị vào mảng để kiểm soát
    for i in range(1, N):
        x = F(dim, x)
        listX.append(x)
        print(x, F.__name__, i, "\n")

        if abs(max(listX[i]) - max(listX[i - 1])) / max(listX[i]) <= 1e-3:

            break

    return listX, i


def saveLog(file: str, N: int, xByhand: list, xMatrix: list, xGauss: list):

    with open(file, "w", newline="") as writefile:
        header = [f"{"k":^11}", f"{"Jacobian By Hand":^49}", f"{"Jacobian Matrix":^49}", f"{"Gaussian Matrix":^49}"]
        writer = csv.DictWriter(writefile, fieldnames=header, delimiter="|")
        writer.writeheader()

        for i in range(1, N + 1):

            writer.writerow(
                {
                    f"{"k":^11}": f"{i:^11}",
                    f"{"Jacobian By Hand":^49}": xByhand[i] if i < len(xByhand) else f"{"":<49}",
                    f"{"Jacobian Matrix":^49}": xMatrix[i] if i < len(xMatrix) else f"{"":<49}",
                    f"{"Gaussian Matrix":^49}": xGauss[i] if i < len(xGauss) else f"{"":<49}",
                }
            )


def main():
    dim = 4
    N = 40
    fileLog = "Jacobi&GaussSeidel.txt"

    xByhand, i = solveLoop(dim, FArr, N)
    xMatrix, i = solveLoop(dim, FArrMatrixJacobian, N)
    xGauss, i = solveLoop(dim, FArrMatrixGauss, N)
    saveLog(fileLog, N, xByhand, xMatrix, xGauss)


if __name__ == "__main__":
    main()
