import numpy as np
from numpy import abs
import csv


AConst = [
    [10, -1, 2, 0],
    [-1, 11, -2, 3],
    [2, -1, 10, -1],
    [0, 3, -1, 8],
]

BConst = [6, 25, -11, 15]


def FArr(dim, x):  ## Đưa mảng vào bằng tay
    F = np.zeros(4)
    F[0] = 1 / 10 * x[1] - 1 / 5 * x[2] + 6 / 10
    F[1] = 1 / 11 * x[0] + 2 / 11 * x[2] - 3 / 11 * x[3] + 25 / 11
    F[2] = -2 / 10 * x[0] + 1 / 10 * x[1] + 1 / 10 * x[3] - 11 / 10
    F[3] = -3 / 8 * x[1] + 1 / 8 * x[2] + 15 / 8

    return F


def FArrMatrixJacobian(dim, x):  ## Đưa mảng vào bằng ma trận sử dụng pp Jacobian
    F = np.zeros(dim)

    for i in range(dim):
        sumAij = 0
        for j in range(dim):
            if j != i:
                sumAij += -AConst[i][j] * x[j]

        F[i] = 1 / AConst[i][i] * (sumAij + BConst[i])

    return F


def FArrMatrixGauss(dim, x):  ## Đưa mảng vào bằng ma trận sử dụng pp Gaussian-Seidel
    F = np.zeros(dim)

    for i in range(dim):

        F[i] = (BConst[i] - sum([AConst[i][j] * x[j] for j in range(i + 1, dim)]) - sum([AConst[i][j] * F[j] for j in range(i)])) / AConst[i][i]

    return F


def Jacobian(dim, F, N):
    x = np.zeros(4)
    listX = [x]  ### Lưu giá trị vào mảng để kiểm soát
    for i in range(1, N):
        x = F(dim, x)
        listX.append(x)
        print(x, F.__name__, i, "\n")

        if abs(max(listX[i]) - max(listX[i - 1])) / max(listX[i]) <= 1e-3:

            break

    return listX, i


def saveLog(file, N, xByhand, xMatrix, xGauss):

    with open(file, "w", newline="") as writefile:
        header = [f"{"k":^11}", f"{"Jacobian By Hand":^49}", f"{"Jacobian Matrix":^49}", f"{"Gaussian Matrix":^49}"]
        writer = csv.DictWriter(writefile, fieldnames=header)
        writer.writeheader()

        for i in range(1, N + 1):

            writer.writerow(
                {
                    f"{"k":^11}": f"{i:^11}",
                    f"{"Jacobian By Hand":^49}": xByhand[i],
                    f"{"Jacobian Matrix":^49}": xMatrix[i],
                    f"{"Gaussian Matrix":^49}": xGauss[i],
                }
            )


def main():
    dim = 4
    fileLog = "Jacobi&GaussSeidel.txt"

    xByhand, i = Jacobian(dim, FArr, 40)
    xMatrix, i = Jacobian(dim, FArrMatrixJacobian, 40)
    xGauss, i = Jacobian(dim, FArrMatrixGauss, 40)
    saveLog(fileLog, 8, xByhand, xMatrix, xGauss)


if __name__ == "__main__":
    main()
