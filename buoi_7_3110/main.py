import numpy as np
from Jacobi import FArrMatrixJacobian
from GaussianSeidel import FArrMatrixGauss
from solveLoop import solveLoop
from math import sqrt
import matplotlib.pyplot as plt
import csv

np.set_printoptions(precision=10, linewidth=120, edgeitems=12)


def createMatrix(subDiagonal: int, dim: int):

    AMatrix = np.zeros([dim, dim])

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

    BMatrix = np.zeros(dim)

    for i in range(0, int(sqrt(dim))):
        BMatrix[i] = -100
        # BMatrix[-i] = 100

    return AMatrix, BMatrix


def MatrixFull(AMatrix, value):

    size = AMatrix.shape[0]

    AMatrixFull = np.zeros([size + 2, size + 2])

    for i in range(1, size + 1):
        for j in range(1, size + 1):
            AMatrixFull[i, j] = AMatrix[i - 1, j - 1]
    for i in range(0, size + 2):
        AMatrixFull[0, i] = value

    return AMatrixFull


def plotSol(grid, vjacobi, vgauss):

    fig = plt.figure()

    ax1 = fig.add_subplot(121, projection="3d")
    ax2 = fig.add_subplot(122, projection="3d")

    steps = grid + 2

    X = np.linspace(0, grid, steps)
    Y = np.linspace(0, grid, steps)
    # print(len(X))
    # print(len(vjacobi))

    X, Y = np.meshgrid(X, Y)
    ax1.plot_wireframe(X, Y, vjacobi)
    ax2.plot_wireframe(X, Y, vgauss)
    plt.show()


def saveLog(file: str, N, Jacobian, Gaussian):
    with open(file, "w", newline="") as writefile:
        header = [f"{"n":^4}", f"{"Jacobian":^18}", f"{"Gaussian":^18}"]
        writer = csv.DictWriter(writefile, fieldnames=header, delimiter="|")
        writer.writeheader()
        for i in range(0, N + 1):
            writer.writerow(
                {
                    f"{"n":^4}": f"{i:^4}",
                    f"{"Jacobian":^18}": f"{Jacobian[i]:^18}" if i < len(Jacobian) else f"{"":<18}",
                    f"{"Gaussian":^18}": f"{Gaussian[i]:^18}" if i < len(Gaussian) else f"{"":<18}",
                }
            )


def main():
    numberLoop = 10
    ic = 100
    unknowPoints = 12
    i, j = unknowPoints, unknowPoints
    dim = unknowPoints**2
    file = "ElectricPotentials.txt"

    AMatrix, BMatrix = createMatrix(unknowPoints, dim)
    eigenFunction = np.zeros(dim)
    print((AMatrix))
    print((BMatrix))
    ujacobi, loops = solveLoop(AMatrix, BMatrix, dim, FArrMatrixJacobian, numberLoop, eigenFunction)
    ugauss, loops = solveLoop(AMatrix, BMatrix, dim, FArrMatrixGauss, numberLoop, eigenFunction)
    # u = np.linalg.solve(AMatrix, BMatrix)
    print(ujacobi)

    vjacobi = np.reshape(ujacobi, [i, j])
    vgauss = np.reshape(ugauss, [i, j])

    vjacobi = MatrixFull(vjacobi, ic)
    vgauss = MatrixFull(vgauss, ic)
    print(vjacobi)

    plotSol(unknowPoints, vjacobi, vgauss)
    saveLog(file, numberLoop, ujacobi, ugauss)


if __name__ == "__main__":
    main()
