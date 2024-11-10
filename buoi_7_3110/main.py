import numpy as np
from Jacobi import FArrMatrixJacobian
from GaussianSeidel import FArrMatrixGauss
from solveLoop import solveLoop
from math import sqrt
import matplotlib.pyplot as plt
import csv
from tqdm import tqdm

np.set_printoptions(precision=10, linewidth=120, edgeitems=12)


def createMatrix(dim: int):

    AMatrix = np.zeros([dim**2, dim**2])

    for i in range(1, dim - 1):
        for j in range(1, dim - 1):
            north = (i - 1) * dim + j
            west = i * dim + j - 1
            index = i * dim + j
            east = i * dim + j + 1
            south = (i + 1) * dim + j

            AMatrix[index, north] = 1
            AMatrix[index, west] = 1
            AMatrix[index, index] = -4
            AMatrix[index, east] = 1
            AMatrix[index, south] = 1

    i = 0
    for j in range(1, dim - 1):
        west = i * dim + j - 1
        index = i * dim + j
        east = i * dim + j + 1
        south = (i + 1) * dim + j

        AMatrix[index, west] = 1
        AMatrix[index, index] = -4
        AMatrix[index, east] = 1
        AMatrix[index, south] = 1

    j = 0
    for i in range(1, dim - 1):
        north = (i - 1) * dim + j
        index = i * dim + j
        east = i * dim + j + 1
        south = (i + 1) * dim + j

        AMatrix[index, north] = 1
        # AMatrix[index,west] =1
        AMatrix[index, index] = -4
        AMatrix[index, east] = 1
        AMatrix[index, south] = 1

    j = dim - 1
    for i in range(1, dim - 1):
        north = (i - 1) * dim + j
        west = i * dim + j - 1
        index = i * dim + j
        # east = i*n+j+1
        south = (i + 1) * dim + j

        AMatrix[index, north] = 1
        AMatrix[index, west] = 1
        AMatrix[index, index] = -4
        # a[index,east] =1
        AMatrix[index, south] = 1

    i = dim - 1
    for j in range(1, dim - 1):
        north = (i - 1) * dim + j
        west = i * dim + j - 1
        index = i * dim + j
        east = i * dim + j + 1
        # south = (i+1)*n+j

        AMatrix[index, north] = 1
        AMatrix[index, west] = 1
        AMatrix[index, index] = -4
        AMatrix[index, east] = 1

    i = 0
    j = 0
    index = i * dim + j
    east = i * dim + j + 1
    south = (i + 1) * dim + j

    AMatrix[index, index] = -4
    AMatrix[index, east] = 1
    AMatrix[index, south] = 1

    i = 0
    j = dim - 1
    west = i * dim + j - 1
    index = i * dim + j
    south = (i + 1) * dim + j

    AMatrix[index, west] = 1
    AMatrix[index, index] = -4
    AMatrix[index, south] = 1

    i = dim - 1
    j = 0
    north = (i - 1) * dim + j
    index = i * dim + j
    east = i * dim + j + 1

    AMatrix[index, north] = 1
    AMatrix[index, index] = -4
    AMatrix[index, east] = 1

    i = dim - 1
    j = dim - 1
    north = (i - 1) * dim + j
    west = i * dim + j - 1
    index = i * dim + j

    AMatrix[index, north] = 1
    AMatrix[index, west] = 1
    AMatrix[index, index] = -4

    # print(AMatrix, "\n")

    BMatrix = np.zeros([dim**2, 1])

    for i in range(0, dim):
        BMatrix[i] = -100
        # BMatrix[-i - 1] = 100

    return AMatrix, BMatrix


def MatrixFull(matrix, value):

    size = matrix.shape[0]

    matrixFull = np.zeros([size + 2, size + 2])

    for i in range(1, size + 1):
        for j in range(1, size + 1):
            matrixFull[i, j] = matrix[i - 1, j - 1]
    for i in range(0, size + 2):
        matrixFull[0, i] = value
        # matrixFull[size + 1, i] = -value

    return matrixFull


def u_values(d):
    u = np.zeros((d, d))
    for i in range(d):
        u[0][i] = 100
        # u[d - 1][i] = -100
    return u


def Gauss(d, N, u):
    for k in tqdm(range(N), desc="Processing", unit="step"):
        for i in range(1, d - 1):
            for j in range(1, d - 1):
                u[i][j] = 0.25 * (u[i + 1][j] + u[i - 1][j] + u[i][j + 1] + u[i][j - 1])
    return u


def plotSol(file, grid, vMatrixJacobi, vMatrixGauss, GaussWOMatrix, vGiaiTich):
    color = "purple"

    fig = plt.figure(figsize=(15, 7))

    ax1 = fig.add_subplot(221, projection="3d")
    ax2 = fig.add_subplot(222, projection="3d")
    ax3 = fig.add_subplot(223, projection="3d")
    ax4 = fig.add_subplot(224, projection="3d")

    steps = grid + 2

    X = np.linspace(0, grid, steps)
    Y = np.linspace(0, grid, steps)
    # print(len(X))
    # print(len(vjacobi))

    X, Y = np.meshgrid(X, Y)
    ax1.plot_wireframe(X, Y, vMatrixJacobi, color=color)
    ax1.set_title(f"Jacobian with matrĩ at {(steps-2)**2} unknow points")

    ax2.plot_wireframe(X, Y, vMatrixGauss, color=color)
    ax2.set_title(f"Gaussian-Seidel with matrix at {(steps-2)**2} unknow points")

    x = np.linspace(0, grid, grid)
    y = np.linspace(0, grid, grid)
    x, y = np.meshgrid(x, y)
    ax3.plot_wireframe(x, y, GaussWOMatrix, color=color)
    ax3.set_title(f"Gaussian-Seidel Iterative at {(steps-2)**2} unknow points")

    ax4.plot_wireframe(X, Y, vGiaiTich, color=color)
    ax4.set_title(f"Nghiem giai tich su dung matrix at {(steps-2)**2} unknow points")

    ax1.view_init(elev=20, azim=45)
    ax2.view_init(elev=20, azim=45)
    ax3.view_init(elev=20, azim=45)
    ax4.view_init(elev=20, azim=45)

    fig.savefig(file, format="png", orientation="landscape")
    plt.show()


def saveLog(file: str, N, Jacobian, Gaussian, i, j):
    with open(file, "w", newline="") as writefile:
        header = [f"{"n":^4}", f"{"i":^4}", f"{"j":^4}", f"{"Jacobian":^18}", f"{"Gaussian":^18}"]
        writer = csv.DictWriter(writefile, fieldnames=header, delimiter="|")
        writer.writeheader()
        for n in range(0, N + 1):
            writer.writerow(
                {
                    f"{"n":^4}": f"{n:^4}",
                    f"{"i":^4}": f"{i:^4}",
                    f"{"j":^4}": f"{j:^4}",
                    f"{"Jacobian":^18}": f"{Jacobian[n]:^18}" if n < len(Jacobian) else f"{"":<18}",
                    f"{"Gaussian":^18}": f"{Gaussian[n]:^18}" if n < len(Gaussian) else f"{"":<18}",
                }
            )


def main():
    numberLoop = 100
    ic = 100
    unknowPoints = 5
    i, j = unknowPoints, unknowPoints
    dim = unknowPoints
    fileSave = "ElectricPotentials.txt"
    filePlot = "ElectricPotentials.png"

    AMatrix, BMatrix = createMatrix(dim)
    eigenFunction = np.zeros([dim**2, 1])
    print(eigenFunction)
    print((AMatrix))
    print((BMatrix))
    ugauss, loops = solveLoop(AMatrix, BMatrix, dim, FArrMatrixGauss, numberLoop, eigenFunction)
    ujacobi, loops = solveLoop(AMatrix, BMatrix, dim, FArrMatrixJacobian, numberLoop, eigenFunction)
    uSolution = np.linalg.solve(AMatrix, BMatrix)
    print(ugauss)
    print(ujacobi)
    print(uSolution)
    #
    vgauss = np.reshape(ugauss, [i, j])
    vjacobi = np.reshape(ujacobi, [i, j])
    vSolution = np.reshape(uSolution, [i, j])

    vjacobi = MatrixFull(vjacobi, ic)
    vgauss = MatrixFull(vgauss, ic)
    vSolution = MatrixFull(vSolution, ic)
    print(vjacobi)

    UMatrix = u_values(unknowPoints)
    GaussWOMatrix = Gauss(unknowPoints, numberLoop, UMatrix)

    plotSol(filePlot, unknowPoints, vjacobi, vgauss, GaussWOMatrix, vSolution)
    saveLog(fileSave, numberLoop, ujacobi, ugauss, i, j)


if __name__ == "__main__":
    main()
