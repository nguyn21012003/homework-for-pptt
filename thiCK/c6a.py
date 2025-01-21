import numpy as np

from math import sqrt
import matplotlib.pyplot as plt
import csv
from tqdm import tqdm
import subprocess


np.set_printoptions(precision=10, linewidth=120, edgeitems=15)


def FArrMatrixJacobian(AMatrix, BMatrix, dim, eigenFunction):  ## Đưa mảng vào bằng ma trận sử dụng pp Jacobian
    F = np.zeros(dim**2)

    for i in range(dim**2):
        sumAij = 0
        for j in range(dim**2):
            if j != i:
                sumAij += -AMatrix[i][j] * (eigenFunction[j])
        F[i] = 1 / AMatrix[i][i] * (sumAij + BMatrix[i])

    return F


def solveLoop(AMatrix, BMatrix, dim, F, N, initInput):

    eigenFunction = initInput
    listEigenFunction = [eigenFunction]  ### Lưu giá trị vào mảng để kiểm soát
    for i in tqdm(range(1, N + 1), desc=f"{F.__qualname__}", unit="step"):

        eigenFunction = F(AMatrix, BMatrix, dim, eigenFunction)
        listEigenFunction.append(eigenFunction)
        # print(eigenFunction, F.__qualname__, f"tại k = {i}", "\n")

    # s = np.array(listEigenFunction)

    return eigenFunction, i


def FArrMatrixGauss(AMatrix, BMatrix, dim, eigenFunction):  ## Đưa mảng vào bằng ma trận sử dụng pp Gaussian-Seidel
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


def createMatrix(dim):

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
    # BMatrix[int(len(BMatrix) / 2)][0] = -10

    for i in range(0, dim):
        BMatrix[i] = -100
        BMatrix[-i - 1] = 100

    return AMatrix, BMatrix


def MatrixFull(matrix, value):

    size = matrix.shape[0]

    matrixFull = np.zeros([size + 2, size + 2])

    for i in range(1, size + 1):
        for j in range(1, size + 1):
            matrixFull[i, j] = matrix[i - 1, j - 1]
    for i in range(0, size + 2):
        matrixFull[0, i] = value
        matrixFull[size + 1, i] = -value

    return matrixFull


def u_values(d):
    u = np.zeros((d, d))
    for i in range(d):
        u[0][i] = 100
        # u[d - 1][i] = -100
    return u


def Gauss(d, N, u):
    for k in tqdm(range(N), desc="Gauss", unit="step"):
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
    ax1.set_title(f"Jacobian with matrix at {(steps-2)**2} unknow points")

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
        header = [f"{'n':^4}", f"{'i':^4}", f"{'j':^4}", f"{'Jacobian':^18}", f"{'Gaussian':^18}"]
        writer = csv.DictWriter(writefile, fieldnames=header, delimiter="\t")
        writer.writeheader()

        n = 0
        for ith in range(0, i + 2):
            for jth in range(0, j + 2):

                writer.writerow(
                    {
                        f"{'n':^4}": f"{n:^4}",
                        f"{'i':^4}": f"{ith:^4}",
                        f"{'j':^4}": f"{jth:^4}",
                        f"{'Jacobian':^18}": f"{Jacobian[ith][jth]:^18}",
                        f"{'Gaussian':^18}": f"{Gaussian[ith][jth]:^18}",
                    }
                )
                # print(ith, jth)
                n += 1
            writer.writerow({})


def gnuPlot():
    with open("gnuPlot.gp", "w") as gnuplot:
        gnuplot.write(
            """
    set ylabel "y"
    set xlabel "x"
    set zlabel "V"


    set multiplot layout 1,2
    set grid
    # set key horiz

    splot "ElectricPotentials.txt" u 2:3:4 with lines tit "Jacobian"
    
    
    set ylabel "y"
    set xlabel "x"
    set zlabel "V"
    splot "ElectricPotentials.txt" u 2:3:5 with lines tit "Gaussian"


    unset multiplot

    pause -1
"""
        )
    subprocess.run(["gnuplot", "gnuPlot.gp"])


def main():
    numberLoop = 100
    ic = 100
    unknowPoints = 4
    i, j = unknowPoints, unknowPoints
    dim = unknowPoints
    fileSave = "ElectricPotentials.txt"
    filePlot = "ElectricPotentials.png"

    AMatrix, BMatrix = createMatrix(dim)
    eigenFunction = np.zeros([dim**2, 1])

    ugauss, loops = solveLoop(AMatrix, BMatrix, dim, FArrMatrixGauss, numberLoop, eigenFunction)
    ujacobi, loops = solveLoop(AMatrix, BMatrix, dim, FArrMatrixJacobian, numberLoop, eigenFunction)
    uSolution = np.linalg.solve(AMatrix, BMatrix)

    vgauss = np.reshape(ugauss, [i, j])
    vjacobi = np.reshape(ujacobi, [i, j])
    vSolution = np.reshape(uSolution, [i, j])

    vjacobi = MatrixFull(vjacobi, ic)
    vgauss = MatrixFull(vgauss, ic)
    vSolution = MatrixFull(vSolution, ic)
    print(vjacobi)
    # (vjacobi[1][2])

    UMatrix = u_values(unknowPoints)
    GaussWOMatrix = Gauss(unknowPoints, numberLoop, UMatrix)

    plotSol(filePlot, unknowPoints, vjacobi, vgauss, GaussWOMatrix, vSolution)
    saveLog(fileSave, numberLoop, vjacobi, vgauss, i, j)
    gnuPlot()


if __name__ == "__main__":
    main()
