import numpy as np
from tqdm import tqdm
import subprocess
import csv
import matplotlib.pyplot as plt
from math import sqrt

L = 1
t = 0.12
rho = 0.01
T = 40
dL = 0.01
dt = 0.0001
c = sqrt(T / rho)
beta = c * dt / dL
# beta = 0.01
print(beta)


def FowardDiff(x: int, t: int, beta):
    U = []
    for _ in range(x + 1):
        row = []
        for _ in range(t + 1):
            row.append(0.0)
        U.append(row)
    for xi in range(1, x):
        U[xi][0] = 100.0

    for ti in range(1, t + 1):
        U[0][ti] = 0.0
        U[x][ti] = 0.0

    for j in range(t):
        for i in range(1, x):
            U[i][j + 1] = 2 * (1 - beta**2) * U[i][j] + beta**2 * (U[i + 1][j] + U[i - 1][j]) - U[i][j - 1]

    return U


def BackwardDiff(x: int, t: int, beta):
    U = []
    for _ in range(x + 1):
        row = []
        for _ in range(t + 1):
            row.append(0.0)
        U.append(row)
    for xi in range(1, x):
        U[xi][0] = 100.0

    for ti in range(1, t + 1):
        U[0][ti] = 0.0
        U[x][ti] = 0.0

    for j in range(t):
        for i in range(1, x):
            U[i][j + 1] = 2 * (1 - beta**2) * U[i][j] + beta**2 * (U[i + 1][j] + U[i - 1][j]) - U[i][j - 1]

    return U


def writeLog(x: int, t: int) -> None:
    """This function responsible for write out data of this code. I need the value of U with respectively x and t, so I also need the value of x and t."""
    file = "HyperbolicEquationData.txt"
    X = np.arange(0, x + dL, dL)
    T = np.arange(0, t + dt, dt)
    U = FowardDiff(len(X) - 1, len(T) - 1, beta)
    with open(file, "w", newline="") as writefile:
        header = [
            f"{'x':^4}",
            f"{'t':^8}",
            f"{'U':^10}",
        ]
        writer = csv.DictWriter(writefile, fieldnames=header, delimiter="\t")
        writer.writeheader()
        prev_xi = None
        for xi in range(len(X)):
            if prev_xi is not None and xi != prev_xi:
                writer.writerow({})
            for ti in range(len(T)):
                writer.writerow(
                    {
                        f"{'x':^4}": f"{xi:^4}",
                        f"{'t':^8}": f"{ti:^8}",
                        f"{'U':^10}": f"{U[xi][ti]:^10}",
                    }
                )
            prev_xi = xi

    Z = np.array(U)
    X, Y = np.meshgrid(T, X)
    fig = plt.figure(figsize=(15, 7))
    ax1 = fig.add_subplot(projection="3d")
    # ax1.plot_wireframe(X, Y, Z, color="purple")
    ax1.plot_surface(X, Y, Z, cmap="inferno")
    ax1.set_xlabel("time (s)")
    ax1.set_zlabel("T (K)")
    ax1.set_ylabel("x")
    ax1.legend([f"time = {t}s , L = {L}m, T = {T} N, rho = {rho}"])

    plt.savefig(f"HyperbolicEquationData & L = {L} & time = {t} & dL = {dL} & dt = {dt}.png")
    plt.show()

    return None


def gnuPlot(file: str) -> None:
    """
    Không xài dòng nào thì # vào dòng đó
    """
    with open("gnuPlot.gp", "w") as gnuplot:
        gnuplot.write(
            f"""

    #set multiplot layout 1,3


    set ylabel "y"
    set xlabel "x"
    set zlabel "z"


    set grid
    #set key horiz

    splot "{file}" u 1:2:3 with lines 

    #set datafile separator '\t'
    

    #unset multiplot

    pause -1

"""
        )
    subprocess.run(["gnuplot", "gnuPlot.gp"])
    return None


def main():

    writeLog(L, t)
    gnuPlot(file="HyperbolicEquationData.txt")


if __name__ == "__main__":
    main()
