import numpy as np
import subprocess
import csv
import matplotlib.pyplot as plt

np.set_printoptions(precision=10, linewidth=1200, edgeitems=24)


def FowardDiff(x, t, eta):

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
            U[i][j + 1] = (1 - 2 * eta) * U[i][j] + eta * (U[i + 1][j] + U[i - 1][j])

    return np.array(U)


def writeLog(x, t, U):
    file = "heatEquationData.txt"
    with open(file, "w", newline="") as writefile:
        header = [
            f"{'x':^4}",
            f"{'t':^8}",
            f"{'U':^10}",
        ]
        writer = csv.DictWriter(writefile, fieldnames=header, delimiter="\t")
        writer.writeheader()
        prev_xi = None
        for xi in range(x + 1):
            if prev_xi is not None and xi != prev_xi:
                writer.writerow({})
            for ti in range(t + 1):
                writer.writerow(
                    {
                        f"{'x':^4}": f"{xi:^4}",
                        f"{'t':^8}": f"{ti:^8}",
                        f"{'U':^10}": f"{U[xi][ti]:^10}",
                    }
                )
            prev_xi = xi
    print(U)

    fig = plt.figure(figsize=(15, 7))
    X, T = np.meshgrid(x, t)
    ax1 = fig.add_subplot(projection="3d")
    ax1.plot_surface(X, T, U, cmap="viridis")
    plt.show()


def gnuPlot(file):
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

    set datafile separator '\t'

    #set ylabel "y"
    #set xlabel "x"
    #set zlabel "z"
    #splot "{file}" u 1:2:4 with lines tit "p_t"


    #set xlabel "Omega"
    #set ylabel "alpha(omega)"
    #plot "{file}" u 2:3 with lines tit "Phân cực toàn phần"
    ##plot "{file}" u 2:4 with lines tit "Mật độ toàn phần"
    #plot "{file}" u 2:5 with lines tit "Phổ hấp thụ"


    #unset multiplot

    pause -1

"""
        )
    subprocess.run(["gnuplot", "gnuPlot.gp"])


def main():
    L = 10  ### Mét
    t = 300

    kappa = 210
    Cv = 900
    rho = 2700
    dL = 0.01
    dt = 0.3
    eta = kappa / (Cv * rho) * dt / (dL**2)

    U = FowardDiff(L, t, eta)
    writeLog(L, t, U)
    gnuPlot(file="heatEquationData.txt")


if __name__ == "__main__":
    main()
