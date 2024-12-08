import numpy as np
from tqdm import tqdm
import csv
import matplotlib.pyplot as plt
from math import sqrt

L = 1
t = 1
rho = 1
Tension = 40
kappa = 1 * 1e-2 * 0
dL = 0.01
dt = 0.001
c = sqrt(Tension / rho)
beta = c * dt / dL
# beta = 0.01
print(beta)
frictionConstant = 2 * kappa / rho
### Thông số tham khảo có thể là L = 0.7m, dL = 0.7 mm , c = 300 m/s , dt = 5*1e-6 s, kappa = 2.6*1e-5 s/m


def FowardDiff(x: int, t: int, beta):

    U = []
    for _ in range(x + 1):
        row = []
        for _ in range(t + 1):
            row.append(0.0)
        U.append(row)

    for xi in range(0, x):
        if x <= 0.8 * L:
            U[xi][0] = 1.25 * xi / L
        elif x > 0.8 * L:
            U[xi][0] = (5 - 5 * xi) / L

    for ti in range(1, t + 1):
        U[0][ti] = 0.0
        U[x][ti] = 0.0
    if beta <= 1:
        for j in range(0, t):
            for i in range(1, x):
                U[i][j + 1] = (
                    1
                    / (1 / (c**2 * dt**2) + 2 * kappa / (dt * rho))
                    * (
                        (U[i + 1][j] - 2 * U[i][j] + U[i - 1][j]) / dL**2
                        - (U[i][j - 1] - 2 * U[i][j]) / (c**2 * dt**2)
                        + 2 * kappa / (dt * rho) * U[i][j - 1]
                    )
                )

    return U


def writeLog(x: int, t: float) -> None:
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
    ax1.plot_wireframe(X, Y, Z, color="purple")
    # ax1.plot_surface(X, Y, Z, cmap="inferno")
    ax1.set_xlabel("time (s)")
    ax1.set_zlabel("Biên độ")
    ax1.set_ylabel("x")
    ax1.legend([f"time = {t}s , L = {L}m, tension = {Tension} N, rho = {rho}"])

    plt.savefig(
        f"HyperbolicEquationData & L = {L} & time = {t} & dL = {dL} & dt = {dt}.png"
    )
    plt.show()

    return None


def main():

    writeLog(L, t)


if __name__ == "__main__":
    main()
