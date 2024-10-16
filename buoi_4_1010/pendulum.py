import numpy as np
import numpy.typing as npt
from numpy import sin, cos, pi
import matplotlib.pyplot as plt
import csv


def fArr1(t: float, y0: float):
    theta, omgega = y0
    g = 9.81
    L = 0.6
    F = np.zeros(2)
    F[0] = omgega
    F[1] = -(g / L) * sin(theta)

    return F


def fArr2(t: float, y0: float):
    theta, omgega = y0
    g = 9.81
    L = 0.6
    F = np.zeros(2)
    F[0] = omgega
    F[1] = -(g / L) * (theta)

    return F


def rk4(fArr: npt.NDArray, tn: npt.NDArray, yn: npt.NDArray, h: float) -> npt.NDArray:

    k1 = fArr(tn, yn)
    k2 = fArr(tn + 0.5 * h, yn + 0.5 * h * k1)
    k3 = fArr(tn + 0.5 * h, yn + 0.5 * h * k2)
    k4 = fArr(tn + h, yn + h * k3)
    return yn + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def solve_sys_ode(fArr1, fArr2: npt.NDArray, a: float, b: float, h: float, solver: npt.NDArray, N: int) -> npt.NDArray:
    big = []
    small = []
    yn = np.zeros(2)
    yn1 = np.zeros(2)
    yn[0] = pi / 6
    yn1[0] = pi / 6
    yn[1] = 0
    yn1[1] = 0
    t = np.arange(a, b + h, h)
    for j in range(N + 1):
        yn = solver(fArr1, t, yn, h)
        big.append(yn)

    for i in range(N + 1):
        yn1 = solver(fArr2, t, yn1, h)
        small.append(yn1)
    return big, small


def plotTheta(t, big_theta, small_theta):
    fig, axs = plt.subplots(figsize=(16, 9))

    plt.rcParams["text.usetex"] = True

    axs.plot(t, big_theta, t, small_theta)

    axs.legend(
        [r"$\displaystyle\frac{d\theta^2}{dt^2} + \frac{g}{L} \sin\theta$", r"$\displaystyle\frac{d\theta^2}{dt^2} + \frac{g}{L}\theta$"],
        loc="best",
    )
    plt.savefig("pendulumTheta.png")

    plt.show()


def saveLog(file, t, big_theta, small_theta):
    with open(file, "w", newline="", encoding="utf8") as writefile:
        header = [
            f"{'t':^4}",
            f"{'sin θ':^25}",
            f"{'θ':^25}",
            f"{'sin θ - θ':^25}",
        ]
        writer = csv.DictWriter(writefile, fieldnames=header, delimiter="|")
        writer.writeheader()
        for i in range(len(t)):
            writer.writerow(
                {
                    f"{'t':^4}": f"{i:^4}",
                    f"{'sin θ':^25}": f"{float(big_theta[i]):^25}",
                    f"{'θ':^25}": f"{float(small_theta[i]):^25}",
                    f"{'sin θ - θ':^25}": f"{float(big_theta[i] - small_theta[i]):^25}",
                }
            )


def main():

    N = 500
    t0 = 0
    tn = 10

    file = "data.txt"

    h = (tn - t0) / N
    t = np.arange(t0, tn + h, h)

    big, small = solve_sys_ode(fArr1, fArr2, t0, tn, h, rk4, N)
    big_theta = []
    small_theta = []

    for i in range(len(big)):
        big_theta.append(big[i][0])
        small_theta.append(small[i][0])

    saveLog(file, t, big_theta, small_theta)
    plotTheta(t, big_theta, small_theta)


if __name__ == "__main__":
    main()
