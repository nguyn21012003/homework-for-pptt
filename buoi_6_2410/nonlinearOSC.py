import numpy as np
from numpy import typing as npt, pi, sqrt, sin
import matplotlib.pyplot as plt
import csv


def FArr(t, initInput):
    x, velocity = initInput

    alphax = 0.01
    k = 2
    m = 1
    omega0 = sqrt(k / m)
    b = 2 * m * omega0
    F = np.zeros(2)

    Fext = 15 * sin(omega0 * t)
    Fviscous = -b * velocity

    F[0] = velocity
    F[1] = -(k * x - k * x * alphax) / m + Fext + Fviscous

    return F


def rk4(fArr: npt.NDArray, tn: npt.NDArray, yn: npt.NDArray, h: float) -> npt.NDArray:

    k1 = fArr(tn, yn)
    k2 = fArr(tn + 0.5 * h, yn + 0.5 * h * k1)
    k3 = fArr(tn + 0.5 * h, yn + 0.5 * h * k2)
    k4 = fArr(tn + h, yn + h * k3)
    return yn + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def savelog(file, N, solution):
    with open(file, "w", newline="") as writefile:
        header = ["n", "x", "velocity"]
        writer = csv.DictWriter(writefile, fieldnames=header)
        writer.writeheader()
        for i in range(N):
            writer.writerow({"n": i})


def solve(N, t0, t1, h, solver):

    initInput = np.zeros(2)
    initInput[0] = 1
    initInput[1] = 0

    t = np.linspace(t0, t1, N)
    solution = []
    # listX = []

    for i in range(N):
        initInput = solver(FArr, t[i], initInput, h)
        solution.append(initInput[1])
        # listX.append(initInput[0])
    plt.plot(t, solution)
    plt.show()

    return solution


def main():
    N = 1000
    t0, t1 = 0, 30
    h = (t1 - t0) / N
    file = "nonlinearOSC.txt"
    solution = solve(N, t0, t1, h, rk4)
    print(solution)
    savelog(file, N, solution)


if __name__ == "__main__":
    main()
