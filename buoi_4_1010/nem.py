import numpy as np
import numpy.typing as npt
from numpy import sin, cos, pi
import matplotlib.pyplot as plt
import csv


def fArrX(t, args):
    x, vx = args
    ax = 0
    F = np.zeros(2)
    F[0] = vx
    F[1] = ax

    return F


def fArrY(t, args):
    g = 9.81
    y, vy = args
    ay = -g
    F = np.zeros(2)
    F[0] = vy
    F[1] = ay

    return F


def fArrFriction(arg):
    pass


def rk4(fArr: npt.NDArray, tn: npt.NDArray, yn: npt.NDArray, h: float) -> npt.NDArray:

    k1 = fArr(tn, yn)
    k2 = fArr(tn + 0.5 * h, yn + 0.5 * h * k1)
    k3 = fArr(tn + 0.5 * h, yn + 0.5 * h * k2)
    k4 = fArr(tn + h, yn + h * k3)
    return yn + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def solve_sys_ode(fArr1, fArr2: npt.NDArray, a: float, b: float, h: float, solver: npt.NDArray, N: int) -> npt.NDArray:
    tn = np.arange(a, b + h, h)
    theta = pi / 6
    argsX = np.zeros(2)
    argsX[0] = 0
    argsX[1] = 50 * cos(theta)

    argsY = np.zeros(2)
    argsY[0] = 10
    argsY[1] = 50 * sin(theta)

    x, y = [], []

    for i in range(len(tn) - 1):
        argsX = rk4(fArr1, tn, argsX, h)
        x.append(argsX)
    for i in range(len(tn) - 1):
        argsY = rk4(fArr2, tn, argsY, h)
        y.append(argsY)

    listX = []
    listY = []

    for i in range(len(x)):
        if float(y[i][0]) > 0:
            listX.append(x[i][0])
            listY.append(y[i][0])

    plt.plot(listX, listY)
    plt.show()

    return listX, listY


def main():
    x, y = solve_sys_ode(fArrX, fArrY, 0, 100, 0.1, rk4, 100)
    print(x)
    print(y)


if __name__ == "__main__":
    main()
