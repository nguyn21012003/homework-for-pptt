import numpy as np
import numpy.typing as npt
from numpy import sin, cos, pi
import matplotlib.pyplot as plt


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
    y = []
    yn = np.zeros(2)
    yn[0] = pi / 6
    yn[1] = 0
    times = np.linspace(a, b, N)
    for j in range(N):
        if yn[0] > 1:
            yn = solver(fArr1, times, yn, h)
            y.append(yn)
        if yn[0] <= 1:
            yn = solver(fArr2, times, yn, h)
            y.append(yn)
    return y


def main():

    N = 100
    t0 = 0
    tn = N / 10

    h = (tn - t0) / N
    times = np.linspace(t0, tn, N)

    sol = solve_sys_ode(fArr1, fArr2, t0, tn, h, rk4, N)
    theta = []
    omega = []

    for i in range(len(sol)):
        theta.append(sol[i][0])
        print(sol[i][0])

    plt.plot(times, theta)

    plt.show()


if __name__ == "__main__":
    main()
