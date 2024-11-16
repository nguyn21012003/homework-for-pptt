import numpy as np
import numpy.typing as npt
from numpy import sin, cos, pi
import matplotlib.pyplot as plt
import csv


def fArr(t, args):
    g = 9.81
    ax, vx, ay, vy = args
    ax = 0
    ay = -g
    F = np.zeros(4)
    F[0] = vx
    F[1] = ax
    F[2] = vy
    F[3] = ay

    return F


def fArrFriction(t, args):
    g = 9.81
    k = 0.8
    ax, vx, ay, vy = args

    n = 1
    ax = -k * (vx**2 + vy**2) ** ((n - 1) / 2) * vx
    ay = -k * (vx**2 + vy**2) ** ((n - 1) / 2) * vy - g

    F = np.zeros(4)
    F[0] = vx
    F[1] = ax
    F[2] = vy
    F[3] = ay

    return F


def rk4(fArr: npt.NDArray, tn: npt.NDArray, yn: npt.NDArray, h: float) -> npt.NDArray:

    k1 = fArr(tn, yn)
    k2 = fArr(tn + 0.5 * h, yn + 0.5 * h * k1)
    k3 = fArr(tn + 0.5 * h, yn + 0.5 * h * k2)
    k4 = fArr(tn + h, yn + h * k3)
    return yn + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def solve_sys_ode(fArr1, fArr2, a: float, b: float, h: float, solver: npt.NDArray, N: int) -> npt.NDArray:
    tn = np.arange(a, b + h, h)
    # theta = pi / 3

    listTheta = {
        0: [pi / 3, r"$\dfrac{\pi}{3}$"],
        1: [pi / 4, r"$\dfrac{\pi}{4}$"],
        2: [pi / 5, r"$\dfrac{\pi}{5}$"],
        3: [pi / 6, r"$\dfrac{\pi}{6}$"],
        4: [pi / 7, r"$\dfrac{\pi}{7}$"],
        5: [pi / 8, r"$\dfrac{\pi}{8}$"],
        6: [pi / 9, r"$\dfrac{\pi}{9}$"],
        7: [pi / 10, r"$\dfrac{\pi}{10}$"],
    }

    v0 = 40

    listXY = {"x": [], "ax": [], "y": [], "ay": []}

    listXYFriction = {"x": [], "ax": [], "y": [], "ay": []}

    for key in listTheta:
        theta = listTheta[key][0]
        dataXY = {"x": [], "ax": [], "y": [], "ay": []}
        dataXYFriction = {"x": [], "ax": [], "y": [], "ay": []}

        args = np.zeros(4)
        args[0] = 0
        args[1] = v0 * cos(theta)
        args[2] = 0
        args[3] = v0 * sin(theta)

        argsFriction = np.zeros(4)
        argsFriction[0] = 0
        argsFriction[1] = v0 * cos(theta)
        argsFriction[2] = 0
        argsFriction[3] = v0 * sin(theta)

        for i in range(len(tn) - 1):
            args = rk4(fArr1, tn, args, h)

            if args[2] > 0:
                dataXY["x"].append(args[0])
                dataXY["ax"].append(args[1])
                dataXY["y"].append(args[2])
                dataXY["ay"].append(args[3])

        for i in range(len(tn) - 1):
            argsFriction = solver(fArr2, tn, argsFriction, h)
            if argsFriction[2] > 0:
                dataXYFriction["x"].append(argsFriction[0])
                dataXYFriction["ax"].append(argsFriction[1])
                dataXYFriction["y"].append(argsFriction[2])
                dataXYFriction["ay"].append(argsFriction[3])

        listXY["x"].append(dataXY["x"])
        listXY["ax"].append(dataXY["ax"])
        listXY["y"].append(dataXY["y"])
        listXY["ay"].append(dataXY["ay"])

        listXYFriction["x"].append(dataXYFriction["x"])
        listXYFriction["ax"].append(dataXYFriction["ax"])
        listXYFriction["y"].append(dataXYFriction["y"])
        listXYFriction["ay"].append(dataXYFriction["ay"])

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 7))

    for i in range(len(listXY["x"])):
        ax1.plot(listXY["x"][i], listXY["y"][i])
        ax1.legend([f"{listTheta[key][1]}" for key in listTheta])

    for i in range(len(listXYFriction["x"])):
        ax2.plot(listXYFriction["x"][i], listXYFriction["y"][i])
        ax2.legend([f"{listTheta[key][1]}" for key in listTheta])

    ax1.axhline(y=0, ls="-.")
    ax2.axhline(y=0, ls="-.")
    ax1.set_title(f"Đồ thị phương trình ném xiên không có ma sát")
    ax2.set_title(f"Đồ thị phương trình ném xiên có ma sát")

    plt.savefig("nemxien.pdf")

    plt.show()

    return listXY["x"], listXY["y"]


def main():
    x, y = solve_sys_ode(fArr, fArrFriction, 0, 100, 0.01, rk4, 100)


if __name__ == "__main__":
    main()
