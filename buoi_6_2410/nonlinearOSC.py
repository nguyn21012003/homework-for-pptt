import numpy as np
from numpy import typing as npt, pi, sqrt, sin
import matplotlib.pyplot as plt
import csv


alphax = 1.1
k = 2
m = 0.5
omega0 = sqrt(k / m)
b = 10 * m * omega0


def FArr(t: npt.NDArray, initInput: npt.NDArray) -> npt.NDArray:
    x, velocity = initInput

    F = np.zeros(2)

    F[0] = velocity
    F[1] = -(k * x + k * x * alphax) / m

    return F


def FArrExt(t: npt.NDArray, initInput: npt.NDArray) -> npt.NDArray:
    x, velocity = initInput

    F = np.zeros(2)

    F_Ext = 15 * sin(omega0 * t)

    F[0] = velocity
    F[1] = -(k * x - k * x * alphax) / m + F_Ext

    return F


def FArrVis(t: npt.NDArray, initInput: npt.NDArray) -> npt.NDArray:
    x, velocity = initInput

    F = np.zeros(2)

    F_Viscous = -b * velocity

    F[0] = velocity
    F[1] = -(k * x) / m + F_Viscous

    return F


def FArrExtVis(t: npt.NDArray, initInput: npt.NDArray) -> npt.NDArray:
    x, velocity = initInput

    F = np.zeros(2)

    F_Ext = 15 * sin(omega0 * t)
    F_Viscous = -b * velocity

    F[0] = velocity
    F[1] = -(k * x) / m + F_Ext + F_Viscous

    return F


def rk4(fArr: npt.NDArray, tn: npt.NDArray, yn: npt.NDArray, h: float) -> npt.NDArray:

    k1 = fArr(tn, yn)
    k2 = fArr(tn + 0.5 * h, yn + 0.5 * h * k1)
    k3 = fArr(tn + 0.5 * h, yn + 0.5 * h * k2)
    k4 = fArr(tn + h, yn + h * k3)
    return yn + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def savelog(file: str, N: int, solution: npt.NDArray):
    with open(file, "w", newline="") as writefile:
        header = ["n", "x", "velocity"]
        writer = csv.DictWriter(writefile, fieldnames=header)
        writer.writeheader()
        for i in range(N):
            writer.writerow({"n": i, "x": solution[i][0], "velocity": solution[i][1]})


def plot(
    file: str,
    t: npt.NDArray,
    solutionWOExtField: npt.NDArray,
    solutionWExtField: npt.NDArray,
    solutionWVis: npt.NDArray,
    solutionWExtFieldVis: npt.NDArray,
):
    x = {"xWOExtField": [], "xWExtField": [], "xWVis": [], "xWExtFieldVis": []}

    velocity = {"vWOExtField": [], "vWExtField": [], "vWVis": [], "vWExtFieldVis": []}

    for n in range(len(t)):
        x["xWOExtField"].append(solutionWOExtField[n][0])
        velocity["vWOExtField"].append(solutionWOExtField[n][1])

        x["xWExtField"].append(solutionWExtField[n][0])
        velocity["vWExtField"].append(solutionWExtField[n][1])

        x["xWVis"].append(solutionWVis[n][0])
        velocity["vWVis"].append(solutionWVis[n][1])

        x["xWExtFieldVis"].append(solutionWExtFieldVis[n][0])
        velocity["vWExtFieldVis"].append(solutionWExtFieldVis[n][1])

    fig, ax = plt.subplots(2, 2, figsize=(15, 7))
    ax[0][0].plot(t, x["xWOExtField"])
    ax[0][0].plot(t, velocity["vWOExtField"])
    ax[0][0].legend([r"$x$", r"v(m/s)"], fontsize=15)
    ax[0][0].set_title(r"Khi chưa có trường ngoài", loc="center")

    ax[0][1].plot(t, x["xWExtField"])
    ax[0][1].plot(t, velocity["vWExtField"])
    ax[0][1].legend([r"$x$", r"v(m/s)"], fontsize=15)
    ax[0][1].set_title(r"Khi có trường ngoài", loc="center")

    ax[1][0].plot(t, x["xWVis"])
    ax[1][0].plot(t, velocity["vWVis"])
    ax[1][0].legend([r"$x$", r"v(m/s)"], fontsize=15)
    ax[1][0].set_title(r"Khi có lực ma sát và $\alpha = 0$", loc="center")

    ax[1][1].plot(t, x["xWExtFieldVis"])
    ax[1][1].plot(t, velocity["vWExtFieldVis"])
    ax[1][1].legend([r"$x$", r"v(m/s)"], fontsize=15)
    ax[1][1].set_title(r"Khi có trường ngoài và lực ma sát", loc="center")

    plt.savefig(file)
    plt.show()


def solve(F: npt.NDArray, N: int, t0: float, t1: float, h: float, solver: npt.NDArray) -> npt.NDArray:

    initInput = np.zeros(2)
    initInput[0] = 5  # Điều kiện đầu cho tại x = 1 so với vị trí cân bằng là x = 0
    initInput[1] = 0  # Điều kiện đầu cho vật tại x = 1 là lúc thả tay ra thì v = 0

    t = np.linspace(t0, t1, N)
    solution = []

    for i in range(N):
        initInput = solver(F, t[i], initInput, h)
        solution.append(initInput)

    return t, solution


def main():
    N = 1000
    t0, t1 = 0, 50
    h = (t1 - t0) / N
    fileWrite = "nonlinearOSC.txt"
    filePlot = "nonlinearOSC.png"

    t, solutionWOExtField = solve(FArr, N, t0, t1, h, rk4)
    t, solutionWExtField = solve(FArrExt, N, t0, t1, h, rk4)
    t, solutionWVis = solve(FArrVis, N, t0, t1, h, rk4)
    t, solutionWExtFieldVis = solve(FArrExtVis, N, t0, t1, h, rk4)

    plot(filePlot, t, solutionWOExtField, solutionWExtField, solutionWVis, solutionWExtFieldVis)

    savelog(fileWrite, N, solutionWOExtField)


if __name__ == "__main__":
    main()
