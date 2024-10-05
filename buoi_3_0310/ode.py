import numpy as np
import matplotlib.pyplot as plt
import csv
from tabulate import tabulate
from math import exp


def f(t, y):

    return y - t**2 + 1


def exact_solution(t0, tn, h):

    t = np.arange(t0, tn + h, h)
    y = np.zeros(len(t))
    for i in range(len(t)):
        y[i] = (t[i] + 1) ** 2 - 0.5 * exp(t[i])

    return y


def euler(f, tn, yn, h):

    return yn + h * f(tn, yn)


def mod_euler(f, tn, yn, h):

    return yn + h / 2 * (f(tn, yn) + f(tn, yn + h * f(tn, yn)))


def midpoint(f, tn, yn, h):

    return yn + h * f(tn + h / 2, yn + h / 2 * f(tn, yn))


def rk3(f, tn, yn, h):

    k1 = f(tn, yn)
    k2 = f(tn + 1 / 3 * h, yn + 1 / 3 * h * k1)
    k3 = f(tn + 2 / 3 * h, yn + 2 / 3 * h * k2)

    return yn + h / 4 * (k1 + 3 * k3)


def rk4(f, tn, yn, h):

    k1 = f(tn, yn)
    k2 = f(tn + 0.5 * h, yn + 0.5 * h * k1)
    k3 = f(tn + 0.5 * h, yn + 0.5 * h * k2)
    k4 = f(tn + h, yn + h * k3)

    return yn + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def solve_ode(f, y0, t0, tn, solver, h):
    t = np.arange(t0, tn + h, h)
    y = np.zeros(len(t))
    y[0] = y0

    for i in range(len(t) - 1):
        y[i + 1] = solver(f, t[i], y[i], h)
    return t, y


def plot_ty(file, t, y, yeuler, yrk4, ymidpoint, ymodeuler, yrk3):
    fig, ax = plt.subplots(figsize=(15, 7))

    ax.plot(t, y, "-.", color="r")
    ax.plot(t, yeuler, color="orange")
    ax.plot(t, ymodeuler, color="gold")
    ax.plot(t, ymidpoint, color="lime")
    ax.plot(t, yrk3, color="navy")
    ax.plot(t, yrk4, color="magenta")

    ax.legend(
        [
            "Nghiệm giải tích",
            "Phương pháp Euler",
            "Phương pháp Euler cải tiến",
            "Phương pháp Midpoint/RK2",
            "Phương pháp RK3",
            "Phương pháp RK4",
        ],
        loc="best",
    )
    ax.grid(True)
    ax.set_title(r"Nghiệm của phương trình vi phân $y'=y-t^2+1$")
    ax.set_xlabel("t")
    ax.set_ylabel("y")

    fig.savefig(file)

    plt.show()


def table_log(file: str, t: list, y: list, yeuler: list, yrk4: list, ymidpoint: list, ymodeuler: list, yrk3: list) -> tabulate:
    with open(file, "w", newline="") as writefile:
        header = [
            "n",
            "t",
            "y",
            "Euler",
            "Euler modified",
            "RK2",
            "RK3",
            "RK4",
        ]
        writer = csv.DictWriter(writefile, fieldnames=header)
        writer.writeheader()
        for i in range(len(t)):
            writer.writerow(
                {
                    "n": i,
                    "t": f"{t[i]:.4f}",
                    "y": y[i],
                    "Euler": yeuler[i],
                    "Euler modified": ymodeuler[i],
                    "RK2": ymidpoint[i],
                    "RK3": yrk3[i],
                    "RK4": yrk4[i],
                }
            )
    with open(file, "r") as readfile:
        rows = []
        reader = csv.DictReader(readfile)
        for row in reader:
            rows.append(row)
        return tabulate(rows, headers="keys", tablefmt="fancy_grid")


def save_log(file: str, t: list, y: list, yeuler: list, yrk4: list, ymidpoint: list, ymodeuler: list, yrk3: list) -> dict:
    with open(file, "w", newline="") as writefile:
        header = [
            f"{'n':^3}",
            f"{'t':^6}",
            f"{'y':^18}",
            f"{'Euler':^18}",
            f"{'Euler modified':^18}",
            f"{'RK2':^18}",
            f"{'RK3':^18}",
            f"{'RK4':^18}",
        ]
        writer = csv.DictWriter(writefile, fieldnames=header, delimiter="|")
        writer.writeheader()
        for i in range(len(t)):
            writer.writerow(
                {
                    f"{'n':^3}": f"{i:^3}",
                    f"{'t':^6}": f"{t[i]:^4.4f}",
                    f"{'y':^18}": f"{y[i]:^18}",
                    f"{'Euler':^18}": f"{yeuler[i]:^18}",
                    f"{'Euler modified':^18}": f"{ymodeuler[i]:^18}",
                    f"{'RK2':^18}": f"{ymidpoint[i]:^18}",
                    f"{'RK3':^18}": f"{yrk3[i]:^18}",
                    f"{'RK4':^18}": f"{yrk4[i]:^18}",
                }
            )


def main():
    N = 3

    tspan = [0, 2]
    t0 = tspan[0]
    tn = tspan[1]
    h = (tn - t0) / N
    y0 = 0.5

    file_log = "data.txt"
    file_data_log = "data_table.txt"
    file_pdf = "data.pdf"

    y = exact_solution(t0, tn, h)
    t1, y1 = solve_ode(f, y0, t0, tn, euler, h)
    t2, y2 = solve_ode(f, y0, t0, tn, rk4, h)
    t3, y3 = solve_ode(f, y0, t0, tn, midpoint, h)
    t4, y4 = solve_ode(f, y0, t0, tn, mod_euler, h)
    t5, y5 = solve_ode(f, y0, t0, tn, rk3, h)

    plot_ty(file_pdf, t1, y, y1, y2, y3, y4, y5)

    table = table_log(file_log, t1, y, y1, y2, y3, y4, y5)
    print(table)
    save_log(file_data_log, t1, y, y1, y2, y3, y4, y5)


if __name__ == "__main__":
    main()
