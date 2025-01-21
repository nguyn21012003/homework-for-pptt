import numpy as np
from numpy import sqrt
from random import uniform
from tqdm import tqdm
import csv
import matplotlib.pyplot as plt
import subprocess
from datetime import datetime
import time
import numba as nb
from numba import njit


def main():
    N = 10
    eps = 0.01
    time_run = datetime.now().strftime("%H%M")
    print(time_run)
    file = f"c3_{time_run}.csv"
    print(file)

    a, b = 1, 5
    start_mc1 = time.time()
    print(start_mc1)

    ans, list_x, list_y = monte(f1, a, b, N, eps, file)
    end_mc1 = time.time()

    print(f"{end_mc1 - start_mc1}s")

    gnuPlot(file)
    plot(list_x, list_y)

    a2 = 0
    b2 = 1

    x1, x3, x5, x7, x9 = [a for i in nb.prange(5)]
    x2, x4, x6, x8, x10 = [b for i in nb.prange(5)]
    start_mc = time.time()
    ans2 = monte10D(f2, a2, b2, N, eps)
    end_mc = time.time()

    start_rm = time.time()
    ans3 = riemann10D(f2, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, a2, b2, N, eps)
    end_rm = time.time()

    saveTime(
        start_mc1,
        end_mc1,
        start_rm,
        end_rm,
        start_mc,
        end_mc,
        end_mc - start_mc,
        end_rm - start_rm,
        end_mc1 - start_mc1,
        f"time_mc_{time_run}.csv",
    )

    print(f"Nghiệm tính theo phương pháp Monter-Carlo là {ans2} với thời gian là {end_mc - start_mc}s \n")
    print(f"Nghiệm tính theo phương pháp tổng Rieman là {ans3} với thời gian là {end_rm - start_rm}s")


def f1(x):
    return x**2


@njit
def f2(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10):
    return (x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10) ** 2


@njit(parallel=True)
def riemann10D(f, x1_0, x2_0, x3_0, x4_0, x5_0, x6_0, x7_0, x8_0, x9_0, x10_0, a, b, N, eps):
    count = 0
    ans = 0

    hx = (b - a) / N

    S = np.zeros(N)
    S[0] = 0

    x1 = np.zeros(N)
    x2 = np.zeros(N)
    x3 = np.zeros(N)
    x4 = np.zeros(N)
    x5 = np.zeros(N)
    x6 = np.zeros(N)
    x7 = np.zeros(N)
    x8 = np.zeros(N)
    x9 = np.zeros(N)
    x10 = np.zeros(N)

    x1[0] = x1_0
    x2[0] = x2_0
    x3[0] = x3_0
    x4[0] = x4_0
    x5[0] = x5_0
    x6[0] = x6_0
    x7[0] = x7_0
    x8[0] = x8_0
    x9[0] = x9_0
    x10[0] = x10_0

    for k1 in nb.prange(1, N):
        x1[k1] = x1_0 + k1 * hx
        for k2 in nb.prange(1, N):
            x2[k2] = x2_0 + k2 * hx
            for k3 in nb.prange(1, N):
                x3[k3] = x3_0 + k3 * hx
                for k4 in nb.prange(1, N):
                    x4[k4] = x4_0 + k4 * hx
                    for k5 in nb.prange(1, N):
                        x5[k5] = x5_0 + k5 * hx
                        for k6 in nb.prange(1, N):
                            x6[k6] = x6_0 + k6 * hx
                            for k7 in nb.prange(1, N):
                                x7[k7] = x7_0 + k7 * hx
                                for k8 in nb.prange(1, N):
                                    x8[k8] = x8_0 + k8 * hx
                                    for k9 in nb.prange(1, N):
                                        x9[k9] = x9_0 + k9 * hx
                                        for k10 in nb.prange(1, N):
                                            x10[k10] = x10_0 + k10 * hx

                                            S[k1] = f(x1[k1], x2[k2], x3[k3], x4[k4], x5[k5], x6[k6], x7[k7], x8[k8], x9[k9], x10[k10])
                                            ans += hx**10 * S[k1]
    return ans


def monte(f, a, b, N, eps, file):
    count = 0
    S = np.zeros(N)
    S[0] = 0
    ans = 0
    list_x = np.zeros(N)
    list_y = np.zeros(N)
    for i in tqdm(nb.prange(N)):
        x = uniform(a, b)
        y = uniform(f(0), f(x))
        S[i] = f(x)
        count += 1
        ans += (b - a) * S[i] / N
        list_x[i] = x
        list_y[i] = y
        print(count)
        if abs(ans - f(x)) < eps:

            break

    with open(file, "w", newline="") as file:
        header = ["x", "y"]
        writer = csv.DictWriter(file, fieldnames=header)
        writer.writeheader()
        for i in nb.prange(N):
            writer.writerow({"x": list_x[i], "y": list_y[i]})

    return ans, list_x, list_y


def monte10D(fx, a, b, N, eps):
    S = np.zeros(N)
    S[0] = 0
    ans = 0
    for i in tqdm(nb.prange(N)):
        x1 = uniform(a, b)
        x2 = uniform(a, b)
        x3 = uniform(a, b)
        x4 = uniform(a, b)
        x5 = uniform(a, b)
        x6 = uniform(a, b)
        x7 = uniform(a, b)
        x8 = uniform(a, b)
        x9 = uniform(a, b)
        x10 = uniform(a, b)

        S[i] = fx(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
        ans += (b - a) ** 10 * S[i] / N
    return ans


def gnuPlot(file):
    with open("gnuPlot.gp", "w") as gnuplot:
        gnuplot.write(
            f"""
set datafile separator ","

plot "{file}" using 1:2 with points title "Monte Carlo" ps 1 pt 7,\
     x**2 with lines title "Nghiệm chính xác"

pause -1
"""
        )
    subprocess.run(["gnuplot", "gnuPlot.gp"])


def saveTime(start_mc1, end_mc1, start_rm, end_rm, start_mc, end_mc, time_mc, time_rm, time_mc1, fileTime):
    with open(fileTime, "w", newline="") as file:
        header = ["start_mc1", "end_mc1", "start_rm", "end_rm", "start_mc", "end_mc", "time_mc", "time_rm", "time_mc1"]
        writer = csv.DictWriter(file, fieldnames=header)
        writer.writeheader()
        writer.writerow(
            {
                "start_mc1": start_mc1,
                "end_mc1": end_mc1,
                "start_rm": start_rm,
                "end_rm": end_rm,
                "start_mc": start_mc,
                "end_mc": end_mc,
                "time_mc": time_mc,
                "time_rm": time_rm,
                "time_mc1": time_mc1,
            }
        )


def plot(list_x, list_y):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(list_x, list_y, "ro")
    plt.show()


if __name__ == "__main__":
    main()
