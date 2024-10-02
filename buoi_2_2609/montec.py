from random import randrange
from random import uniform
import numpy as np
import csv
import matplotlib.pyplot as plt
from time import time
from tabulate import tabulate


def monte(N):
    count_o = 0
    pi = []
    count = []
    xin, yin = [], []
    xout, yout = [], []
    for i in range(1, N):
        xi = uniform(-1, 1)
        yi = uniform(-1, 1)

        if xi**2 + yi**2 < 1:
            count_o += 1
            pi.append(4 * count_o / i)
            count.append(count_o)
            xin.append(xi)
            yin.append(yi)
        else:
            xout.append(xi)
            yout.append(yi)
    return pi, count, xin, yin, xout, yout


def save_log(file, pi, x, y, count, N):

    with open(file, "w", newline="") as writefile:
        header = [f"{'n':<4}", f"{'pi':^18}", f"{'x':^25}", f"{'y':^25}", f"{'x^2+y^2':^25}"]
        writer = csv.DictWriter(writefile, fieldnames=header)
        writer.writeheader()
        for i in range(len(pi)):
            writer.writerow(
                {
                    f"{'n':<4}": f"{count[i]:<4}",
                    f"{'pi':^18}": f"{pi[i]:^18}",
                    f"{'x':^25}": f"{x[i]:^25}",
                    f"{'y':^25}": f"{y[i]:^25}",
                    f"{'x^2+y^2':^25}": f"{x[i]**2+y[i]**2:^25}",
                },
            )


def read_log(file, pi, x, y, count, N):
    with open(file, "r") as file:
        reader = csv.DictReader(file)
        rows = []
        for row in reader:
            rows.append(row)

        return tabulate(rows, headers="keys", tablefmt="fancy-grid")


def plot_fig(count, pi, xin, yin, xout, yout):

    plt.figure(figsize=(15, 7))
    plt.subplot(1, 2, 1)

    plt.plot(count, pi)

    plt.subplot(1, 2, 2)
    plt.scatter(xin, yin, color="red")
    plt.scatter(xout, yout, color="blue")
    plt.grid(True, axis="both")

    plt.show()


def main():
    N = 1000
    start = time()
    pi, count, xin, yin, xout, yout = monte(N)
    end = time()
    print(end - start)
    ########################### save file
    file = "log_draw_circle.txt"
    save_log(file, pi, xin, yin, count, N)
    read = read_log(file, pi, xin, yin, count, N)
    print(read)
    ########################### setting plot
    plot_fig(count, pi, xin, yin, xout, yout)


if __name__ == "__main__":
    main()
