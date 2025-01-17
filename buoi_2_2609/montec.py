import csv
import subprocess
from random import randrange, uniform
from time import time

import matplotlib.pyplot as plt
import numpy as np
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
    print(len(xin))
    print(len(xout))

    return pi, count, xin, yin, xout, yout


def save_log(fileWrite_IN, fileWrite_Out, pi, xin, yin, xout, yout, count, N):

    with open(fileWrite_IN, "w", newline="") as writefile:
        header = [f"{'n':<4}", f"{'pi':^18}", f"{'x':^25}", f"{'y':^25}", f"{'x^2+y^2':^25}"]
        writer = csv.DictWriter(writefile, fieldnames=header)
        writer.writeheader()
        for i in range(len(xin)):
            writer.writerow(
                {
                    f"{'n':<4}": f"{count[i]:<4}",
                    f"{'pi':^18}": f"{pi[i]:^18}",
                    f"{'x':^25}": f"{xin[i]:^25}",
                    f"{'y':^25}": f"{yin[i]:^25}",
                    f"{'x^2+y^2':^25}": f"{xin[i]**2+yin[i]**2:^25}",
                },
            )
    with open(fileWrite_Out, "w", newline="") as writefile:
        header = [f"{'xout':^25}", f"{'yout':^25}"]
        writer = csv.DictWriter(writefile, fieldnames=header)
        writer.writeheader()
        for i in range(len(xout)):
            writer.writerow(
                {
                    f"{'xout':^25}": f"{xout[i]:^25}",
                    f"{'yout':^25}": f"{yout[i]:^25}",
                }
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
    N = 100000
    start = time()
    pi, count, xin, yin, xout, yout = monte(N)
    end = time()
    print(end - start)
    ########################### save file
    file = "dataCircle.csv"
    fileWriteOut = "dataCircle_out.csv"
    save_log(file, fileWriteOut, pi, xin, yin, xout, yout, count, N)
    # read = read_log(file, pi, xin, yin, count, N)
    # print(read)
    ########################### setting plot
    # plot_fig(count, pi, xin, yin, xout, yout)

    gnuplot(file, fileWriteOut)


def gnuplot(fileCircle, fileOut):
    with open("gnuPlot.gp", "w") as writefile:
        writefile.write(
            f"""
    set xlabel "x"
    set ylabel "y"
    set grid 
    set datafile separator ","

    plot "{fileCircle}" u 3:4 with points ps 1 pt 7,\
    "{fileOut}" u 1:2 with points ps 1 pt 7
    pause -1

            """
        )
    subprocess.run(["gnuplot", "gnuPlot.gp"])


if __name__ == "__main__":
    main()
