import numpy as np
from numpy import sin, cos, exp
from random import uniform
from tqdm import tqdm
import subprocess
import csv


def monte(N, fileIN, fileOUT, filePi):
    xin, yin = [], []
    xout, yout = [], []
    listPi = []
    count = []
    count_o = 0
    ##ném 1000 lần thì sẽ có ví dụ là 731 lần là vô bên trong vậy thì số lần ném bên ngoài sẽ là 1000 - 731 lần

    for i in tqdm(range(1, N)):
        x = uniform(-1, 1)
        y = uniform(-1, 1)
        if x**2 + y**2 < 1:  ### (x - a)^2 + (y - b)^2 = R^2 ma ban kinh R = 1
            xin.append(x)
            yin.append(y)
            count_o += 1
            listPi.append(4 * count_o / i)
            count.append(count_o)
        else:
            xout.append(x)
            yout.append(y)

    with open(fileIN, "w", newline="") as fileIN:
        headerIN = ["xin", "yin"]
        IN = csv.DictWriter(fileIN, fieldnames=headerIN)
        IN.writeheader()
        for i in range(len(xin)):
            IN.writerow({"xin": xin[i], "yin": yin[i]})

    with open(fileOUT, "w", newline="") as fileOUT:
        headerOUT = ["xout", "yout"]
        OUT = csv.DictWriter(fileOUT, fieldnames=headerOUT)
        OUT.writeheader()
        for i in range(len(xout)):
            OUT.writerow({"xout": xout[i], "yout": yout[i]})

    with open(filePi, "w") as filePi:
        headerPi = ["count", "pi"]
        Pi = csv.DictWriter(filePi, fieldnames=headerPi)
        Pi.writeheader()
        for i in range(len(listPi)):
            Pi.writerow({"count": i, "pi": listPi[i]})
    return xout, yout, xin, yin


def gnuPlot(fileIN, fileOUT, filePi):
    with open("gnuPlot.gp", "w") as gnuplot:
        gnuplot.write(
            f"""
set multiplot layout 1,2
set datafile separator ","
plot "{fileIN}" using 1:2 with points title "Circle OUT", "{fileOUT}" using 1:2 with points title "Circle IN"
plot "{filePi}" using 1:2 with lines title "Pi"


pause -1

"""
        )
    subprocess.run(["gnuplot", "gnuPlot.gp"])


def main():
    N = 10000
    fileOut = "circleOUT.csv"
    fileIn = "circleIN.csv"
    filePi = "pi.csv"
    monte(N, fileIn, fileOut, filePi)
    gnuPlot(fileIn, fileOut, filePi)


if __name__ == "__main__":
    main()
