from numpy import pi
from random import uniform
import csv
from tqdm import tqdm


def monte(N):
    xin, yin = [], []
    xout, yout = [], []

    ##ném 1000 lần thì sẽ có ví dụ là 731 lần là vô bên trong vậy thì số lần ném bên ngoài sẽ là 1000 - 731 lần

    for i in tqdm(range(N)):
        x = uniform(-1, 1)
        y = uniform(-1, 1)
        if x**2 + y**2 < 1:  ### (x - a)^2 + (y - b)^2 = R^2 ma ban kinh R = 1
            xin.append(x)
            yin.append(y)
        else:
            xout.append(x)
            yout.append(y)

    return xout, yout, xin, yin


xout, yout, xin, yin = monte(10000)

with open("circleOUT.csv", "w", newline="") as fileIN:
    headerOUT = ["xout", "yout"]
    OUT = csv.DictWriter(fileIN, fieldnames=headerOUT)
    OUT.writeheader()
    for i in range(len(xout)):
        OUT.writerow({"xout": xout[i], "yout": yout[i]})

with open("circleIN.csv", "w", newline="") as fileOUT:
    headerIN = ["xin", "yin"]
    IN = csv.DictWriter(fileOUT, fieldnames=headerIN)
    IN.writeheader()
    for i in range(len(xin)):
        IN.writerow({"xin": xin[i], "yin": yin[i]})
