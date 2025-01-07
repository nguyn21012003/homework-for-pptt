import math
from math import exp as e, cos, sin, sqrt
import numpy as np
import csv


def fx1(x):
    return x - 2 ** (-x)


def fx2(x):
    return e**x - 2 - cos(e**x - 2)


def bisection(fx, a, b, N, eps, file):
    count = 0  ## dung de dem so lan lap

    listA = np.zeros(N)
    listB = np.zeros(N)
    c = np.zeros(N)
    if a > b:
        a = b
        b = a

    listA[0] = a
    listB[0] = b

    for i in range(N):
        if i + 1 < N:
            c[i] = (listA[i] + listB[i]) / 2
            if abs(fx(c[i])) < eps:
                break
            if fx(listA[i]) * fx(c[i]) < 0:
                listA[i + 1] = listA[i]
                listB[i + 1] = c[i]

                count += 1
            elif fx(c[i]) * fx(listB[i]):
                listA[i + 1] = c[i]
                listB[i + 1] = listB[i]

                count += 1
    ################ Ghi file

    with open(file, "w", newline="") as writefile:
        header = [
            f"n{'':^2}",
            f"listA{'':^18}",
            f"listB{'':^18}",
            f"listC{'':^18}",
            "f(c)",
        ]

        writer = csv.DictWriter(writefile, fieldnames=header)
        writer.writeheader()
        for i in range(count):
            writer.writerow(
                {
                    f"n{'':^2}": f"{i:^2}",
                    f"listA{'':^18}": f"{listA[i]:^18}",
                    f"listB{'':^18}": f"{listB[i]:^18}",
                    f"listC{'':^18}": f"{c[i]:^18}",
                    "f(c)": f"{fx(c[i])}",
                }
            )


def main():
    N = 100  ### so lan lap
    eps = 1e-6
    a, b = 0, 1
    fileWriteFile = "dataBisection.txt"
    bisection(fx1, a, b, N, eps, fileWriteFile)


if __name__ == "__main__":
    main()
