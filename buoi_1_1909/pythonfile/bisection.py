import math
from math import exp as e
from math import cos, sin

import csv
import numpy as np


def fx1(x):
    return float(x - 2 ** (-x))


def fx2(x):
    return float(e(x) - 2 - cos(e(x) - 2))


def bisection(fx, a, b, file):
    eps = 1e-6
    count = 0

    N = 100
    c = np.zeros(N)
    na = np.zeros(N)
    nb = np.zeros(N)
    if a > b:
        a = b
        b = a

    na[0] = a
    nb[0] = b
    for i in range(N):
        if i + 1 < N:
            c[i] = (na[i] + nb[i]) / 2
            if fx(c[i]) == 0 or abs(fx(c[i])) < eps:
                break

            if fx(na[i]) * fx(c[i]) < 0:
                nb[i + 1] = c[i]
                na[i + 1] = na[i]
                count += 1
            elif fx(c[i]) * fx(nb[i]) < 0:
                na[i + 1] = c[i]
                nb[i + 1] = nb[i]
                count += 1

    with open(file, "w", newline="") as writefile:
        header = [f"n{'':^2}", f"an{'':^18}", f"bn{'':^18}", f"cn{'':^18}", "f(cn)"]

        writer = csv.DictWriter(writefile, fieldnames=header)
        writer.writeheader()
        for i in range(count):
            writer.writerow(
                {
                    f"n{'':^2}": f"{i:^2}",
                    f"an{'':^18}": f"{na[i]:^18}",
                    f"bn{'':^18}": f"{nb[i]:^18}",
                    f"cn{'':^18}": f"{c[i]:^18}",
                    "f(cn)": f"{fx(c[i])}",
                }
            )


def main():
    file = "nghiem.txt"
    # bisection(fx1, 0, 1, file)
    bisection(fx2, 0.5, 1.5, file)


if __name__ == "__main__":
    main()
