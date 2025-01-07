from random import uniform
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm


def monte(radius, N):
    count = 0
    pi = []
    xin_n, yin_n = [], []
    xout_n, yout_n = [], []
    for i in tqdm(range(1, N)):
        xi = uniform(-radius, radius)
        yi = uniform(-radius, radius)

        if xi**2 + yi**2 <= 1:
            count += 1
            pi.append(4 * count / i)
            xin_n.append(xi)
            yin_n.append(yi)
        else:
            xout_n.append(xi)
            yout_n.append(yi)

    return xin_n, yin_n, xout_n, yout_n


def plot(xin, yin, xout, yout):
    plt.figure(figsize=(15, 7))

    plt.scatter(xin, yin, c="purple")
    plt.scatter(xout, yout, c="red")
    plt.grid(True)
    plt.show()


def main():

    N = 100000
    xin, yin, xout, yout = monte(1, N)
    plot(xin, yin, xout, yout)


if __name__ == "__main__":
    main()
