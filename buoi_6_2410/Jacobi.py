import matplotlib.pyplot as plt
import numpy as np


def a_values(d):
    a = np.zeros((d, d))
    a[0][0] = 10
    a[0][1] = -1
    a[0][2] = 2
    a[0][3] = 0
    a[1][0] = -1
    a[1][1] = 11
    a[1][2] = -2
    a[1][3] = 3
    a[2][0] = 2
    a[2][1] = -1
    a[2][2] = 10
    a[2][3] = -1
    a[3][0] = 0
    a[3][1] = 3
    a[3][2] = -1
    a[3][3] = 8
    return a


def b_values(d=4):
    b = np.zeros(d)
    b[0] = 6
    b[1] = 25
    b[2] = -11
    b[3] = 15
    return b


def main():
    N = 1000
    d = 4
    eps = 10e-20

    a = a_values(d)
    b = b_values(d)

    x = np.zeros((4, N + 1))
    for i in range(4):
        x[i][0] = 0

    k = 1
    while k <= N:
        for i in range(d):
            x_next = x[i]
            x[i][k] = (1 / a[i][i]) * (b[i])
            for j in range(d):
                if j != i:
                    x[i][k] += -(1 / a[i][i]) * (a[i][j] * x[j][k - 1])

                if abs(x[i][k] - x[i][k - 1]) < eps:
                    print(f"Nghiệm {[i+1]} hội tụ tại {x[i][k]} với {k} lần lặp.")
                    break
        k += 1
        if k + 1 > N:
            print(f"Nghiệm {[i+1]} không hội tụ sau {k} lần lặp.")


main()
