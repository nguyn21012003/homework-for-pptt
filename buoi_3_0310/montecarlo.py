import time
from math import pi, sqrt
from random import random, uniform

import matplotlib.pyplot as plt
import numba as nb
import numpy as np
import sympy as sym
from numba import njit
from numpy import cos, sin
from sympy.printing import latex

nb.set_num_threads(6)


def promp_user(id: int):
    N = int(input("Nhập số vòng lặp N(nhập N từ 10->100 sẽ cho ra kết quả nhanh hơn, từ 1000->10000 kết quả sẽ 'smooth' hơn): "))
    eps = sqrt(N)
    match id:
        case 1:
            ##################################################################################
            # Câu 1: tính tích phân x = [0,5]  của x^2 và vẽ minh họa số sỏi ném trong miền tạo bởi f(x) và trục tung
            a, b, c, d = 0, 5, 0, 1
            x0 = a
            ans, list_x, list_y = monte(f, a, b, N, eps)
            ans1, list_x1, h1 = riemann(f, x0, a, b, N, eps)
            print(f"Nghiệm tính theo phương pháp Monter-Carlo là {ans} \n")
            print(f"Nghiệm tính theo phương pháp tổng Rieman là {ans1}")
            #########################################
            plot(f, a, b, list_x, list_y, list_x1, h1, N)
        case 2:
            ##################################################################################
            # Câu 2: Tính tích phân x = [4,6] và y = [0,1] của cos(x^4) + 3y^2
            a, b, c, d = 4, 6, 0, 1
            x0 = a
            y0 = c
            ans2, list_x2, list_y2 = monte2D(f1, a, b, c, d, N, eps)
            ans3, list_x3, list_y3, hx3, hy3 = riemann2D(f1, x0, y0, a, b, c, d, N, eps)
            print(f"Nghiệm tính theo phương pháp Monter-Carlo là {ans2} \n")
            print(f"Nghiệm tính theo phương pháp tổng Rieman là {ans3}")
            ##################################################################################
        case 3:
            ##################################################################################
            # Câu 3: Tính tích phân x = [0,pi/2] và y = [0,pi/2] của xsin(x+y)
            a, b, c, d = 0, pi / 2, 0, pi / 2
            x0 = a
            y0 = c

            start_mc = time.time()
            ans2, list_x2, list_y2 = monte2D(f2, a, b, c, d, N, eps)
            end_mc = time.time()

            start_rm = time.time()
            ans3, list_x3, list_y3, hx3, hy3 = riemann2D(f2, x0, y0, a, b, c, d, N, eps)
            end_rm = time.time()
            print(f"Nghiệm tính theo phương pháp Monter-Carlo là {ans2} với thời gian là {end_mc - start_mc}s \n")
            print(f"Nghiệm tính theo phương pháp tổng Rieman là {ans3} với thời gian là {end_rm - start_rm}s")
        case 4:
            a, b, c, d, g, h = 0, 1, 0, 1, 0, 1
            x0 = a
            y0 = c
            z0 = g

            start_mc = time.time()
            ans2, list_x2, list_y2, listz2 = monte3D(f3, a, b, c, d, g, h, N, eps)
            end_mc = time.time()

            start_rm = time.time()
            ans3, list_x3, list_y3, list_z3, hx3, hy3, hz3 = riemann3D(f3, x0, y0, z0, a, b, c, d, g, h, N, eps)
            end_rm = time.time()
            print(f"Nghiệm tính theo phương pháp Monter-Carlo là {ans2} với thời gian là {end_mc - start_mc}s \n")
            print(f"Nghiệm tính theo phương pháp tổng Rieman là {ans3} với thời gian là {end_rm - start_rm}s")
        case 5:
            a = 0
            b = 1
            x1, x3, x5, x7, x9 = [a for i in range(5)]
            x2, x4, x6, x8, x10 = [b for i in range(5)]
            start_mc = time.time()
            ans2 = monte10D(f4, a, b, N, eps)
            end_mc = time.time()

            start_rm = time.time()
            ans3 = riemann10D(f4, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, a, b, N, eps)
            end_rm = time.time()
            print(f"Nghiệm tính theo phương pháp Monter-Carlo là {ans2} với thời gian là {end_mc - start_mc}s \n")
            # print(f"Nghiệm tính theo phương pháp tổng Rieman là {ans3} với thời gian là {end_rm - start_rm}s")


def f(x):
    return x**2


def f1(x, y):
    return cos(x**4) + 3 * y**2


@njit
def f2(x, y):
    return x * sin(x + y)


@njit
def f3(x, y, z):
    return x**2 + y**2 + z**2


@njit
def f4(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10):
    return (x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10) ** 2


def monte(f, a, b, N, eps):
    count = 0
    S = np.zeros(N)
    S[0] = 0
    ans = 0
    list_x = np.zeros(N)
    list_y = np.zeros(N)
    for i in range(N):
        x = uniform(a, b)
        y = uniform(f(0), f(x))
        S[i] = f(x)
        count += 1
        ans += (b - a) * S[i] / N
        list_x[i] = x
        list_y[i] = y
    return ans, list_x, list_y


def riemann(f, x0, a, b, N, eps):
    h = (b - a) / N
    S = np.zeros(N)
    S[0] = 0
    ans = 0
    x = np.zeros(N)
    x[0] = x0
    for i in range(1, N):
        x[i] = a + i * h
        S[i] = f(x[i])
        ans += h * S[i]
    return ans, x, h


@njit(parallel=True)
def monte2D(f, a, b, c, d, N, eps):
    S = np.zeros(N)
    count = 0
    S[0] = 0
    ans = 0
    list_x = np.zeros(N)
    list_y = np.zeros(N)
    for i in range(N):
        x = uniform(a, b)
        y = uniform(c, d)
        S[i] = f(x, y)
        ans += (b - a) * (d - c) * S[i] / N
        list_x[i] = x
        list_y[i] = y
        count += i
    return ans, list_x, list_y


@njit(parallel=True)
def monte3D(f, a, b, c, d, g, h, N, eps):
    S = np.zeros(N)
    count = 0
    S[0] = 0
    ans = 0
    list_x = np.zeros(N)
    list_y = np.zeros(N)
    list_z = np.zeros(N)
    for i in nb.prange(N):
        x = uniform(a, b)
        y = uniform(c, d)
        z = uniform(g, h)
        S[i] = f(x, y, z)
        ans += (b - a) * (d - c) * (h - g) * S[i] / N
        list_x[i] = x
        list_y[i] = y
        list_z[i] = z
        count += i
    return ans, list_x, list_y, list_z


@njit(parallel=True)
def monte10D(fx, a, b, N, eps):
    S = np.zeros(N)
    S[0] = 0
    ans = 0
    for i in range(N):
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


@njit(parallel=True)
def riemann2D(f, x0, y0, a, b, c, d, N, eps):
    ans = 0

    hx = (b - a) / N
    hy = (d - c) / N

    S = np.zeros(N)
    S[0] = 0

    x = np.zeros(N)
    y = np.zeros(N)
    x[0] = x0
    y[0] = y0

    for i in range(1, N):
        x[i] = x0 + i * hx
        for j in range(1, N):
            y[j] = y0 + j * hy

            S[i] = f(x[i], y[j])
            ans += hx * hy * S[i]
    return ans, x, y, hx, hy


@njit(parallel=True)
def riemann3D(f, x0, y0, z0, a, b, c, d, g, h, N, eps):
    ans = 0

    hx = (b - a) / N
    hy = (d - c) / N
    hz = (h - g) / N

    S = np.zeros(N)
    S[0] = 0

    x = np.zeros(N)
    y = np.zeros(N)
    z = np.zeros(N)
    x[0] = x0
    y[0] = y0
    z[0] = z0

    for i in nb.prange(1, N):
        x[i] = x0 + i * hx
        for j in range(1, N):
            y[j] = y0 + j * hy
            for k in range(1, N):
                z[k] = z0 + k * hz
                S[i] = f(x[i], y[j], z[k])
                ans += hx * hy * hz * S[i]
    return ans, x, y, z, hx, hy, hz


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


def plot(f, a, b, xmonte1D, ymonte1D, xrieman1D, h, N):
    #################################################################
    # parameters
    transparent = 0.5
    x = sym.symbols("x")
    fx = sym.Integral(x**2)
    ##################################################################
    # cài đặt 1 số thứ cho plot
    fig, axs = plt.subplots(2, 2, figsize=(18, 7))
    ##################################################################
    # nghiệm của phương trình f(x) dc biểu diễn trên đoạn [0,5]
    X = np.linspace(a, b, N)
    Y = f(X)
    ####################################################################################################################################
    axs[0, 0].plot(X, Y, color="r", label=r"y=$x^2$", lw=2)
    axs[0, 0].scatter(xmonte1D, ymonte1D, marker="o", s=10)
    axs[0, 0].fill_between(X, Y, color="lightgreen", alpha=transparent)

    axs[0, 0].legend([r"y=$x^2$", "sỏi", f"${latex(fx)}$"])
    axs[0, 0].set_title(r"Monte-Carlo cho phương trình $y = x^2$")
    axs[0, 0].set_xlabel(r"x")
    axs[0, 0].set_ylabel(r"$f(x)$", rotation=0)
    ####################################################################################################################################
    axs[0, 1].plot(X, Y, color="r", label=r"y=$x^2$", lw=2)
    axs[0, 1].plot(xrieman1D[:-1], f(xrieman1D[:-1]), "b.")
    axs[0, 1].bar(xrieman1D[:-1], f(xrieman1D[:-1]), width=h, alpha=transparent, align="edge", edgecolor="b")

    axs[0, 1].legend([r"y=$x^2$"])
    axs[0, 1].set_title(r"Phuwong pháp tổng Riemann(trái) cho phương trình $y = x^2$")
    axs[0, 1].set_xlabel(r"x")
    axs[0, 1].set_ylabel(r"$f(x)$", rotation=0)
    ####################################################################################################################################
    xrieman_mid = (xrieman1D[1:] + xrieman1D[:-1]) * 0.5
    yrieman_mid = f(xrieman_mid)
    axs[1, 0].plot(X, Y, color="r", label=r"y=$x^2$", lw=2)
    axs[1, 0].plot(xrieman_mid, yrieman_mid, "b.")
    axs[1, 0].bar(xrieman_mid, yrieman_mid, width=h, alpha=transparent, edgecolor="b")

    axs[1, 0].legend([r"y=$x^2$"])
    axs[1, 0].set_title(r"Phương pháp tổng Riemann(giữa) cho phương trình $y = x^2$")
    axs[1, 0].set_xlabel(r"x")
    axs[1, 0].set_ylabel(r"$f(x)$", rotation=0)
    ####################################################################################################################################
    axs[1, 1].plot(X, Y, color="r", label=r"y=$x^2$", lw=2)
    axs[1, 1].plot(xrieman1D[1:], f(xrieman1D[1:]), "b.")
    axs[1, 1].bar(xrieman1D[1:], f(xrieman1D[1:]), width=-h, alpha=transparent, align="edge", edgecolor="b")

    axs[1, 1].legend([r"y=$x^2$"])
    axs[1, 1].set_title(r"Phương pháp tổng Riemann(phải) cho phương trình $y = x^2$")
    axs[1, 1].set_xlabel(r"x")
    axs[1, 1].set_ylabel(r"$f(x)$", rotation=0)
    ####################################################################################################################################
    fig.savefig(f"Cau1.Monte-Carlo-Riemann1D.pdf")

    plt.show()


def main():
    id = int(input("Nhập câu thứ ... để in ra kết quả của câu đó(có 5 câu): "))

    promp_user(id)
    # test print
    # print(list_x2)


if __name__ == "__main__":
    main()
