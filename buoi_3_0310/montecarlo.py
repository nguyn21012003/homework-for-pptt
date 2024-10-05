import numpy as np
from random import random, uniform
from math import sqrt, pi
from numpy import sin, cos
import sympy as sym
from sympy.printing import latex
import matplotlib.pyplot as plt
import time


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
            plotc2(f1, a, b, c, d, list_x2, list_y2, list_x3, N)

            ##################################################################################
        case 3:
            ##################################################################################
            # Câu 3: Tính tích phân x = [0,pi/2] và y = [0,pi/2] của xsin(x6y)
            a, b, c, d = 0, pi / 2, 0, pi / 2
            x0 = a
            y0 = c
            ans2, list_x2, list_y2 = monte2D(f2, a, b, c, d, N, eps)
            ans3, list_x3, list_y3, hx3, hy3 = riemann2D(f2, x0, y0, a, b, c, d, N, eps)
            print(f"Nghiệm tính theo phương pháp Monter-Carlo là {ans2} \n")
            print(f"Nghiệm tính theo phương pháp tổng Rieman là {ans3}")
        case 4:
            a, b, c, d, g, h = 0, 1, 0, 1, 0, 1
            x0 = a
            y0 = c
            z0 = c
            ans2, list_x2, list_y2, listz2 = monte3D(f3, a, b, c, d, g, h, N, eps)
            print(f"Nghiệm tính theo phương pháp Monter-Carlo là {ans2} \n")


def f(x):
    return x**2


def f1(x, y):
    return cos(x**4) + 3 * y**2


def f2(x, y):
    return x * sin(x + y)


def f3(x, y, z):
    return x**2 + y**2 + z**2


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
    count = 0
    x = np.zeros(N)
    x[0] = x0
    for i in range(1, N):
        x[i] = a + i * h
        count += i
        S[i] = f(x[i])
        ans += (b - a) * S[i] / N
    return ans, x, h


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


def monte3D(f, a, b, c, d, g, h, N, eps):
    S = np.zeros(N)
    count = 0
    S[0] = 0
    ans = 0
    list_x = np.zeros(N)
    list_y = np.zeros(N)
    list_z = np.zeros(N)
    for i in range(N):
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


def riemann2D(f, x0, y0, a, b, c, d, N, eps):
    count = 0
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
        y[i] = y0 + i * hy
        count += i
        S[i] = f(x[i], y[i])
        ans += (b - a) * (d - c) * S[i] / N
    return ans, x, y, hx, hy


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


def plotc2(f1, a, b, c, d, xmonte2D, ymonted2D, xrieman2D, N):
    pass


def main():
    id = int(input("Nhập câu thứ ... để in ra kết quả của câu đó(có 3 câu): "))
    promp_user(id)
    # test print
    # print(list_x2)


if __name__ == "__main__":
    main()
