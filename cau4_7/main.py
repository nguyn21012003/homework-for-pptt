import numpy as np
import matplotlib.pyplot as plt

import re, sys, csv


# Caau 4
def function_x(x):
    f = 2 * np.cos(x) - x
    return f


def bisection_algorith(f, a, b):
    if a > b:
        a = b
        b = a
    count = 0
    while True:
        try:
            c = (a + b) / 2
            if f(c) == 0:
                return f"nghiệm là {c} với số lần lặp là {count}"
            elif f(a) * f(c) < 0:
                a = a
                b = c
                count += 1
            elif f(c) * f(b) < 0:
                a = c
                b = b
                count += 1
        except EOFError:
            return c


############################################################################################################


# câu 7
def ivpfunction(t, y) -> float:  # chỗ return có thể thay bằng bất kì hàm nào, ở đây là đề bài: 1 + (x+y)^2 trong đó x thay bằng t

    return 1 + (t - y) ** 2


def euler_method(f, tn, yn, h):
    return yn + h * f(tn, yn)


def t_arr(arg: str) -> tuple:
    pattern = re.sub(r"[\([{})\]]", "", arg)
    t0, tend = pattern.split(",")
    return (float(t0), float(tend))


def y0_init():
    y0 = input("Input y0: ")
    return float(y0)


def hstep():
    h = input("Input h: ")
    return float(h)


def solve_odefunction(f, t, y0, h, solver):
    t = np.arange(t[0], t[1] + h, h)
    y = np.zeros(len(t))  # Create an array from t0 to tend with step h
    y[0] = y0
    for n in range(len(t) - 1):
        y[n + 1] = solver(f, t[n], y[n], h)
    return t, y


def plot_solution(name):
    with open(name, "r") as file:
        reader = csv.DictReader(file)
        rows = []
        for row in reader:
            rows.append(row)
    tpoints = []
    ypoints = []
    for i in range(len(rows)):
        t = rows[i]["t"]
        y = rows[i]["y"]
        tpoints.append(t)
        ypoints.append(y)

    fig, ax = plt.subplots()
    plt.plot(tpoints, ypoints)
    plt.xlabel("$t")
    plt.ylabel("$y")
    plt.show()


def main():
    """
    File gồm 2 câu chính là câu 4 và câu 7, để kiểm tra câu 4 thì mình xóa dấu hashtag "#", kiểm tra câu 7 thì giữ nguyên file main.py. Câu 7 sẽ tạo ra một file tên là euler.txt. File này sẽ chưas lời giải số cho bất kì phương trình giá trị đầu nào
    """
    # câu 4
    # print(bisection_algorith(function_x, 0, 10))
    #########################################################################################################
    # câu 7
    t, y = solve_odefunction(ivpfunction, t_arr(input("Input t range: ")), y0_init(), hstep(), euler_method)
    with open("euler.txt", "w", newline="") as writefile:
        header = ["t", "y"]
        writer = csv.DictWriter(writefile, fieldnames=header)
        writer.writeheader()
        for i in range(len(t)):
            writer.writerow({"t": f"{t[i]:.2f}", "y": f"{y[i]:.4f}"})
    plot_solution("euler.txt")
    #########################################################################################################


if __name__ == "__main__":
    main()
