from math import sqrt
import numpy as np
import time


def fx1(xn):
    return sqrt(10 - xn**3) / 2


def fx2(xn):
    return xn - (xn**3 + 4 * xn**2 - 10)


def fx3(xn):
    return sqrt(10 / xn - 4 * xn)


def fx4(xn):
    return sqrt(10 / (xn + 4))


def fx5(xn):
    return xn - (xn**3 + 4 * xn**2 - 10) / (3 * xn**2 + 8 * xn)


def fx6(xn):
    return 1 / (xn - 1)


def fixed_point(fx, p0, a, b, N, eps):
    while True:
        try:
            if a <= p0 <= b:
                x = np.zeros(N)
                x[0] = p0
                count = 0
                for i in range(N):
                    if i + 1 < N:

                        try:
                            x[i + 1] = fx(x[i])
                            print(f"x[{i}] = {x[i]}")
                            if abs(x[i] - x[i - 1]) <= eps:
                                return f"nghiem hoi tu {fx(x[count])} voi {count} vong lap."
                        except ZeroDivisionError as Z:
                            return f"Loi chia cho khong :{Z}"
                        except ValueError as ve:
                            return f"Loi chia cho khong :{ve}"

                        count += 1

                    if i + 1 >= N:
                        return f"nghiem khong hoi tu voi {count} lap."
            else:
                p0 = float(input("Nhập lại p0: "))
        except KeyboardInterrupt:
            raise KeyboardInterrupt


def main():
    p0 = float(input("Nhập p0: "))

    no1 = fixed_point(fx1, p0, 1, 2, 10, 1e-5)
    print(no1)

    no2 = fixed_point(fx2, p0, 1, 2, 100, 1e-5)
    print(no2)

    no3 = fixed_point(fx3, p0, 1, 2, 100, 1e-5)
    print(no3)

    no4 = fixed_point(fx4, p0, 1, 2, 100, 1e-5)
    print(no4)

    no5 = fixed_point(fx5, p0, 1, 2, 100, 1e-5)
    print(no5)


if __name__ == "__main__":
    main()
