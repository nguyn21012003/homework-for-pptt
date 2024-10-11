import numpy as np


def f(z):
    return np.tan(z) - np.sqrt((z0 / z) ** 2 - 1)


def df(z):
    return 1 / (np.cos(z)) ** 2 + z0**2 / (z**2 * np.sqrt((z0) ** 2 - z**2))


def newton_raphson(f, df, p0, e, N):
    n = 1
    while n < N:
        p = p0 - f(p0) / df(p0)
        if abs(p - p0) < e:
            break
        else:
            p0 = p
            n = n + 1
    return p


hbar = 1.055e-34
m = 9.11e-31
a = 1 * 0.529e-9
V0 = 10 * 1.6e-19
z0 = a * np.sqrt(2 * m * V0) / hbar
e = 1e-15
N = 100

while True:
    p0 = float(input("Nhập giá trị khỏi tạo = "))
    if ((z0 / p0) ** 2 - 1) < 0:
        print("căn bậc hai không thể bị âm. Vui lòng nhập lại.")
    else:
        break

giai = newton_raphson(f, df, p0, e, N)
print(f"Giá trị nghiệm: {giai}")
