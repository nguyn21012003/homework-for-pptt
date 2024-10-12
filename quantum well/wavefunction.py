from numpy import sin, exp, cos, tan, sqrt
import numpy as np
import matplotlib.pyplot as plt
import random


def F(z):
    alpha = (1 / 2 + sin(2 * z) / (4 * z)) * a * (exp(-2 * z * tan(z)) / (cos(z) ** 2)) + a * exp(-2 * z * tan(z)) / (2 * z * tan(z))
    return sqrt(1 / (2 * alpha))


def D(z):
    return F(z) * exp(-z * tan(z)) / cos(z)


def Psi1(x):
    wave = []
    for i in range(len(z)):
        f = F(z)[i]
        Z = z[i]
        wave.append(f * exp(-z * tan(Z) * x / a))
    return np.array(wave)


def Psi2(x):
    wave = []
    for i in range(len(z)):
        d = D(z)[i]
        Z = z[i]
        wave.append(d * cos(Z * x / a))
    return np.array(wave)


a = 0.529e-10
z = np.array(
    [
        1.4707963267948965,
        4.5123889803846895,
        7.521329284916163,
        10.467110630045587,
        13.48433559328026,
        16.515365075041686,
        19.586333260392166,
        22.134885568292177,
        24.9,
    ]
)


x2 = np.linspace(0, 0.529e-10, 100)
x1 = np.linspace(0.529e-10, 1e-10, 100)
print(F(z))
print(D(z))
# print(Psi2(x2))


def random_color():
    return (random.random(), random.random(), random.random())


for i in range(len(Psi2(x2))):
    plt.plot(x2, Psi2(x2)[i], color=random_color(), label=str(i))
plt.xlim(0, 1e-10)  # Giới hạn trục x từ 0 đến 5
plt.ylim(-0.05, 0.05)
plt.legend()
plt.show()
