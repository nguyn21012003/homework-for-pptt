from numpy import exp, cos, tan, sqrt
import numpy as np
import matplotlib.pyplot as plt
import random


# Hệ số chuẩn hóa
# ---------------------------------
def F(z):
    return exp(z * tan(z)) * cos(z) / sqrt((z * tan(z) + 1) / (z * tan(z)))


def D(z):
    return 1 / sqrt((z * tan(z) + 1) / (z * tan(z)))


# ---------------------------------


# Hàm sóng trong giếng (-a < x < a)
# ---------------------------------
def Psi1(x):
    wave = []
    for i in range(len(z)):
        d = D(z)[i]
        Z = z[i]
        wave.append(d * cos(Z * x))
    return np.array(wave)


# Hàm sóng ngoài giếng ( x > a )
# ---------------------------------
def Psi2B(x):
    wave = []
    for i in range(len(z)):
        f = F(z)[i]
        Z = z[i]
        wave.append(f * exp(-Z * tan(Z) * x))
    return np.array(wave)


# ---------------------------------


# Hàm sóng ngoài giếng ( x < a )
# ---------------------------------
def Psi2A(x):
    wave = []
    for i in range(len(z)):
        f = F(z)[i]
        Z = z[i]
        wave.append(f * exp(Z * tan(Z) * x))
    return np.array(wave)


# ---------------------------------


# Giá trị a và z
# ---------------------------------
a = 1e-9
z = np.array([1.3637840514290067, 4.054915182534266, 6.491617892999432])
# ---------------------------------


# Vẽ đồ thị
# ---------------------------------
# Khởi tạo x
# -------------
x1 = np.linspace(-1, 1, 1000)
x2A = np.linspace(-10, -1, 1000)
x2B = np.linspace(1, 10, 1000)
# -------------


for i in range(len(z)):
    m = (random.random(), random.random(), random.random())
    plt.plot(x1, Psi1(x1)[i], color=m, label="φ" + str(i+1) + " ,z" + str(i+1) + " = " + str(z[i]))
    plt.plot(x2A, Psi2A(x2A)[i], color=m)
    plt.plot(x2B, Psi2B(x2B)[i], color=m)
plt.title("Hàm riêng")
plt.ylabel("φ(x) ( $a^{-1/2}$)")
plt.xlabel("x ($a$)")
plt.grid()
plt.legend(loc="best")
plt.savefig("wavefunctions.pdf")
plt.show()
