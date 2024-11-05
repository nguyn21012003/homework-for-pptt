import numpy as np
from math import pi, sqrt, log, exp
import matplotlib.pyplot as plt
import time


# Define constants
N = 100
hbar = 658.5  # meV.fs
khi_o = 0.5
E_R = 4.2  # meV
delta_t = 100  # femtosecond
delta_0 = 100  # meV
dt = 2  # femtosecond
t_max = 500  # femtosecond
t_0 = -3 * delta_t
e_max = 300  # meV
delta_e = e_max / N
T2 = 200  # femtosecond

Y = np.zeros((2, N + 1), dtype="complex")


# Hàm tính toán


def g(n, n1):
    G = (1 / sqrt((n + 1) * delta_e)) * log(abs((sqrt(n + 1) + sqrt(n1 + 1)) / (sqrt(n + 1) - sqrt(n1 + 1))))
    return G


def En(Y, n):
    summ = 0
    for i in range(1, N + 1):
        if i != n:
            summ += g(n, i) * ((Y[0][i]).real + (Y[0][i]).imag)
    E = (sqrt(E_R) / pi) * delta_e * summ
    return E


def Omega(Y, n, t):
    summ = 0
    for i in range(1, N + 1):
        if i != n:
            summ += g(n, i) * Y[1][i]
    Omega1 = (sqrt(pi) / (2 * delta_t)) * khi_o * exp(-((t / delta_t) ** 2))
    Omega2 = (sqrt(E_R) / (hbar * pi)) * delta_e * summ
    return Omega1 + Omega2


def Ft(t, Y):
    F = np.zeros((2, N + 1), dtype="complex")
    print(np.real(Y[0]))
    for i in range(N):
        F[0][i] = -2 * (Omega(Y, i, t) * np.conj(Y[1][i])).imag + 1j * -2 * (Omega(Y, i, t) * np.conj(Y[1][i])).imag
        F[1][i] = (
            (-1j / hbar) * (i * delta_e - delta_0 - En(Y, i)) * (Y[1][i])
            + 1j * (1 - np.real(Y[0][i]) - np.imag(Y[0][i])) * Omega(Y, i, t)
            - Y[1][i] / T2
        )
    return F


def RK4(Ft, t, Y):
    k1 = dt * Ft(t, Y)
    k2 = dt * Ft(t + dt / 2, Y + k1 / 2)
    k3 = dt * Ft(t + dt / 2, Y + k2 / 2)
    k4 = dt * Ft(t + dt, Y + k3)
    Y = Y + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return Y


# Lưu trữ kết quả theo thời gian
epsilon = [(i + 1) * delta_e for i in range(N)]
Y1_e_collection = []
Y2_p_collection = []
time_steps = np.linspace(t_0, t_max, 100)
print(time_steps)
# Vòng lặp tính toán theo các bước thời gian
with open("SBE.txt", "w") as file:
    file.write(f"{'Time t':^20} {'Epsilon':^20} {'Y1_e':^20} {'Y2_p':^40} \n")
    file.write("=" * 100 + "\n")
    start = time.time()
    # Vòng lặp qua các bước thời gian
    for t_step in time_steps:
        Y = RK4(Ft, t_step, Y)
        Y1_e = []
        Y2_p = []

        # Vòng lặp qua các mức năng lượng
        for i in range(N):
            Y1_e.append(float((Y[0][i]).real))  # Phần thực của Y[0]
            Y2_p.append(abs(Y[1][i]))  # Độ lớn của Y[1]

        # Lưu kết quả cho biểu đồ
        Y1_e_collection.append(Y1_e)
        Y2_p_collection.append(Y2_p)
        # Ghi kết quả vào file
        for i in range(len(Y1_e)):
            file.write(f"{t_step:^20.5e} {epsilon[i]:^20.5e} {Y1_e[i]:^20.5e} {Y2_p[i]:^40} \n")
        file.write("=" * 100 + "\n")
    end = time.time()
    print(end - start)


# Chuyển đổi kết quả thành mảng numpy để vẽ biểu đồ
Y1_e_collection = np.array(Y1_e_collection).T  # Hình dạng: (N, số bước thời gian)
Y2_p_collection = np.array(Y2_p_collection).T  # Hình dạng: (N, số bước thời gian)

# Tạo lưới để vẽ biểu đồ
T, Eps = np.meshgrid(time_steps, epsilon)

# Tạo biểu đồ 3D cho Y1_e_collection (Hình 1)
fig1 = plt.figure(figsize=(10, 6))
ax1 = fig1.add_subplot(111, projection="3d")
surf1 = ax1.plot_surface(T, Eps, Y1_e_collection, cmap="viridis")
ax1.set_xlabel("Time (fs)")
ax1.set_ylabel("Energy (meV)")
ax1.set_zlabel("Y1_e")
ax1.set_title("Y1_e theo thời gian và mức năng lượng")

# Thêm thanh màu cho đồ thị thứ nhất
fig1.colorbar(surf1, ax=ax1, shrink=0.5, aspect=5)

# Tạo biểu đồ 3D cho Y2_p_collection (Hình 2)
fig2 = plt.figure(figsize=(10, 6))
ax2 = fig2.add_subplot(111, projection="3d")
surf2 = ax2.plot_surface(T, Eps, Y2_p_collection, cmap="plasma")
ax2.set_xlabel("Time (fs)")
ax2.set_ylabel("Energy (meV)")
ax2.set_zlabel("Y2_p")
ax2.set_title("Y2_p theo thời gian và mức năng lượng")

# Thêm thanh màu cho đồ thị thứ hai
fig2.colorbar(surf2, ax=ax2, shrink=0.5, aspect=5)

# Hiển thị cả hai biểu đồ
plt.show()
