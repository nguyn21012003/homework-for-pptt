import numpy as np
from math import pi, sqrt, log, exp
import matplotlib.pyplot as plt
import time

# Định nghĩa các hằng số
N = 100
hbar = 658.5  # meV.fs
khi_o = 0.5
E_R = 4.2  # meV
delta_t = 20  # femtosecond
delta_0 = 100  # meV
dt = 2  # femtosecond
t_max = 500  # femtosecond
t_0 = -3 * delta_t
e_max = 300  # meV
delta_e = e_max / N
T2 = 200  # femtosecond


# Hàm tính toán
def g(n, n1):
    d = (sqrt(n + 1) + sqrt(n1 + 1)) / (sqrt(n + 1) - sqrt(n1 + 1))
    G = (1 / sqrt((n + 1) * delta_e)) * log(abs(d))
    return G


def En(Y, n):
    sum = 0
    for i in range(1, N + 1):
        if i != n:
            sum += g(n, i) * ((Y[0][i]).real + (Y[0][i]).imag)
    E = (sqrt(E_R) / pi) * delta_e * sum
    return E


def Omega(Y, n, t):
    sum = 0
    for i in range(1, N + 1):
        if i != n:
            sum += g(n, i) * Y[1][i]
    Omega1 = (sqrt(pi) / (2 * delta_t)) * khi_o * exp(-((t / delta_t) ** 2))
    Omega2 = (sqrt(E_R) / (hbar * pi)) * delta_e * sum
    return Omega1 + Omega2


def F(t, Y):
    f = np.zeros((2, N + 1), dtype="complex")
    for i in range(N):
        f[0][i] = -2 * (Omega(Y, i, t) * np.conj(Y[1][i])).imag
        f[1][i] = (
            (-1j / hbar) * (i * delta_e - delta_0 - En(Y, i)) * (Y[1][i])
            + 1j * (1 - np.real(Y[0][i]) - np.imag(Y[0][i])) * Omega(Y, i, t)
            - Y[1][i] / T2
        )
    return f


# Hàm RK4
def RK4(Ft, t, Y):
    k1 = dt * Ft(t, Y)
    k2 = dt * Ft(t + dt / 2, Y + k1 / 2)
    k3 = dt * Ft(t + dt / 2, Y + k2 / 2)
    k4 = dt * Ft(t + dt, Y + k3)
    Y = Y + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return Y


# Hàm vẽ kết quả
def plot(time_steps, epsilon, Y0_real, Y1_abs):
    X, Y = np.meshgrid(time_steps, epsilon)  # Tạo lưới cho t và epsilon

    fig = plt.figure(figsize=(12, 6))

    # Đồ thị 3D cho Y[0][n].real
    ax1 = fig.add_subplot(121, projection="3d")
    ax1.plot_surface(X, Y, Y0_real, cmap="viridis")
    ax1.set_xlabel("t (femtosecond)")
    ax1.set_ylabel("epsilon (meV)")
    ax1.set_zlabel("Y[0][n].real")
    ax1.set_title("Đồ thị 3D của Y[0][n].real theo t và epsilon")

    # Đồ thị 3D cho |Y[1][n]|
    ax2 = fig.add_subplot(122, projection="3d")
    ax2.plot_surface(X, Y, Y1_abs, cmap="plasma")
    ax2.set_xlabel("t (femtosecond)")
    ax2.set_ylabel("epsilon (meV)")
    ax2.set_zlabel("|Y[1][n]|")
    ax2.set_title("Đồ thị 3D của |Y[1][n]| theo t và epsilon")

    plt.tight_layout()
    plt.show()


def main():
    # 2 Điều kiện ban đầu
    t = t_0
    Y = np.zeros((2, N + 1), dtype="complex")
    Y[0][0] = 0 + 1j * 0
    Y[1][0] = 0 + 1j * 0

    time_steps = np.linspace(t_0, t_max, 100)
    epsilon = [(j + 1) * delta_e for j in range(N)]
    Y0_real = np.zeros((N, len(time_steps)))
    Y1_abs = np.zeros((N, len(time_steps)))

    start_time = time.time()  # Bắt đầu đo thời gian
    with open("sbe_f(n,t).txt", "w") as f:
        # Tiêu đề của file
        s0 = "{:^20}|{:^20}|{:^20}|{:^20}\n".format("t", "epsilon", "Y1", "Y2")
        f.write(s0)
        f.write("-" * 84 + "\n")  # Dòng gạch ngang chia cách tiêu đề và dữ liệu

        index = 0  # Khởi tạo chỉ số cho time_steps
        for t in time_steps:
            Y = RK4(F, t, Y)

            # Cập nhật kết quả vào file và cập nhật giá trị để vẽ
            for n in range(N):
                Y0_real[n, index] = Y[0][n].real
                Y1_abs[n, index] = abs(Y[1][n])
                s1 = "{:^20.5e}|{:^20.5e}|{:^20.5e}|{:^20.5e}\n".format(t, epsilon[n], Y[0][n].real, abs(Y[1][n]))
                f.write(s1)

            # Cập nhật thời gian ước tính sau mỗi bước
            if index % 1 == 0 and index != 0:
                elapsed_time = time.time() - start_time
                estimated_total_time = (elapsed_time / (index + 1)) * len(time_steps)
                remaining_time = estimated_total_time - elapsed_time
                print(f"Bước {index}/{len(time_steps)}. Thời gian còn lại ước tính: {remaining_time:.2f} giây.")

            index += 1  # Tăng chỉ số

    end_time = time.time()  # Kết thúc đo thời gian
    total_time = end_time - start_time
    print(f"Chương trình đã hoàn tất! Thời gian chạy: {total_time:.2f} giây.")

    # Gọi hàm vẽ kết quả
    plot(time_steps, epsilon, Y0_real, Y1_abs)


if __name__ == "__main__":
    main()
