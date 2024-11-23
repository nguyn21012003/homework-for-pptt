import numpy as np
import matplotlib.pyplot as plt
from numpy import pi, sqrt, log, exp

import time
from tqdm import tqdm

# Định nghĩa các hằng số
N = 100  # Số phương trình
hbar = 658.5  # meV.fs
E_R = 4.2  # meV
delta_t = 25  # femtosecond
delta_0 = 100  # meV (Độ bơm để electron đi sâu vào dải -> delta_0 càng lớn thì electron càng vào sâu bên trong dải dẫn!) = hbar*omega_0 - E_g
t_max = 1000  # femtosecond
t_0 = -3 * delta_t  # femtosecond
e_max = 300  # meV
delta_e = e_max / N  # meV
T20 = 210  # femtosecond
dt = 2  # femtosecond
khi_o = 0.001  # Tham số cường độ xung (0.1 - 2)
a0 = 125  # Bán kính Bohr (Amstrong)
gamma = 6.5 * 1e-20  # cm^{3}.fs^{-1}
M = int((t_max - t_0) / dt)  # Số bước lặp RK4
steps = np.arange(t_0, t_max, dt)  # Mảng thời gian

# Tính C_0
C_0 = delta_e * sqrt(delta_e) * (1e24) / (2 * pi**2 * E_R ** (3 / 2) * a0**3)

# Tính trước giá trị của g(n, n1)
n_values = np.arange(1, N + 1)
sqrt_n = np.sqrt(n_values)
g = np.zeros((N, N))
for n in range(N):
    for n1 in range(N):
        if n1 != n:
            g[n, n1] = (1 / sqrt((n + 1) * delta_e)) * log(
                abs((sqrt_n[n] + sqrt_n[n1]) / (sqrt_n[n] - sqrt_n[n1]))
            )  # Hệ số tương tác giữa mức năng lượng n với n1?


def En(Y, n):  # Năng lượng rời rạc
    n1 = np.arange(N) != n
    tong = np.dot(g[n, n1], (Y[0][n1].real + Y[0][n1].imag))
    return (sqrt(E_R) / pi) * delta_e * tong


def Omega(Y, n, t):  # Tần só tái chuẩn hóa Rabi
    n1 = np.arange(N) != n
    tong = np.dot(g[n, n1], Y[1][n1])
    Omega1 = (sqrt(pi) / (2 * delta_t)) * khi_o * exp(-((t / delta_t) ** 2))
    Omega2 = (sqrt(E_R) / (hbar * pi)) * delta_e * tong
    return Omega1 + Omega2


def F(t, Y, T2) -> float:
    f = np.zeros((2, N), dtype="complex")
    for i in range(N):
        f[0][i] = -2 * (Omega(Y, i, t) * np.conj(Y[1][i])).imag + 1j * (-2) * (Omega(Y, i, t) * np.conj(Y[1][i])).imag
        f[1][i] = (
            (-1j / hbar) * ((i + 1) * delta_e - delta_0 - En(Y, i)) * Y[1][i]
            + 1j * (1 - np.real(Y[0][i]) - np.imag(Y[0][i])) * Omega(Y, i, t)
            - Y[1][i] / T2
        )
    return f


# Hàm E(t)
def Et(t):
    Et = ((hbar * sqrt(pi)) / delta_t) * khi_o * np.exp(-(t**2) / (delta_t**2))
    return Et


# Hàm RK4 tính f_e, p, N(t), P(t)
def RK4(Ft, t, Y, T2):
    k1 = dt * Ft(t, Y, T2)
    k2 = dt * Ft(t + dt / 2, Y + k1 / 2, T2)
    k3 = dt * Ft(t + dt / 2, Y + k2 / 2, T2)
    k4 = dt * Ft(t + dt, Y + k3, T2)

    Y[0] = Y[0] + (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6

    Nt = C_0 * np.sum(sqrt_n[:N] * Y[0].real)

    T2 = 1 / (1 / T20 + gamma * Nt)

    Y[1] = Y[1] + (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6

    Pt = delta_e * sqrt(delta_e) * np.sum(sqrt_n[:N] * Y[1])
    PT = delta_e * sqrt(delta_e) * np.sum(sqrt_n[:N] * np.abs(Y[1]))

    return Y, Nt, T2, Pt, PT


# Fourier Transform
def FT(t, Pt):
    a = -200
    b = 200
    S = 1000
    delta = (b - a) / S
    H = np.linspace(a, b, S)

    P = np.zeros([S], dtype="complex")
    E = np.zeros([S], dtype="complex")
    alpha = np.zeros((S))

    for n in range(S):
        homega = a + n * delta
        P[n] = dt * np.sum(Pt * exp(1j * homega * t / hbar)) * (1 / hbar)
        E[n] = dt * np.sum(Et(t) * exp(1j * homega * t / hbar)) * (1 / hbar)
        alpha[n] = np.imag(P[n] / E[n])

    write_FT("Absorption Spectra Coulomb.txt", H, alpha, np.real(P), np.imag(P), np.real(E), np.imag(E))


# Hàm ghi kết quả vào file
def write_result(filename, steps, epsilon, Y0_real, Y0_imag, Y1_abs):
    with open(filename, "w") as f:
        f.write(f"#{'t':^20}{'epsilon':^20}{'Y0_real':^20}{'Y0_imag':^20}{'Y1_abs':^20}\n")
        for i, t in enumerate(steps):
            for j in range(N):
                f.write(f"{t:^20.5f}{epsilon[j]:^20.5f}{Y0_real[j, i]:^20.5e}{Y0_imag[j, i]:^20.5e}{Y1_abs[j, i]:^20.5e}\n")
            f.write("\n")


def write_Nt_Pt(filename, steps, Nt, Pt):
    with open(filename, "w") as f:
        f.write(f"#{'t':^20}{'Nt':^20}{'Pt':^20}\n")
        for i, t in enumerate(steps):
            f.write(f"{t:^20.5f}{Nt[i]:^20.5e}{Pt[i]:^20.5e}\n")


def write_FT(filename, H, alpha, RE_P, IM_P, RE_E, IM_E):
    with open(filename, "w") as f:
        f.write(f"#{'omega':^20}{'alpha':^20}{'RE_P':^20}{'IM_P':^20}{'RE_E':^20}{'IM_E':^20}\n")
        for i in range(len(H)):
            f.write(f"{H[i]:^20.5f}{alpha[i]:^20.5e}{RE_P[i]:^20.5e}{IM_P[i]:^20.5e}{RE_E[i]:^20.5e}{IM_E[i]:^20.5e}\n")


# Hàm vẽ kết quả
def plot(steps, epsilon, Y0_real, Y1_abs, Nt, Pt):
    X, Y = np.meshgrid(steps, epsilon)
    fig1 = plt.figure(figsize=(15, 10))

    ax1 = fig1.add_subplot(121, projection="3d")
    ax1.plot_surface(X, Y, Y0_real, cmap="viridis")
    ax1.set_xlabel("t (femtosecond)")
    ax1.set_ylabel("epsilon (meV)")
    ax1.set_zlabel(r"$f_e$")
    ax1.set_title(r"Đồ thị của $f_e$ theo t và epsilon")

    ax2 = fig1.add_subplot(122, projection="3d")
    ax2.plot_surface(X, Y, Y1_abs, cmap="plasma")
    ax2.set_xlabel("t (femtosecond)")
    ax2.set_ylabel("epsilon (meV)")
    ax2.set_zlabel(r"$|p_n|$")
    ax2.set_title(r"Đồ thị của $|p_n|$ theo t và epsilon")

    fig2 = plt.figure(figsize=(15, 6))
    # Vẽ đồ thị 2D cho N(t) theo t
    ax3 = fig2.add_subplot(121)
    ax3.plot(steps, Nt, color="blue")
    ax3.set_xlabel("t (femtosecond)")
    ax3.set_ylabel("N(t)")
    ax3.set_title("Đồ thị của N(t) theo t")

    # Vẽ đồ thị 2D cho |P(t)| theo t
    ax4 = fig2.add_subplot(122)
    ax4.plot(steps, Pt, color="red")
    ax4.set_xlabel("t (femtosecond)")
    ax4.set_ylabel("|P(t)|")
    ax4.set_title("Đồ thị của |P(t)| theo t")

    # fig3 = plt.figure(figsize=(15, 6))
    # Vẽ đồ thị alpha theo M
    # M = np.linspace(-100, 100, len(alpha))  # Giả sử M thay đổi từ -100 đến 100
    # plt.plot(M, alpha, color='purple')
    # plt.xlabel('M')
    # plt.ylabel(r'$\alpha$')
    # plt.title('Đồ thị của $\\alpha$ theo M')

    plt.tight_layout()
    plt.show()


# Hàm chính
def main():
    Y = np.zeros((2, N), dtype="complex")

    epsilon = [(j + 1) * delta_e for j in range(N)]

    Y0_real = np.zeros((N, len(steps)))
    Y0_imag = np.zeros((N, len(steps)))
    Y1_abs = np.zeros((N, len(steps)))
    Nt = np.zeros((len(steps)))
    Pt = np.zeros((len(steps)), dtype="complex")
    PT = np.zeros((len(steps)))
    T2_values = np.zeros_like(Nt)

    start_time = time.time()

    # Thực hiện tính
    for i, t in enumerate(tqdm(steps, desc="Processing", unit="steps")):
        Y, Nt[i], T2_values, Pt[i], PT[i] = RK4(F, t, Y, T20)
        print(Y)

        Y0_real[:, i] = Y[0, :N].real
        Y0_imag[:, i] = Y[0, :N].imag
        Y1_abs[:, i] = abs(Y[1, :N])

    FT(steps, Pt)

    # Ghi dữ liệu vào file
    write_result("Data Coulomb.txt", steps, epsilon, Y0_real, Y0_imag, Y1_abs)
    write_Nt_Pt("Nt Pt Coulomb.txt", steps, Nt, PT)

    end_time = time.time()
    print(f"Thời gian chạy: {end_time - start_time:.2f} giây!!!")
    plot(steps, epsilon, Y0_real, Y1_abs, Nt, PT)


if __name__ == "__main__":
    main()
