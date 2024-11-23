import numpy as np
from math import pi, sqrt, log, exp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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
dt = 2  # femtosecond
khi_o = 0.0001  # Tham số cường độ xung (0.1 - 2)
a0 = 125 * 1e-8  # Bán kính Bohr (cm)
gamma = 6.5 * 1e-20  # cm^{3}.fs^{-1}

# Tính hệ số
C0 = (delta_e * sqrt(delta_e)) / (2 * (pi**2) * E_R ** (3 / 2) * a0**3)
E0 = 1.42 * 1e3  # meV

# Tạo mảng cho hàm phân bố fe và p
Y = np.zeros((2, N), dtype="complex")
Y[0][0] = 0 + 1j * 0  # Khởi tạo giá trị cho Y[0][0]
Y[1][0] = 0 + 1j * 0  # Khởi tạo giá trị cho Y[1][0]

# Hàm tính toán
sqrt_n = np.sqrt(np.arange(1, N + 1))
g = np.zeros((N, N))
for n in range(N):
    for n1 in range(N):
        if n1 != n:
            g[n, n1] = (1 / sqrt((n + 1) * delta_e)) * log(abs(((sqrt_n[n]) + sqrt_n[n1]) / ((sqrt_n[n] - sqrt_n[n1]))))


def En(Y, n):  # Năng lượng rời rạc
    n1 = np.arange(N) != n
    tong = np.dot(g[n, n1], (Y[0][n1].real + Y[0][n1].imag))
    return (sqrt(E_R) / pi) * delta_e * tong


def Omega(Y, n, t):  # Tần số tái chuẩn hóa Rabi
    n1 = np.arange(N) != n
    tong = np.dot(g[n, n1], Y[1][n1])
    Omega1 = (sqrt(pi) / (2 * delta_t)) * khi_o * exp(-((t / delta_t) ** 2))
    Omega2 = (sqrt(E_R) / (hbar * pi)) * delta_e * tong
    return Omega1 + Omega2


def Ft(t, Y, T2):
    F = np.zeros((2, N), dtype="complex")
    for i in range(N):
        F[0][i] = -2 * (Omega(Y, i, t) * np.conj(Y[1][i])).imag + 1j * -2 * (Omega(Y, i, t) * np.conj(Y[1][i])).imag
        F[1][i] = (
            (-1j / hbar) * ((i + 1) * delta_e - delta_0 - En(Y, i)) * (Y[1][i])
            + 1j * (1 - np.real(Y[0][i]) - np.imag(Y[0][i])) * Omega(Y, i, t)
            - (Y[1][i] / T2)
        )
    return F


def RK4(Ft, t, Y, T2):
    k1 = dt * Ft(t, Y, T2)
    k2 = dt * Ft(t + dt / 2, Y + k1 / 2, T2)
    k3 = dt * Ft(t + dt / 2, Y + k2 / 2, T2)
    k4 = dt * Ft(t + dt, Y + k3, T2)

    Y[0] = Y[0] + (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6

    Nt = C0 * np.sum(sqrt_n[:N] * Y[0].real)

    T2 = 1 / (1 / T20 + gamma * Nt)

    Y[1] = Y[1] + (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6

    Pt = C0 * np.sum(sqrt_n[:N] * abs(Y[1]))

    return Y, Nt, T2, Pt


# Lưu trữ kết quả theo thời gian
epsilon = [(i + 1) * delta_e for i in range(N)]
Y1_e_collection1 = []
Y2_p_collection1 = []
Y2_p_collection2 = []
Y2_p_collection3 = []
total_density_arr = []
total_polarz1_arr = []

time_steps = np.arange(t_0, t_max, dt)
omega = np.linspace(-250, 250, 200)  # meV

start = time.time()

################################################# WRITE FILE ###########################################################
with open("SBE_Tdiff.txt", "w") as file:
    # Ghi tiêu đề các cột vào file
    file.write(f"{'#Time t':^20} {'Epsilon':^20} {'Omega':^20} {'Y1_e':^20} {'Y2_p':^40} {'Total Density':^20} {'Total Polarization':^20}\n")

    # Vòng lặp qua các bước thời gian
    for t_step in tqdm(time_steps):
        T20 = 210
        Y, Nt, T2, Pt = RK4(Ft, t_step, Y, T20)
        Y1_e = []
        Y2_p = []

        # Vòng lặp qua các mức năng lượng
        for i in range(N):
            Y1_e.append(float((Y[0][i]).real))  # Phần thực của Y[0]
            Y2_p.append(abs(Y[1][i]))  # Độ lớn của Y[1]

        Y1_e_collection1.append(Y1_e)
        Y2_p_collection1.append(Y2_p)

        # Tính tổng mật độ (total_density)
        total_density_arr.append(Nt)

        # Tính tổng phân cực (total_polarz1)
        total_polarz1_arr.append(Pt)

        # Ghi các giá trị vào file
        for i in range(len(Y1_e)):
            file.write(f"{t_step:^20.5e} {epsilon[i]:^20.5e} {omega[i]:^20.5e} {Y1_e[i]:^20.5e} {Y2_p[i]:^40} {Nt:^20.5e} {Pt:^20.5e}\n")
        file.write("\n")

##################################### FOURIER TRANSFORMATION #####################################
with open("SBE_2Tdiff.txt", "w") as file:
    Et = np.zeros(len(time_steps))
    energy = np.zeros(200, dtype="complex")
    Re_energy = np.zeros(200)
    p_omega = np.zeros(200, dtype="complex")
    a_omega = np.zeros(200)
    file.write(f"{'#Omega':^20} {'E(omega)':^30} {'P(omega)':^30} {'a(omega)':^30}\n")

    for t in range(len(time_steps)):
        Et[t] = E0 * exp(-time_steps[t] ** 2 / delta_t**2)

    for o in range(200):
        summ1 = 0 + 0j
        summ2 = 0 + 0j
        for t in range(len(time_steps)):
            summ1 += dt * E0 * exp(-(time_steps[t] * time_steps[t]) / (delta_t * delta_t)) * np.exp(1j * omega[o] * time_steps[t] / hbar)
            summ2 += dt * total_polarz1_arr[t] * np.exp(1j * omega[o] * time_steps[t] / hbar)

        energy[o] = summ1
        Re_energy[o] = np.real(summ1)
        p_omega[o] = summ2
        a_omega[o] = np.imag(p_omega[o] / energy[o])

        file.write(f"{omega[o]:^20.5e} {energy[o]:^30e} {p_omega[o]:^30e} {a_omega[o]:^30e}\n")

########################################### CALCULATING WITH DIFFERENT T2 ########################################
"""total_polarz2 = []
for t in tqdm(range(len(time_steps))):
    T20 = 500
    Y, Nt, T2, Pt = RK4(Ft, time_steps[t], Y, T20)
    Y1_e = []
    Y2_p = []
    
    # Vòng lặp qua các mức năng lượng
    for i in range(N):
        Y1_e.append(float((Y[0][i]).real))  
        Y2_p.append(abs(Y[1][i])) 
    
    Y2_p_collection2.append(Y2_p)
    total_polarz2.append(Pt)

total_polarz3 = []
for t in tqdm(range(len(time_steps))):
    T20 = 700
    Y, Nt, T2, Pt = RK4(Ft, time_steps[t], Y, T20)
    Y1_e = []
    Y2_p = []

    for i in range(N):
        Y1_e.append(float((Y[0][i]).real))  
        Y2_p.append(abs(Y[1][i]))  
  
    Y2_p_collection3.append(Y2_p)
    total_polarz3.append(Pt)"""

####################################################################################################################

Y1_e_collection = np.array(Y1_e_collection1).T
Y2_p_collection1 = np.array(Y2_p_collection1).T
Y2_p_collection2 = np.array(Y2_p_collection2).T
Y2_p_collection3 = np.array(Y2_p_collection3).T

# Tạo lưới để vẽ biểu đồ
T, Eps = np.meshgrid(time_steps, epsilon)

# 3D Surface plots
fig = plt.figure(figsize=(16, 12), constrained_layout=True)

# Create 3D subplot for ABSORTION COEFFICIENT ALPHA(W)
ax1 = fig.add_subplot(221)
ax1.plot(omega, a_omega, label=r"$\alpha(\omega)$", color="purple", linewidth=2)
ax1.set_xlabel(r"$\hbar\omega(eV)$", fontsize=10)
ax1.set_ylabel(r"$\alpha(\omega)$", fontsize=10)
ax1.set_title("Absorption spectrum", fontsize=12)
ax1.grid(True)
ax1.legend()
ax1.figure.savefig("absorption_spectrum.png", dpi=300, bbox_inches="tight")

# Create 3D subplot
ax2 = fig.add_subplot(222)
ax2.plot(omega, energy, label=r"$E(\omega)$", color="blue", linewidth=2)
ax2.set_xlabel(r"$\omega(f^{-1})$", fontsize=10)
ax2.set_ylabel("E(w)", fontsize=10)
ax2.set_title("Electric Fourier transform", fontsize=12)
ax2.grid(True)
ax2.legend()

# Create 2D subplot for Total Density (N(t)) vs Time
ax3 = fig.add_subplot(223)
ax3.plot(time_steps, Et, label=r"$E(t)$", color="blue", linewidth=2)
ax3.set_xlabel("Time(fs)", fontsize=10)
ax3.set_ylabel("E(t)", fontsize=10)
ax3.set_title("Electric field(V/m)", fontsize=12)
ax3.grid(True)
ax3.legend()

# Create 2D subplot for Total Polarization (P(t)) vs Time
ax4 = fig.add_subplot(224)
ax4.plot(time_steps, total_polarz1_arr, label="|P1(t)| with T2 = 200", color="red", linewidth=2)
"""ax4.plot(time_steps, total_polarz2, label='|P2(t)| with T2 = 500', color='green', linewidth=2)
ax4.plot(time_steps, total_polarz3, label='|P3(t)| with T2 = 700', color='blue', linewidth=2)"""
ax4.set_xlabel("Time (fs)", fontsize=10)
ax4.set_ylabel("Total Polarization P(t)", fontsize=10)
ax4.set_title("Total Polarization vs Time", fontsize=12)
ax4.grid(True)
ax4.legend()

plt.subplots_adjust(wspace=0.3, hspace=0.3)
plt.show()

end = time.time()
interval = end - start
print("Quá trình chạy kết thúc trong ", interval)
