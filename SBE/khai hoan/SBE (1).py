import numpy as np
from math import pi, sqrt, log, exp
import matplotlib.pyplot as plt
import time
from tqdm import tqdm

# Định nghĩa các hằng số
N = 300                     # Số phương trình
hbar = 658.5                # meV.fs
E_R = 4.2                   # meV
delta_t = 30                # femtosecond
delta_0 = 50                # meV (Độ bơm để electron đi sâu vào dải -> delta_0 càng lớn thì electron càng vào sâu bên trong dải dẫn!) = hbar*omega_0 - E_g
t_max = 300                 # femtosecond
t_0 = -3 * delta_t          # femtosecond
e_max = 300                 # meV
delta_e = e_max / N         # meV
T2 = 200                    # femtosecond
dt = 2                      # femtosecond
khi_o = 2
a0 = 125                    # Bán kính Bohr

# Tính C_0
#C_0 = delta_e * sqrt(delta_e) * (1e24) / (pi**2 * E_R**(3/2)* a0**3)

# Tính trước giá trị của g(n, n1)
n_values = np.arange(1, N + 1)
sqrt_n = np.sqrt(n_values)
g = np.zeros((N, N))
for n in range(N):
    for n1 in range(N):
        if n1 != n:
            g[n, n1] = (1 / sqrt((n + 1) * delta_e)) * log(abs((sqrt_n[n] + sqrt_n[n1]) / (sqrt_n[n] - sqrt_n[n1])))

def En(Y, n):
    n1 = np.arange(N) != n
    tong = np.dot(g[n, n1], (Y[0][n1].real + Y[0][n1].imag))
    return (sqrt(E_R) / pi) * delta_e * tong

def Omega(Y, n, t):
    n1 = np.arange(N) != n
    tong = np.dot(g[n, n1], Y[1][n1])
    Omega1 = (sqrt(pi) / (2 * delta_t)) * khi_o * exp(-(t / delta_t) ** 2)
    Omega2 = (sqrt(E_R) / (hbar * pi)) * delta_e * tong
    return Omega1 + Omega2

def F(t, Y):
    f = np.zeros((2, N), dtype='complex')
    for i in range(N):
        f[0][i] = -2 * (Omega(Y, i, t) * np.conj(Y[1][i])).imag + 1j* (-2) * (Omega(Y, i, t) * np.conj(Y[1][i])).imag
        f[1][i] = (-1j / hbar) * ((i + 1) * delta_e - delta_0 - En(Y, i)) * Y[1][i] \
                  + 1j * (1 - np.real(Y[0][i]) - np.imag(Y[0][i])) * Omega(Y, i, t) \
                  - Y[1][i] / T2
    return f

# Hàm RK4
def RK4(Ft, t, Y):
    k1 = dt * Ft(t, Y)
    k2 = dt * Ft(t + dt / 2, Y + k1 / 2)
    k3 = dt * Ft(t + dt / 2, Y + k2 / 2)
    k4 = dt * Ft(t + dt, Y + k3)
    return Y + (k1 + 2 * k2 + 2 * k3 + k4) / 6

# Hàm ghi kết quả vào file
def write_result(filename, steps, epsilon, Y0_real, Y0_imag ,Y1_abs, Nt, Pt):
    with open(filename, 'w') as f:
        f.write(f"#{'t':^20}{'epsilon':^20}{'Y0_real':^20}{'Y0_imag':^20}{'Y1_abs':^20}{'Nt':^20}{'Pt':^20}\n")
        for i, t in enumerate(steps):
            for j in range(N):
                f.write(f"{t:^20.5f}{epsilon[j]:^20.5f}{Y0_real[j, i]:^20.5e}{Y0_imag[j, i]:^20.5e}{Y1_abs[j, i]:^20.5e}{Nt[j, i]:^20.5e}{Pt[j, i]:^20.5e}\n")
            f.write('\n')

# Hàm vẽ kết quả
def plot(steps, epsilon, Y0_real, Y1_abs, Nt, Pt):
    X, Y = np.meshgrid(steps, epsilon)
    fig = plt.figure(figsize=(15, 10))

    ax1 = fig.add_subplot(221, projection='3d')
    ax1.plot_surface(X, Y, Y0_real, cmap='viridis')
    ax1.set_xlabel('t (femtosecond)')
    ax1.set_ylabel('epsilon (meV)')
    ax1.set_zlabel(r'$f_e$')
    ax1.set_title(r"Đồ thị 3D của $f_e$ theo t và epsilon")

    ax2 = fig.add_subplot(222, projection='3d')
    ax2.plot_surface(X, Y, Y1_abs, cmap='plasma')
    ax2.set_xlabel('t (femtosecond)')
    ax2.set_ylabel('epsilon (meV)')
    ax2.set_zlabel(r'$|p_n|$')
    ax2.set_title(r"Đồ thị 3D của $|p_n|$ theo t và epsilon")

    # Vẽ đồ thị 2D cho N(t) theo t
    ax3 = fig.add_subplot(223)
    ax3.plot(steps, Nt.mean(axis = 0), color='blue')
    ax3.set_xlabel('t (femtosecond)')
    ax3.set_ylabel('N(t)')
    ax3.set_title('Đồ thị 2D của N(t) theo t')

    # Vẽ đồ thị 2D cho |P(t)| theo t
    ax4 = fig.add_subplot(224)
    ax4.plot(steps, Pt.mean(axis = 0), color='red')
    ax4.set_xlabel('t (femtosecond)')
    ax4.set_ylabel('|P(t)|')
    ax4.set_title('Đồ thị 2D của |P(t)| theo t')


    plt.tight_layout()
    plt.show()

# Hàm chính
def main():
    Y = np.zeros((2, N), dtype='complex')

    steps = np.arange(t_0, t_max, dt)
    epsilon = [(j + 1) * delta_e for j in range(N)]
    Y0_real = np.zeros((N, len(steps)))
    Y0_imag = np.zeros((N, len(steps)))
    Y1_abs = np.zeros((N, len(steps)))
    Nt = np.zeros((N, len(steps)))
    Pt = np.zeros((N, len(steps)))

    start_time = time.time()

    # Cái hàm này Khôi Nguyên chỉ và nhìn cũng khá hay :)))
    for i, t in enumerate(tqdm(steps, desc="Processing", unit="step")):
        Y = RK4(F, t, Y)
        Y0_real[:, i] = Y[0, :N].real
        Y0_imag[:, i] = Y[0, :N].imag
        Y1_abs[:, i] = abs(Y[1, :N])
        Nt[:, i] =  np.sum(sqrt_n[:N] * Y[0, :N].real, axis=0) 
        Pt[:, i] =  abs(np.sum(delta_e * sqrt(delta_e) * sqrt_n[:N] *Y[1, :N], axis = 0))


    # Ghi dữ liệu vào file
    write_result('sbe1.txt', steps, epsilon, Y0_real, Y0_imag, Y1_abs, Nt, Pt)

    end_time = time.time()
    print(f"Thời gian chạy: {end_time - start_time:.2f} giây!!!")
    plot(steps, epsilon, Y0_real, Y1_abs, Nt, Pt)

if __name__ == "__main__":
    main()
