import numpy as np
from math import pi, sqrt, log, exp
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from tqdm import tqdm
import csv

# Input variables
N = 200
a0 = 125  # armstrong
hbar = 658.5  # meV fs
E_max = 300  # meV
dE = E_max / N  # meV
d0 = 10  # Độ bơm để electron đi vào dải dẫn, d0 càng lớn electron đi vào càng sâu

delta_t = 25  # fs
dt = 2  # fs
tmax = 1000  # fs
t0 = -3 * delta_t

# T2_0 = 210
E_R = 4.2
chi_0 = 0.01  # Tham số xung đưa vào
gamma = 6.5 * 1e-20  # cm^{3}.fs^{-1}

n_values = np.arange(1, N + 1)
sqrt_n = np.sqrt(np.arange(1, N + 1))
C0 = (dE * sqrt(dE) * 1e24) / (2 * (pi**2) * E_R ** (3 / 2) * a0**3)


t_values = np.arange(t0, tmax, dt)
E_values = np.arange(dE, E_max + dE, dE)

E_g = 1420
omega_0 = -200
omega_max = 200


# Def fj and fn equation
def fp_init(f_e, f_h, p_n):
    fp = np.zeros((2, N), dtype="complex")
    fp[0, :] = f_e + 1j * f_h
    fp[1, :] = p_n
    return fp


# Define g function
def g_func(N):
    g = np.zeros((N, N))
    for n in range(N):
        for n1 in range(N):
            if n1 != n:
                g[n, n1] = log(abs((sqrt_n[n] + sqrt_n[n1]) / (sqrt_n[n] - sqrt_n[n1]))) / sqrt((n + 1) * dE)
    return g


g = g_func(N)


# Define energy function
def E_n(fp, n):
    n1 = np.arange(N) != n
    sum = np.dot(g[n, n1], np.real(fp[0, n1]) + np.imag(fp[0, n1]))
    return sqrt(E_R) * dE * sum / pi


# Define omega_R function
def omega_R(fp, n, t):
    n1 = np.arange(N) != n
    sum = np.dot(g[n, n1], fp[1, n1])
    return (1 / hbar) * (0.5 * (hbar * sqrt(np.pi) / delta_t) * chi_0 * exp(-(t**2) / delta_t**2)) + (1 / hbar) * ((sqrt(E_R) / np.pi) * dE * sum)


# Define sbe function
def sbe(t, fp, T2):
    F = np.zeros((2, N), dtype="complex")
    for n in range(N):
        F[0, n] = -2 * np.imag(omega_R(fp, n, t) * np.conjugate(fp[1, n]))
        F[1, n] = (
            -(1j / hbar) * (E_values[n] - d0 - E_n(fp, n)) * fp[1, n]
            + 1j * (1 - np.real(fp[0, n]) - np.imag(fp[0, n])) * omega_R(fp, n, t)
            - (fp[1, n] / T2)
        )
    return F


def RK4(F, t, fp, T2):
    k1 = dt * F(t, fp, T2)
    k2 = dt * F(t + dt / 2, fp + k1 * dt / 2, T2)
    k3 = dt * F(t + dt / 2, fp + k2 * dt / 2, T2)
    k4 = dt * F(t + dt, fp + k3 * dt / 2, T2)
    fp = fp + (k1 + 2 * k2 + 2 * k3 + k4) / 6

    # N_t = C0 * np.sum(sqrt_n[:] * np.real(fp[0, :]))

    # T2 = T2_0 / (1 + T2_0 * gamma * N_t)

    # P_t = np.sum(dE * np.sqrt(dE) * sqrt_n[:] * fp[1, :])

    return fp


# Define function
def alpha_o(omega, P_t):
    E_omega = dt * np.sum(chi_0 * (hbar * np.pi / delta_t) * np.exp(-(t_values**2) / delta_t**2) * np.exp(1j * omega * t_values / hbar))
    P_omega = dt * np.sum(P_t * np.exp(1j * omega * t_values / hbar))
    alpha = np.imag(P_omega / E_omega)
    return E_omega, P_omega, alpha


def plot(E_values, t_values, f_e_real, p_values, N_t, P_t, omega_values, alpha_omega):
    E, t = np.meshgrid(E_values, t_values)

    fig1 = plt.figure(figsize=(12, 6))
    fig2 = plt.figure(figsize=(12, 6))

    ax1 = fig1.add_subplot(121, projection="3d")
    ax2 = fig1.add_subplot(122, projection="3d")
    ax3 = fig2.add_subplot(121)
    ax4 = fig2.add_subplot(122)

    # 3D plot for f_e
    ax1.plot_surface(E, t, f_e_real, cmap="plasma")
    ax1.set_title(f"f_e plot with d0 = {d0}")
    ax1.set_xlabel("Energy (meV)")
    ax1.set_ylabel("Time (fs)")
    ax1.set_zlabel("f_e")

    # 3D plot for p_n
    ax2.plot_surface(E, t, p_values, cmap="plasma")
    ax2.set_title(f"p_n plot with d0 = {d0}")
    ax2.set_xlabel("Energy (meV)")
    ax2.set_ylabel("Time (fs)")
    ax2.set_zlabel("p_n")

    # 2D plot for N_t
    ax3.plot(t, N_t)
    ax3.set_title(f"N_t plot with d0 = {d0}")
    ax3.set_xlabel("Time (fs)")
    ax3.set_ylabel("N(t)")

    # 2D plot for P_t
    ax4.plot(t, abs(P_t))
    ax4.set_title(f"P_t plot with d0 = {d0}")
    ax4.set_xlabel("Time (fs)")
    ax4.set_ylabel("P(t)")

    plt.figure(figsize=(8, 6))
    plt.plot(omega_values, alpha_omega, label=f"d0 = {d0}, chi_0 = {chi_0}")
    plt.xlabel("Omega (rad/s)")
    plt.ylabel("Alpha (a.u.)")
    plt.title("Phổ hấp thụ")

    plt.legend()
    plt.tight_layout()
    plt.show()


def main():
    # Initial conditions
    f0_e, f0_h, p0 = 0, 0, 0
    fp = fp_init(f0_e, f0_h, p0)

    filedata = f"sbe with d0 = {d0} and chi_0 = {chi_0} T2thaydoi.txt"
    results = []
    N_t = np.zeros(len(t_values))
    P_t = np.zeros(len(t_values), dtype=complex)

    T2_0 = 210
    with open(filedata, "w") as sbe_w:
        # header = f"{'E':^10} {'t':^10} {'f_e':^20} {'f_h':^20} {'p_n':^20} {'N_t':^20} {'P_t':^20}\n" ### Chỗ này mình sửa lại một xíu
        header = ["t", "Nt_t"]
        writer = csv.DictWriter(sbe_w, fieldnames=header, delimiter="\t")
        writer.writeheader()
        # sbe_w.write(header)
        for ti in tqdm(range(len(t_values))):
            fp = RK4(sbe, t_values[ti], fp, T2_0)
            Nt_sum = 0
            Pt_sum = 0
            for n in range(N):
                fe = np.real(fp[0][n])
                #    print(fe)
                fh = np.imag(fp[0][n])
                pt = fp[1][n]
                Nt_sum += C0 * sqrt(n) * fe
                Pt_sum += dE * sqrt(dE) * sqrt(n) * pt

            N_t[ti] = Nt_sum
            P_t[ti] = Pt_sum
            T2_0 = 1 / (1 / 210 + gamma * N_t[ti])
        for i in range(len(t_values)):
            writer.writerow({"t": t_values[i], "Nt_t": N_t[i]})

            ############ Từ sau khúc này là mình ko có đụng


#            results.append(fp)
#            for n in range(N):
#                line = f"{E_values[n]:<10.5f} {t_values[ti]:<10.5f} {np.real(fp[0][n]):<20.15f} {np.imag(fp[0][n]):<20.15f} {abs(fp[1][n]):<20.15f} {N_t[ti]:<20.15f} {P_t[ti]:<20.15f} \n"
#                sbe_w.write(line)
#            sbe_w.write(" ")
#
#        E_omega = np.zeros(1000, dtype="complex")
#        P_omega = np.zeros(1000, dtype="complex")
#        alpha_omega = np.zeros(1000)
#        omega_values = np.linspace(omega_0, omega_max, 1000)


#    with open(f"absorption with d0 = {d0} and chi_0 = {chi_0} T2thaydoi.txt", "w") as alpha_w:
#        alpha_w.write(f"{'omega':^20} {'E_omega':^20} {'P_omega':^20} {'alpha':^20} \n")
#
#        for i in range(1000):
#            E_omega[i], P_omega[i], alpha_omega[i] = alpha_o(omega_values[i], P_t)
#
#            alpha_w.write(f"{omega_values[i]:<20.5f} {E_omega[i]:<20.5f} {P_omega[i]:<20.5f} {alpha_omega[i]:<20.5f} \n")
#
#    f_e_real = np.array([np.real(fp[0]) for fp in results])
#    p_values = np.array([abs(fp[1]) for fp in results])
#
# plot (E_values, t_values, f_e_real, p_values, N_t, P_t, omega_values, alpha_omega)


main()
