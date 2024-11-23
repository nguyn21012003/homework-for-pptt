import numpy as np
from math import pi, sqrt, log, exp
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from tqdm import tqdm

# Input variables
N = 100
a0 = 125  # armstrong
hbar = 658.5  # meV fs
E_max = 300  # meV
dE = E_max / N  # meV
d0 = 50  # Độ bơm để electron đi vào dải dẫn, d0 càng lớn electron đi vào càng sâu

delta_t = 50  # fs
dt = 2  # fs
tmax = 500  # fs
t0 = -3 * delta_t

T2 = 200  # 200, 500, 600 fs
E_R = 4.2
chi_0 = 2  # Tham số xung đưa vào [0.01, 0.1, 0.5, 1, 2]

E_g = 1420
omega_0 = -100
omega_max = 100


def RK4(F, t, fp, dt):
    k1 = dt * F(t, fp)
    k2 = dt * F(t + dt / 2, fp + k1 * dt / 2)
    k3 = dt * F(t + dt / 2, fp + k2 * dt / 2)
    k4 = dt * F(t + dt, fp + k3 * dt / 2)
    fp = fp + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return fp


# Def fj and fn equation
def fp_init(f_e, f_h, p_n):
    fp = np.zeros((2, N), dtype="complex")
    for n in range(N):
        fp[0][n] = f_e + 1j * f_h
        fp[1][n] = p_n
    return fp


# Define g function
def g(n, n1):
    g1 = (sqrt(n + 1) + sqrt(n1 + 1)) / (sqrt(n + 1) - sqrt(n1 + 1))
    return log(abs(g1)) / sqrt((n + 1) * dE)


# Define energy function
def E_n(n, fp):
    sum = 0
    for n1 in range(N):
        if n1 != n:
            sum += g(n, n1) * (np.real(fp[0][n1]) + np.imag(fp[0][n1]))
    return (sqrt(E_R) / pi) * dE * sum


# Define omega_R function
def omega_R(n, t, fp):
    sum = 0
    for n1 in range(N):
        if n1 != n:
            sum += g(n, n1) * fp[1][n1]
    return (1 / hbar) * (0.5 * (hbar * sqrt(np.pi) / delta_t) * chi_0 * exp(-(t**2) / delta_t**2)) + (1 / hbar) * ((sqrt(E_R) / np.pi) * dE * sum)


# Define N(t) and P(t) functions
def NP_func(fp):
    for n in range(N):
        N_t = (1 / 2 * np.pi**2) * (dE ** (3 / 2) / (E_R ** (3 / 2) * a0**3)) * np.sum(np.sqrt(n) * np.real(fp[0, :]))
        P_t = np.sum(dE * np.sqrt(dE) * np.sqrt(n) * fp[1, :])
    return N_t, P_t


# Define sbe function
def sbe(t, fp):
    F = np.zeros((2, N), dtype="complex")
    for n in range(N):  # Iterate over range N to stay within bounds
        F[0][n] = -2 * np.imag(omega_R(n, t, fp) * np.conjugate(fp[1][n]))
        F[1][n] = (
            -(1j / hbar) * (n * dE - d0 - E_n(n, fp)) * fp[1][n]
            + 1j * (1 - np.real(fp[0][n]) - np.imag(fp[0][n])) * omega_R(n, t, fp)
            - (fp[1][n] / T2)
        )
    return F


# Define function
def alpha_o(omega, t, P_t):
    E_omega = dt * np.sum(chi_0 * exp(-(t**2) / delta_t**2) * exp(1j * omega * t / hbar))
    P_omega = dt * np.sum(P_t * exp(1j * omega * t / hbar))
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
    plt.plot(omega_values, alpha_omega)
    plt.xlabel("Omega")
    plt.ylabel("alpha_omega")
    plt.title("Phổ hấp thụ")

    plt.tight_layout()
    plt.show()


def main():
    # Initial conditions
    f0_e, f0_h, p0 = 0, 0, 0
    fp = fp_init(f0_e, f0_h, p0)

    t_values = np.arange(t0, tmax, dt)
    E_values = np.arange(dE, E_max + dE, dE)
    results = []
    N_t = np.zeros(len(t_values))
    P_t = np.zeros(len(t_values))

    with open("sbe.txt", "w") as sbe_w:
        header = f"{'E':^10} {'t':^10} {'f_e':^20} {'f_h':^20} {'p_n':^20} {'N_t':^20} {'P_t':^20} {'alpha':^20} \n"
        sbe_w.write(header)

        for i, t in enumerate(tqdm(t_values)):

            fp = RK4(sbe, t, fp, dt)
            results.append(fp)

            N_t[i], P_t[i] = NP_func(fp)

            for n in range(N):
                line = f"{E_values[n]:<10.1f} {t:<10.1f} {np.real(fp[0][n]):<20.15f} {np.imag(fp[0][n]):<20.15f} {abs(fp[1][n]):<20.15f} {N_t[i]:<20.15f} {abs(P_t[i]):<20.15f} \n"
                sbe_w.write(line)
            sbe_w.write(" ")

    E_omega = np.zeros(1000, dtype="complex")
    P_omega = np.zeros(1000, dtype="complex")
    alpha_omega = np.zeros(1000)
    omega_values = np.linspace(omega_0, omega_max, 1000)

    with open("alpha.txt", "w") as alpha_w:
        alpha_w.write(f"{'omega':^20} {'E_omega':^20} {'P_omega':^20} {'alpha':^20} \n")

        for i in range(1000):
            E_omega[i], P_omega[i], alpha_omega[i] = alpha_o(omega_values[i], t_values, P_t)

            alpha_w.write(f"{omega_values[i]:<20.15f} {E_omega[i]:<20.15f} {P_omega[i]:<20.15f} {alpha_omega[i]:<20.15f} \n")

    f_e_real = np.array([np.real(fp[0]) for fp in results])
    p_values = np.array([abs(fp[1]) for fp in results])

    plot(E_values, t_values, f_e_real, p_values, N_t, P_t, omega_values, alpha_omega)


main()
