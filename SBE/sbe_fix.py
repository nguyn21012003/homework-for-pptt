import numpy as np
from math import pi, sqrt, log, exp
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from tqdm import tqdm

# Input variables
a0 = 125  # armstrong
hbar = 658.5  # meV fs
N = 100
E_max = 300  # meV
dE = E_max / N
d0 = 100  # meV

delta_t = 20  # fs
dt = 2  # fs
tmax = 500  # fs
t0 = -3 * delta_t

T2 = 200  # fs
E_R = 4.2
chi_0 = 0.5




def RK4(F, t, fp, dt, omega_R, E_n):
    k1 = dt * F(t, omega_R, E_n, fp)
    k2 = dt * F(t + dt / 2, omega_R, E_n, fp + k1 * dt / 2)
    k3 = dt * F(t + dt / 2, omega_R, E_n, fp + k2 * dt / 2)
    k4 = dt * F(t + dt, omega_R, E_n, fp + k3 * dt / 2)
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
    return (1 / sqrt((n + 1) * dE)) * log(abs(g1))


# Define energy function
def E_n(n, g, fp):
    E_n = 0
    for n1 in range(N):
        if n1 != n:
            E_n += (sqrt(E_R) / pi) * dE * g(n, n1) * (np.real(fp[0][n1]) + np.imag(fp[0][n1]))
    return E_n


# Define omega_R function
def omega_R(n, t, g, fp):
    omega_R = (1 / hbar) * (0.5 * (hbar * sqrt(pi) / delta_t) * chi_0 * exp(-(t**2) / delta_t**2))
    for n1 in range(N):
        if n1 != n:
            omega_R += (1 / hbar) * ((sqrt(E_R) / pi) * dE * g(n, n1) * fp[1][n1])
    return omega_R


# Define sbe function
def sbe(t, omega_R, E_n, fp):
    F = np.zeros((2, N), dtype="complex")
    for n in range(N):  # Iterate over range N to stay within bounds
        F[0][n] = -2 * np.imag(omega_R * np.conjugate(fp[1][n]))
        F[1][n] = -(1j / hbar) * (n * dE - d0 - E_n) * fp[1][n] + 1j * (1 - np.real(fp[0][n]) - np.imag(fp[0][n])) * omega_R - (fp[1][n] / T2)
    return F


def main():
    # Initial conditions
    f0_e, f0_h, p0 = 0, 0, 0
    fp = fp_init(f0_e, f0_h, p0)

    results = []
    nt = int(tmax - t0 / dt)
    print(nt)

    with open("sbe.txt", "w") as sbe_w:
        s0 = "{0:^10} {delim} {1:^10} {delim} {2:^20} {delim} {3:^10} {delim} {4:^20} \n"
        sbe_w.write(s0.format("E", "t", "f_e", "f_h", "p_n", delim="|"))

        for i in tqdm(range(nt)):
            t = t0 + i * dt

            E_n_values = np.zeros(N)
            omega_R_values = np.zeros(N, dtype="complex")
            E_values = np.zeros(N)

            # Calculate E_n and omega_R at each time step
            for n in range(N):
                E_n_values[n] += E_n(n, g, fp)
                omega_R_values[n] += omega_R(n, t, g, fp)
                E_values[n] = (n + 1) * dE

            # Update fp using RK4
            fp = RK4(sbe, t, fp, dt, omega_R_values[n], E_n_values[n])

            results.append(fp)
            for n in range(N):
                s1 = "{0:^10} {delim} {1:^10} {delim} {2:^20.15f} {delim} {3:^10} {delim} {4:^20.15f} \n"
                sbe_w.write(s1.format(E_values[n], t0 + i * dt, np.real(fp[0][n]), np.imag(fp[0][n]), fp[1][n], delim="|"))

    # Generate the 2D grids for E and t using meshgrid
    E_values = np.array([(n + 1) * dE for n in range(N)])  # energy levels
    t_values = np.linspace(t0, tmax, nt)  # time steps
    f_e_real = np.array([np.real(fp[0]) for fp in results])
    p_values = np.array([np.real(fp[1]) for fp in results])

    E_grid, t_grid = np.meshgrid(E_values, t_values)

    fig = plt.figure(figsize=(12, 6))
    ax1 = fig.add_subplot(121, projection="3d")
    ax2 = fig.add_subplot(122, projection="3d")
    ax1.grid(), ax2.grid()

    ax1.plot_surface(E_grid, t_grid, f_e_real, cmap="viridis")
    ax1.set_title("f_e plot")

    ax1.set_xlabel("Energy (meV)", labelpad=20)
    ax1.set_ylabel("Time (fs)", labelpad=20)
    ax1.set_zlabel("f_e", labelpad=20)

    ax2.plot_surface(E_grid, t_grid, p_values, cmap="viridis")
    ax2.set_title("p_n plot")

    ax2.set_xlabel("Energy (meV)", labelpad=20)
    ax2.set_ylabel("Time (fs)", labelpad=20)
    ax2.set_zlabel("p_n", labelpad=20)

    plt.show()


main()
