import numpy as np
from numpy import exp
import matplotlib.pyplot as plt
from math import sqrt, log, pi
from numpy import real as RE, imag as IM

import csv
import time
from tqdm import tqdm

# Định nghĩa các hằng số
N = 1000  # Số phương trình
hbar = 658.5  # meV.fs
E_R = 4.2  # meV
delta_t = 25  # femtosecond
delta_0 = 10  # meV (Độ bơm để electron đi sâu vào dải -> delta_0 càng lớn thì electron càng vào sâu bên trong dải dẫn!) = hbar*omega_0 - E_g  #### Để từ 10 - 50
t_max = 1000  # femtosecond
t_0 = -3 * delta_t  # femtosecond
e_max = 300  # meV
delta_e = e_max / N  # meV
T20 = 210  # femtosecond
dt = 2  # femtosecond
khi_o = 0.01  # Tham số cường độ xung (0.1 - 2)
a0 = 125  # Bán kính Bohr (Amstrong)
gamma = 6.5 * 1e-20  # cm^{3}.fs^{-1}
M = int((t_max - t_0) / dt)  # Số bước lặp RK4


# Tính C_0
C0 = (delta_e * sqrt(delta_e) * 1e24) / (2 * (pi**2) * E_R ** (3 / 2) * a0**3)
E0 = khi_o  ## biên độ của trường điện

#################################################################################

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


def F(t, Y, T2):
    f = np.zeros((2, N), dtype="complex")  ## Mảng 2 chiều
    for i in range(N):
        f[0][i] = -2 * (Omega(Y, i, t) * np.conj(Y[1][i])).imag

        f[1][i] = (
            (-1j / hbar) * ((i + 1) * delta_e - delta_0 - En(Y, i)) * Y[1][i]
            + 1j * (1 - np.real(Y[0][i]) - np.imag(Y[0][i])) * Omega(Y, i, t)
            - Y[1][i] / T2
        )
    return f


# Hàm RK4 tính f_e, p, N(t), P(t)
def rk4(dF, tn, yn, h, T2):
    # tn là giá trị trong steps thứ n
    k1 = dF(tn, yn, T2)
    k2 = dF(tn + 0.5 * h, yn + 0.5 * h * k1, T2)
    k3 = dF(tn + 0.5 * h, yn + 0.5 * h * k2, T2)
    k4 = dF(tn + h, yn + h * k3, T2)

    return yn + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


# Hàm chính
def main():
    Y = np.zeros((2, N), dtype="complex")  ### Để chứa nghiệm từ giải lặp RK4
    #### Y[0][N] sẽ là fe,n + i fh,n
    #### Y[1][N] sẽ là pn

    steps = np.arange(t_0, t_max, dt)  # Mảng thời gian
    omega = np.linspace(-100, 100, len(steps))  ### Chia đoạn omega để Fourier Tranform

    # print(steps)

    Mat_do_toan_phan_t = np.zeros(len(steps))  #### Để lưu giá trị hàm mật độ toàn phần
    Ham_phan_cuc_t = np.zeros(len(steps), dtype="complex")  #### Để lưu giá trị của hàm phân cực toàn phần chỗ này phải là REAL vì ta lấy abs của Pt
    E_dientruong = np.zeros(len(steps), dtype="complex")
    Ham_phan_cuc_omega = np.zeros(len(steps), dtype="complex")
    # PT = np.zeros((len(steps)))
    # T2_values = np.zeros_like(Nt)

    epsilon = [(n + 1) * delta_e for n in range(N)]  ### e = n * delta e
    T, E = np.meshgrid(steps, epsilon)
    # Thực hiện tính
    listFE_t = []
    listFH_t = []
    listPt_t = []

    T2 = 210
    for ti in tqdm(range(len(steps)), desc="Processing", unit="step ", ascii=" #"):
        ## Tại thời gian t = ti thì sẽ có 1 lần lặp RK4 và trả về 100 pt của Y. Thì listFE,FH,Pt có tác dụng là sẽ ghi thành phần thứ n trong 100 pt của Y
        Y = rk4(F, steps[ti], Y, dt, T2)
        # print(Y)
        listFE = []
        listFH = []
        listPt = []
        Mat_do_toan_phan = 0  ### Để tính tổng
        Ham_phan_cuc = 0  ### Để tính tổng
        for n in range(N):
            #### lấy ra thành phần thứ n của phương trình fe,fh,pt để tính toán mật độ toàn phần. Mật độ toàn phần (Mat_do_toan_phan) , Hàm phân cực (Ham_phan_cuc) được lấy tổng qua N = 100 phương trình.
            fe = RE(Y[0][n])
            fh = IM(Y[0][n])
            pt = Y[1][n]
            listFE.append(fe)
            listFH.append(fh)
            listPt.append(pt)
            Mat_do_toan_phan += C0 * sqrt(n) * fe  ### Cộng dồn chỗ Mật độ toàn phần = 0
            Ham_phan_cuc += delta_e * sqrt(delta_e) * sqrt(n) * pt  ### Cộng dồn chỗ Phân cực toàn phần = 0

        listFE_t.append(listFE)
        listFH_t.append(listFH)
        listPt_t.append(listPt)  ### Lưu dữ liệu để vẽ

        #### Sau khi lấy tổng xong thì phải lưu giá trị của Mật độ toàn phần (Nt), Hàm phân cực (Pt) tại thời điểm t = ti. Mật độ toàn phần sẽ mang giá trị thực, Hàm phân cực sau khi lấy suất sẽ là mang giá trị thực.
        Mat_do_toan_phan_t[ti] = Mat_do_toan_phan
        Ham_phan_cuc_t[ti] = Ham_phan_cuc
        T2 = 1 / (1 / 210 + gamma * Mat_do_toan_phan)
        # print(T2)

    listAlpha = []  ### Để chứa giá trị của hệ số hấp thụ

    for omega_i in range(len(omega)):  ### Lấy giá trị thứ i trong mảng omega
        ham_phan_cuc_iomega = 0  ### Để lấy tổng
        truong_dien_iomega = 0  ### Để lấy tổng
        for i in range(len(steps)):
            ham_phan_cuc_iomega += dt * Ham_phan_cuc_t[i] * np.exp(1j * omega[omega_i] * steps[i] / hbar)
            ### CỘng dồn giá trị của ham_phan_cuc_iomega

            ### Tại vì đang cần là P(omega) nên lúc này ta phải chia hbar xuống dưới trong thành phần exp để thu được exp(i omega t).
            truong_dien_iomega += E0 * dt * np.exp(1j * omega[omega_i] * steps[i] / hbar) * np.exp(-(steps[i] ** 2) / (delta_t**2))
            ### CỘng dồn giá trị của truong_dien_iomega, Biên độ điện trường
            #### Lúc này là có được E và P(ham_phan_cuc_iomega lúc này sẽ là số phức bởi vì hệ số exp(i omgega t))
        Alpha = IM(ham_phan_cuc_iomega / truong_dien_iomega)
        E_dientruong[omega_i] = truong_dien_iomega
        Ham_phan_cuc_omega[omega_i] = ham_phan_cuc_iomega
        ### Alpha sẽ lấy tổng của ham_phan_cuc_iomega truong_dien_iomega và chỉ lấy thành phần phức.
        listAlpha.append(Alpha)

    with open("sbe.dat", "w", newline="") as writefile:
        header = [
            "t",
            "RE_E",
            "IM_E",
            "RE_P",
            "IM_P",
            "Alpha",
        ]
        writer = csv.DictWriter(writefile, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for i in range(len(steps)):
            writer.writerow(
                {
                    "t": i,
                    "RE_E": RE(E_dientruong[i]),
                    "IM_E": IM(E_dientruong[i]),
                    "RE_P": RE(Ham_phan_cuc_omega[i]),
                    "IM_P": IM(Ham_phan_cuc_omega[i]),
                    "Alpha": listAlpha[i],
                }
            )
    with open("sbeOMEGA.dat", "w", newline="") as writefile:
        header = [
            "omega",
            "RE_E",
            "IM_E",
            "RE_P",
            "IM_P",
            "Alpha",
        ]
        writer = csv.DictWriter(writefile, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for omega_ii in range(len(omega)):
            writer.writerow(
                {
                    "omega": omega[omega_ii],
                    "RE_E": RE(E_dientruong[omega_ii]),
                    "IM_E": IM(E_dientruong[omega_ii]),
                    "RE_P": RE(Ham_phan_cuc_omega[omega_ii]),
                    "IM_P": IM(Ham_phan_cuc_omega[omega_ii]),
                    "Alpha": listAlpha[omega_ii],
                }
            )

    # np.array(listAlpha)  ### Cchuyển sang mảng của np để vẽ.
    plt.plot(omega, listAlpha)  ### Vẽ alpha
    # plt.show()

    # plot(steps, epsilon, Y0_real, Y1_abs, Nt, PT)


main()
