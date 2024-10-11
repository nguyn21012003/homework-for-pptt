import numpy as np


def f(z):
    return np.tan(z) - np.sqrt((z0 / z) ** 2 - 1)


def df(z):
    return 1 / (np.cos(z)) ** 2 + z0 / (z**3 * np.sqrt((z0 / z) ** 2 - 1))


def newton_raphson(f, df, p0, e, N, file_name, step_limit=0.05):
    with open(file_name, "w") as file:
        file.write(f"n \t\t\t z\n")
        n = 1
        while n < N:
            p = p0 - f(p0) / df(p0)
            # Giới hạn bước nhảy ổn định chương trình nhưng cần nhiều thời gian hơn
            if abs(p - p0) > step_limit:
                p = p0 + np.sign(p - p0) * step_limit
            file.write(f"{n}\t\t\t{p:.15f}\n")
            if abs(p - p0) < e:
                break
            else:
                p0 = p
                n += 1
        return p


# Giá trị hằng số và tham số đầu vào
hbar = 1.055e-34  # hằng số Planck rút gọn
m = 9.11e-31  # khối lượng electron
a0 = 0.529e-10  # bán kính Bohr
a = 50 * a0
V0 = 100 * 1.6e-19
z0 = a * np.sqrt(2 * m * V0) / hbar
e = 1e-15
N = 100

# Nhập giá trị ban đầu
while True:
    p0 = float(input("Nhập giá trị khởi tạo = "))
    if ((z0 / p0) ** 2 - 1) < 0:
        print("Căn bậc hai không thể bị âm. Vui lòng nhập lại.")
    else:
        break

# Tên file xuất kết quả
file_name = "ketqua.txt"

# Tìm nghiệm
giai = newton_raphson(f, df, p0, e, N, file_name)
print(f"Giá trị nghiệm: {giai}")
