import matplotlib.pyplot as plt
import numpy as np

# Đọc dữ liệu từ file dataFT.txt
file = np.loadtxt("dataFT.txt")

# Lấy các cột dữ liệu từ file
time = file[:, 0]
alpha = file[:, 1]
dist = file[:, 2]
dist1 = file[:, 3]
polari = file[:, 4]
polari1 = file[:, 5]

# Tạo figure với 3 subplot
plt.figure(figsize=(10, 12))

# Subplot 1: Vẽ dist và dist1
plt.subplot(3, 1, 1)
plt.plot(time, dist, color='blue', linewidth=1, linestyle='solid', label="Real")
plt.plot(time, dist1, color='red', linewidth=1, linestyle='solid', label="Imag")
plt.xlabel("$\omega$")
plt.ylabel("$P(\omega)$")
plt.legend()
plt.grid()

# Subplot 2: Vẽ polari và polari1
plt.subplot(3, 1, 2)
plt.plot(time, polari, color='red', linewidth=2, linestyle='solid', label="Real")
plt.plot(time, polari1, color='green', linewidth=2, linestyle='solid', label="Imag")
plt.xlabel("$\omega$")
plt.ylabel("E($\omega$)")
plt.legend()
plt.grid()

# Subplot 3: Vẽ alpha
plt.subplot(3, 1, 3)
plt.plot(time, alpha, color='purple', linewidth=2, linestyle='solid')
plt.xlabel("$\omega$")
plt.ylabel("$\\alpha(\omega)$")
plt.grid()

# Lưu hình ảnh vào file
plt.tight_layout()
plt.savefig("combined_plot.png", dpi=200)
plt.show()
