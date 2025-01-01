from sklearn.datasets import load_iris
import numpy as np
import matplotlib.pyplot as plt

# Gọi ra bộ dữ liệu
data = load_iris()
X = data.data
Y = data.target
target_names = data.target_names

# Bước 1: Chuẩn hóa dữ liệu (mean = 0)
X_mean = np.mean(X, axis=0)
X_centered = X - X_mean

# Bước 2: Tính ma trận hiệp phương sai
S = np.cov(X_centered.T)  # = (X_centered.T @ X_centered)/150

# Bước 3: Tính giá trị riêng và vector riêng
eigenvalues, eigenvectors = np.linalg.eig(S)

# Bước 4: Sắp xếp giá trị riêng theo thứ tự giảm dần
sorted_indices = np.argsort(eigenvalues)[::-1]  # Sắp xếp lại chỉ số theo giá trị giảm dần của trị riêng
eigenvalues = eigenvalues[sorted_indices]
eigenvectors = eigenvectors[:, sorted_indices]  # vc1 = vc[i,1]
X_centered = X_centered[:, sorted_indices]  # xep lai cho dac tinh bo dw lieu

# Bước 5: Chọn 3 thành phần chính
eigenvectors_reduced = eigenvectors[:, :3]  # Lấy 3 vector riêng đầu tiên

# Bước 6: Chuyển dữ liệu sang không gian mới
X_pca = X_centered @ eigenvectors_reduced  # 150x4 x 4x3 -> biểu diễn thông qua cơ sở là 3 vector riêng

# Vẽ biểu đồ 3D
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection="3d")

# Gán màu theo nhãn
colors = ["r", "g", "b"]
for color, i, target_name in zip(colors, [0, 1, 2], target_names):
    ax.scatter(X_pca[Y == i, 0], X_pca[Y == i, 1], X_pca[Y == i, 2], color=color, label=target_name, edgecolor="k", s=50)

ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
ax.set_zlabel("PC3")
ax.legend(loc="best", title="Loài hoa")
ax.set_title("PCA giảm từ 4 chiều xuống 3 chiều")

plt.show()
