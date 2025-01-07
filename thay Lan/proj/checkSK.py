#### Kiểm tra bằng scikitlearn
from sklearn import datasets, decomposition
import matplotlib.pyplot as plt
import numpy as np


iris = datasets.load_iris() 
X = iris.data
y = iris.target

pca = decomposition.PCA(n_components=3)
pca.fit(X)
X = pca.transform(X)

fig = plt.figure(figsize=(10, 8))
ax1 = fig.add_subplot(121, projection="3d")
ax2 = fig.add_subplot(122)
colors = np.array(['r', 'g', 'b'])

y = np.choose(y, [1, 2, 0]).astype(float)
ax1.scatter(X[:, 0], -1 * X[:, 1], X[:, 2], c=colors[(y - 1).astype(int)], s=50, edgecolor="k")
ax1.set_xlabel("Principal Component 1")
ax1.set_ylabel("Principal Component 2")
ax1.set_title("PCA of Iris Dataset with Sklearn") #### Giảm từ 4 chiều còn 3 chiều
ax2.scatter(X[:, 0], -1 * X[:, 1], c=colors[(y - 1).astype(int)], s=50, edgecolor="k" )
ax2.set_xlabel("Principal Component 1")
ax2.set_ylabel("Principal Component 2")
ax2.set_title("PCA of Iris Dataset with Sklearn")  #### Giảm từ 4 chiều còn 2 chiều

############# REFERENCES 
# https://scikit-learn.org/1.5/auto_examples/decomposition/plot_pca_iris.html

# Show the plot
plt.tight_layout()
plt.show()
