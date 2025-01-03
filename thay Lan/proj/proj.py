import kagglehub
import pandas as pd
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from sklearn import datasets, decomposition


def mean(array):
    summ = 0
    for element in array:
        if type(element) == type(np.array([1])):
            summ += mean(element)
        else:
            summ += element

    mean = summ / len(array)

    return mean


def substrac_mean(data):

    mean_vector1 = mean(data["sepal_length"])
    mean_vector2 = mean(data["sepal_width"])
    mean_vector3 = mean(data["petal_length"])
    mean_vector4 = mean(data["petal_width"])

    substrac_meanarr = []

    for i in range(len(data)):
        substrac_meanarr.append(
            np.array(
                [
                    data["sepal_length"][i] - mean_vector1,
                    data["sepal_width"][i] - mean_vector2,
                    data["petal_length"][i] - mean_vector3,
                    data["petal_width"][i] - mean_vector4,
                ]
            )
        )
    return np.array(substrac_meanarr)


def covariant_matrix(N, data):
    substrac_meanarr = substrac_mean(data)

    S = (substrac_meanarr.T @ substrac_meanarr) / (N)

    return S


def DiagonalizeMatrx(matrix):
    eigenValues, eigenVecs = LA.eig(matrix)

    return eigenValues, eigenVecs


def pickK(eigenValues, eigenVecs):

    dim = len(eigenValues) - 1
    eigenValues = sorted(eigenValues, reverse=True)

    sorted_eigenValues = np.array(eigenValues[:dim])

    # ki hieu ":" trong numpy la lay het cac hang, cot, nhng ma "0:dim" la lay tu 0 den dim-1
    sorted_eigenVecs = np.array(eigenVecs[:, :dim])

    return sorted_eigenVecs


def projectedData(data, eigenVecs):
    return data @ eigenVecs


def main():
    iris = datasets.load_iris()  ### Load data from sklearn

    path = kagglehub.dataset_download("arshid/iris-flower-dataset")
    df = pd.read_csv(path + "/Iris.csv")

    sub_mearr = substrac_mean(df)

    S = covariant_matrix(len(sub_mearr), df)
    eigenValues, eigenVecs = DiagonalizeMatrx(S)
    new_eigenVecs = pickK(eigenValues, eigenVecs)
    projected_data = projectedData(sub_mearr, new_eigenVecs)

    fig = plt.figure(figsize=(10, 8))
    ax1 = fig.add_subplot(121, projection="3d")
    ax2 = fig.add_subplot(122, projection="3d")

    color_map = {"Iris-setosa": 0, "Iris-versicolor": 1, "Iris-virginica": 2}
    colors = []

    for species in df["species"]:
        colors.append(color_map[species])
    ax2.scatter(
        projected_data[:, 0],
        projected_data[:, 1],
        projected_data[:, 2],
        c=colors,
        s=50,
        edgecolors="k",
    )
    ax2.set_xlabel("Principal Component 1")
    ax2.set_ylabel("Principal Component 2")
    ax2.set_title("PCA of Iris Dataset")

    ############ check with Sklearn ############

    X = iris.data
    y = iris.target

    pca = decomposition.PCA(n_components=3)
    pca.fit(X)
    X = pca.transform(X)

    y = np.choose(y, [1, 2, 0]).astype(float)
    ax1.scatter(X[:, 0], -1 * X[:, 1], X[:, 2], c=colors, s=50, edgecolor="k")
    ax1.set_xlabel("Principal Component 1")
    ax1.set_ylabel("Principal Component 2")
    ax1.set_title("PCA of Iris Dataset with Sklearn")

    print(projected_data, X)

    # Show the plot
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
