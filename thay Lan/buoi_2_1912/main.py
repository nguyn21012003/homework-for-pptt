from sklearn.datasets import load_iris
from sklearn.naive_bayes import GaussianNB as GNB


model = GNB()

data_iris = load_iris()
data_iris.target[[10, 25, 50]]
list_iris_name = list(data_iris.target_names)

