from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import accuracy_score

# loading data from iris
iris = load_iris()
X = iris.data   # feature   
Y = iris.target # class

# spliting data into training and testing 
X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.2, random_state=42)

# calling the model
clf = GaussianNB()

# training dataset
clf.fit(X_train, y_train)

# predicting the class
y_predict = clf.predict(X_test)

# evaluating the accuracy 
accuracy = accuracy_score(y_test, y_predict)
print(f"Độ chính xác của mô hình: {accuracy:.2f}")

# predicting with new data
new_data = [[5.1, 3.5, 1.4, 0.2]]  # Một mẫu mới để thử
prediction = clf.predict(new_data)
print(f"Dự đoán cho mẫu dữ liệu mới: {iris.target_names[prediction][0]}")