# Thuật toán Bisection

## Giả sử f(x) liên tục trên [a,b] có nghĩa là f(x) có đạo hàm tại mọi điểm trong a,b
### Bước 1 tính c = (a+b) / 2
### Bước 2 Tính coi f(a) * f(c) < 0, f(b) * f(c) < 0 ### Điều kiện trái dấu của f(a) và f(c)
### Nếu f(a) * f(c) < 0 thì a(mới) = a(cũ) , b(mới) = c(cũ) hoặc f(b) * f(c) < 0 thì a(mới) = c(cũ) , b(mới) = b(cũ)
### Bước 3 lặp lại bước 1 nếu thoả điều kiện hội tụ tức có nghĩa là |f(c)| < epsilon thì dừng vòng lặp

### Code

from numpy import sin, cos, sqrt, exp
import csv
import numpy as np

a, b = 0.5, 5


def f1(x):
    return x - 2 ** (-x)


def f2(x):
    return exp(x) - 2 - cos(exp(x) - 2)


def bisection(f, a, b):
    N = 1000
    epsilon = 10 ** (-6)  # hoặc ghi như này 1e-6

    listC = np.zeros(N)
    listA = np.zeros(N)
    listB = np.zeros(N)

    listA[0] = a
    listB[0] = b

    for i in range(N):
        if i + 1 < N:
            listC[i] = (listA[i] + listB[i]) / 2  ## Step 1

            fa = f(listA[i])
            fb = f(listB[i])
            fc = f(listC[i])

            if fa * fc < 0:  ### Step 2
                listA[i + 1] = listA[i]
                listB[i + 1] = listC[i]
            elif fb * fc < 0:
                listA[i + 1] = listC[i]
                listB[i + 1] = listB[i]

            if abs(fc) < epsilon:
                break  ### Điều kiện hội tụ

    with open("file_ghi.txt", "w", newline="") as writefile:
        header = ["a", "b", "c"]
        writer = csv.DictWriter(writefile, fieldnames=header)
        writer.writeheader()
        for i in range(N):
            writer.writerow(
                {
                    "a": listA[i],
                    "b": listB[i],
                    "c": listC[i],
                }
            )

    return listC


# x1 = bisection(f1, a, b)
x2 = bisection(f2, a, b)
print(x2)
