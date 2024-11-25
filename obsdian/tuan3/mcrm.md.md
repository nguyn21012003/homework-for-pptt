- [[#Source code|Source code]]
		- [[#Các hàm $f(...)$ cần tính|Các hàm $f(...)$ cần tính]]
		- [[#Hàm để giải bằng Monte - Carlo 1D (MC 1D):|Hàm để giải bằng Monte - Carlo 1D (MC 1D):]]
		- [[#Hàm để giải bằng phương pháp lấy tổng Riemann 1 lớp (RM 1D)|Hàm để giải bằng phương pháp lấy tổng Riemann 1 lớp (RM 1D)]]
		- [[#Hàm để giải bằng phương pháp MC 2D và RM 2D|Hàm để giải bằng phương pháp MC 2D và RM 2D]]
		- [[#Hàm để giải bằng phương pháp MC 3D, RM 3D và MC 10D và RM 10D|Hàm để giải bằng phương pháp MC 3D, RM 3D và MC 10D và RM 10D]]
- [[#Tính toán song song|Tính toán song song]]

# Các bài toán sử dụng phương pháp Monte - Carlo và tổng Riemann

- Tích phân 1 lớp của hàm $f(x) = x^2$

- Tích phân 2 lớp của hàm $f(x,y) = \cos(x^4) + 3 * y^2$

- Tích phân 3 lớp của hàm $f(x,y,z) = x^2 + y^2 + z^2$

- Tích phân 10 lớp của hàm $f(x_1,x_2,...x_{10}) = (x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8 + x_9 + x_{10}) ^ 2$

# Các phương pháp sử dụng

- Monte - Carlo 1 chiều, 2 chiều, 3 chiều

- Tổng Riemann 1 chiều, 2 chiều , 3 chiều

## Source code
#### Các hàm $f(...)$ cần tính

$f(x) = x^2$

```python
def f(x):
    return x**2
```

$f(x,y) = \cos(x^4) + 3 * y^2$

```python
def f1(x, y):
    return cos(x**4) + 3 * y**2
```

$f(x,y,z) = x^2 + y^2 + z^2$

```python
def f3(x, y, z):
    return x**2 + y**2 + z**2
```

$f(x_1,x_2,...x_{10}) = (x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8 + x_9 + x_{10}) ^ 2$

```python
def f4(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10):
    return (x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10) ** 2
```

#### Hàm để giải bằng Monte - Carlo 1D (MC 1D):

```python
def monte(f, a, b, N, eps):
    S = np.zeros(N)
    S[0] = 0
    ans = 0
    list_x = np.zeros(N)
    list_y = np.zeros(N)
    for i in range(N):
        x = uniform(a, b)
        y = uniform(f(0), f(x))
        S[i] = f(x)
        ans += (b - a) * S[i] / N
        list_x[i] = x
        list_y[i] = y
    return ans, list_x, list_y
```

Mô tả: 

Ta có $\int_a^b f(x)dx = S$, chia $[a,b]$ thành $N$ đoạn $(h \equiv \frac{b-a}{N})$. Ta có thể viết lại $S$ thành: 
$$
S = \frac{b-a}{N} \sum_n^N f(x_n)
$$
trong đó $x$ được gieo ngẫu nhiên trong $[a,b]$, $y$ được lấy ngẫu nhiên trong đoạn $[f(0),f(x)]$.
Ta  phải để `y = uniform(f(0), f(x))` ngẫu nhiên là vì: ta đang cần tính toán phần diện tích bên dưới đường cong $f(x)$, nếu đặt `y = uniform(a, b)` thì phần điểm gieo ngẫu nhiên, thì "viên sỏi" có thể rơi ra phía trên đường cong.
Kết quả trả về sẽ là 1 `tuple` chứa nghiệm của tích phân dưới dạng `float` , 2 list của x và y để vẽ ra hình.

#### Hàm để giải bằng phương pháp lấy tổng Riemann 1 lớp (RM 1D)

```python 
def riemann(f:float, x0:float, a:float, b:float, N:int, eps:float)->tuple:
    h = (b - a) / N
    S = np.zeros(N)
    x = np.zeros(N)
    S[0] = 0
    ans = 0
    for i in range(0, N):
        x[i] = a + i * h
        S[i] = f(x[i])
        ans += h * S[i]
    return ans, x, h
```

Mô tả:
Đặt $\triangle x \equiv h = \frac{b-a}{N}$. Công thức tính tổng Riemann cho tích phân 1 lớp là:
$$
S = \triangle x \sum_{i=1}^{N} f(x_i) \equiv h \sum_{i=1}^{N} f(x_i)
$$
Kết quả trả về sẽ là 1 `tuple` chứa `ans:float`, `x:dict` , `h:float`. Ta không cần `return y`  , thay vào đó ta `return h` để khi plot thì thỏa $x_n - x_{n-1} = h$.

#### Hàm để giải bằng phương pháp MC 2D và RM 2D

Monte - Carlo:

```python
def monte2D(f, a, b, c, d, N, eps):
    S = np.zeros(N)
    S[0] = 0
    ans = 0
    list_x = np.zeros(N)
    list_y = np.zeros(N)
    for i in range(N):
        x = uniform(a, b)
        y = uniform(c, d)
        S[i] = f(x, y)
        ans += (b - a) * (d - c) * S[i] / N
        list_x[i] = x
        list_y[i] = y
    return ans, list_x, list_y
```

Ở Riemann, công thức tổng 2 lớp tích phân:
$$
S = \triangle x \triangle y \sum_{i=1}^{N}\sum_{i=j}^{N} f(x_i,y_i) \equiv hx\, hy \sum_{i=1}^{N} \sum_{j=1}^{N} f(x_i,y_i)
$$
Riemann: 

```python
def riemann2D(f, a, b, c, d, N, eps):
    ans = 0
    hx = (b - a) / N
    hy = (d - c) / N
    S = np.zeros(N)
    S[0] = 0
    x = np.zeros(N)
    y = np.zeros(N)
    x[0] = a
    y[0] = c
    for i in range(0, N):
        x[i] = a + i * hx
        for j in range(1, N):
            y[j] = c + j * hy
            S[i] = f(x[i], y[j])
            ans += hx * hy * S[i]
    return ans, x, y, hx, hy
```

Do $i,j$ đang là 2 tổng chạy độc lập cho nên vòng lặp `for j` cần phải lồng bên trong vòng lặp `for i`. 

Tương tự ta có cho MC 3D, RM 3D và MC 10D và RM 10D

#### Hàm để giải bằng phương pháp MC 3D, RM 3D và MC 10D và RM 10D

```python
def monte3D(f, a, b, c, d, g, h, N, eps):
    S = np.zeros(N)
    S[0] = 0
    ans = 0
    list_x = np.zeros(N)
    list_y = np.zeros(N)
    list_z = np.zeros(N)
    for i in range(N):
        x = uniform(a, b)
        y = uniform(c, d)
        z = uniform(g, h)
        S[i] = f(x, y, z)
        ans += (b - a) * (d - c) * (h - g) * S[i] / N
        list_x[i] = x
        list_y[i] = y
        list_z[i] = z
    return ans, list_x, list_y, list_z
```

```python
def monte10D(fx, a, b, N, eps):
    S = 0
    ans = 0
    for i in range(N):
        x1 = uniform(a, b)
        x2 = uniform(a, b)
        x3 = uniform(a, b)
        x4 = uniform(a, b)
        x5 = uniform(a, b)
        x6 = uniform(a, b)
        x7 = uniform(a, b)
        x8 = uniform(a, b)
        x9 = uniform(a, b)
        x10 = uniform(a, b)
        np.random.rand(100)
        S = fx(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
        ans += (b - a) ** 10 * S / N
    return ans
```

```python 
def riemann3D(f, x0, y0, z0, a, b, c, d, g, h, N, eps):
    ans = 0
    hx = (b - a) / N
    hy = (d - c) / N
    hz = (h - g) / N
    S = np.zeros(N)
    S[0] = 0
    x = np.zeros(N)
    y = np.zeros(N)
    z = np.zeros(N)
    x[0] = x0
    y[0] = y0
    z[0] = z0
    for i in range(1, N):
        x[i] = x0 + i * hx
        for j in range(1, N):
            y[j] = y0 + j * hy
            for k in range(1, N):
                z[k] = z0 + k * hz
                S[i] = f(x[i], y[j], z[k])
                ans += hx * hy * hz * S[i]
    return ans, x, y, z, hx, hy, hz
```

```python
def riemann10D(f, x1_0, x2_0, x3_0, x4_0, x5_0, x6_0, x7_0, x8_0, x9_0, x10_0, a, b, N, eps):
    ans = 0
    h = (b - a) / N
    S = np.zeros(N)
    S[0] = 0
    x1 = np.zeros(N)
    x2 = np.zeros(N)
    x3 = np.zeros(N)
    x4 = np.zeros(N)
    x5 = np.zeros(N)
    x6 = np.zeros(N)
    x7 = np.zeros(N)
    x8 = np.zeros(N)
    x9 = np.zeros(N)
    x10 = np.zeros(N)
    x1[0] = x1_0
    x2[0] = x2_0
    x3[0] = x3_0
    x4[0] = x4_0
    x5[0] = x5_0
    x6[0] = x6_0
    x7[0] = x7_0
    x8[0] = x8_0
    x9[0] = x9_0
    x10[0] = x10_0
    for k1 in range(1, N):
        x1[k1] = x1_0 + k1 * hx
        for k2 in range(1, N):
            x2[k2] = x2_0 + k2 * hx
            for k3 in range(1, N):
                x3[k3] = x3_0 + k3 * hx
                for k4 in range(1, N):
                    x4[k4] = x4_0 + k4 * hx
                    for k5 in range(1, N):
                        x5[k5] = x5_0 + k5 * hx
                        for k6 in range(1, N):
                            x6[k6] = x6_0 + k6 * hx
                            for k7 in range(1, N):
                                x7[k7] = x7_0 + k7 * hx
                                for k8 in range(1, N):
                                    x8[k8] = x8_0 + k8 * hx
                                    for k9 in range(1, N):
                                        x9[k9] = x9_0 + k9 * hx
                                        for k10 in range(1, N):
                                            x10[k10] = x10_0 + k10 * hx
                                            S[k1] = f(x1[k1], x2[k2], x3[k3], x4[k4], x5[k5], x6[k6], x7[k7], x8[k8], x9[k9], x10[k10])
                                            ans += hx**10 * S[k1]
    return ans
```

Có thể thấy tích phân càng nhiều lớp thì phương pháp lấy tổng Riemann sẽ tốn rất nhiều thời gian, công sức hơn là phương pháp Monte - Carlo.

### Hàm vẽ đồ thị cho MC 1D và RM 1D
```python
def plot(f, a, b, xmonte1D, ymonte1D, xrieman1D, h, N):
    # parameters
    transparent = 0.5
    x = sym.symbols("x")
    fx = sym.Integral(x**2)
    ##################################################################
    # cài đặt 1 số thứ cho plot
    fig, axs = plt.subplots(2, 2, figsize=(18, 7))
    ##################################################################
    # nghiệm của phương trình f(x) dc biểu diễn trên đoạn [0,5]
    X = np.linspace(a, b, N)
    Y = f(X)
    ##################################################################
    axs[0, 0].plot(X, Y, color="r", label=r"y=$x^2$", lw=2)
    axs[0, 0].scatter(xmonte1D, ymonte1D, marker="o", s=10)
    axs[0, 0].fill_between(X, Y, color="lightgreen", alpha=transparent)

    axs[0, 0].legend([r"y=$x^2$", "sỏi", f"${latex(fx)}$"])
    axs[0, 0].set_title(r"Monte-Carlo cho phương trình $y = x^2$")
    axs[0, 0].set_xlabel(r"x")
    axs[0, 0].set_ylabel(r"$f(x)$", rotation=0)
    ##################################################################
    axs[0, 1].plot(X, Y, color="r", label=r"y=$x^2$", lw=2)
    axs[0, 1].plot(xrieman1D[:-1], f(xrieman1D[:-1]), "b.")
    axs[0, 1].bar(xrieman1D[:-1], f(xrieman1D[:-1]), width=h, alpha=transparent, align="edge", edgecolor="b")

    axs[0, 1].legend([r"y=$x^2$"])
    axs[0, 1].set_title(r"Phuwong pháp tổng Riemann(trái) cho phương trình $y = x^2$")
    axs[0, 1].set_xlabel(r"x")
    axs[0, 1].set_ylabel(r"$f(x)$", rotation=0)
    ##################################################################
    xrieman_mid = (xrieman1D[1:] + xrieman1D[:-1]) * 0.5
    yrieman_mid = f(xrieman_mid)
    axs[1, 0].plot(X, Y, color="r", label=r"y=$x^2$", lw=2)
    axs[1, 0].plot(xrieman_mid, yrieman_mid, "b.")
    axs[1, 0].bar(xrieman_mid, yrieman_mid, width=h, alpha=transparent, edgecolor="b")
    
    axs[1, 0].legend([r"y=$x^2$"])
    axs[1, 0].set_title(r"Phương pháp tổng Riemann(giữa) cho phương trình $y = x^2$")
    axs[1, 0].set_xlabel(r"x")
    axs[1, 0].set_ylabel(r"$f(x)$", rotation=0)
    ##################################################################
    axs[1, 1].plot(X, Y, color="r", label=r"y=$x^2$", lw=2)
    axs[1, 1].plot(xrieman1D[1:], f(xrieman1D[1:]), "b.")
    axs[1, 1].bar(xrieman1D[1:], f(xrieman1D[1:]), width=-h, alpha=transparent, align="edge", edgecolor="b")
    
    axs[1, 1].legend([r"y=$x^2$"])
    axs[1, 1].set_title(r"Phương pháp tổng Riemann(phải) cho phương trình $y = x^2$")
    axs[1, 1].set_xlabel(r"x")
    axs[1, 1].set_ylabel(r"$f(x)$", rotation=0)
    ##################################################################
    fig.savefig(f"Cau1.Monte-Carlo-Riemann1D.pdf")
    plt.show()
```

Mô tả:


## Tính toán song song

Để tiết kiệm được thời gian cho phương pháp Riemann nhất có thể, ta sẽ dùng thêm thuật toán tính toán song song ( Multi processing hoặc là Parallel processing ). Có nhiều thư viện để tính toán song song, như là Cython, Numba, Pythran,...

Ở đây mình dùng Numba và function decorator có trong python.

- Ở mỗi trước hàm $f(...)$ cần tính toán song song, ta thêm decorator `@njit` từ `numba`. Ví dụ như tính toán 3 lớp tích phân bằng phương pháp tích phân Monte - Carlo và phương pháp lấy tổng Riemann.

 ```python
@njit
def f3(x, y, z):
    return x**2 + y**2 + z**2
    
@njit
def f4(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10):
    return (x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10) ** 2    
```

```python
@njit(parallel=True)
def monte3D(f, a, b, c, d, g, h, N, eps):
    S = np.zeros(N)
    S[0] = 0
    ans = 0
    list_x = np.zeros(N)
    list_y = np.zeros(N)
    list_z = np.zeros(N)
    for i in nb.prange(N):
        x = uniform(a, b)
        y = uniform(c, d)
        z = uniform(g, h)
        S[i] = f(x, y, z)
        ans += (b - a) * (d - c) * (h - g) * S[i] / N
        list_x[i] = x
        list_y[i] = y
        list_z[i] = z
    return ans, list_x, list_y, list_z
```

```python
@njit(parallel=True)
def riemann3D(f, x0, y0, z0, a, b, c, d, g, h, N, eps):
    ans = 0
    hx = (b - a) / N
    hy = (d - c) / N
    hz = (h - g) / N
    S = np.zeros(N)
    S[0] = 0
    x = np.zeros(N)
    y = np.zeros(N)
    z = np.zeros(N)
    x[0] = x0
    y[0] = y0
    z[0] = z0
    for i in nb.prange(1, N):
        x[i] = x0 + i * hx
        for j in range(1, N):
            y[j] = y0 + j * hy
            for k in range(1, N):
                z[k] = z0 + k * hz
                S[i] = f(x[i], y[j], z[k])
                ans += hx * hy * hz * S[i]
    return ans, x, y, z, hx, hy, hz
```

Bên cạnh đó, ở vòng lặp for i đầu tiên, ta cũng phải thay `range(1, N)` thành `nb.prange(1, N)`. Mục đích là tính `range` theo `Numba`.

## Kết quả
Ở câu 1, câu 2, thuật toán chưa yêu cầu tính toán quá nặng nên thời gian tính toán của 2 phương pháp gần như bằng nhau.
### Câu 1:

Để nhìn rõ đồ thị hơn, ta set up `N = 10`

```text
Nhập số vòng lặp N(nhập N từ 10->100 sẽ cho ra kết quả nhanh hơn, từ 1000->10000 kết quả sẽ 'smooth' hơn): 10

Nghiệm tính theo phương pháp Monter-Carlo là 43.10278365811478 

Nghiệm tính theo phương pháp tổng Rieman là 35.625
```

![N=10](c1.mcrm1dN10.png)

Để kết quả chính xác hơn, ta set up `N = 100000`

```text
Nhập số vòng lặp N(nhập N từ 10->100 sẽ cho ra kết quả nhanh hơn, từ 1000->10000 kết quả sẽ 'smooth' hơn): 100000  

Nghiệm tính theo phương pháp Monter-Carlo là 41.65246168054039 

Nghiệm tính theo phương pháp tổng Rieman là 41.66604166875005
```

![N=100000](c1.mcrm1dN100000.png)

### Câu 2,3:

```text
Nhập số vòng lặp N(nhập N từ 10->100 sẽ cho ra kết quả nhanh hơn, từ 1000->10000 kết quả sẽ 'smooth' hơn): 1000

Nghiệm tính theo phương pháp Monter-Carlo là 2.0682555700910847

Nghiệm tính theo phương pháp tổng Rieman là 2.0017167065262487
```

### Câu 4,5: Tính toán song song

```text
Nhập số vòng lặp N(nhập N từ 10->100 sẽ cho ra kết quả nhanh hơn, từ 1000->10000 kết quả sẽ 'smooth' hơn): 1000

Nghiệm tính theo phương pháp Monter-Carlo là 0.9979397147369866 với thời gian là 4.437032222747803s

Nghiệm tính theo phương pháp tổng Rieman là 0.9965044975003308 với thời gian là 2.1190381050109863s
```

Có thể thấy, khi dùng tính toán song song, thời gian tính toán nghiệm gần như nhanh gấp đôi với `N =1000`

# Toàn bộ Source code
```python
import numpy as np

from random import random, uniform

from math import sqrt, pi

from numpy import sin, cos

import sympy as sym

from sympy.printing import latex

import matplotlib.pyplot as plt

import time

from numba import njit

import numba as nb

  
  

nb.set_num_threads(6)

  
  

def promp_user(id: int):

    N = int(input("Nhập số vòng lặp N(nhập N từ 10->100 sẽ cho ra kết quả nhanh hơn, từ 1000->10000 kết quả sẽ 'smooth' hơn): "))

    eps = sqrt(N)

    match id:

        case 1:

            ##################################################################################

            # Câu 1: tính tích phân x = [0,5]  của x^2 và vẽ minh họa số sỏi ném trong miền tạo bởi f(x) và trục tung

            a, b, c, d = 0, 5, 0, 1

            ans, list_x, list_y = monte(f, a, b, N, eps)

            ans1, list_x1, h1 = riemann(f, a, b, N, eps)

            print(f"Nghiệm tính theo phương pháp Monter-Carlo là {ans} \n")

            print(f"Nghiệm tính theo phương pháp tổng Rieman là {ans1}")

            #########################################

            plot(f, a, b, list_x, list_y, list_x1, h1, N)

        case 2:

            ##################################################################################

            # Câu 2: Tính tích phân x = [4,6] và y = [0,1] của cos(x^4) + 3y^2

            a, b, c, d = 4, 6, 0, 1

            ans2, list_x2, list_y2 = monte2D(f1, a, b, c, d, N, eps)

            ans3, list_x3, list_y3, hx3, hy3 = riemann2D(f1, a, b, c, d, N, eps)

            print(f"Nghiệm tính theo phương pháp Monter-Carlo là {ans2} \n")

            print(f"Nghiệm tính theo phương pháp tổng Rieman là {ans3}")

            ##################################################################################

        case 3:

            ##################################################################################

            # Câu 3: Tính tích phân x = [0,pi/2] và y = [0,pi/2] của xsin(x+y)

            a, b, c, d = 0, pi / 2, 0, pi / 2

            x0 = a

            y0 = c

  

            start_mc = time.time()

            ans2, list_x2, list_y2 = monte2D(f2D, a, b, c, d, N, eps)

            end_mc = time.time()

  

            start_rm = time.time()

            ans3, list_x3, list_y3, hx3, hy3 = riemann2D(f2D, x0, y0, a, b, c, d, N, eps)

            end_rm = time.time()

            print(f"Nghiệm tính theo phương pháp Monter-Carlo là {ans2} với thời gian là {end_mc - start_mc}s \n")

            print(f"Nghiệm tính theo phương pháp tổng Rieman là {ans3} với thời gian là {end_rm - start_rm}s")

        case 4:

            a, b, c, d, g, h = 0, 1, 0, 1, 0, 1

            x0 = a

            y0 = c

            z0 = g

  

            start_mc = time.time()

            ans2, list_x2, list_y2, listz2 = monte3D(f3, a, b, c, d, g, h, N, eps)

            end_mc = time.time()

  

            start_rm = time.time()

            ans3, list_x3, list_y3, list_z3, hx3, hy3, hz3 = riemann3D(f3, x0, y0, z0, a, b, c, d, g, h, N, eps)

            end_rm = time.time()

            print(f"Nghiệm tính theo phương pháp Monter-Carlo là {ans2} với thời gian là {end_mc - start_mc}s \n")

            print(f"Nghiệm tính theo phương pháp tổng Rieman là {ans3} với thời gian là {end_rm - start_rm}s")

        case 5:

            a = 0

            b = 1

            x1, x3, x5, x7, x9 = [a for i in range(5)]

            x2, x4, x6, x8, x10 = [b for i in range(5)]

            start_mc = time.time()

            ans2 = monte10D(f4, a, b, N, eps)

            end_mc = time.time()

  

            start_rm = time.time()

            ans3 = riemann10D(f4, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, a, b, N, eps)

            end_rm = time.time()

            print(f"Nghiệm tính theo phương pháp Monter-Carlo là {ans2} với thời gian là {end_mc - start_mc}s \n")

            print(f"Nghiệm tính theo phương pháp tổng Rieman là {ans3} với thời gian là {end_rm - start_rm}s")

  
  

def f(x):

    return x**2

  
  

def f1(x, y):

    return cos(x**4) + 3 * y**2

  
  

def f2D(x, y):

    return x**2 + y**2

  
  

@njit

def f3(x, y, z):

    return x**2 + y**2 + z**2

  
  

@njit

def f4(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10):

    return (x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10) ** 2

  
  

def monte(f, a, b, N, eps):

    S = np.zeros(N)

    S[0] = 0

    ans = 0

    list_x = np.zeros(N)

    list_y = np.zeros(N)

    for i in range(N):

        x = uniform(a, b)

        y = uniform(f(0), f(x))

        S[i] = f(x)

        ans += (b - a) * S[i] / N

        list_x[i] = x

        list_y[i] = y

  

    return ans, list_x, list_y

  
  

def riemann(f, a, b, N, eps):

    h = (b - a) / N

    S = np.zeros(N)

    x = np.zeros(N)

    S[0] = 0

    ans = 0

    for i in range(0, N):

        x[i] = a + i * h

        S[i] = f(x[i])

        ans += h * S[i]

    return ans, x, h

  
  

def monte2D(f, a, b, c, d, N, eps):

    S = np.zeros(N)

    S[0] = 0

    ans = 0

    list_x = np.zeros(N)

    list_y = np.zeros(N)

    for i in range(N):

        x = uniform(a, b)

        y = uniform(c, d)

        S[i] = f(x, y)

        ans += (b - a) * (d - c) * S[i] / N

        list_x[i] = x

        list_y[i] = y

    return ans, list_x, list_y

  
  

@njit(parallel=True)

def monte3D(f, a, b, c, d, g, h, N, eps):

    S = np.zeros(N)

    count = 0

    S[0] = 0

    ans = 0

    list_x = np.zeros(N)

    list_y = np.zeros(N)

    list_z = np.zeros(N)

    for i in nb.prange(N):

        x = uniform(a, b)

        y = uniform(c, d)

        z = uniform(g, h)

        S[i] = f(x, y, z)

        ans += (b - a) * (d - c) * (h - g) * S[i] / N

        list_x[i] = x

        list_y[i] = y

        list_z[i] = z

        count += i

    return ans, list_x, list_y, list_z

  
  

@njit(parallel=True)

def monte10D(fx, a, b, N, eps):

    S = 0

    ans = 0

    for i in nb.prange(N):

        x1 = uniform(a, b)

        x2 = uniform(a, b)

        x3 = uniform(a, b)

        x4 = uniform(a, b)

        x5 = uniform(a, b)

        x6 = uniform(a, b)

        x7 = uniform(a, b)

        x8 = uniform(a, b)

        x9 = uniform(a, b)

        x10 = uniform(a, b)

        np.random.rand(100)

        S = fx(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)

        ans += (b - a) ** 10 * S / N

    return ans

  
  

def riemann2D(f, a, b, c, d, N, eps):

    ans = 0

    hx = (b - a) / N

    hy = (d - c) / N

    S = np.zeros(N)

    S[0] = 0

    x = np.zeros(N)

    y = np.zeros(N)

    x[0] = a

    y[0] = c

    for i in range(0, N):

        x[i] = a + i * hx

        for j in range(1, N):

            y[j] = c + j * hy

            S[i] = f(x[i], y[j])

            ans += hx * hy * S[i]

    return ans, x, y, hx, hy

  
  

@njit(parallel=True)

def riemann3D(f, x0, y0, z0, a, b, c, d, g, h, N, eps):

    ans = 0

  

    hx = (b - a) / N

    hy = (d - c) / N

    hz = (h - g) / N

  

    S = np.zeros(N)

    S[0] = 0

  

    x = np.zeros(N)

    y = np.zeros(N)

    z = np.zeros(N)

    x[0] = x0

    y[0] = y0

    z[0] = z0

  

    for i in nb.prange(1, N):

        x[i] = x0 + i * hx

        for j in range(1, N):

            y[j] = y0 + j * hy

            for k in range(1, N):

                z[k] = z0 + k * hz

                S[i] = f(x[i], y[j], z[k])

                ans += hx * hy * hz * S[i]

    return ans, x, y, z, hx, hy, hz

  
  

@njit(parallel=True)

def riemann10D(f, x1_0, x2_0, x3_0, x4_0, x5_0, x6_0, x7_0, x8_0, x9_0, x10_0, a, b, N, eps):

    ans = 0

  

    hx = (b - a) / N

  

    S = np.zeros(N)

    S[0] = 0

  

    x1 = np.zeros(N)

    x2 = np.zeros(N)

    x3 = np.zeros(N)

    x4 = np.zeros(N)

    x5 = np.zeros(N)

    x6 = np.zeros(N)

    x7 = np.zeros(N)

    x8 = np.zeros(N)

    x9 = np.zeros(N)

    x10 = np.zeros(N)

  

    x1[0] = x1_0

    x2[0] = x2_0

    x3[0] = x3_0

    x4[0] = x4_0

    x5[0] = x5_0

    x6[0] = x6_0

    x7[0] = x7_0

    x8[0] = x8_0

    x9[0] = x9_0

    x10[0] = x10_0

  

    for k1 in nb.prange(1, N):

        x1[k1] = x1_0 + k1 * hx

        for k2 in nb.prange(1, N):

            x2[k2] = x2_0 + k2 * hx

            for k3 in nb.prange(1, N):

                x3[k3] = x3_0 + k3 * hx

                for k4 in nb.prange(1, N):

                    x4[k4] = x4_0 + k4 * hx

                    for k5 in nb.prange(1, N):

                        x5[k5] = x5_0 + k5 * hx

                        for k6 in nb.prange(1, N):

                            x6[k6] = x6_0 + k6 * hx

                            for k7 in nb.prange(1, N):

                                x7[k7] = x7_0 + k7 * hx

                                for k8 in nb.prange(1, N):

                                    x8[k8] = x8_0 + k8 * hx

                                    for k9 in nb.prange(1, N):

                                        x9[k9] = x9_0 + k9 * hx

                                        for k10 in nb.prange(1, N):

                                            x10[k10] = x10_0 + k10 * hx

  

                                            S[k1] = f(x1[k1], x2[k2], x3[k3], x4[k4], x5[k5], x6[k6], x7[k7], x8[k8], x9[k9], x10[k10])

                                            ans += hx**10 * S[k1]

    return ans

  
  

def plot(f, a, b, xmonte1D, ymonte1D, xrieman1D, h, N):

    #################################################################

    # parameters

    transparent = 0.5

    x = sym.symbols("x")

    fx = sym.Integral(x**2)

    ##################################################################

    # cài đặt 1 số thứ cho plot

    fig, axs = plt.subplots(2, 2, figsize=(18, 7))

    ##################################################################

    # nghiệm của phương trình f(x) dc biểu diễn trên đoạn [0,5]

    X = np.linspace(a, b, N)

    Y = f(X)

    ####################################################################################################################################

    axs[0, 0].plot(X, Y, color="r", label=r"y=$x^2$", lw=2)

    axs[0, 0].scatter(xmonte1D, ymonte1D, marker="o", s=10)

    axs[0, 0].fill_between(X, Y, color="lightgreen", alpha=transparent)

  

    axs[0, 0].legend([r"y=$x^2$", "sỏi", f"${latex(fx)}$"])

    axs[0, 0].set_title(r"Monte-Carlo cho phương trình $y = x^2$")

    axs[0, 0].set_xlabel(r"x")

    axs[0, 0].set_ylabel(r"$f(x)$", rotation=0)

    ####################################################################################################################################

    axs[0, 1].plot(X, Y, color="r", label=r"y=$x^2$", lw=2)

    axs[0, 1].plot(xrieman1D[:-1], f(xrieman1D[:-1]), "b.")

    axs[0, 1].bar(xrieman1D[:-1], f(xrieman1D[:-1]), width=h, alpha=transparent, align="edge", edgecolor="b")

  

    axs[0, 1].legend([r"y=$x^2$"])

    axs[0, 1].set_title(r"Phuwong pháp tổng Riemann(trái) cho phương trình $y = x^2$")

    axs[0, 1].set_xlabel(r"x")

    axs[0, 1].set_ylabel(r"$f(x)$", rotation=0)

    ####################################################################################################################################

    xrieman_mid = (xrieman1D[1:] + xrieman1D[:-1]) * 0.5

    yrieman_mid = f(xrieman_mid)

    axs[1, 0].plot(X, Y, color="r", label=r"y=$x^2$", lw=2)

    axs[1, 0].plot(xrieman_mid, yrieman_mid, "b.")

    axs[1, 0].bar(xrieman_mid, yrieman_mid, width=h, alpha=transparent, edgecolor="b")

  

    axs[1, 0].legend([r"y=$x^2$"])

    axs[1, 0].set_title(r"Phương pháp tổng Riemann(giữa) cho phương trình $y = x^2$")

    axs[1, 0].set_xlabel(r"x")

    axs[1, 0].set_ylabel(r"$f(x)$", rotation=0)

    ####################################################################################################################################

    axs[1, 1].plot(X, Y, color="r", label=r"y=$x^2$", lw=2)

    axs[1, 1].plot(xrieman1D[1:], f(xrieman1D[1:]), "b.")

    axs[1, 1].bar(xrieman1D[1:], f(xrieman1D[1:]), width=-h, alpha=transparent, align="edge", edgecolor="b")

  

    axs[1, 1].legend([r"y=$x^2$"])

    axs[1, 1].set_title(r"Phương pháp tổng Riemann(phải) cho phương trình $y = x^2$")

    axs[1, 1].set_xlabel(r"x")

    axs[1, 1].set_ylabel(r"$f(x)$", rotation=0)

    ####################################################################################################################################

    fig.savefig(f"Cau1.Monte-Carlo-Riemann1D.png")

  

    plt.show()

  
  

def main():

    id = int(input("Nhập câu thứ ... để in ra kết quả của câu đó(có 5 câu): "))

  

    promp_user(id)

    # test print

    # print(list_x2)

  
  

if __name__ == "__main__":

    main()
```

