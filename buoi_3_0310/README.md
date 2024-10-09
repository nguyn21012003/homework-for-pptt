# Bài toàn giá trị đầu
## Hàm f(t,y)

```python 
def f(t:list, y:list)->float:
    return y - t**2 + 1
```

f phụ thuộc theo 2 biến t,y hoặc x,y kiểu float và trả về kiểu float
## Nghiệm giải tích

```python 
def exact_solution(t0:float, tn:float, h:float)->float:
    t = np.arange(t0, tn + h, h)
    y = np.zeros(len(t))
    for i in range(len(t)):
        y[i] = (t[i] + 1) ** 2 - 0.5 * exp(t[i])
    return y
```

## Phương pháp Euler 

```python
def euler(f:float, tn:float, yn:float, h:float)->float:
    return yn + h * f(tn, yn)
```

## Phương pháp Euler - cải tiến(Modified Euler)

```python 
def mod_euler(f:float, tn:float, yn:float, h:float)->float:
    return yn + h / 2 * (f(tn, yn) + f(tn, yn + h * f(tn, yn)))
```

## Phương pháp Runge-Kutta 2(phương pháp midpoints)

```python
def midpoint(f:float, tn:float, yn:float, h:float)->float:
    return yn + h * f(tn + h / 2, yn + h / 2 * f(tn, yn))
```

## Phương pháp Runge-Kutta 3

```python
def rk3(f:float, tn:float, yn:float, h:float)->float:
    k1 = f(tn, yn)
    k2 = f(tn + 1 / 3 * h, yn + 1 / 3 * h * k1)
    k3 = f(tn + 2 / 3 * h, yn + 2 / 3 * h * k2)
    return yn + h / 4 * (k1 + 3 * k3)
```

## Function để giải hàm f(t,y)

```python
def solve_ode(f:float, y0:float, t0:float, tn:float, solver:float, h:float)->tuple:
    t = np.arange(t0, tn + h, h)
    y = np.zeros(len(t))
    y[0] = y0
    for i in range(len(t) - 1):
        y[i + 1] = solver(f, t[i], y[i], h)
    return t, y
```

###### Mô tả:

Function này sẽ lấy tên hàm được định nghĩa từ đầu và $y0,t0,tn$ là các thông số ban đầu và tạo một mảng chứa các giá trị của $y$ và gán cho toàn bộ phần tử của $y$ bên trong mảng đều bằng 0.
Đặt phần tử đầu tiên của $y$ = $y0$ và sử dụng vòng lặp `for` để gán lại các giá trị bên trong mảng $y$. Sử dụng `np.zeros(len(t))` sẽ tốn bộ hơn so với tạo một list trống, bù lại ta sẽ kiểm soát được các biến bên trong $y$ dễ hơn. Vòng lặp for chạy i từ 0 cho tới `len(t) - 1`. `-1` ở đây có ý nghĩa là để tránh bị lỗi `out of index range của mảng` $y$. 
###### Trong đó:
- $f$ lấy từ hàm $f(t,y)$.
- cho $t \in [a,b]$.
- $y0$ là giá trị ban đầu của y tại $t=t0$  vd : $y(t0) = y(a) = num$.
- t$0,tn = a,b$.
- solver chính là các phương pháp để giải.
- Hàm này cuối cùng sẽ return ra 1 `tuple` chứa 1 list $t$ và 1 list $y$.

## Hàm để plot

```python
def plot_ty(file:str, t:list, y:list, yeuler:list, yrk4:list, ymidpoint:list, ymodeuler:list, yrk3:list):
    fig, ax = plt.subplots(figsize=(15, 7))
  
    ax.plot(t, y, "-.", color="r")
    ax.plot(t, yeuler, color="orange")
    ax.plot(t, ymodeuler, color="gold")
    ax.plot(t, ymidpoint, color="lime")
    ax.plot(t, yrk3, color="navy")
    ax.plot(t, yrk4, color="magenta")

    ax.legend(
        ["Nghiệm giải tích", "Phương pháp Euler", "Phương pháp Euler cải tiến", "Phương pháp Midpoint/RK2", "Phương pháp RK3", "Phương pháp RK4"],
        loc="best")
    ax.grid(True)
    ax.set_title(r"Nghiệm của phương trình vi phân $y'=y-t^2+1$")
    ax.set_xlabel("t")
    ax.set_ylabel("y")
  
    fig.savefig(file)

    plt.show()  
```
###### Mô tả:

Hàm này sẽ lấy giá trị $y$ của tất cả các phương pháp, và $t$. Do các phương pháp được lấy cùng trên 1 khoảng $t$ và được chia cùng cho `N` khoảng,  nên `return t` thì ta không cần phân biệt $t$.

Trong đó:
- Ta cần phải lấy giá trị của t `tspan = [0,2]` được định nghĩa bằng list chứa 2 thành phần gồm `t0 = 0`  `tn =2` 

## Hàm để ghi file
```python
def save_log(file: str, t: list, y: list, yeuler: list, yrk4: list, ymidpoint: list, ymodeuler: list, yrk3: list) -> dict:
    with open(file, "w", newline="") as writefile:
        header = [f"{'n':^3}", f"{'t':^6}", f"{'y':^18}", f"{'Euler':^18}", f"{'Euler modified':^18}", f"{'RK2':^18}", f"{'RK3':^18}", f"{'RK4':^18}"]
        writer = csv.DictWriter(writefile, fieldnames=header, delimiter="|")
        writer.writeheader()
        for i in range(len(t)):
            writer.writerow(
                {
                    f"{'n':^3}": f"{i:^3}",
                    f"{'t':^6}": f"{t[i]:^4.4f}",
                    f"{'y':^18}": f"{y[i]:^18}",
                    f"{'Euler':^18}": f"{yeuler[i]:^18}",
                    f"{'Euler modified':^18}": f"{ymodeuler[i]:^18}",
                    f"{'RK2':^18}": f"{ymidpoint[i]:^18}",
                    f"{'RK3':^18}": f"{yrk3[i]:^18}",
                    f"{'RK4':^18}": f"{yrk4[i]:^18}",
                }
            )
```

## Kết quả
Chi tiết kết quả thì trong code.

![Result](res/data.pdf)


|  n  |  t  |         y          |     Euler     |  Euler modified   |
| :-: | :-: | :----------------: | :-----------: | :---------------: |
|  0  | 0.0 |        0.5         |      0.5      |        0.5        |
|  1  | 0.2 | 0.829298620919915  |      0.8      |       0.83        |
|  2  | 0.4 | 1.2140876511793646 |     1.152     |      1.2238       |
|  3  | 0.6 | 1.648940599804746  |    1.5504     |     1.677836      |
|  4  | 0.8 | 2.1272295357537665 |    1.98848    |    2.18775992     |
|  5  | 1.0 | 2.6408590857704777 |   2.458176    |   2.7482671024    |
|  6  | 1.2 | 3.1799415386317267 |   2.9498112   |  3.352885864928   |
|  7  | 1.4 | 3.732400016577664  |  3.45177344   | 3.99372075521216  |
|  8  | 1.6 | 4.283483787802443  |  3.950128128  | 4.661139321358835 |
|  9  | 1.8 | 4.815176267793525  | 4.4281537536  | 5.343389972057778 |
| 10  | 2.0 | 5.305471950534675  | 4.86578450432 | 6.02613576591049  |


|  n  |  t  |         y          |        RK2         |        RK3         |        RK4         |
| :-: | :-: | :----------------: | :----------------: | :----------------: | :----------------: |
|  0  | 0.0 |        0.5         |        0.5         |        0.5         |        0.5         |
|  1  | 0.2 | 0.829298620919915  |       0.828        | 0.8292444444444445 | 0.8292933333333334 |
|  2  | 0.4 | 1.2140876511793646 |      1.21136       | 1.2139749925925927 | 1.2140762106666667 |
|  3  | 0.6 | 1.648940599804746  |     1.6446592      | 1.6487659020641976 | 1.6489220170416001 |
|  4  | 0.8 | 2.1272295357537665 |    2.121284224     | 2.1269905328321843 | 2.1272026849479437 |
|  5  | 1.0 | 2.6408590857704777 |   2.63316675328    | 2.6405555485434853 | 2.6408226927287517 |
|  6  | 1.2 | 3.1799415386317267 |  3.1704634390016   | 3.179576287732221  | 3.1798941702322305 |
|  7  | 1.4 | 3.732400016577664  | 3.721165395581952  | 3.731980283861397  | 3.7323400728549796 |
|  8  | 1.6 | 4.283483787802443  | 4.2706217826099815 | 4.283023031133831  | 4.283409498318405  |
|  9  | 1.8 | 4.815176267793525  | 4.800958574784177  | 4.814696573135896  | 4.815085694579433  |
| 10  | 2.0 | 5.305471950534675  | 5.290369461236696  | 5.305007192434418  | 5.305363000692653  |
## Source code 

```python
import numpy as np
import matplotlib.pyplot as plt
import csv
from tabulate import tabulate
from math import exp

def f(t, y):
    return y - t**2 + 1

def exact_solution(t0, tn, h):
    t = np.arange(t0, tn + h, h)
    y = np.zeros(len(t))
    for i in range(len(t)):
        y[i] = (t[i] + 1) ** 2 - 0.5 * exp(t[i])
    return y

def euler(f, tn, yn, h):
    return yn + h * f(tn, yn)

def mod_euler(f, tn, yn, h):
    return yn + h / 2 * (f(tn, yn) + f(tn, yn + h * f(tn, yn)))
    
def midpoint(f, tn, yn, h):
    return yn + h * f(tn + h / 2, yn + h / 2 * f(tn, yn))

def rk3(f, tn, yn, h):  
    k1 = f(tn, yn)
    k2 = f(tn + 1 / 3 * h, yn + 1 / 3 * h * k1)
    k3 = f(tn + 2 / 3 * h, yn + 2 / 3 * h * k2)
    return yn + h / 4 * (k1 + 3 * k3)

def rk4(f, tn, yn, h):
    k1 = f(tn, yn)
    k2 = f(tn + 0.5 * h, yn + 0.5 * h * k1)
    k3 = f(tn + 0.5 * h, yn + 0.5 * h * k2)
    k4 = f(tn + h, yn + h * k3)
    return yn + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

def solve_ode(f, y0, t0, tn, solver, h):
    t = np.arange(t0, tn + h, h)
    y = np.zeros(len(t))
    y[0] = y0
    for i in range(len(t) - 1):
        y[i + 1] = solver(f, t[i], y[i], h)
    return t, y
    
def plot_ty(file, t, y, yeuler, yrk4, ymidpoint, ymodeuler, yrk3):
    fig, ax = plt.subplots(figsize=(15, 7))
  
    ax.plot(t, y, "-.", color="r")
    ax.plot(t, yeuler, color="orange")
    ax.plot(t, ymodeuler, color="gold")
    ax.plot(t, ymidpoint, color="lime")
    ax.plot(t, yrk3, color="navy")
    ax.plot(t, yrk4, color="magenta")

    ax.legend(
        ["Nghiệm giải tích", "Phương pháp Euler", "Phương pháp Euler cải tiến", "Phương pháp Midpoint/RK2", "Phương pháp RK3", "Phương pháp RK4"],
        loc="best")
    ax.grid(True)
    ax.set_title(r"Nghiệm của phương trình vi phân $y'=y-t^2+1$")
    ax.set_xlabel("t")
    ax.set_ylabel("y")
  
    fig.savefig(file)

    plt.show()  

def table_log(file: str, t: list, y: list, yeuler: list, yrk4: list, ymidpoint: list, ymodeuler: list, yrk3: list) -> tabulate:
    with open(file, "w", newline="") as writefile:
        header = ["n", "t", "y", "Euler", "Euler modified", "RK2", "RK3", "RK4"]
        writer = csv.DictWriter(writefile, fieldnames=header)
        writer.writeheader()
        for i in range(len(t)):
            writer.writerow(
                {
                    "n": i,
                    "t": f"{t[i]:.4f}",
                    "y": y[i],
                    "Euler": yeuler[i],
                    "Euler modified": ymodeuler[i],
                    "RK2": ymidpoint[i],
                    "RK3": yrk3[i],
                    "RK4": yrk4[i],
                }
            )

    with open(file, "r") as readfile:
        rows = []
        reader = csv.DictReader(readfile)
        for row in reader:
            rows.append(row)
        return tabulate(rows, headers="keys", tablefmt="fancy_grid")
  
def save_log(file: str, t: list, y: list, yeuler: list, yrk4: list, ymidpoint: list, ymodeuler: list, yrk3: list) -> dict:
    with open(file, "w", newline="") as writefile:
        header = [f"{'n':^3}", f"{'t':^6}", f"{'y':^18}", f"{'Euler':^18}", f"{'Euler modified':^18}", f"{'RK2':^18}", f"{'RK3':^18}", f"{'RK4':^18}"]
        writer = csv.DictWriter(writefile, fieldnames=header, delimiter="|")
        writer.writeheader()
        for i in range(len(t)):
            writer.writerow(
                {
                    f"{'n':^3}": f"{i:^3}",
                    f"{'t':^6}": f"{t[i]:^4.4f}",
                    f"{'y':^18}": f"{y[i]:^18}",
                    f"{'Euler':^18}": f"{yeuler[i]:^18}",
                    f"{'Euler modified':^18}": f"{ymodeuler[i]:^18}",
                    f"{'RK2':^18}": f"{ymidpoint[i]:^18}",
                    f"{'RK3':^18}": f"{yrk3[i]:^18}",
                    f"{'RK4':^18}": f"{yrk4[i]:^18}",
                }
            )

def main():
    N = 3
  
    tspan = [0, 2]
    t0 = tspan[0]
    tn = tspan[1]
    h = (tn - t0) / N
    y0 = 0.5

    file_log = "data.txt"
    file_data_log = "data_table.txt"
    file_pdf = "data.pdf"

    y = exact_solution(t0, tn, h)
    t1, y1 = solve_ode(f, y0, t0, tn, euler, h)
    t2, y2 = solve_ode(f, y0, t0, tn, rk4, h)
    t3, y3 = solve_ode(f, y0, t0, tn, midpoint, h)
    t4, y4 = solve_ode(f, y0, t0, tn, mod_euler, h)
    t5, y5 = solve_ode(f, y0, t0, tn, rk3, h)

    plot_ty(file_pdf, t1, y, y1, y2, y3, y4, y5)
    table = table_log(file_log, t1, y, y1, y2, y3, y4, y5)
    print(table)
    save_log(file_data_log, t1, y, y1, y2, y3, y4, y5)

if __name__ == "__main__":
    main()
```
