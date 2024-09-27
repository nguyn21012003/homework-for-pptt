from math import sqrt
import numpy as np
import time, csv


def fx1(xn):
    return sqrt(10 - xn**3) / 2


def fx2(xn):
    return xn - (xn**3 + 4 * xn**2 - 10)


def fx3(xn):
    return sqrt(10 / xn - 4 * xn)


def fx4(xn):
    return sqrt(10 / (xn + 4))


def fx5(xn):
    return xn - (xn**3 + 4 * xn**2 - 10) / (3 * xn**2 + 8 * xn)


def fx6(xn):
    return 1 / (xn - 1)

def fx7(xn):
    return 1 / (xn - 1)    


def fixed_point(fx, p0, a, b, N, eps):
    """Description


    Find the root of a function by using Fixed Point method

    Return the root and the number of loops of a function.
    Args:
        fx (callable function): function address
        p0 (float): the init value in from a to b 
        a (float): the boundary condition
        b (float): _description_
        N (int): the maximum loop
        eps (float): aka TOL, the value error is allowed

    Returns:
        n (int): the number of loops 
        x (float): the value of x at n_th loops  


    """
    if a <= p0 <= b:
        x = np.zeros(N)
        x[0] = p0
        count = 0
        for n in range(N):
            if n + 1 < N:
                try:
                    x[n + 1] = fx(x[n])
                    print(f"x[{n}] = {x[n]}")
                    if x[n]>= b*10:
                        print(f"nghiem khong hoi tu voi {count} vong lap.")
                        return n, x
                    elif abs(x[n] - x[n - 1]) <= eps and x[n]<= b*10:
                        print(f"nghiem hoi tu {x[count]} voi {count} vong lap.")
                        return n, x
                except (ZeroDivisionError, ValueError) as Z:
                    print(f"{Z}")
                    return n, x

                count += 1
            if n + 1 >= N:
                print(f"nghiem khong hoi tu voi {count} vong lap.")
                return n, x
    else:
        p0 = float(input("Nhập lại p0: "))


def main():
    """Co the mo rong them bang cach so sanh thuat toan nao toi uu hon, bang cach so sanh cac gia tri n_i
    """
    p0 = float(input("Nhập p0: "))
    N = 20
    eps = 1e-5
    a,b = 1,2
    res1 = []
    res2 = []
    res3 = []
    res4 = []
    res5 = []
    file = "table.txt"
    
    n1, x1 = fixed_point(fx1, p0, a, b, N, eps)
    for n in range(len(x1)):
        res1.append(float(x1[n]))
    
    n2, x2 = fixed_point(fx2, p0, a, b, N, eps)
    for n in range(len(x2)):
        res2.append(float(x2[n]))
    
    n3, x3 = fixed_point(fx3, p0, a, b, N, eps)
    for n in range(len(x3)):
        res3.append(float(x3[n]))
    
    n4, x4 = fixed_point(fx4, p0, a, b, N, eps)
    for n in range(len(x4)):
        res4.append(float(x4[n]))
    
    n5, x5 = fixed_point(fx5, p0, a, b, N, eps)
    for n in range(len(x5)):
        res5.append(float(x5[n]))
    
    with open(file, "w", newline="") as writefile:
        header = [f"{"n":<2}" ,f"{"fx1":^18}",f"{"fx2":^25}",f"{"fx3":^18}",f"{"fx4":^18}",f"{"fx5":^18}"]
        writer = csv.DictWriter(writefile, fieldnames=header, delimiter="|")
        writer.writeheader()
        for i in range(20):
            writer.writerow({f"{"n":<2}": f"{i:<2}", f"{"fx1":^18}": f"{res1[i]:<18}",f"{"fx2":^25}":f"{res2[i]:<25}",f"{"fx3":^18}":f"{res3[i]:<18}",f"{"fx4":^18}":f"{res4[i]:<18}",f"{"fx5":^18}":f"{res5[i]:<18}"})


if __name__ == "__main__":
    main()
