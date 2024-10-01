import numpy as np
import csv
from numpy import exp, cos
from math import log as ln
import time


def f(x):
    return exp(x) - 4 * x - 5


def df_newton(x):
    return exp(x) - 4


def df_secant(f, pn1, pn2):
    return (f(pn1) - f(pn2)) / (pn1 - pn2)


##################################################################
def bisection(f, a, b, N, eps):
    a = float(a)
    b = float(b)
    if a > b:
        a = b
        b = a

    na = np.zeros(N)
    nb = np.zeros(N)
    nc = np.zeros(N)
    na[0] = a
    nb[0] = b

    count = 0

    for i in range(N):
        if i + 1 < N:
            nc[i] = (na[i] + nb[i]) / 2

            if f(nc[i]) == 0:
                break

            if abs(f(nc[i])) < eps:
                break

            if f(na[i]) * f(nc[i]) < 0:
                nb[i + 1] = nc[i]
                na[i + 1] = na[i]
                count += 1

            elif f(nc[i]) * f(nb[i]) < 0:
                nb[i + 1] = nb[i]
                na[i + 1] = nc[i]
                count += 1

    return (na, nb, nc), count


##########################################################################
def newtonr(f, df, p0, N, eps):
    p = np.zeros(N)
    p[0] = p0
    count = 0
    listp = []
    for i in range(1, N + 1):
        if i + 1 < N:
            p[i] = p[i - 1] - f(p[i - 1]) / df(p[i - 1])
            count += 1
            if abs(p[i] - p[i - 1]) <= eps:
                break
        else:
            print(f"Nghiem khong hoi tu voi {count} vong lap !")
            break

    # for i in range(len(p)):
    #    if p[i] != 0:
    #        listp.append(float(p[i]))

    return p, count


################################################################################# fixed point
def g1(x):
    return (exp(x) - 5) / 4


def g2(x):
    return ln(4 * x + 5)


def g3(x):
    return ln(5 + exp(x)) / ln(4)


def g4(x):
    return ln(exp(x) / 4) + 5 / 4


def fixed_point(fx, p0, a, b, N, eps):
    if a <= p0 <= b:
        x = np.zeros(N)
        x[0] = p0
        count = 0
        for n in range(1, N + 1):
            if n + 1 > N:
                # print(f"nghiem khong hoi tu voi {n} vong lap voi {fx.__name__}")
                break
            else:
                try:
                    x[n] = fx(x[n - 1])
                    if x[n] >= b * 10:
                        # print(f"nghiem khong hoi tu voi {count} vong lap.")
                        break
                    elif abs(x[n] - x[n - 1]) <= eps and x[n] <= b * 10:
                        # print(f"nghiem hoi tu {x[count]} voi {count} vong lap.")
                        break
                except (ZeroDivisionError, ValueError) as Z:
                    print(f"{Z}")
                    break

                count += 1

        return x, n, fx.__name__
    else:
        p0 = float(input("Nhập lại p0: "))


def find_optimal(count_arr, N):
    optimal_fucn = min(count_arr, key=count_arr.get)
    # print(optimal_fucn, count_arr[optimal_fucn])
    return optimal_fucn, count_arr[optimal_fucn]


##############################################################################################
def secant(f, df, p0, p1, N, eps):
    p = np.zeros(N)
    p[0] = p0
    p[1] = p1
    count = 0
    listp = []
    for i in range(2, N + 2):
        if i + 1 < N:
            p[i] = p[i - 1] - f(p[i - 1]) / df(f, p[i - 1], p[i - 2])
            count += 1
            if abs(p[i] - p[i - 1]) <= eps:
                break
        else:
            print(f"Nghiem khong hoi tu voi {count} vong lap !")
            break

    # for i in range(len(p)):
    #    if p[i] != 0:
    #        listp.append(float(p[i]))

    return p, count


def save_log(file, sol_secant, sol_newton, sol_bisection, sol_fixedpoint, N):
    na, nb, nc = sol_bisection

    with open(file, "w", newline="") as writefile:
        header = [
            f"{'n':<3}",
            f"{'Secant':^28}",
            f"{'Newton Rasphson':^28}",
            f"{'Bisection':^96}",
            f"{'Fixed Point':^28}",
        ]
        writer = csv.DictWriter(writefile, fieldnames=header, delimiter="|")
        writer.writeheader()
        for i in range(N):
            writer.writerow(
                {
                    f"{'n':<3}": f"{i:^3}",
                    f"{'Secant':^28}": f"{sol_secant[i]:^28}",
                    f"{'Newton Rasphson':^28}": f"{sol_newton[i]:^28}",
                    f"{'Bisection':^96}": [
                        f"{na[i]:^28}",
                        f"{nb[i]:^28}",
                        f"{nc[i]:^28}",
                    ],
                    f"{'Fixed Point':^28}": f"{sol_fixedpoint[i]:^28}",
                }
            )


def main():
    ######################################################### parameters
    a, b = 0, 10
    N = 100
    p0 = 2
    p1 = 3
    p = 2.4
    eps1 = 10e-8
    eps2 = 10e-15
    file = "table.txt"
    if not a <= p0 <= b:
        p0 = input("Nhap lai p0: ")
    #################################################################### call function

    solution_secant, count1 = secant(f, df_secant, p0, p1, N, eps1)
    solution_newton, count2 = newtonr(f, df_newton, p0, N, eps1)
    n, count3 = bisection(f, a, b, N, eps1)
    solution_fixedpoint, count, g = fixed_point(g2, p, a, b, N, eps1)

    ##########################################################################################
    ## phan nay de sau nay update them cai gi nua thi lam
    list_g = [g1, g2, g3, g4]
    solution = []
    counts = {}
    for i in list_g:
        solution, count, g = fixed_point(i, p, a, b, N, eps1)
        # solution_fixedpoint.append(solution)
        counts[f"{g}"] = count

    gfunc, loops = find_optimal(counts, N)
    print(f"voi phuong phap fixed point thi {gfunc} la toi uu nhat voi {loops} vong lap")

    ##########################################################################################
    # test print
    # print(solution_secant, solution_newton, solution_fixedpoint, na, nb)
    # print()
    ##########################################################################################

    # save file
    save_log(
        file,
        solution_secant,
        solution_newton,
        n,
        solution_fixedpoint,
        N,
    )


if __name__ == "__main__":
    main()
