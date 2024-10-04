import numpy as np
import csv
from numpy import exp, cos
from math import log as ln
from tabulate import tabulate
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy.optimize import fsolve


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

    solution = []

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

    for i in range(len(nc)):
        if nc[i] == 0:
            continue
        solution.append(float(nc[i]))

    return (na, nb, solution), count


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

    for i in range(len(p)):
        if p[i] != 0:
            listp.append(float(p[i]))

    return listp, count


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
    solution = []
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

        for i in range(len(x)):
            if x[i] == 0:
                continue
            solution.append(float(x[i]))

        return solution, n, fx.__name__

    else:
        p0 = float(input("Nhập lại p0: "))


def find_optimal(count_arr, N):
    optimal_fucn = min(count_arr, key=count_arr.get)
    # print(optimal_fucn, count_arr[optimal_fucn])
    return optimal_fucn, count_arr[optimal_fucn]


########################################################################################
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

    for i in range(len(p)):
        if p[i] != 0:
            listp.append(float(p[i]))

    return listp, count


def save_log(file, sol_secant, sol_newton, sol_bisection, sol_fixedpoint, maxcount):
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

        for i in range(maxcount):
            bisection_values = f"{na[i]:^28} | {nb[i]:^28} | {nc[i]:^28}"

            if i < len(sol_secant):
                secant_vals = f"{sol_secant[i]:^28}"
            else:
                secant_vals = "0".center(28)

            if i < len(sol_newton):
                newton_vals = f"{sol_newton[i]:^28}"
            else:
                newton_vals = "0".center(28)

            if i < len(sol_fixedpoint):
                fixed_point_vals = f"{sol_fixedpoint[i]:^28}"
            else:
                fixed_point_vals = "0".center(28)

            writer.writerow(
                {
                    f"{'n':<3}": f"{i:^3}",
                    f"{'Secant':^28}": f"{secant_vals}",
                    f"{'Newton Rasphson':^28}": f"{newton_vals}",
                    f"{'Bisection':^96}": bisection_values,
                    f"{'Fixed Point':^28}": f"{fixed_point_vals}",
                }
            )


def read_log(file, sol_secant, sol_newton, sol_bisection, sol_fixedpoint, maxcount):
    na, nb, nc = sol_bisection

    with open(file, "w", newline="") as writefile:
        header = ["n", "Secant", "Newton Rasphson", "Bisection", "Fixed Point"]
        writer = csv.DictWriter(writefile, fieldnames=header)
        writer.writeheader()
        for i in range(maxcount):
            bisection_values = nc[i]
            if i < len(sol_secant):
                secant_vals = sol_secant[i]
            else:
                secant_vals = 0

            if i < len(sol_newton):
                newton_vals = sol_newton[i]
            else:
                newton_vals = 0

            if i < len(sol_fixedpoint):
                fixed_point_vals = sol_fixedpoint[i]
            else:
                fixed_point_vals = 0

            writer.writerow(
                {
                    "n": i,
                    "Secant": secant_vals,
                    "Newton Rasphson": newton_vals,
                    "Bisection": bisection_values,
                    "Fixed Point": fixed_point_vals,
                }
            )
    with open(file, "r") as readfile:
        rows = []
        reader = csv.DictReader(readfile)
        for row in reader:
            rows.append(row)

        return tabulate(rows, headers="keys", tablefmt="fancy_grid")


def plot_solution(file, sol_secant, sol_newton, sol_bisection, sol_fixedpoint, ex_sol, N):
    an, bn, cn = sol_bisection

    fig, ax = plt.subplots()

    ax.plot(sol_secant, lw=1)
    ax.plot(sol_newton, lw=1)
    ax.plot(cn, lw=1)
    ax.plot(sol_fixedpoint, lw=1)
    ax.axhline(y=ex_sol, ls="-.", lw=1, color="r")

    ax.legend(["Secant", "Newton Raphson", "Bisection", "Fixed point"])
    ax.set_title("Table compare data")
    ax.set_xlabel("Number of loops")
    ax.set_ylabel("Data")

    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(MultipleLocator(5))

    fig.savefig("data1.pdf")
    plt.show()


def main():
    ######################################################################################## parameters
    a, b = 0, 10
    N = 50
    p0 = 2
    p1 = 4
    p = 2
    eps1 = 10e-8
    eps2 = 10e-8
    file = "table.txt"
    file1 = "table1.txt"
    if not a <= p0 <= b:
        p0 = input("Nhap lai p0: ")

    ######################################################################################### call function

    solution_secant, count1 = secant(f, df_secant, p0, p1, N, eps1)
    solution_newton, count2 = newtonr(f, df_newton, p0, N, eps1)
    n, count3 = bisection(f, a, b, N, eps1)
    solution_fixedpoint, count4, g = fixed_point(g2, p0, a, b, N, eps1)

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
    # print(f"voi phuong phap fixed point thi {gfunc} la toi uu nhat voi {loops} vong lap")

    ##########################################################################################

    ################################## ddeems sos vong lap
    maxcount = count1
    if maxcount < count2:
        maxcount = count2

    if maxcount < count3:
        maxcount = count3

    if maxcount < count4:
        maxcount = count4
    #########################################
    # test print
    # print(solution_secant)
    # print()
    ##########################################################################################

    ####### save file
    save_log(file, solution_secant, solution_newton, n, solution_fixedpoint, maxcount)
    read = read_log(file1, solution_secant, solution_newton, n, solution_fixedpoint, maxcount)
    print(read)
    ###### plot file

    exact_solution = fsolve(f, p0)

    # plot_solution(file, solution_secant, solution_newton, n, solution_fixedpoint, exact_solution, maxcount)


if __name__ == "__main__":
    main()
