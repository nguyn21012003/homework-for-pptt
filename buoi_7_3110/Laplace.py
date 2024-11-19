import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm


def u_values(d):
    u = np.zeros((d, d))
    for i in range(d):
        u[0][i] = 100
        u[d - 1][i] = -100
    print(u)
    return u


def Gauss(d, N, u):
    for k in tqdm(range(N)):
        for i in range(1, d - 1):
            for j in range(1, d - 1):
                u[i][j] = 0.25 * (u[i + 1][j] + u[i - 1][j] + u[i][j + 1] + u[i][j - 1])
    return u



def main():
    N = 10000
    d = 40
    h = d / N

    u = u_values(d)

    u_results = Gauss(d, N, u)
    print(u_results)

    x = np.linspace(0, d, d)
    y = np.linspace(0, d, d)
    x_grid, y_grid = np.meshgrid(x, y)

    with open("laplace.txt", "w") as laplace:
        s0 = "{0:^10} {delim} {1:^11} {delim} {2:^10} \n"
        # laplace.write(s0.format("x", "y ", "u(x,y)", delim="|"))
        for i in range(d):
            for j in range(d):
                s1 = "{0:10.6f} {delim} {1:11.6f} {delim} {2:10.6f} \n"
                # laplace.write(s1.format(x[i], y[i], u_results[i][j], delim="|"))

    fig = plt.figure(figsize=(12, 6))
    ax = plt.axes(projection="3d")
    ax.grid()

    ax.plot_wireframe(x_grid, y_grid, u_results)
    ax.set_title("Electrostatic Potentials")

    ax.set_xlabel("x", labelpad=20)
    ax.set_ylabel("y", labelpad=20)
    ax.set_zlabel("u(x,y)", labelpad=20)

    plt.show()


main()
