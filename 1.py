# KN va Khai Hoan


import torch
import numpy as np
import matplotlib.pyplot as plt
import torch.nn as nn
import torch.optim as optim


class PINN(nn.Module):
    def __init__(self):
        super(PINN, self).__init__()
        self.fc1 = nn.Linear(2, 50)
        self.fc2 = nn.Linear(50, 50)
        self.fc3 = nn.Linear(50, 50)
        self.fc4 = nn.Linear(50, 1)

    def f_nn(self, t, x):
        u = self.u_nn(t, x)
        u_t = get_derivative(u, t, 1)
        u_x = get_derivative(u, x, 1)
        u_xx = get_derivative(u, x, 2)
        k = 0.01 * u + 7
        k_u = 0.01
        c = 0.0005 * u**2 + 500
        s = self.source_term(t, x)
        f = c * u_t - k_u * u_x * u_x - k * u_xx - s
        return f

    def get_derivative(y, x, n):
        if n == 0:
            return y
        else:
            dy_dx = grad(y, x, torch.ones_like(y), create_graph=True, retain_graph=True, allow_unused=True)[0]

            return get_derivative(dy_dx, x, n - 1)

    def source_term(self, t, x):
        t_max = 0.5
        sigma = 0.02
        u_max = 1
        p = 0.25 * torch.cos(2 * np.pi * t / t_max) + 0.5
        p_t = -0.5 * torch.sin(2 * np.pi * t / t_max) * np.pi / t_max
        u_sol = u_max * torch.exp(-((x - p) ** 2) / (2 * sigma**2))
        k_sol = 0.01 * u_sol + 7
        k_u_sol = 0.01
        c_sol = 0.0005 * u_sol**2 + 500
        factor = 1 / (sigma**2)
        s = factor * k_sol * u_sol + u_sol * (x - p) * factor * (c_sol * p_t - (x - p) * factor * (k_sol + u_sol * k_u_sol))
        return s

    def cost_function(self, t0, x0, t_lb, x_lb, t_ub, x_ub, t_f, x_f, u0):
        u0_pred = self.u_nn(t0, x0)
        u_lb_pred = self.u_nn(t_lb, x_lb)
        u_x_lb_pred = get_derivative(u_lb_pred, x_lb, 1)
        u_ub_pred = self.u_nn(t_ub, x_ub)
        u_x_ub_pred = get_derivative(u_ub_pred, x_ub, 1)
        f_pred = self.f_nn(t_f, x_f)
        mse_0 = torch.mean((u0 - u0_pred) ** 2)
        mse_b = torch.mean(u_x_lb_pred**2) + torch.mean(u_x_ub_pred**2)
        mse_f = torch.mean((f_pred) ** 2)
        return mse_0, mse_b, mse_f

    def train_adam(self, t0, x0, t_lb, x_lb, t_ub, x_ub, t_f, x_f, u0, epochs=3000):
        optimizer = torch.optim.Adam(self.parameters(), lr=1e-3)
        mse_0_vals = []
        mse_b_vals = []
        mse_f_vals = []

        for epoch in range(epochs):
            optimizer.zero_grad()

            # Calculate the cost function
            mse_0, mse_b, mse_f = self.cost_function(t0, x0, t_lb, x_lb, t_ub, x_ub, t_f, x_f, u0)
            loss = mse_0 + mse_b + mse_f
            loss.backward()
            optimizer.step()

            # Track the cost function values
            mse_0_vals.append(mse_0.item())
            mse_b_vals.append(mse_b.item())
            mse_f_vals.append(mse_f.item())

            if epoch % 500 == 0:
                print(f"Epoch {epoch}/{epochs}, Loss: {loss.item()}")

        # Plot the cost function over epochs
        plt.figure(figsize=(10, 6))
        plt.plot(mse_0_vals, label="Initial condition MSE")
        plt.plot(mse_b_vals, label="Boundary condition MSE")
        plt.plot(mse_f_vals, label="Residual MSE")
        plt.xlabel("Epochs")
        plt.ylabel("Cost")
        plt.legend()
        plt.title("Cost Function over Epochs (Adam Optimization)")
        plt.show()

    def train_lbfgs(self, t0, x0, t_lb, x_lb, t_ub, x_ub, t_f, x_f, u0, epochs=1000):
        optimizer = torch.optim.LBFGS(self.parameters(), lr=1e-1)
        mse_0_vals = []
        mse_b_vals = []
        mse_f_vals = []

        def closure():
            optimizer.zero_grad()
            mse_0, mse_b, mse_f = self.cost_function(t0, x0, t_lb, x_lb, t_ub, x_ub, t_f, x_f, u0)
            loss = mse_0 + mse_b + mse_f
            loss.backward()
            return loss

        for epoch in range(epochs):
            optimizer.step(closure)

            mse_0, mse_b, mse_f = self.cost_function(t0, x0, t_lb, x_lb, t_ub, x_ub, t_f, x_f, u0)
            mse_0_vals.append(mse_0.item())
            mse_b_vals.append(mse_b.item())
            mse_f_vals.append(mse_f.item())

            if epoch % 100 == 0:
                print(f"Epoch {epoch}/{epochs}, Loss: {mse_0 + mse_b + mse_f}")

        # Plot the cost function over epochs
        plt.figure(figsize=(10, 6))
        plt.plot(mse_0_vals, label="Initial condition MSE")
        plt.plot(mse_b_vals, label="Boundary condition MSE")
        plt.plot(mse_f_vals, label="Residual MSE")
        plt.xlabel("Epochs")
        plt.ylabel("Cost")
        plt.legend()
        plt.title("Cost Function over Epochs (L-BFGS Optimization)")
        plt.show()

    def plot_u_nn(self, t_vals, x_vals):
        T, X = torch.meshgrid(t_vals, x_vals)
        U = self.u_nn(T, X)

        plt.figure(figsize=(10, 6))
        plt.contourf(X, T, U.detach().numpy(), 100, cmap="viridis")
        plt.colorbar(label="u(t, x)")
        plt.xlabel("x")
        plt.ylabel("t")
        plt.title("2D plot of u(t, x)")
        plt.show()

    def plot_u_comparison(self, exact_solution, t_vals, x_vals):
        for t in [0.125, 0.25, 0.375]:
            X = torch.tensor(x_vals)
            u_pred = self.u_nn(torch.tensor([t]).repeat(len(X)), X).detach().numpy()
            u_exact = exact_solution(t, X)  # Assuming exact_solution is defined

            plt.figure(figsize=(10, 6))
            plt.plot(X, u_pred, label="Predicted")
            plt.plot(X, u_exact, label="Exact", linestyle="dashed")
            plt.xlabel("x")
            plt.ylabel("u(t, x)")
            plt.title(f"u(t, x) comparison at t={t}")
            plt.legend()
            plt.show()

    def plot_phi_comparison(self, exact_phi_solution, t_vals, x_vals):
        for t in [0.125, 0.25, 0.375]:
            X = torch.tensor(x_vals)
            phi_pred = self.phi_nn(torch.tensor([t]).repeat(len(X)), X).detach().numpy()
            phi_exact = exact_phi_solution(t, X)  # Assuming exact_phi_solution is defined

            plt.figure(figsize=(10, 6))
            plt.plot(X, phi_pred, label="Predicted")
            plt.plot(X, phi_exact, label="Exact", linestyle="dashed")
            plt.xlabel("x")
            plt.ylabel("phi(t, x)")
            plt.title(f"phi(t, x) comparison at t={t}")
            plt.legend()
            plt.show()
