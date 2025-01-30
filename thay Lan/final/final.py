import torch
import torch.nn as nn
import numpy as np
import matplotlib.pyplot as plt


def buildModel(input_dim, output_dim, hidden_dim):
    """
    ### nn.Sequential is a container for modules. Modules will be added to it in the order they are passed in the constructor.

    ### nn.Linear(input_dim, hidden_dim) is a linear transformation, y = xA^T + b, where x is the input, A is the weight matrix, and b is the bias.

    ### nn.Tanh() is the hyperbolic tangent non-linear activation function, y = tanh(x).

    ### There is no activation function after the last linear transformation, so the output is a linear combination of the input.
    Args:
        input_dim: _description_
        output_dim: _description_
        hidden_dim: _description_

    Returns:
        _description_
    """
    model = nn.Sequential(
        nn.Linear(input_dim, hidden_dim),
        nn.Tanh(),
        nn.Linear(hidden_dim, output_dim),
    )

    nn.Linear(input_dim, output_dim)
    return model


def fArr(model, x, EA, p):
    u = model(x)
    udx = derivative(u, x)
    EAudx = derivative(EA(x) * udx, x)
    f = EAudx + p(x)
    return f


def derivative(y, x):
    dydx = torch.autograd.grad(
        y,
        x,
        torch.ones(x.size()[0], 1),
        create_graph=True,
        retain_graph=True,
    )[0]

    return dydx


def main():
    model = buildModel(1, 1, 10)
    x = torch.linspace(0, 1, 100).view(-1, 1).requires_grad_()

    EA = lambda x: 1 + 0 * x
    p = lambda x: 4 * np.pi**2 * torch.sin(2 * np.pi * x)

    f = fArr(model, x, EA, p)

    u0 = 0
    u1 = 0
    loss_values = []

    optimizer = torch.optim.Adam(model.parameters(), lr=0.01)

    for epoch in range(1000):
        optimizer.zero_grad()
        f = fArr(model, x, EA, p)
        MSE_f = torch.mean(f**2)

        u0_pred = model(torch.tensor([[0.0]]))
        u1_pred = model(torch.tensor([[1.0]]))

        MSE_b = (u0_pred - u0) ** 2 + (u1_pred - u1) ** 2

        loss = MSE_f + MSE_b
        loss.backward()

        optimizer.step()
        loss_values.append(loss.item())

        # print(f"Epoch: {epoch}, Loss: {loss.item()}")

    # Plot model prediction
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax1.plot(x.detach().numpy(), model(x).detach().numpy(), label="Model Prediction")

    # Define and plot analytic solution
    analytic_solution = lambda x: torch.sin(2 * np.pi * x)
    ax1.plot(x.detach().numpy(), analytic_solution(x).detach().numpy(), label="Analytic Solution", linestyle="dashed")

    # Plot loss
    ax2 = fig.add_subplot(122)
    ax2.plot(loss_values, label="Loss")
    plt.xlabel("Epoch")
    plt.ylabel("Loss")

    plt.legend()
    plt.show()
    

if __name__ == "__main__":
    main()
