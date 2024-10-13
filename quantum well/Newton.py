from numpy import sqrt, cos


df = lambda z, z0: 1 / (cos(z) ** 2) + z0**2 / (z**3 * sqrt(z0**2 / z**2 - 1))


def newton(f: float, p0: float, eps: float, N: int, z0: float) -> float:
    """_summary_

    Args:
        f (float): _description_
        p0 (float): _description_
        eps (float): _description_
        N (int): _description_
        z0 (float): _description_

    Returns:
        float: _description_
    """
    p0_n = p0
    for n in range(0, N + 1):

        p0_n = p0_n - f(p0_n, z0) / df(p0_n, z0)
        if df(p0_n, z0) == 0:
            break
        if abs(f(p0_n, z0)) < eps:
            break

    return p0_n
