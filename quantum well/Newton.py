from numpy import sqrt, cos


df = lambda z, z0: 1 / (cos(z) ** 2) + z0**2 / (z**3 * sqrt(z0**2 / z**2 - 1))


def newton(f: float, p0: float, eps: float, N: int, z0: float) -> float:
    """Description

    Sử dụng thuật toán Newton Raphson cho f(x) = 0 cho ra nghiệm xấp xỉ quanh p0

    Parameters:
        f (function): function
                    Hàm cần tính toán.
        p0 (float): Giá trị đầu nằm trong khoảng [a,b].
        eps (float): Sai số
        N (int): Số vòng lặp
        z0 (float): Tham số.

    Returns:
        p_n: Nghiệm hay điểm cắt nhau hay điểm làm cho f(x) = 0.
    """
    p0_n = p0
    for n in range(0, N + 1):

        p0_n = p0_n - f(p0_n, z0) / df(p0_n, z0)
        if df(p0_n, z0) == 0:
            break
        if abs(f(p0_n, z0)) < eps:
            break

    return p0_n
