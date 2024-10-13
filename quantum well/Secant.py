def secant(f: float, p0: float, p1: float, eps: float, N: int, z0: float) -> float:
    """Description

    Sử dụng thuật toán Secant cho f(x) = 0 cho ra nghiệm xấp xỉ trên đoạn [p0,p1].

    Parameters:
        f (function): function
                    Hàm cần tính toán.
        p0,p1 (float): Giá trị đầu nằm trong khoảng [a,b].
        eps (float): Sai số
        N (int): Số vòng lặp
        z0 (float): Tham số.

    Returns:
        p_n: Nghiệm hay điểm cắt nhau hay điểm làm cho f(x) = 0. Nếu như hàm f(x) tại các điểm p0,p1,p_n cùng dấu nhau thì chúng ta được nghiệm không hội tụ với n vòng lặp
    """
    if f(p0, z0) * f(p1, z0) >= 0:
        return None
    p0_n = p0
    p1_n = p1
    for n in range(N + 1):
        p_n = p0_n - f(p0_n, z0) * (p1_n - p0_n) / (f(p1_n, z0) - f(p0_n, z0))

        if f(p0_n, z0) * f(p_n, z0) < 0:
            p0_n = p0_n
            p1_n = p_n
        elif f(p1_n, z0) * f(p_n, z0) < 0:
            p0_n = p_n
            p1_n = p1_n

        if abs(f(p_n, z0)) < eps:
            break

        if abs(p1_n - p0_n) < eps:
            break

        if f(p_n, z0) == 0:
            break
    return p_n
