def bisection(f: float, a: float, b: float, eps: float, N: int, z0: float) -> float:
    """Description

    Sử dụng thuật toán Bisection để tìm nghiêm xấp xỉ của phương trình f(x) = 0 nằm trong [a,b].

    Parameters:
        f (function): function
                    Hàm cần tính toán.
        a,b (float): Giá trị đầu nằm trong khoảng [a,b].
        eps (float): Sai số
        N (int): Số vòng lặp
        z0 (float): Tham số.

    Returns:
        c: Nghiệm hay điểm cắt nhau hay điểm làm cho f(x) = 0. Nếu như hàm f(x) tại các điểm p0,p1,p_n cùng dấu nhau thì chúng ta được nghiệm không hội tụ với n vòng lặp
    """
    for i in range(N):
        c = (a + b) / 2
        if abs(f(c, z0)) < eps:
            break
        if f(a, z0) * f(c, z0) < 0:
            b = c
        elif f(c, z0) * f(b, z0) < 0:
            a = c
        else:
            break
    return c
