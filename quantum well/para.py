from numpy import sqrt


def para(_):
    """Description

    Args:
        _ (int): Nhập số trường hợp để tính z0 từ bé tới lớn tương đương với thế nông, chiều dài giếng hẹp cho tới thế sâu, chiều dài giếng rộng,

    Returns:
        tuple: Chứa các parameter cần dùng để giải bài toán ví dụ như a,m,V0,z0,....
    """
    a0 = 1e-9  # m
    m = 0.067 * 9.11e-31  # kg
    Ce = 1.6e-19  # eV to Joule(J)
    hbar_ev = 6.582119569e-16  # eV*s

    hbar_Js = 1.054571817e-34  # J*s
    match _:
        case 1:  # z0 ~ 6.707666865904062
            V0 = 1
            a = 10 * a0
            z0 = sqrt(m * V0 * a * a / (2 * (hbar_ev) ** 2 * Ce))
            return a, m, Ce, V0, z0, hbar_ev, hbar_Js
        case 2:  # z0 ~ 12.779821793191344
            V0 = 3
            a = 11 * a0
            z0 = sqrt(m * V0 * a * a / (2 * (hbar_ev) ** 2 * Ce))
            return a, m, Ce, V0, z0, hbar_ev, hbar_Js
        case 3:  # z0 ~ 26.830667463616248
            V0 = 4
            a = 20 * a0
            z0 = sqrt(m * V0 * a * a / (2 * (hbar_ev) ** 2 * Ce))
            return a, m, Ce, V0, z0, hbar_ev, hbar_Js
        case 4:  # z0 ~ 42.42301016380012
            V0 = 10
            a = 20 * a0
            z0 = sqrt(m * V0 * a * a / (2 * (hbar_ev) ** 2 * Ce))
            return a, m, Ce, V0, z0, hbar_ev, hbar_Js
        case 5:  # z0 ~ 237.15183634105392
            V0 = 50
            a = 50 * a0
            z0 = sqrt(m * V0 * a * a / (2 * (hbar_ev) ** 2 * Ce))
            return a, m, Ce, V0, z0, hbar_ev, hbar_Js
        case 6:  # z0 ~ 1341.5333731808123
            V0 = 100
            a = 200 * a0
            z0 = sqrt(m * V0 * a * a / (2 * (hbar_ev) ** 2 * Ce))
            return a, m, Ce, V0, z0, hbar_ev, hbar_Js
        case 7:  # z0 ~ 42423.01016380013
            V0 = 0.1
            a = 0.2 * a0
            z0 = sqrt(m * V0 * a * a / (2 * (hbar_ev) ** 2 * Ce))
            return a, m, Ce, V0, z0, hbar_ev, hbar_Js
        case 8:  # z0 ~ 42423.01016380013
            V0 = 83e6
            a = 2.0e-15
            z0 = 0.094
            return a, m, Ce, V0, z0, hbar_ev, hbar_Js
