from numpy import sqrt


def para(_):
    """_summary_

    Args:
        _ (_type_): _description_

    Returns:
        _type_: _description_
    """
    a0 = 1e-9
    m = 9.31e-31  # kg
    Ce = 1.6e-19  # eV to Joule(J)
    hbar_ev = 6.582119569e-16  # eV*s

    hbar_Js = 1.054571817e-34  # J*s
    match _:
        case 1:  # z0 ~ 6.707666865904062
            V0 = 1
            a = 10 * a0
            z0 = sqrt(0.067 * m * V0 * a * a / (2 * (hbar_ev) ** 2 * Ce))
            return a, m, Ce, V0, z0, hbar_ev, hbar_Js
        case 2:  # z0 ~ 12.779821793191344
            V0 = 3
            a = 11 * a0
            z0 = sqrt(0.067 * m * V0 * a * a / (2 * (hbar_ev) ** 2 * Ce))
            return a, m, Ce, V0, z0, hbar_ev, hbar_Js
        case 3:  # z0 ~ 26.830667463616248
            V0 = 4
            a = 20 * a0
            z0 = sqrt(0.067 * m * V0 * a * a / (2 * (hbar_ev) ** 2 * Ce))
            return a, m, Ce, V0, z0, hbar_ev, hbar_Js
        case 4:  # z0 ~ 42.42301016380012
            V0 = 10
            a = 20 * a0
            z0 = sqrt(0.067 * m * V0 * a * a / (2 * (hbar_ev) ** 2 * Ce))
            return a, m, Ce, V0, z0, hbar_ev, hbar_Js
        case 5:  # z0 ~ 42.42301016380012
            V0 = 50
            a = 50 * a0
            z0 = sqrt(0.067 * m * V0 * a * a / (2 * (hbar_ev) ** 2 * Ce))
            return a, m, Ce, V0, z0, hbar_ev, hbar_Js
