from numpy import sqrt

a0 = 5.29e-11  ## Hang so Bohr
m = 9.31e-31  # kg
Ce = 1.6e-19  # eV to Joule(J)
hbar_ev = 6.582119569e-16  # eV*s

hbar_Js = 1.054571817e-34  # J*s


def case_(_):
    match _:
        case 1:  # z0 ~ 3.468005224411644
            V0 = 100
            a = 50 * a0
            z0 = sqrt(0.064 * m * V0 * a * a / (2 * (hbar_ev) ** 2 * Ce))
            return a, m, Ce, V0, z0, hbar_ev, hbar_Js
        case 2:  # z0 ~ 109.6679544650417
            V0 = 40
            a = 50 * a0
            z0 = sqrt(0.064 * m * V0 * a * a / (2 * (hbar_ev) ** 2 * Ce))
            return a, m, Ce, V0, z0, hbar_ev, hbar_Js
        case 3:  # z0 ~ 109.6679544650417
            V0 = 1
            a = 15 * a0
            z0 = sqrt(0.064 * m * V0 * a * a / (2 * (hbar_ev) ** 2 * Ce))
            return a, m, Ce, V0, z0, hbar_ev, hbar_Js
