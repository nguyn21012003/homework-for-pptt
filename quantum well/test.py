import numpy as np
from numpy import sqrt, sin, cos, tan, exp
from math import pi
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider
import csv
from tabulate import tabulate

L = 10  # well width in nm
Vo = 1  # Well depth in eV
m = 0.067

theta0 = sqrt(m * 9.31e-31 * Vo * L * L / (2 * (6.58e-16) ** 2 * 1e18 * 1.6e-19))

print(theta0)
