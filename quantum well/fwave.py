from numpy import exp, cos, tan, sqrt
import numpy as np
import matplotlib.pyplot as plt
import random

# def F(z):
#    alpha = ((1/2 + sin(2*z)/(4*z))*(exp(-2*z*tan(z))/(cos(z)**2)) + exp(-2*z*tan(z))/(2*z*tan(z)))*a
#    return sqrt(1/(2*alpha))


# def D(z):
#    return F(z)*exp(-z*tan(z))/cos(z)
def F(z):
    return exp(z * tan(z)) * cos(z) / sqrt((z * tan(z) + 1) / (z * tan(z) / a))


def D(z):
    return 1 / sqrt((z * tan(z) + 1) / (z * tan(z) / a))


def Psi2B(x):
    wave = []
    for i in range(len(z)):
        f = F(z)[i]
        Z = z[i]
        wave.append(f * exp(-Z * tan(Z) * x / a))
    return np.array(wave)


def Psi2A(x):
    wave = []
    for i in range(len(z)):
        f = F(z)[i]
        Z = z[i]
        wave.append(f * exp(Z * tan(Z) * x / a))
    return np.array(wave)


def Psi1(x):
    wave2 = []
    for i in range(len(z)):
        d = D(z)[i]
        Z = z[i]
        wave2.append(d * cos(Z * x / a))
    return np.array(wave2)


a = 0.529e-10
z = np.array(
    [
        1.5103754495117065,
        3.2590322977622193,
        6.342461412835801,
        9.462989360316994,
        12.593270514088786,
        15.727410465230284,
        18.863296250475514,
        21.999727035830396,
    ]
)


x1 = np.linspace(-a, a, 1000)
x2A = np.linspace(-10 * a, -a, 1000)
x2B = np.linspace(a, 10 * a, 1000)
v = np.sum(Psi1(x1), axis=0)
c = np.sum(Psi2A(x2A), axis=0)
l = np.sum(Psi2B(x2B), axis=0)
cc = (random.random(), random.random(), random.random())
plt.plot(x1, v, color=cc)
plt.plot(x2A, c, color=cc)
plt.plot(x2B, l, color=cc)
# for i in range(8):
#    cc=(random.random(), random.random(), random.random())
#    plt.plot(x1,Psi1(x1)[i],color=cc,label='z'+str(i)+' = '+str(z[i]))
#    plt.plot(x2A,Psi2A(x2A)[i],color=cc)
#    plt.plot(x2B,Psi2B(x2B)[i],color=cc)
plt.grid()
# plt.legend(loc='best')
plt.show()
