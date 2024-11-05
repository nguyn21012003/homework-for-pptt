import matplotlib.pyplot as plt
import numpy as np

file = np.loadtxt("dataFT.txt")

time = file[:,0]

dist = file[:,1]



plt.plot(time,dist, color = 'red', linewidth = 2, linestyle = 'solid')
plt.xlabel("$\omega$")
plt.ylabel("$\\alpha(\omega)$")
plt.grid()
plt.savefig("alpha.png",dpi = 200)

