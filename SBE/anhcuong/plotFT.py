import matplotlib.pyplot as plt
import numpy as np

file = np.loadtxt("dataFT.txt")
time = file[:,0]

dist = file[:,2]

dist1 = file[:,3]

polari = file[:,4]

polari1 = file[:,5]

plt.subplot(2,1,1)
plt.plot(time,dist, color = 'blue', linewidth = 1, linestyle = 'solid',label = "Real")
plt.plot(time,dist1, color = 'red', linewidth = 1, linestyle = 'solid',label = "Imag")
plt.xlabel("$\omega$")
plt.ylabel("$P(\omega)$")
plt.legend()
plt.grid()

plt.subplot(2,1,2)
plt.plot(time,polari,color = 'red', linewidth = 2, linestyle = 'solid',label = "Real")
plt.plot(time,polari1,color = 'green', linewidth = 2, linestyle = 'solid',label = "Imag")
plt.xlabel("$\omega$")
plt.ylabel("E($\omega$)")
plt.legend()
plt.grid()
plt.savefig("PEw.png",dpi = 200)

