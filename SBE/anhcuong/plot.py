import matplotlib.pyplot as plt
import numpy as np

file = np.loadtxt("data.txt")

time = file[:,0]

dist = file[:,1]

polari = file[:,2]

plt.subplot(2,1,1)
plt.plot(time,dist, color = 'red', linewidth = 2, linestyle = 'solid')
plt.xlabel("time")
plt.ylabel("N(t)")
plt.grid()

plt.subplot(2,1,2)
plt.plot(time,polari,color = 'red', linewidth = 2, linestyle = 'solid')
plt.xlabel("time")
plt.ylabel("P(t)")
plt.grid()
plt.savefig("fig_polari.png",dpi = 100)

