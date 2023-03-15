import numpy as np
import matplotlib.pyplot as plt

x = np.array([0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000])
y = np.array([0, 50, 55, 60, 65, 72, 78, 83, 85, 85, 82, 77, 72])
p = np.polyfit(x,y,9)
plt.plot(x,y,'ks-')
plt.grid()
plt.show()
print(p)