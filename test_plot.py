import matplotlib.pyplot as plt
import numpy as np

a=np.zeros((10,10))

for i in range(10):
	for j in range(10):
		a[i][j]	=	i*j

plt.figure()
for i in range(10):
	plt.plot(a[i])
plt.show()