import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def magnitude(a):
	mag	=	np.sqrt((a[0]**2)+(a[1]**2)+(a[2]**2))
	return mag

a	=	np.ones(3)
s	=	np.ones(3)
s	=	s/magnitude(s)
a[2]	=	0
a	=	-1*a/magnitude(a)
T		=	np.cross(a,s)
T		=	np.cross(s,T)
T		=	T/magnitude(T)
print(T)
print(np.dot(s,[-0.57735027,-0.57735027,0.57735027]))
origin = [0],[0],[0]
V = np.array([T,s,a])
fig=plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.quiver(origin[0],origin[1],origin[2], V[:,1], V[:,2], V[:,0],color=['r','b','g'])
ax.set_xlim([-1, 1])
ax.set_ylim([-1, 1])
ax.set_zlim([-1, 1])

plt.show()

