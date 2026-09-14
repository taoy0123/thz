import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

VecStart_x = [0,1,3,5]
VecStart_y = [2,2,5,5]
VecStart_z = [0,1,1,5]
VecEnd_x = [1,2,-1,6]
VecEnd_y = [3,1,-2,7]
VecEnd_z  =[1,0,4,9]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for i in range(4):
    ax.plot([VecStart_x[i], VecEnd_x[i]], [VecStart_y[i],VecEnd_y[i]],zs=[VecStart_z[i],VecEnd_z[i]])

plt.show()