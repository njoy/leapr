
import numpy as np

x = [1,0,3, 0,3,4,56,8,5,3,7,8]
y = [2,1,4,0,1,3,4,6,73,4,5,6,2,6,3,4,5,6,7]

Lx = 0
Rx = 1
Ly = 0
Ry = 1

x_onRightSide = False
y_onRightSide = False
x_onLeftSide = True
y_onLeftSide = True 

z = []
for i in range(-len(x)+1,len(y)):
    xSlice = x[Lx:Rx]
    ySlice = y[Ly:Ry][::-1]
    z.append(np.dot(xSlice,ySlice))
    if (Rx != len(x)):
        Rx += 1
        Ry += 1
    elif (Ry != len(y)):
        Ry += 1
        Ly += 1
    else:
        Lx += 1
        Ly += 1



print(np.convolve(x,y))
print(np.array(z))
print(np.convolve(x,y)-(np.array(z)))

