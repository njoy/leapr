import matplotlib.pyplot as plt
import numpy as np



def interpolate(x,y,val):
    if val < x[0] or val > x[-1]:
        return 0.0
    for i in range(len(x)-1):
        if x[i] <= val <= x[i+1]:
            m = (y[i+1]-y[i])/(x[i+1]-x[i])
            b = y[i]
            return m*(val-x[i])+b


beta_vec_orig = [-2, 0, 2]
beta_vec_orig = [ 0, 2, 4]
beta_vec_orig = [ -1, 0, 1]

y1 = [1.0, 2.0, 3.0]
y2 = [2.0, 4.0, 1.0]


numZeros = 4
# find max spaing
maxSpacing = beta_vec_orig[1]-beta_vec_orig[0]
beta_vec = beta_vec_orig
for i in range(numZeros):
    beta_vec = [beta_vec[0]-maxSpacing] + beta_vec + [beta_vec[-1]+maxSpacing]
y1 = [0.0]*numZeros + y1 + [0.0]*numZeros
y2 = [0.0]*numZeros + y2 + [0.0]*numZeros
y3 = [0.0]*len(y1)
y4 = [0.0]*len(y1)

print()

for b,beta in enumerate(beta_vec):
    integrand = [0.0]*len(y1)
    for b_p,beta_p in enumerate(beta_vec):
        y1_piece = interpolate(beta_vec,y1,beta_p)
        y2_piece = interpolate(beta_vec,y2,beta-beta_p)
        integrand[b_p] = y1_piece*y2_piece
        #print("Beta'",beta_p,"      Beta-Beta'",beta-beta_p)
        #print("     ",interpolate(beta_vec,y1,beta_p),"      Beta-Beta'",interpolate(beta_vec,y2,beta-beta_p))
        #print()

    y3[b] = np.trapz(integrand,x=beta_vec)




print()
print(y3)
print(y4)
print()
print(np.convolve(y1,y2))
print(np.convolve(y1,y3))

#plt.plot(beta_vec,y1)
#plt.plot(beta_vec,y2)
#plt.plot(beta_vec,y3)
#plt.show()


print()






