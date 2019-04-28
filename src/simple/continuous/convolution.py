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

T1 = [1.0, 2.0, 3.0]
T2 = [2.0, 4.0, 1.0]


num0s = 4
# find max spaing
maxSpacing = beta_vec_orig[1]-beta_vec_orig[0]
beta_vec = beta_vec_orig
for i in range(num0s):
    beta_vec = [beta_vec[0]-maxSpacing] + beta_vec + [beta_vec[-1]+maxSpacing]
T1 = [0.0]*num0s + T1 + [0.0]*num0s
T2 = [0.0]*num0s + T2 + [0.0]*num0s

print()

def getTnext(beta_vec,T1,Tn):
    Tnext = [0.0]*len(T1)
    for b,beta in enumerate(beta_vec):
        integrand = [0.0]*len(T1)
        # Since T1(b') will only be nonzero for b' on the original grid, I know
        # that b' only has to traverse beta_vec_orig. 
        for b_p,beta_p in enumerate(beta_vec_orig):
            T1_piece = interpolate(beta_vec,T1,beta_p)
            Tn_piece = interpolate(beta_vec,Tn,beta-beta_p)
            integrand[num0s+b_p] = T1_piece*Tn_piece
        Tnext[b] = np.trapz(integrand,x=beta_vec)
    return Tnext

T3 = getTnext(beta_vec,T1,T2)
T4 = getTnext(beta_vec,T1,T3)


print()
print(T3)
print(T4)
print()
print(np.convolve(T1,T2))
T3 = np.convolve(T1,T2)
print(np.convolve(T1,T3))

#plt.plot(beta_vec,y1)
#plt.plot(beta_vec,y2)
#plt.plot(beta_vec,y3)
#plt.show()


print()






