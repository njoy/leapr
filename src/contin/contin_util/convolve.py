import numpy as np


def convol(t1,t2,delta):
    def getTnVal(tn,i,delta):
        if abs(i) >= len(tn):
            return 0.0
        if i < 0:
            return tn[abs(i)]*np.exp(i*delta)
        return tn[i]

    i_range = range(-len(t1)+1,len(t1))
    t3 = [0.0]*(2*len(t1)-1)

    for j in range(len(t3)):
        for i in i_range:
            if (j-i >= len(t2)): 
                break
            value = getTnVal(t1,i,delta)*getTnVal(t2,j-i,delta)
            t3[j] += value * delta * 0.5 if abs(i) == len(t1)-1 else value*delta
    return t3


#t1 = [0.1, 0.6, 0.4]
#t2 = [0.1, 0.6, 0.4, 0, 0]
#delta = 1.0
#t3 = convol(t1,t2,delta)
#print(np.array(t3))
"""
t1 = [0.1, 0.2, 0.2, 0.2, 0.4]
t2 = [0.1, 0.2, 0.2, 0.2, 0.4, 0, 0, 0, 0]
delta = 1.0
t3 = convol(t1,t2,delta)
print(np.array(t3))
"""



##############################################################################


def convol_under_construction(t1,t2,betaVals):
    def getTnVal_2(tn,i,delta):
        return 0.0 if abs(i) >= len(tn) else \
               tn[abs(i)]*np.exp(i*delta) if i < 0 else \
               tn[i]

    deltaVals = [betaVals[i+1]-betaVals[i] for i in range(len(betaVals)-1)]

    i_range = range(-len(t1)+1,len(t1))

    t3 = [0.0]*(2*len(t1)-1)

    for j in range(len(t3)):
        for i in i_range[:-1]:
            if (j-i >= len(t2)): break
            value_left = getTnVal_2(t1,i,deltaVals[i])*getTnVal_2(t2,j-i,deltaVals[i])
            value_right = getTnVal_2(t1,i+1,deltaVals[i])*getTnVal_2(t2,j-i-1,deltaVals[i])
            value = (value_left + value_right)*0.5*deltaVals[i]
            if (j == 0): print(value_left,value_right,value,t3[j])
            t3[j] += value 
    return t3

#t1 = [0.1, 0.6, 0.4]
#t2 = [0.1, 0.6, 0.4, 0, 0]
#betaVals = [0.0,1.0,2.0]
#t3 = convol_under_construction(t1,t2,betaVals)
#print(np.array(t3))

"""
t1 = [0.1, 0.2, 0.2, 0.2, 0.4]
t2 = [0.1, 0.2, 0.2, 0.2, 0.4, 0, 0, 0, 0]
betaVals = [0.0,1.0,2.0,3.0,4.0]
t3 = convol_under_construction(t1,t2,betaVals)
print(np.array(t3))
"""



t1 = [0.2, 0.6, 0.8, 2.0, 6.0, 8.0]
t2 = [0.2, 0.6, 0.8, 2.0, 6.0, 8.0,\
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
print()
delta = 0.03
t3 = convol(t1,t2,delta)
print(np.array(t3))

print()
betaVals = [i*delta for i in range(2*len(t1)-1)]
t3 = convol_under_construction(t1,t2,betaVals)
print(np.array(t3))






