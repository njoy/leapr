import numpy as np


def getTnVal(tn,i,delta):
    if abs(i) >= len(tn):
        return 0.0
    if i < 0:
        return tn[abs(i)]*np.exp(i*delta)
    return tn[i]

t1 = [0.1, 0.6, 0.4]
t2 = [0.1, 0.6, 0.4, 0, 0]
delta = 1.0

i_range = [-2,-1,0,1,2]
t3 = [0.0]*(2*len(t1)-1)

for j in range(len(t3)):
    for i in i_range:
        if (j-i >= len(t2)): 
            break
        value = getTnVal(t1,i,delta)*getTnVal(t2,j-i,delta)
        t3[j] += value * delta * 0.5 if abs(i) == len(t1)-1 else value*delta



print(np.array(t3))

##############################################################################


def getTnVal_2(tn,i,delta):
    if abs(i) >= len(tn):
        return 0.0
    if i < 0:
        return tn[abs(i)]*np.exp(i*delta)
    return tn[i]



t1 = [0.1, 0.6, 0.4]
t2 = [0.1, 0.6, 0.4, 0, 0]
betaVals = [0.0,1.0,2.0]
deltaVals = [betaVals[i+1]-betaVals[i] for i in range(len(betaVals)-1)]

i_range = [-2,-1,0,1,2]

t3 = [0.0]*(2*len(t1)-1)

for j in range(len(t3)):
    for i in i_range[:-1]:
        if (j-i >= len(t2)): 
            break
        value_left = getTnVal_2(t1,i,deltaVals[i])*getTnVal_2(t2,j-i,deltaVals[i])
        value_right = getTnVal_2(t1,i+1,deltaVals[i])*getTnVal_2(t2,j-i-1,deltaVals[i])
        value = (value_left + value_right)*0.5*deltaVals[i]
        t3[j] += value 

print(np.array(t3))









