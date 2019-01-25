import matplotlib.pyplot as plt
import numpy as np



def interp(Tn,betas,beta):
    isNeg = beta < 0
    beta = abs(beta)
    if (beta < betas[0] or betas[-1] < beta): return 0.0
    for i in range(len(betas)-1):
        if betas[i] <= beta <= betas[i+1]:
            y1, y2 = Tn[i], Tn[i+1]
            x1, x2 = betas[i], betas[i+1]
            value = ((y2-y1)/(x2-x1)*beta + y1 - (y2-y1)/(x2-x1)*x1)
            return value*np.exp(beta) if isNeg else value
            

def getNextTn(betas,T1,Tn):
    def getT1Tn(b,beta): return T1[b] * interp(Tn,betas,(beta-betas[b]))
    T_next = []
    for beta in betas:
        T_next_val = 0.0
        for b in range(len(betas)-1):
            T_next_val += ( getT1Tn(b,beta)+getT1Tn(b+1,beta) ) * 0.5 * (betas[b+1]-betas[b])
        T_next.append(T_next_val)
    return T_next




# What is given by the user
phononGrid = [0.0, 0.1, 0.4, 0.6]
T1_orig    = [1.0, 2.0, 1.5, 2.5] 

# What is requested by the user
betaGrid = list(np.linspace(0,0.59,3))
T1_newGrid = [(interp(T1_orig,phononGrid,beta)) for beta in betaGrid] 

# Reflect so that T1 is defined for positive and negative beta
betas = [-x for x in betaGrid[::-1]]+betaGrid[1:]
negPart = T1_newGrid[None:0:-1]
T1 = [ negPart[i]*np.exp(betas[i+1]-betas[i]) for i in range(len(negPart))] + T1_newGrid[:]

T2 = getNextTn(betas,T1,T1)
T3 = getNextTn(betas,T1,T2)
T4 = getNextTn(betas,T1,T3)
T5 = getNextTn(betas,T1,T4)


print(T2)
#plt.plot(phononGrid,T1_orig,'go')
#plt.plot(betaGrid,T1)
plt.plot(betas,T1,'ro')
plt.plot(betas,T2,'yo')
plt.plot(betas,T3,'go')
plt.plot(betas,T4,'bo')
plt.plot(betas,T5,'ko')
plt.show()

























