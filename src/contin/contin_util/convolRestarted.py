import matplotlib.pyplot as plt
import numpy as np



def interp(Tn,betas,beta):
    if (beta < betas[0] or betas[-1] < beta): return 0.0
    for i in range(len(betas)-1):
        if betas[i] <= beta <= betas[i+1]:
            y1, y2 = Tn[i], Tn[i+1]
            x1, x2 = betas[i], betas[i+1]
            value = ((y2-y1)/(x2-x1)*beta + y1 - (y2-y1)/(x2-x1)*x1)
            return value
            

def getNextTn(betas,T1,Tn):
    def getT1Tn(b,beta): return T1[b] * interp(Tn,betas,(beta-betas[b]))
    T_next = []
    for beta in betas:
        T_next_val = 0.0
        for b in range(len(betas)-1):
            T_next_val += ( getT1Tn(b,beta)+getT1Tn(b+1,beta) ) * 0.5 * (betas[b+1]-betas[b])
        T_next.append(T_next_val)
    return T_next



T1_orig = [0.01,0.04,0.09,0.11,0.16,0.21]
phononGrid = [0,0.5,1.0,1.5,2.0,2.5]
betas = [-x for x in phononGrid[::-1]]+phononGrid[1:]
negPart = T1_orig[None:0:-1]
T1 = [ negPart[i]*np.exp(betas[i]) for i in range(len(negPart))] + T1_orig[:]
print(betas)
print(T1)
print("T1\n","%.5f"%T1[5],"%.5f"%T1[6],"%.5f"%T1[7],"%.5f"%T1[8],"%.5f"%T1[9],"%.5f"%T1[10])
print("\n\n")
T2 = getNextTn(betas,T1,T1)
print("T2\n",T2,"\n")
print("T2\n","%.5f"%T2[5],"%.5f"%T2[6],"%.5f"%T2[7],"%.5f"%T2[8],"%.5f"%T2[9],"%.5f"%T2[10])
print()


"""
alphaGrid = [0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28]
betaGrid = [0.0,0.15,0.30,0.60,1.20]
T1_newGrid = [(interp(T1_orig,phononGrid,beta)) for beta in betaGrid] 


# Reflect so that T1 is defined for positive and negative beta
betas = [-x for x in betaGrid[::-1]]+betaGrid[1:]
negPart = T1_newGrid[None:0:-1]
T1 = [ negPart[i]*np.exp(betas[i+1]-betas[i]) for i in range(len(negPart))] + T1_newGrid[:]


Tlast = T1[:]
Tnext = T1[:]
denom = 1.0

sab = T1[:]
max_n = 3

for n in range(max_n+1):

    if (n>0):
        Tnext = getNextTn(betas,T1,Tlast)
    Tlast = Tnext
    for b in range(len(betas)):
        sab[b] += (Tnext[b]/np.math.factorial(n))*(alphaGrid[0]*lambda_s)**n
    print(sab[0])
print(sab)


"""
























