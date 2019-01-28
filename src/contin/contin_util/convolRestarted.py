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


# Requested by the users
alphaGrid = [0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28]
betaGrid = [0.0,0.15,0.30,0.60,1.20]
scaling = 1.46796
tev = 0.01723477
nphon = 3

# Calculated earlier in contin
T1_orig = [0.0010065, 0.00585766, 0.0146, 0.0194665, 0.0729995, 0.116799]
phononGrid = [0,0.1,0.2,0.3,0.4,0.5]
lambda_s = 0.0438473

# Scale phonon grid
phononGrid = [x/tev for x in phononGrid]

# Make new grid
T1_newGrid = [(interp(T1_orig,phononGrid,beta)) for beta in betaGrid] 

# Reflect beta values and T1 values
betas = [-x for x in betaGrid[::-1]]+betaGrid[1:]
negPart = T1_newGrid[None:0:-1]
T1 = [ negPart[i]*np.exp(betas[i]) for i in range(len(negPart))] + T1_newGrid[:]
print(T1)

T_now = T1[:]
T_next = T1[:]


sab = [[0.0]*len(betaGrid) for a in range(len(alphaGrid))]
factorialTerm = 1.0

for n in range(1,nphon+1):
    if n > 1:
        factorialTerm /= n
        T_next = getNextTn(betas,T1,T_now)
    for a in range(len(alphaGrid)):
        alpha_term = np.exp(-alphaGrid[a]*scaling*lambda_s)
        for b in range(len(betaGrid)):
            sab[a][b] += alpha_term*factorialTerm*(alphaGrid[a]*scaling*lambda_s)**n*T_next[b+int(len(T_next)/2)]
    T_now = T_next


print()
print(sab[0])
print()
print(sab[1])
print()



































