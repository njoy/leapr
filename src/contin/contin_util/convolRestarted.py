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
alpha_scaling = 1.46796
tev = 0.01723477

# Calculated earlier in contin
T1_orig = [0.0010065, 0.00585766, 0.0146, 0.0194665, 0.0729995, 0.116799]
phononGrid = [0,0.1,0.2,0.3,0.4,0.5]
lambda_s = 0.0438473

# Reflect beta values and T1 values
betas = [-x for x in phononGrid[::-1]]+phononGrid[1:]
negPart = T1_orig[None:0:-1]
T1 = [ negPart[i]*np.exp(betas[i]) for i in range(len(negPart))] + T1_orig[:]

# Scale phonon grid
phononGrid = [x/tev for x in phononGrid]



print()
print("T1\n","%.7E"%T1[5],"%.7E"%T1[6],"%.7E"%T1[7],"%.7E"%T1[8],"%.7E"%T1[9],"%.7E"%T1[10])
print()
T2 = getNextTn(betas,T1,T1)
#print("T2\n",T2,"\n")
print("T2\n","%.7E"%T2[5],"%.7E"%T2[6],"%.7E"%T2[7],"%.7E"%T2[8],"%.7E"%T2[9],"%.7E"%T2[10])
print()
T3 = getNextTn(betas,T1,T2)
#print("T3\n",T3,"\n")
print("T3\n","%.7E"%T3[5],"%.7E"%T3[6],"%.7E"%T3[7],"%.7E"%T3[8],"%.7E"%T3[9],"%.7E"%T3[10])
print()

T1_interp = [(interp(T1_orig,phononGrid,beta)) for beta in betaGrid] 
T2_interp = [(interp(T2,phononGrid,beta)) for beta in betaGrid] 
T3_interp = [(interp(T3,phononGrid,beta)) for beta in betaGrid] 

sab_alpha0 = []
alpha_term = np.exp(-alphaGrid[0]*alpha_scaling*lambda_s)
for b in range(len(betaGrid)):
    T1_term = alpha_term * (1/1) * (alphaGrid[0]*alpha_scaling*lambda_s)**1 * T1_interp[b]
    T2_term = alpha_term * (1/2) * (alphaGrid[0]*alpha_scaling*lambda_s)**2 * T2_interp[b]
    T3_term = alpha_term * (1/6) * (alphaGrid[0]*alpha_scaling*lambda_s)**3 * T3_interp[b]
    sab_alpha0.append(T1_term+T2_term+T3_term)

print(sab_alpha0)
print()

sab_alpha1 = []
alpha_term = np.exp(-alphaGrid[1]*alpha_scaling*lambda_s)
for b in range(len(betaGrid)):
    T1_term = alpha_term * (1/1) * (alphaGrid[1]*alpha_scaling*lambda_s)**1 * T1_interp[b]
    T2_term = alpha_term * (1/2) * (alphaGrid[1]*alpha_scaling*lambda_s)**2 * T2_interp[b]
    T3_term = alpha_term * (1/6) * (alphaGrid[1]*alpha_scaling*lambda_s)**3 * T3_interp[b]
    sab_alpha1.append(T1_term+T2_term+T3_term)

print(sab_alpha1)
print()


"""


T1_orig = [0.0010065, 0.00585766, 0.0146, 0.0194665, 0.0729995, 0.116799]
phononGrid = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]

alphaGrid = [0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28]
betaGrid = [0.0,0.15,0.30,0.60,1.20]
T1_newGrid = [(interp(T1_orig,phononGrid,beta)) for beta in betaGrid] 



# Reflect so that T1 is defined for positive and negative beta
betas = [-x for x in betaGrid[::-1]]+betaGrid[1:]
negPart = T1_newGrid[None:0:-1]
T1 = [ negPart[i]*np.exp(betas[i+1]-betas[i]) for i in range(len(negPart))] + T1_newGrid[:]

print(T1)
T2 = getNextTn(betas,T1_orig,T1_orig)
print(T2)

Tlast = T1[:]
Tnext = T1[:]
denom = 1.0

sab = T1[:]
max_n = 3
"""

"""
for n in range(max_n+1):

    if (n>0):
        Tnext = getNextTn(betas,T1,Tlast)
    Tlast = Tnext
    for b in range(len(betas)):
        sab[b] += (Tnext[b]/np.math.factorial(n))*(alphaGrid[0]*lambda_s)**n
    print(sab[0])
print(sab)
"""


























