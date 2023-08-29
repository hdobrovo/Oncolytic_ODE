import numpy as np
from scipy.integrate import odeint

#define all the variables within the five ODE's
kc = 979.525
knc = 979.525
lamc = 0.130252
lamnc = 0.130252
betc = 4.80261e-06
betnc = 4.80261e-06
delc = 0.0999653
delnc = 0.0999653
pc = 2.24909e+06
pnc = 2.24909e+06
c = 0.803196

#setting up the system of equations with the variables
def oncol(y,t):
    uc = y[0]
    ic = y[1]
    unc = y[2]
    inc = y[3]
    v = y[4]
    ducdt = lamc*uc*(1-(uc/kc)) - betc*uc*v
    dicdt = betc*uc*v - delc*ic
    duncdt = lamnc*unc*(1-(unc/knc)) - betnc*unc*v
    dincdt = betnc*unc*v - delnc*inc
    dvdt = pc*ic + pnc*inc - c*v
    return[ducdt, dicdt, duncdt, dincdt, dvdt]


#define your initial conditions
x0 = [217.957,0.0,950.,0.0,1e-4]

#time points
tp = np.linspace(0,250,10000)

#solve ODE
x = odeint(oncol,x0,tp)

#lists
ucl = x[:,0]
icl = x[:,1]
uncl = x[:,2]
incl = x[:,3]
vl = x[:,4]

datasetbvl = np.empty([10000,7])

count=0
for rl in np.logspace(0,4,100):
    kc = rl*kc
    for rb in np.logspace(-6,0,100):
        print(rl,rb)
        betnc = rb*betc
        x = odeint(oncol, x0, tp)
        ucl = x[:, 0]
        icl = x[:, 1]
        uncl = x[:, 2]
        incl = x[:, 3]
        vl = x[:, 4]
        dur, = np.where(vl>1)
        datasetbvl[count,:] =[rl,rb,max(vl),(dur[-1]-dur[0])/10000,ucl[-1],min(uncl),vl[-1]]
        count=count+1

np.savetxt("bvske-4.dat", datasetbvl)
