from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import copy 



def init_lattice(L):

    
    lattice = np.random.choice([1, -1], size=(L, L))
    return lattice

def deltaE(S0, Sn, J, H):

    
    return 2 * S0 * (H + J * Sn)




def ising(H,J,T,lattice):
    

    
    for step in range(n):

        i = np.random.randint(L)
        j = np.random.randint(L)
        
        # Periodic Boundary Condition
        Sn = lattice[(i - 1) % L, j] + lattice[(i + 1) % L, j] + \
             lattice[i, (j - 1) % L] + lattice[i, (j + 1) % L]

        dE = deltaE(lattice[i, j], Sn, J, H)

        if dE <= 0:
            lattice[i, j] = -lattice[i, j]
    
    
    return lattice

def isingCTMC(L, lattice, t, H, J):  
    
    lattice2=np.zeros((L,L))
    V = np.zeros(5)
    W = np.array([1, 0, 1, 0, 1])*(1/n)
    
    for i in range(0,L,1):
        
        for j in range(0,L,1):

            S0 = lattice[i, j]
            
            
            Sn = np.array([lattice[(i - 1) % L, j], lattice[(i + 1) % L, j],
                lattice[i, (j - 1) % L], lattice[i, (j + 1) % L]])
            
            
            B = np.count_nonzero(Sn == S0)
            
            if B == 0.:
                
                lattice2[i,j] = 0
                V[0] += 1
                 
            if B == 4.:
                
                lattice2[i,j] = 1
                V[1] += 1
            
            if B == 1.:
                
                lattice2[i,j] = 2
                V[2] +=1

            if B == 3.:
                
                lattice2[i,j] = 3
                V[3] += 1

                
            if B == 2.:
                
                lattice2[i,j] = 4
                V[4] += 1
    
    E = np.sum(V*W)    
    F = 1-E
    
    if (F == 1 or F ==0):
        
        dt = 1
        t += dt

    
    
    else:
        
        H = 1-np.random.rand()
        dt = 1+int(np.log(H)/np.log(F))
    
        t += dt
        op = V*W  
        
        a = np.random.rand()*np.sum(op)  
        
        if a <= op[0]:
            
            d = 0
        
        if op[0] < a <= op[0]+op[2]:
            
            d = 2
            
        if op[0]+op[2] < a:
            
            d=4
            
        i = np.random.randint(0,10)
        j = np.random.randint(0,10)
        
        while lattice2[i,j] != d:  
            
            i = np.random.randint(0,10)
            j = np.random.randint(0,10) 
        
        lattice[i,j] = -lattice[i,j]

    return t, lattice



L=10
MTCsteps=1000
n=L*L
H=0
J=1
T=1
a=10
t=0
temps1=np.linspace(0,MTCsteps,MTCsteps)
tplot=np.zeros(MTCsteps) 


lattice = np.zeros((L,L))
MTOT=[]
MPROM=[]
MTOT2=[]
MPROM2=[]
A=[]
Att=[]

for i in range(a):
    
    lattice = init_lattice(L)
    M_list=[]
    M_list.append(1/n*np.sum(lattice))
    MTOT.append(1/n*np.sum(lattice))
    
    
    for k in range(1,MTCsteps,1):
        
        P=ising(H,J,T,lattice)
        suma=np.sum(P)
        M=1/n*suma
        Mabs=np.abs(M)
        MTOT.append(Mabs)
        M_list.append(M)
        
    
        
    res=[abs(ele) for ele in M_list]
    plt.plot(temps1,res, color='red', label='FSMC')
    plt.xlabel("t")
    plt.ylabel("M")
    plt.xscale('log')
    plt.yscale('log')
    
    


    
for jj in range(a):
        
    lattice = init_lattice(L)
    M_list2=[]
    M_list2.append(1/n*np.sum(lattice))
    tplot=np.zeros(MTCsteps)
    A.append(0)
    t=0
    
    for l in range(1,MTCsteps,1):
        
        t, lattice = isingCTMC(L, lattice, t,H,J)
        M2 = (1/n)*np.sum(lattice)
        M2abs=np.abs(M2)
        MTOT2.append(M2abs)
        tplot[l]=t
        A.append(t)
        M_list2.append(M2)
        
    res2=[abs(ele) for ele in M_list2]
    plt.plot(tplot, res2, color='green', label='CMTC')
    plt.xlabel("t")
    plt.ylabel("M")
    plt.xscale('log')
    plt.yscale('log')
    
plt.show()

for ii in range(0,1000,1):
    
    MPROMT=1/a*(MTOT[ii]+MTOT[1000+ii]+MTOT[2000+ii]+MTOT[3000+ii]+MTOT[4000+ii]+MTOT[5000+ii]+MTOT[6000+ii]+MTOT[7000+ii]+MTOT[8000+ii]+MTOT[9000+ii])
    MPROM.append(MPROMT)
    
plt.plot(temps1,MPROM, label='FTMC')
plt.xscale('log')
plt.yscale('log')


for ii in range(0,1000,1):
    
    MPROMT2=1/a*(MTOT2[ii]+MTOT2[1000+ii]+MTOT2[2000+ii]+MTOT2[3000+ii]+MTOT2[4000+ii]+MTOT2[5000+ii]+MTOT2[6000+ii]+MTOT2[7000+ii]+MTOT2[8000+ii]+MTOT[9000+ii])
    MPROM2.append(MPROMT2)
    At=1/a*(A[ii]+A[1000+ii]+A[2000+ii]+A[3000+ii]+A[4000+ii]+A[5000+ii]+A[6000+ii]+A[7000+ii]+A[8000+ii]+A[9000+ii])
    Att.append(At)
      
plt.plot(Att,MPROM2, label='CTMC')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.xlabel('t')
plt.ylabel('<M>')
plt.show()
    
    

    
    
    
