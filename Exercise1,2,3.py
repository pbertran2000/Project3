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



L=50
lattice=init_lattice(L)
lattice1 = copy.deepcopy(lattice)


Xup = []
Yup = []

Xdown = []
Ydown = []

for i in range(L):
    
    for j in range(L):
        
        if lattice[i,j] == 1:
            Xup.append(i)
            Yup.append(j)
            
        if lattice[i,j]== -1:
            Xdown.append(i)
            Ydown.append(j)

plt.plot(Xup, Yup, 'o', marker = r'$\uparrow$', 
         markersize=5)            
plt.plot(Xdown, Ydown, 'o', marker = r'$\downarrow$', 
         markersize=5)
plt.show()




L=10
lattice=init_lattice(L)
lattice1 = copy.deepcopy(lattice)
MTCsteps=1000
n=L*L
H=0
J=1
T=1
M_list=[]
M_list.append(1/n*np.sum(lattice)) 
temps1=np.linspace(0,MTCsteps,MTCsteps)


for i in range(1,MTCsteps,1):
    
    P=ising(H,J,T,lattice)
    suma=np.sum(P)
    M=1/n*suma
    
    M_list.append(M)
   

res=[abs(ele) for ele in M_list]
plt.plot(temps1,res)
plt.xlabel("t")
plt.ylabel("M")
plt.xscale('log')
plt.yscale('log')
plt.show()

res=[abs(ele) for ele in M_list]
plt.plot(temps1,res)
plt.xlabel("t")
plt.ylabel("M")
plt.show()



#------------------------------------------------------------------------------

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





t=0
t_CTMC=np.zeros(MTCsteps) 
M_list2=[]
M_list2.append(1/n*np.sum(lattice))




for k in range(1, MTCsteps):
    
    t, lattice = isingCTMC(L, lattice1, t,H,J)
    
    M = (1/n)*np.sum(lattice)
    t_CTMC[k]=t
    M_list2.append(M)
    


res=[abs(ele) for ele in M_list]
plt.plot(temps1,res, label='FSMC')
plt.xlabel("t")
plt.ylabel("M")


res2=[abs(ele) for ele in M_list2]
plt.plot(t_CTMC, res2, label='CTMC')
plt.xlabel("t")
plt.ylabel("M")
plt.xscale('log')
plt.yscale('log')

plt.legend()
plt.show()


    

    
    
    
    
