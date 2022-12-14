from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import copy
import os 


#Import Data

HOME='/Users/montserrattrullsjuanola/Documents'

    
TEST=os.path.join(HOME, 'TEST.txt')
MISTERI=os.path.join(HOME, 'MISTERI1.txt')
MISTERI2=os.path.join(HOME, 'MISTERI2.txt')
MISTERI3=os.path.join(HOME, 'MISTERI3.txt')


TEST_array=np.loadtxt(TEST)
MISTERI_array=np.loadtxt(MISTERI)
MISTERI2_array=np.loadtxt(MISTERI2)
MISTERI3_array=np.loadtxt(MISTERI3)


TEST_spins=TEST_array[:,1]
MISTERI_spins=MISTERI_array[:,1]
MISTERI2_spins=MISTERI2_array[:,1]
MISTERI3_spins=MISTERI3_array[:,1]
    

#Defining functions------------------------------------------------------------

def GenerateLattice(L, spins):
    
    LATTICE=np.zeros([L,L])
    
    i=0

    for k in range(0,L,1):
        
        for l in range(0,L,1):
            
            LATTICE[k,l]=spins[i]
            
            i=i+1  
            
    return LATTICE
            

def CalculateMandE(L, LATTICE):
    
    N=L*L
    SUM=0

    for i in range(0,L,1):
    
        for j in range(0,L,1):
        
            Sn = LATTICE[(i - 1) % L, j] + LATTICE[(i + 1) % L, j] + \
                LATTICE[i, (j - 1) % L] + LATTICE[i, (j + 1) % L]
             
            SUM = SUM + LATTICE[i,j]*Sn
             
    M=1/N*np.sum(LATTICE)
    E=-1/2*SUM


    return M,E
    
#Results-----------------------------------------------------------------------

L=10

lattice=GenerateLattice(L, TEST_spins)
lattice2=GenerateLattice(L, MISTERI_spins)
lattice3=GenerateLattice(L, MISTERI2_spins)
lattice4=GenerateLattice(L, MISTERI3_spins)

M_test, E_test = CalculateMandE(L, lattice)
M_mystery, E_mystery = CalculateMandE(L, lattice2)
M_mystery2, E_mystery2 = CalculateMandE(L, lattice3)
M_mystery3, E_mystery3 = CalculateMandE(L, lattice4)

print(M_test, M_mystery, M_mystery2, M_mystery3, E_test, E_mystery, E_mystery2, E_mystery3)



    

             
        

        


        