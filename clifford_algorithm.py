# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 15:26:13 2014

@author: Tan Hao Qin
"""

#import gaussian_elimination
import numpy as np


def check_clifford_equivalance(matA, matB):
    assert matA.shape[0] == matA.shape[1],"The first matrix must be a square matrix. Please insert a square matrix."    
    assert matB.shape[0] == matB.shape[1],"The first matrix must be a square matrix. Please insert a square matrix."
    assert matA.shape[0] == matB.shape[0],"Comparison must be done between matrices of the same size"    
   
    clifford_array = generate_linear_systems(matA,matB)
    gaussian_elimination(clifford_array)       
   
    return True

def generate_linear_systems(matA,matB):
    '''
    linear system generated is generated with the following sequence:
    j = 1, k = 1 : a1 a2 [...] an-1 an b1 b2 [...] bn-1 bn c1 [...] cn d1 [...] dn
    j = 1, k = 2: a1 a2 [...] an-1 an b1 b2 [...] bn-1 bn c1 [...] cn d1 [...] dn    
    [...]
    j = 1, k = n : a1 a2 [...] an-1 an b1 b2 [...] bn-1 bn c1 [...] cn d1 [...] dn    
    j = 2, k = 1 : a1 a2 [...] an-1 an b1 b2 [...] bn-1 bn c1 [...] cn d1 [...] dn
    [...]
    j = n, k = n : a1 a2 [...] an-1 an b1 b2 [...] bn-1 bn c1 [...] cn d1 [...] dn
    
    Assuming matA and matB to be n*n matrices, an output matrix generated will be of size 
    n^2(number of equations) by 
    4*n(number of coefficients)    
    '''    
    n =  matA.shape[0]    
    #Generating a zero matrix of n^2 by 4n dimensions    
    output = np.zeros((n**2,4*n),dtype = float)    
    for j in range(n):
        for k in range(n):
            equationNumber = j*n+k
            #Generating a-coefficients
            output[equationNumber,k] = matA[j,k]
            #Generating d-coefficients
            output[equationNumber,3*n+j] = matB[j,k]
            #Generating c-coefficients
            for i in range(n):
                output[equationNumber,2*n+i] = matA[i,j]*matB[i,k]
            #Generating b-coefficients
            if j==k:
                output[equationNumber,n+j] = 1                
    return output
    
def parse_result(matrix):
    m,n = matrix.shape
    k = m/4
    for i in range(m):
        matrix[i]


            
'''        
matA = np.array([[0,1,0],[1,0,1],[0,1,0]])
matB = np.array([[0,1,1],[1,0,1],[1,1,0]])
matC = generate_linear_systems(matA,matB)
print matC
'''

