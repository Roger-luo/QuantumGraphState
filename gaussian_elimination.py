# -*- coding: utf-8 -*-
 """
Created on Tue Oct 28 14:52:12 2014

@author: Tan Hao Qin
"""
#this class is written to perform gaussian elimination on matrices

import Quantum_270914 as Quantum
import numpy as np
import itertools

def check_clifford_equivalance(matA, matB):
    assert matA.shape[0] == matA.shape[1],"The first matrix must be a square matrix. Please insert a square matrix."    
    assert matB.shape[0] == matB.shape[1],"The first matrix must be a square matrix. Please insert a square matrix."
    assert matA.shape[0] == matB.shape[0],"Comparison must be done between matrices of the same size"    
   
    clifford_array = generate_linear_systems(matA,matB)
    row_echelon_form(clifford_array)       
    return resolve_ref(clifford_array)

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
    output = np.zeros((n**2,4*n),dtype = int)    
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



def swap_rows(array,row1,row2):
    #swap two rows in a matrix
    temp = array[row1].copy()
    array[row1] = array[row2]
    array[row2] = temp

def multiply_row(array,row1,constant):
    #multiplies a row in a matrix by a constant
    array[row1]*=np.array(constant)  

def row_addition(array,row1,row2,constant):
    #adds the second row multiplied by a constant into the first row
    array[row1] = array[row1]+constant*array[row2]

def find_non_zero_entry_in_column(array,column_index,startIndex = 0):
    #finds the first row in an array that has a non-zero entry at column j starting from row startIndex
    #returns row number
    m,n = array.shape
    for i in range(startIndex,m):
        if array[i][column_index] != 0:
            return i
    else:
        return -1
        
def ensure_all_other_row_entries_zero(array,row_index,column_index):
    x = np.transpose(array)[column_index].copy()
    for i in range(array.shape[0]):
        if i != row_index:  
            row_addition(array,i,row_index,-1*(float(x[i])/x[row_index]))            

def row_echelon_form(array):
    m,n = array.shape
    current_row = 0
    for current_column in range(n):
        entry = find_non_zero_entry_in_column(array,current_column,current_row)
        print "current_column : "+str(current_column)
        print "current_row : "+str(current_row)
        print "first_non_zero_entry : " + str(find_non_zero_entry_in_column(array,current_column,current_row)) 
        if  entry != -1:
            if entry != current_row:
                swap_rows(array,find_non_zero_entry_in_column(array,current_column,current_row),current_row)
                entry = current_row                
                print "swap"                
                print array
            if array[current_row][current_column] != 1:
                multiply_row(array,current_row,float(1)/array[current_row][current_column])
                print "divide"
                print array
            ensure_all_other_row_entries_zero(array,current_row,current_column)
            print "minus"            
            print array            
            current_row+=1
            
def identify(identifier,k):
    var = {0:"a",1:"b",2:"c",3:"d"}
    return var[int(identifier)/k]+str(int(identifier)%k)

def parse_ref(array):
    arraycopy = np.copy(array)
    m,n = array.shape    
    k = n/4
    leading_ones = []
    floating_vars = []
    for i in range(m):
        found_leading_one = False
        for j in range(n):
            arraycopy[i][j]=array[i][j]%2
            if array[i][j] == 1 and not found_leading_one:
                leading_ones.append(i)
                found_leading_one = True
                #arraycopy[i][j] = 0
    print arraycopy
    for j in range(n):
        if not j in leading_ones:
            floating_vars.append(j)
    print "There are "+str(len(leading_ones))+" leading_ones"+str([identify(i,k) for i in leading_ones])
    print "There are "+str(len(floating_vars))+" floating variables"+str([identify(i,k) for i in floating_vars])
    floating_vars_count = len(floating_vars)
    solution_set_matrix = np.zeros((n,floating_vars_count),dtype = int)
    print leading_ones
    print floating_vars
    for i in range(m):
        selected_var = -1
        for j in range(n):
            if arraycopy[i][j] == 1:
                if not j in floating_vars:
                    selected_var = j
                else:
                    solution_set_matrix[selected_var][floating_vars.index(j)] = 1
    for i in floating_vars:
        solution_set_matrix[i][floating_vars.index(i)] = 1
    print solution_set_matrix
    for a,b,c,d in itertools.product([0,1],repeat = 4):
        solution = np.transpose(np.dot(solution_set_matrix,np.array([[a],[b],[c],[d]])))
        solution_is_valid = True
        for i in range(k):
            if (solution[0][i]*solution[0][i+3*k] + solution[0][i+1*k]*solution[0][i+2*k])%2:
                print "Valid Solution Found"
                print a,b,c,d
                print solution
                return True
    return False



    
    count = [np.count_nonzero(array[i]) for i in range(m)]
    for i in range(m):
        if count[i] == 1:
            break            
    
    '''
    count = [np.count_nonzero(array[i]) for i in range(m)]
    for i in range(m):
        if count[i] == 1:
            print var[int(np.nonzero(array[i])[0]/k)]+str(int(np.nonzero(array[i])[0]%k))
            #print "Variable"+var_name+" = 0"
        elif count[i] != 0:
            output_string = ""
            for j in range(count[i]):
                identifier = np.nonzero(array[i])[0][j]
                output_string+=var[int(identifier)/k]+str(int(identifier)%k)
                if j != count[i]-1:
                    output_string+=" + "
                else:
                    output_string+=" = 0"
            print output_string
    '''
    
    
def resolve_ref(array):
    m,n = array.shape    
    k = n/4
    quadcnsts = np.ones(shape=(k,4))
    for i in range(m):
        for j in range(n):
            if array[i][j] == 1:
                quadcnsts[j%k][j/k] = 0
    for equation in quadcnsts:
        if equation[0]*equation[3] + equation[1]*equation[2] == 0:
            return False
    return True




        
matA = np.array([[1,2],[3,4]]).astype(float)
swap_rows(matA,0,1)
assert (matA == np.array([[3,4],[1,2]])).all()
multiply_row(matA,0,2)
assert (matA == np.array([[6,8],[1,2]])).all()
row_addition(matA,0,1,-6)
assert (matA == np.array([[0,-4],[1,2]])).all()
assert find_non_zero_entry_in_column(matA,0)==1
swap_rows(matA,0,1)
assert find_non_zero_entry_in_column(matA,0,1)==-1
matA = np.array([[6,8],[1,2]]).astype(float)
ensure_all_other_row_entries_zero(matA,1,0)
assert (matA == np.array([[0,-4],[1,2]])).all()
print "Test Success"

matA = np.array([[0,1,0],[1,0,1],[0,1,0]])
matB = np.array([[0,1,1],[1,0,1],[1,1,0]])

qA = Quantum.Quantum()
qB = Quantum.Quantum()
qA.AdMatrix(matA)
qB.AdMatrix(matA)
qA.local_complement(1)
print qA.adMatrix
print qB.adMatrix
matC = generate_linear_systems(qA.adMatrix,qB.adMatrix)
row_echelon_form(matC)
print parse_ref(matC)
'''
matD = np.array([[1,2,3,4],[2,4,6,8],[5,6,7,8],[4,5,7,12]],dtype = float)
row_echelon_form(matD)
parse_ref(matD)
'''