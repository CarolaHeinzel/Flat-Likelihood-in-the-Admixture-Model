import itertools
import numpy as np
import pandas as pd
from scipy.optimize import linprog
import math
from sympy import symbols, Matrix, Eq, solve
#%%
# Load Data
file_path = 'C:\\Users\\carol\\Downloads\\q_migtration7_mutation1'
data = pd.read_csv(file_path, sep=' ', header=None)
print(data)
#%%
# In this case: 5 populations
data = data.values.tolist()
liste = data
q1_vec = [float(subliste[0]) for subliste in liste]
q2_vec = [float(subliste[1]) for subliste in liste]
q3_vec = [float(subliste[2]) for subliste in liste]
#q4_vec = [float(subliste[3]) for subliste in liste]
#q5_vec = [float(subliste[4]) for subliste in liste]
#%%
# Output of STUCTURE
q_alle = [q1_vec, q2_vec, q3_vec] #, q4_vec, q5_vec]
#%%
# Daten für p einlesen
file_path_p = 'C:\\Users\\carol\\Downloads\\p_migtration7_mutation1'
data_p = pd.read_csv(file_path_p, sep=' ', header=None)
print(data_p)
data_p = data_p.values.tolist()
liste = data_p
p1 = [float(subliste[0]) for subliste in liste]
p2 = [float(subliste[1]) for subliste in liste]
p3 = [float(subliste[2]) for subliste in liste]
#p4 = [float(subliste[3]) for subliste in liste]
#p5 = [float(subliste[4]) for subliste in liste]

p_alle = [p1,p2,p3]
#%%
def list_to_column_matrix(lst):
    return Matrix([[x] for x in lst])
#%%
# Evaluates whether every value in vector falls between 0 and 1
def check_vector_bounds(vector):
    for i in vector:
        if(i < 0 or i > 1):
            return 1
    return 0
#%%
def create_b(q_vectors):
    '''
    
    Calculates the vector b in Ax \leq b.
    Parameters
    ----------
    q_vectors : List
        IA of every Individual for every population.

    Returns
    -------
    parameters : List
        Vector b.

    '''
    N = len(q_vectors[0])
    K = len(q_vectors)
    parameters = []
    for i in range(N):# Individuals
        for k in range(K):
            parameters.extend([1 - q_vectors[k][i], q_vectors[k][i]])
    if(K == 2):
        parameters.append(1)
    return parameters
print(create_b([[1,2,3],[4,5,6]]))
#%%
def extact_coeff(vector, K):
    sym = symbols('a0:%d' % (K * (K - 1) ))
    coefficients = [vector[0].coeff(symbol) for symbol in sym]
    return coefficients
#%%
# Calculates the matrix P_K
def create_matrix(K):
    symbols_list = symbols('a0:%d' % (K * (K - 1)))
    M = Matrix.zeros(K)
    symbol_index = 0
    for i in range(K):
        for j in range(K):
            if i != j:
                M[i, j] = symbols_list[symbol_index]
                symbol_index += 1
        M[i, i] = 1 - sum(M.row(i))
    return M
K = 3
matrix = create_matrix(K)
print(matrix)
#%%
def create_row(K,v, P_K, i): 
    prod = np.transpose(v) * P_K 
    temp = extact_coeff([prod[i]], K)
    return temp
v = Matrix([[5], [6], [7], [8]])
K = 4
P_K = create_matrix(4)
print(create_row(K,v,P_K,3))
#%%
def list_to_column_matrix(lst):
    return Matrix([[x] for x in lst])
#%%

def create_matrix_p(K, parameters):
    matrix = np.zeros((K, K))
    remaining_params = parameters
    for i in range(K):
        matrix[i, i] = 1 - sum(parameters[i * (K - 1): (i + 1) * (K - 1)])
    for i in range(K):
        for j in range(K):
            if j != i:
                matrix[i, j] = remaining_params.pop(0)
    return np.linalg.inv(matrix)

# Vektor with zeile_0 *x should be maximized or minimized
def create_cons(A):
    zeile_0 = A[0,]
    return zeile_0
#%%
def create_A(q_alle):
    '''
    Create Matrix A for Linear Optimization, i.e. for the Ax \leq b

    Parameters
    ----------
    q_alle : List
        All IA.

    Returns
    -------
    A : List
        Matrix for Optimization.

    '''
    N = len(q_alle[0])
    K = len(q_alle)
    if(K == 2):
        A = np.zeros((2*N*K+1, K*(K-1)))
        A[2*N*K,] = [1,1]
    else:
        A = np.zeros((2*N*K, K*(K-1)))
    ind = 0
    P_K = create_matrix(K)
    for i in range(N):
        lst = []
        for l in range(K):
            lst.append(q_alle[l][i])
        temp = list_to_column_matrix(lst)
        for k in range(K):
            zeile = create_row(K,temp, P_K, k)
            A[ind, ] = zeile
            A[ind+1, ] = np.multiply(zeile, -1)
            ind += 2
    return A
print(create_A([[1,2,3],[4,5,6]]))
#%%
def algorithm_max(cons, A, b_vek,p_alle, K):
    '''
    
    Calculates the optimal parameter for one case, e.g. find the parameters that
    maximize cons.
    
    Parameters
    ----------
    cons : List
        Vecotor that should be maximized.
    A : List
        Matrix for Ax \leq b.
    b_vek : List
        Vecotr for Ax \leq b.
    p_alle : List
        All allele frequencies in the STRUCTURE output. Format: 
        [[Allele Frequecies Population 1, Marker 1,...,M], ....,
         [Allele Frequecies Population 1, Marker 1,...,M]]
    K : Int
        Number of ancestral populations.

    Returns
    -------
    List
        Parameters that maximize cons and make sure that the constraint 
        P_K^{-1}p \in [0,1] is fulfilled.

    '''

    res_ok = 10
    value_range = [i * 0.01 for i in range(0,30)]
    #cons = replace_first_k_minus_1_with_zero(cons,K)
    zero = find_zero_positions(cons)
    rep = K*(K-1) - len(zero) -( K-1)
    combinations = list(itertools.product(value_range, repeat= rep))
    param = [0]*K*(K-1)
    combi_0 = add_zeros_to_tuples(combinations, zero,K)
    #print(zero, combinations, combi_0)
    for combi in combi_0:
        x_bounds = []
        for i in range(K*(K-1)):
            if(i < K-1):
                x_bounds.append((-2,0))
            else:
                x_bounds.append((-2,combi[i-K+1]))
        #print(x_bounds)
        result = linprog(cons, A_ub=A, b_ub=b_vek, bounds=x_bounds, method='highs')
        resi = []
        for i in result.x:
            resi.append(i)
        test = p_korrekt(p_alle, K, resi) # = 1: Fehler, = 0 alles ok!
        if(result.fun < res_ok and test == 0):
            res_ok = result.fun
            param = result.x            
    return param# , res_ok

#%%
def algorithm_min(cons, A, b_vek,p_alle, K):
    '''
    
    Calculates the optimal parameter for one case, e.g. find the parameters that
    maximize cons.
    
    Parameters
    ----------
    cons : List
        Vecotor that should be maximized.
    A : List
        Matrix for Ax \leq b.
    b_vek : List
        Vecotr for Ax \leq b.
    p_alle : List
        All allele frequencies in the STRUCTURE output. Format: 
        [[Allele Frequecies Population 1, Marker 1,...,M], ....,
         [Allele Frequecies Population 1, Marker 1,...,M]]
    K : Int
        Number of ancestral populations.

    Returns
    -------
    List
        Parameters that maximize cons and make sure that the constraint 
        P_K^{-1}p \in [0,1] is fulfilled.

    '''
    res_ok = 10
    value_range = [i * 0.01 for i in range(0,30)]
    rep =  K-1
    combinations = list(itertools.product(value_range, repeat= rep))
    param = [0]*K*(K-1)
    for combi in combinations:
        x_bounds = []
        for i in range(K*(K-1)):
            if(i < K-1):
                x_bounds.append((-2,combi[i]))
            else:
                x_bounds.append((-2,0))
        result = linprog(cons, A_ub=A, b_ub=b_vek, bounds=x_bounds, method='highs')
        resi = []
        for i in result.x:
            resi.append(i)
        test = p_korrekt(p_alle, K, resi) # = 1: Fehler, = 0 alles ok!
        if(result.fun < res_ok and test == 0):
            res_ok = result.fun
            param = result.x
    return param #, res_ok
  
#cons = [-0.358, -0.358, 0.164,0, 0.478, 0]
cons = [0.478, 0.478, -0.164,0, -0.358, 0]
b = create_b(q_alle)
#print(algorithm_max(cons, res_matrix, b, p_alle, 3))
#%%
def p_korrekt(p_alle, K, parameters):
    '''  
    Evaluate whether P_K^{-1}p falls within the range [0,1]
    Parameters
    ----------
    p_alle : List
        All estimated allele frequencies.
    K : Int
        Number of ancestral populations.
    parameters : List
        Values for the variables in P_K.

    Returns
    -------
    Int
        Either 0 or 1. 0 means that everything is ok and 1 means that an error 
        occured.

    '''
    M = len(p_alle[0])
    matrix = create_matrix_p(K,parameters)
    print(matrix)
    z = 0
    for i in range(M):
        temp = []
        for k in range(K):
            temp.append(p_alle[k][i])
        #print(temp)
        vector = matrix@temp
        #print(vector)
        test = check_vector_bounds(np.array(vector))
        if(test == 1):
            return 1
        else:
            z += 1
    if(z == M):
        return 0 # ALles ok
# para= [-0.02128378, -0.05764499,  0.8,         0.,          0.8,         0.        ]
para= [ 0.        ,  0.        ,  0.17718047, -0.16487617, -0.37136588,
       -0.10510951]
para = [ 0.        ,  0.        ,  0.17718047, -0.16487617, -0.37136588,
       -0.10510951]
 #[  1.02463399  ]
 #[    0.68543411]
 # matrix [[ 1.          0.          0.        ]
 #[-0.1390536   1.02463399  0.11441961]
 #[ 0.24162274  0.07294316  0.68543411]]
p_alle_anders = [p1,p2,p3]

print(p_korrekt(p_alle_anders, 3, para))
#%%
# Finds the positions of every 0 in lst
def find_zero_positions(lst):
    zero_positions = [i for i, elem in enumerate(lst) if elem == 0]
    return zero_positions

print(find_zero_positions([0,1,4,5,0]))
#%%
def add_zeros_to_tuples(tuples_list, lst,K):
    t = 0
    result = tuples_list.copy()
    for k in lst:
        k = k-(K-1)
        if(t == 0):
            result_1 = [(elem[:k] + (0,) + elem[k:]) for elem in tuples_list]
            t = 1
        else:
            result = [(elem[:k+1] + (0,) + elem[k+1:]) for elem in result_1]
    return result
lst = [3, 5]
value_range = [i * 0.1 for i in range(3)]
num_variables = 3
combinations = [(0.7000000000000001, 0.7000000000000001, 0, 0), (0.7000000000000001, 0.8, 0, 0), (0.7000000000000001, 0.9, 0, 0), (0.7000000000000001, 1.0, 0, 0), (0.7000000000000001, 1.1, 0, 0), (0.8, 0.7000000000000001, 0, 0), (0.8, 0.8, 0, 0), (0.8, 0.9, 0, 0), (0.8, 1.0, 0, 0), (0.8, 1.1, 0, 0), (0.9, 0.7000000000000001, 0, 0), (0.9, 0.8, 0, 0), (0.9, 0.9, 0, 0), (0.9, 1.0, 0, 0), (0.9, 1.1, 0, 0), (1.0, 0.7000000000000001, 0, 0), (1.0, 0.8, 0, 0), (1.0, 0.9, 0, 0), (1.0, 1.0, 0, 0), (1.0, 1.1, 0, 0), (1.1, 0.7000000000000001, 0, 0), (1.1, 0.8, 0, 0), (1.1, 0.9, 0, 0), (1.1, 1.0, 0, 0), (1.1, 1.1, 0, 0)]
[(-2, 0), (-2, 0), (-2, 0.7000000000000001), (-2, 0.7000000000000001), (-2, 0), (-2, 0)]
print(add_zeros_to_tuples(combinations, lst,3))
#%%
def replace_first_k_minus_1_with_zero(lst, K):
    lst[:K-1] = [0] * (K-1)
    return lst
# Beispielaufruf:
print(replace_first_k_minus_1_with_zero(cons,3))
#%%
##############################################################################################
# Save the values for the IA and the allele frequencies
# Therefore: Maximize and minimize every IA of the first individual
##############################################################################################
def tausche_erste_stelle(liste1):
    res = [liste1]
    liste_temp = liste1[0]
    k = len(liste1)  
    for i in range(1, k):
        liste = [elem.copy() for elem in liste1]  
        liste[0] = liste1[i]
        liste[i] = liste_temp
        res.append(liste)
    return res

liste1 = [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]]
ergebnis = tausche_erste_stelle(liste1)
#%%
# Try every K IA of individual that should be maximized and minimized
def repeat_algo(q_vectors, p_alle):
    '''
    Parameters
    ----------
    q_vectors : List
        Estimated IAs.
    p_alle : List
        Estimated allele frequencies.

    Returns
    -------
    res_a: List
        Parameters for minimization and maximization of cons.

    '''
    res_a = []
    K = len(q_vectors)
    # Gibt K Möglichkeiten, die IA aus welcher Population zu
    # maximieren/minimieren
    q_tausch = tausche_erste_stelle(q_vectors)
    p_tausch = tausche_erste_stelle(p_alle)
    for i in range(K):
        q_vec = q_tausch[i]
        p_vec = p_tausch[i]
        A = create_A(q_vec)
        conse = create_cons(A)
        b_vek = create_b(q_vec)
        # Minimieren
        res = algorithm_min(conse, A, b_vek, p_vec, K)
        res_a.append(res)
        # Maximieren
        conse = np.multiply(conse, -1)
        res = algorithm_max(conse, A, b_vek, p_vec,K)
        res_a.append(res)
        print("res_a", res_a)
    return res_a

#%%

def create_matrix_bestimmt(K, param):
    '''
    
    Calculates P_K for concrete values of the parameters.
    
    Parameters
    ----------
    K : Int
        Number of ancestral populations.
    param : List
        Parameters for the matrix P_K.

    Returns
    -------
    M : List
        P_K.

    '''
    M = Matrix.zeros(K)
    param_index = 0
    for i in range(K):
        for j in range(K):
            if i != j:
                M[i, j] = param[param_index]
                param_index += 1
        M[i, i] = 1 - sum(M.row(i))
    return M
K = 6
param = np.linspace(1,7,10)
param = np.array([-0.13956835, -0.01377282, -0.1192623 , -0.06928571,  0.03      ,
        0.        ,  0.        ,  0.        ,  0.03      ,  0.        ,
        0.        ,  0.        ,  0.03      ,  0.        ,  0.        ,
        0.        ,  0.03      ,  0.        ,  0.        ,  0.        ])

print(create_matrix_bestimmt(5, param))
#%%

def rearrange_matrix(input_matrix):
    '''
    Brings input_matrix to the correct form, i.e. from 
    [[a_1,...,a_N], [b_1,..., b_N],...]
    to 
    [[a_1,b_1,...], [a_2,b_2,...], ....]
    Here, a_1,...,a_N can be the allele frequencies of population 1 with markers 1,...
    ,M or the IA of Individual 1,...,N from population 1.

    Parameters
    ----------
    input_matrix : List
        Allele Frequencies or IA.

    Returns
    -------
    output_matrix : List
        Same elements like input_matrix, but with different order of the elements.

    '''
    rows, cols = len(input_matrix), len(input_matrix[0])
    output_matrix = Matrix.zeros(cols, rows)

    for i in range(rows):
        for j in range(cols):
            output_matrix[j, i] = input_matrix[i][j]

    return output_matrix

input_matrix = [[-2, 1, 2], [3, -6, 4]]
result_matrix = rearrange_matrix(input_matrix)
print(result_matrix)

#%%
# Jetzt: Wende jeweils die gefundenen Matrizen auf die dazugehörigen Daten an
def create_daten(q_alle, p_alle, K, param):
    '''
    
    Apply the optimal matric P_k to the estimated allele frequencies and IAs.
    
    Parameters
    ----------
    q_alle : List
        Estimated IAs.
    p_alle : List
        Estimated Allele Frequencies.
    K : Int
        Number of ancestral populations.
    param : List
        Optimal Paprameters for the matrix.

    Returns
    -------
    res_q: List
       IAs for every individual and for every population, 
        calculated with P_K. 
        1. Row: Individiual 1, Population 1,...,K
        N Rows
        
    res_p: List
        Allele Frequencies for every marker and for every population, 
        calculated with P_K^{-1}.
        1. Row: Marker 1, Population 1,...,K
        M Rows

    '''
    res_q = []
    res_p = []        
    matrix = create_matrix_bestimmt(K, param)
    matrix = np.array(matrix, dtype=np.float64)
    #print(matrix)
    matrix_inv = np.linalg.inv(matrix)    
    N = len(q_alle[0])
    M = len(p_alle[0])
    q_anders = np.array(rearrange_matrix(q_alle))
    #print("matrix",matrix_inv)
    p_anders = np.array(rearrange_matrix(p_alle))
    #print(p_anders[1])
    for i in range(N):
        #print(p_anders[i])
        q_temp = q_anders[i]@matrix
        #print(q_temp)
        res_q.append(q_temp)
    for m in range(M):
        #print(p_anders[m])
        p_temp = matrix_inv@p_anders[m]
        #print(p_temp)
        res_p.append(p_temp) 
    return res_p, res_q
K = 3
para= [ 0.        ,  0.        ,  0.17718047, -0.16487617, -0.37136588,
       -0.10510951]

#[-0.02128378, -0.05764499,  0.8,         0.,          0.8,         0.        ]

#p_alle = q_alle
#print(create_daten(q_alle, p_alle, K, para))
#%%
m = np.matrix([[1,0,0],[0.17718047, 1 - 0.17718047 + 0.16487617, -0.16487617],[-0.37136588,-0.10510951,1+0.37136588+0.10510951]])
m1 = np.linalg.inv(m)  
#%%
print(m1@[0.941, 0.987, 0.962])
#%%
# Transpos the matrix input_matrix
def transpose_matrix(input_matrix):
    return Matrix(input_matrix).transpose()
#%%
def speichere_werte(temp, name, i):
    '''
    
    Saves the output as txt.-file
    
    Parameters
    ----------
    temp : List
        List that should be saved, i.e. allele frequencies and IAs.
    name : String
        Name of the output file.
    i : Int
        Number of output file.

    Returns
    -------
    None.

    '''
    dateipfad = f"{name}_{i}.txt"  # Füge i zum Dateinamen hinzu
    temp = transpose_matrix(temp)
    temp = temp.tolist()
    # Öffnen der Datei im Schreibmodus
    with open(dateipfad, 'w') as datei:
        # Schreibe die Werte in die Datei
        for row in temp:
            zeilen_text = "\t".join(map(str, row)) + "\n"
            datei.write(zeilen_text)
            #%%
temp = Matrix([[0.899856501448725, 0.921055088348454],[1,2],[3,4]])
temp = transpose_matrix(temp)
print(speichere_werte(temp, "test", 1))
#%%
def repeat_create_daten(q_alle, p_alle, K, parameters):
    '''
    Calculates the maximal and minimal IA for every Individual at the first 
    value of the q_alle-list.    

    Parameters
    ----------
    q_alle : List
        Estimated IAs.
    p_alle : List
        Estimated allele frequencies.
    K : Int
        Number of ancestral populations.
    parameters : List
        Parameters of the matrix P_K to maximize or minimize the IA
        from one population.

    Returns
    -------
    None.

    '''
    l = len(parameters)
    q_tausch = tausche_erste_stelle(q_alle)
    p_tausch = tausche_erste_stelle(p_alle)
    t = 0
    #for i in range(l):
    for i in range(l):
        if(i % 2 == 0):
            q_temp = q_tausch[t]
            p_temp = p_tausch[t]
            t += 1
        p,q = create_daten(q_temp, p_temp, K, parameters[i])
        p = transpose_matrix(p)
        q = transpose_matrix(q)
        speichere_werte(q, "test_q1_K3", i)
        speichere_werte(p, "test_p1_K3", i)
    return
#%%
def algo_final(q_vectors, p_alle, K):
    '''
    Saves all values of q and p in .txt-file. For every calulated matrix
    P_K exisits one file that saves the corresponding qP_K and an other file
    with the P_K^{-1}p values.

    Parameters
    ----------
    q_vectors : List
        All IAs.
    p_alle : List
        All allele frequencies.

    Returns
    -------
    None.

    '''
    temp = repeat_algo(q_vectors, p_alle)
    print("temp", temp)
    repeat_create_daten(q_alle, p_alle, K, temp)
    return
#%%
print(algo_final(q_alle, p_alle, 3))
                            
#%%
print(0.688269544556121 *0.959050599823302+ 	0.112137432188065 * 0.987000000000000+	0.199593023255814 *0.941000000000000
      
)

print(0.110424479395030 * 0.987000000000000 + 	0.304022361785890  * 0.948381386610691+ 	0.958519105144749*0.585553158819079)

print(p_alle[0][0] * q_alle[0][0] + p_alle[1][0] *q_alle[1][0] + p_alle[2][0] * q_alle[2][0])
#%%
matrix = [[1,0,0],[0.17718047, 1-0.17718047+ 0.16487617 ,-0.16487617],[-0.37136588, -0.10510951 , 1+0.37136588+0.10510951]]
print(np.linalg.inv(matrix))  
print(np.linalg.inv(matrix)@[p_alle[0][2], p_alle[1][2], p_alle[2][2]])