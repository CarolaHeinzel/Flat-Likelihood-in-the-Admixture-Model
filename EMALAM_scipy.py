# To do:
# Implementiere den Algorithmus so, dass Label Switching ausgeschlossen wird.
# D.h. ändere die Grenzen dementsprechend!!!


import numpy as np
import pandas as pd
from sympy import symbols, Matrix
from scipy.optimize import minimize
#%%
# EMALAM
# The file 'p_values' has to contain K columns with M rows, i.e.
# at row m, column k stand the allele frequency for marker m in population k
# If we have J > 3, the rows contain p_{1,j,m}, ..., p_{K,j,m} for all j,m
# The file 'q_values' has to contain K columns with N rows, i.e.
# at row n, column k stand the allele frequency for individual n in population k
# See example files for both
# The correct format can be obtained by applying the Code 
# "Extract_q_p.R"
# to the Output of STRUCTURE
# Output: 
# Saves the values for q and p for the maximal and minimal IA
# in test_q1_K3_1.txt, ..., test_q1_K3_(2K).txt
# Application
# 1) Load Data, i.e. insert here the correct names for q and p
file_path_p = 'C:\\Users\\carol\\Downloads\\p_migtration7_mutation1'
file_path_q = 'C:\\Users\\carol\\Downloads\\q_migtration7_mutation1'
# 2) Execute the functions
#%%
data_q = pd.read_csv(file_path_q, sep=' ', header=None)
data_p = pd.read_csv(file_path_p, sep=' ', header=None)
#%%
# Correct format of q and p for later
def correct_format(data_q, data_p):
    data_q = data_q.values.tolist()
    data_p = data_p.values.tolist()
    K = len(data_q[0]) # Number of populations
    q_alle = []
    p_alle = []
    for i in range(K):
        q1_vec = [float(subliste[i]) for subliste in data_q]
        q_alle.append(q1_vec)
        liste = data_p
        p1 = [float(subliste[i]) for subliste in liste]
        p_alle.append(p1)
    return q_alle, p_alle

def list_to_column_matrix(lst):
    return Matrix([[x] for x in lst])

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

def extact_coeff(vector, K):
    sym = symbols('a0:%d' % (K * (K - 1) ))
    coefficients = [vector[0].coeff(symbol) for symbol in sym]
    return coefficients

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

def create_row(K,v, P_K, i): 
    prod = np.transpose(v) * P_K 
    temp = extact_coeff([prod[i]], K)
    return temp

def create_matrix_p(K, parameters):
    matrix = np.zeros((K, K))
    remaining_params = parameters.tolist()
    parameters = parameters.tolist()
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

def objective(x, conse):
    #conse = [0.1, 0.2, 0.4, 0.1, 0.2, 0.3]
    return - np.dot(conse, x)

def objective_min(x, conse):
    return np.dot(conse, x)

def constraint1(x, A, b):
    return b - np.dot(A, x)

def constraint2(x, p, K):
    
    M_inv = create_matrix_p(K, x)    
    C = np.dot(M_inv, p)
    C = C.flatten()
    d = [1]*len(C)
    return d - C

def constraint3(x, p, K):
    
    M_inv = create_matrix_p(K, x)    
    C = np.dot(M_inv, p)
    C = C.flatten()
    return C 

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

    x_bounds = []
    for i in range(K*(K-1)):
        x_bounds.append((-2,2))
    
    constr = [
    {'type': 'ineq', 'fun': lambda x: constraint1(x, A, b_vek)},
    {'type': 'ineq', 'fun': lambda x: constraint2(x, p_alle, K)},
    {'type': 'ineq', 'fun': lambda x: constraint3(x, p_alle, K)}
    ]
    x0 = [0]*len(x_bounds)
    result = minimize(objective, x0, args=(cons,), constraints=constr, bounds = x_bounds)
             
    return result.x

def algorithm_min(cons, A, b_vek,p_alle, K):
    '''
    
    Calculates the optimal parameter for one case, e.g. find the parameters that
    minimize cons.
    
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
    x_bounds = []
    for i in range(K*(K-1)):
        x_bounds.append((-2,2))
    
    constr = [
    {'type': 'ineq', 'fun': lambda x: constraint1(x, A, b_vek)},
    {'type': 'ineq', 'fun': lambda x: constraint2(x, p_alle, K)},
    {'type': 'ineq', 'fun': lambda x: constraint3(x, p_alle, K)}
    ]
    x0 = [0]*len(x_bounds)
    result = minimize(objective_min, x0, args=(cons,), constraints=constr, bounds = x_bounds)
    return result.x

# Finds the positions of every 0 in lst

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
    for i in range(1):
        q_vec = q_tausch[i]
        p_vec = p_tausch[i]
        A = create_A(q_vec)
        conse = create_cons(A)
        b_vek = create_b(q_vec)
        # Minimieren
        res = algorithm_min(conse, A, b_vek, p_vec, K)
        res_a.append(res)
        # Maximieren
        res = algorithm_max(conse, A, b_vek, p_vec,K)
        res_a.append(res)
        print("res_a", res_a)
    return res_a

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


# Apply Matrices to Data
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

# Transpose the matrix input_matrix
def transpose_matrix(input_matrix):
    return Matrix(input_matrix).transpose()

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

q_alle, p_alle = correct_format(data_q, data_p)
# Final Result
print(algo_final(q_alle, p_alle, 3))
