import pandas as pd
from scipy.optimize import linprog
import numpy as np
import math
from sympy import symbols, Matrix, Eq, solve
import itertools
#%%
# Load Data
file_path = 'r1(1)'
data = pd.read_csv(file_path, sep='\t', header=None)
print(data)
#%%
# In this case: 5 populations
data = data.values.tolist()
liste = data
q1_vec = [float(subliste[0]) for subliste in liste]
q2_vec = [float(subliste[1]) for subliste in liste]
q3_vec = [float(subliste[2]) for subliste in liste]
q4_vec = [float(subliste[3]) for subliste in liste]
q5_vec = [float(subliste[4]) for subliste in liste]
#%%
# Output of STUCTURE
q_alle = [q1_vec, q2_vec, q3_vec, q4_vec, q5_vec]
#%%
# Daten für p einlesen
file_path_p = 'q_551.txt'
data_p = pd.read_csv(file_path_p, sep='\t', header=None)
print(data_p)
data_p = data_p.values.tolist()
liste = data_p
p1 = [float(subliste[0]) for subliste in liste]
p2 = [float(subliste[1]) for subliste in liste]
p3 = [float(subliste[2]) for subliste in liste]
p4 = [float(subliste[3]) for subliste in liste]
p5 = [float(subliste[4]) for subliste in liste]

p_alle = [p1,p2,p3,p4,p5]

#%%
def check_vector_bounds(vector):
    for i in vector:
        if(i < 0 or i > 1):
            return 1
    return 0
#%%
# Koeffizienten extrahieren
def extact_coeff(vector, K):
    sym = symbols('a0:%d' % (K * (K - 1) ))
    coefficients = [vector[0].coeff(symbol) for symbol in sym]
    return coefficients
#%%
def create_matrix(K):
    # Symbole definieren
    symbols_list = symbols('a0:%d' % (K * (K - 1)))
    # Matrix initialisieren
    M = Matrix.zeros(K)
    # Matrixelemente setzen
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
# Matrix ausgeben
print(matrix)
#%%
# Jetzt: 
def create_row(K,v, P_K, i): 
    prod = np.transpose(v) * P_K 
    #print(prod[0])
    temp = extact_coeff([prod[i]], K)
    return temp
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

def create_cons(A):
    zeile_0 = A[0,]
    return zeile_0

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
#%%

cons = [0.056, 0.062, 0.014, 0.747, 0.122]
A = create_A(q_alle)
#%%
conse = create_cons(A)
#%%
b = create_b(q_alle)
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
    value_range = [i * 0.01 for i in range(0,10)]
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
        print(x_bounds)
        result = linprog(cons, A_ub=A, b_ub=b_vek, bounds=x_bounds, method='highs')
        resi = []
        for i in result.x:
            resi.append(i)
        test = p_korrekt(p_alle, K, resi) # = 1: Fehler, = 0 alles ok!
        if(result.fun < res_ok and test == 0):
            res_ok = result.fun
            param = result.x
    return param
#%%
def find_zero_positions(lst):
    zero_positions = [i for i, elem in enumerate(lst) if elem == 0]
    return zero_positions
#%%
def add_zeros_to_tuples(tuples_list, lst,K):
    result = tuples_list.copy()
    for k in lst:           
        result = [(elem[:k] + (0,) + elem[k:]) for elem in result]
    return result
lst = [1, 2, 4]

combinations = [(0.7000000000000001, 0.7000000000000001, 0.1, 0.3),(0.7000000000000001, 0.9, 0, 0)]
print(add_zeros_to_tuples(combinations, lst,3))
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
    value_range = [i * 0.01 for i in range(0,10)]
    #cons = replace_first_k_minus_1_with_zero(cons,K)
    zero = find_zero_positions(cons)
    zero_final  = []
    for i in zero:
        zero_final.append(i-K+1)
    rep = K-1
    combinations = list(itertools.product(value_range, repeat= rep))
    param = [0]*K*(K-1)
    combi_0 = add_zeros_to_tuples(combinations, zero_final,K)
    #print(combi_0)
    for combi in combi_0:
        x_bounds = []
        for i in range(K*(K-1)):
            if(i < K-1):
                x_bounds.append((-3,0))
            else:
                
                x_bounds.append((-3,combi[i-K+1]))
        print(x_bounds)
        result = linprog(cons, A_ub=A, b_ub=b_vek, bounds=x_bounds, method='highs')
        resi = []
        for i in result.x:
            resi.append(i)
        test = p_korrekt(p_alle, K, resi) # = 1: Fehler, = 0 alles ok!
        if(result.fun < res_ok and test == 0):
            res_ok = result.fun
            param = result.x
            
    return param
#%%
A = create_A(q_alle)
b = create_b(q_alle)
A_temp = A[1:5,:]
b_temp = b[1:5]
#%%
conse = create_cons(A)
#%%
print(algorithm_max(conse, A_temp, b_temp,p_alle, 5))
#%%
#%%
      
res_a =  [np.array([ 0.02      ,  0.        ,  0.01      ,  0.02      ,  0.        ,
       -0.00458333,  0.        , -0.34766063, -0.00532807, -0.00855587,
       -0.00643338,  0.        , -0.03743651, -0.08536053, -0.00627303,
       -0.08020191, -0.04083378, -0.10573175, -0.01018605, -0.14993084]), 
          np.array([-0.05157953, -0.00755784, -0.03683544, -0.05058482,  0.03      ,
        0.        ,  0.        , -0.0182356 ,  0.02      ,  0.        ,
        0.        ,  0.        ,  0.03      ,  0.        ,  0.        ,
        0.        ,  0.03      , -0.06675806,  0.        ,  0.        ]), 
        np.array([ 0.01      ,  0.        ,  0.01      ,  0.01      , -0.00635561,
       -0.00614683, -0.03533828, -0.06096499, -0.00088745, -0.00627369,
       -0.0061539 ,  0.        , -0.06955843, -0.04534856, -0.00680408,
       -0.10323866, -0.06402884, -0.04859521, -0.00941899, -0.1733958 ]),np.array([-0.04990028, -0.0094742 , -0.04273079, -0.26263538,  0.03      ,
       -0.00588511, -0.03234313,  0.        ,  0.03      , -0.00562289,
       -0.00573596,  0.        ,  0.03      , -0.03722799, -0.00610597,
        0.        ,  0.03      ,  0.        ,  0.        ,  0.        ]), np.array([ 0.02      ,  0.        ,  0.01      ,  0.02      ,  0.        ,
       -0.03944308,  0.        , -0.34892121, -0.00576276, -0.05968404,
       -0.03434477, -0.01162443, -0.00637894, -0.01715206, -0.03734696,
       -0.10292416, -0.00883046, -0.17322417, -0.04426264, -0.13845081]), np.array([-0.00827178, -0.00597431, -0.00579359, -0.00665436,  0.        ,
       -0.03917229,  0.        , -0.34857841,  0.        , -0.05818656,
       -0.03316156, -0.01226031,  0.        , -0.0177404 , -0.03701055,
       -0.07524977,  0.        , -0.17171088, -0.04377176, -0.13401152]), np.array([ 0.02      ,  0.        ,  0.02      ,  0.02      ,  0.        ,
       -0.00342324, -0.03197369, -0.37900364, -0.00481944, -0.0083733 ,
       -0.00621641, -0.00743094, -0.02758571, -0.05901481, -0.00652134,
       -0.01721691, -0.11147854, -0.21587414, -0.01130828, -0.07191942]), np.array([-0.09698633, -0.00815126, -0.04543811, -0.06012184,  0.02      ,
        0.        ,  0.        , -0.25717635,  0.03      ,  0.        ,
        0.        ,  0.        ,  0.01      ,  0.        ,  0.        ,
        0.        ,  0.03      , -0.01631191,  0.        ,  0.        ]), np.array([ 0.02      ,  0.        ,  0.02      ,  0.02      , -0.24218013,
       -0.00986878, -0.05792506, -0.06021517,  0.        , -0.00924689,
       -0.00691093, -0.00607229, -0.05706705, -0.10022131, -0.00636027,
       -0.04136578, -0.00840568, -0.02058811, -0.00613021, -0.03889469]), np.array([-0.13956835, -0.01377282, -0.1192623 , -0.06928571,  0.03      ,
        0.        ,  0.        ,  0.        ,  0.03      ,  0.        ,
        0.        ,  0.        ,  0.03      ,  0.        ,  0.        ,
        0.        ,  0.03      ,  0.        ,  0.        ,  0.        ])]
            #%%
                                                                                    
repeat_create_daten(q_alle, p_alle, 5, res_a)                                                                                   
                                                                                 