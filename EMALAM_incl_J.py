# To do:
# Implementiere den Algorithmus so, dass Label Switching ausgeschlossen wird.
# D.h. ändere die Grenzen dementsprechend!!!
import numpy as np
import pandas as pd
from sympy import symbols, Matrix
from scipy.optimize import minimize
import itertools
import sympy as sp
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
data_q = pd.read_csv(file_path_q, sep=' ', header=None)
data_p = pd.read_csv(file_path_p, sep=' ', header=None)
# Name the vector J_m
M = len(data_p)
J_m = [1]*M # Default (for only biallelic marker)
J_m[0] = 2
poss = "P3"
simi = 0 # Does not take C4 into account
# Alternatily, you can also define the function that should be maximized or
# minimized in the function create_conse by yourself.
# Names of the output file
names = ["test_q1_K3_4_J", "test_p1_K3_4_J"]
# 2) Execute the functions
#%%
def change_format(array):
    transponiert = np.transpose(array)
    umgeschrieben = transponiert.tolist()
    return umgeschrieben

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



def create_S(K, parameters):
    '''
    Creates the Matrix S_K in the paper

    Parameters
    ----------
    K : Int
        Number of ancestral populations.
    parameters : List
        Parmaters that should be optimized.

    Returns
    -------
    matrix : List
        S_K.

    '''
    matrix = np.zeros((K, K))
    remaining_params = parameters.tolist()
    parameters = parameters.tolist()
    for i in range(K):
        matrix[i, i] = 1 - sum(parameters[i * (K - 1): (i + 1) * (K - 1)])
    for i in range(K):
        for j in range(K):
            if j != i:
                matrix[i, j] = remaining_params.pop(0)
    return matrix

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
def create_cons(A, poss, q_vectors):
    if(poss == "P1"):
        conse = A[0,]
    else:
        conse = 0
    return conse

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

def get_indices_of_med_values_from_multiple_lists(*matrices):
    result = []
    K = len(matrices[0])
    num_rows = len(matrices[0])
    num_cols = len(matrices[0][0])
    for row in range(num_rows):
        combined_list = [min(np.abs(matrix[row][col] - 1/K) for matrix in matrices) for col in range(num_cols)]
        print(combined_list)
        max_index = combined_list.index(min(combined_list))
        result.append(max_index)
    
    return result

def get_indices_of_max_values_from_multiple_lists(*matrices):
    result = []
    num_rows = len(matrices[0])
    num_cols = len(matrices[0][0])
    for row in range(num_rows):
        combined_list = [max(matrix[row][col] for matrix in matrices) for col in range(num_cols)]
        max_index = combined_list.index(max(combined_list))
        result.append(max_index)
    return result

def permute_and_sort_lists(q):
    '''
    Example:
    q = [[1, 2, 3], [2, 3, 4], [1, 2, 5]]
    sorted_permutations = permute_and_sort_lists(q)

    [[1, 2, 3], [2, 3, 4], [1, 2, 5]]
    [[1, 2, 3], [1, 2, 5], [2, 3, 4]]
    [[2, 3, 4], [1, 2, 3], [1, 2, 5]]
    [[2, 3, 4], [1, 2, 5], [1, 2, 3]]
    [[1, 2, 5], [1, 2, 3], [2, 3, 4]]
    [[1, 2, 5], [2, 3, 4], [1, 2, 3]]

    Parameters
    ----------
    q : List
        DESCRIPTION.

    Returns
    -------
    sorted_permuted_lists : List
        DESCRIPTION.

    '''
    # Find all permutations of the list q
    permuted_lists = list(itertools.permutations(q))
    sorted_permuted_lists = [list(map(sorted, perm)) for perm in permuted_lists]
    return sorted_permuted_lists


def norm_similarity(q, q1):
    q = np.array(q)
    q1 = np.array(q1)
    temp = np.sum(np.sum((q - q1)**2, axis=0))
    sum_squared = sp.sqrt(temp)
    return sum_squared

def format_q(arrays):
    lists = [[arr.tolist() for arr in arrays]]
    return lists


def constraint4(x, q_alle, p_alle, K):
    p_data, q_data = create_daten(q_alle, p_alle, K, x)
    q_hier = format_q(q_data)
    lists = q_hier[0]
    transposed_array = np.array(lists).T
    q_data = transposed_array.tolist()
    q_permutiert = permute_and_sort_lists(q_alle)
    norm_comparisons = []
    norm_target = norm_similarity(q_alle, q_data)
    for i in range(1,len(q_permutiert)):
        temp =  norm_similarity(q_permutiert[i], q_data)
        norm_comparisons.append(temp)
    constraints_similarity = min(norm_comparisons) - norm_target
    return constraints_similarity

# Function that should be minimized with respect to x
def objective(x, conse):
    return - np.dot(conse, x)

# Function that should be minimized with respect to x
# Have an other minimization problem as above
def objective_min(x, conse):
    return np.dot(conse, x)

def entropy(x, q_alle):
    q_hier = change_format(q_alle)
    S = create_S(len(q_alle), np.array(x))
    N = len(q_alle[0])
    ent = 0
    for i in range(N):
        p = np.dot(q_hier[0], S)
        if(p.all() > 0):
            ent += -sum(p * np.log(p))
    return ent

def entropy_max(x, q_alle):
    q_hier = change_format(q_alle)
    S = create_S(len(q_alle), np.array(x))
    N = len(q_alle[0])
    ent = 0
    for i in range(N):
        p = np.dot(q_hier[0], S)
        if(p.all() > 0):
            ent += sum(p * np.log(p))
    return ent

# Define teh constraints C1 - C4
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

# Important for the function constraint_j
def split_list(original_list, split_sizes):
    result = []
    current_index = 0
    for size in split_sizes:
        result.append(original_list[current_index:current_index + size])
        current_index += size
    return result

# Important for the function constraint_j
def split_all_lists(list_of_lists, split_sizes):
    return [split_list(sublist, split_sizes) for sublist in list_of_lists]

# Important for the function constraint_j
def sum_sublists(list_of_lists, split_sizes):
    split_result = split_all_lists(list_of_lists, split_sizes)
    list_of_lists = split_result
    return [[sum(sublist) for sublist in sublist_group] for sublist_group in list_of_lists]


def constraint_j(x, p, K, J_m, M):
    M_inv = create_matrix_p(K, x)    
    C = [0]*sum(J_m)
    p_aufgeteilt = sum_sublists(p, J_m)
    for m in range(M):
        p_temp = []
        for k in range(K):
            p_temp.append(p_aufgeteilt[k][m])
        temp = np.dot(M_inv, p_temp)
        C.append(sum(temp))
    C = np.array(C)
    d = [1]*len(C)
    return d - C


# Constraints for the minimization problem
def constraints_all(A, b_vek, K, p_alle, q_alle, J_m, M, simi):
    if(sum(J_m) > M and simi == 1):

        constr = [
            {'type': 'ineq', 'fun': lambda x: constraint1(x, A, b_vek)},
            {'type': 'ineq', 'fun': lambda x: constraint2(x, p_alle, K)},
            {'type': 'ineq', 'fun': lambda x: constraint3(x, p_alle, K)},
            {'type': 'ineq', 'fun': lambda x: constraint4(x, q_alle, p_alle, K)},
            {'type': 'ineq', 'fun': lambda x: constraint_j(x, p_alle, K, J_m, M)}]
            
    elif(sum(J_m) == M and simi == 1):

        constr = [
            {'type': 'ineq', 'fun': lambda x: constraint1(x, A, b_vek)},
            {'type': 'ineq', 'fun': lambda x: constraint2(x, p_alle, K)},
            {'type': 'ineq', 'fun': lambda x: constraint3(x, p_alle, K)},
            {'type': 'ineq', 'fun': lambda x: constraint4(x, q_alle, p_alle, K)}
            ]
    
    elif(sum(J_m) > M and simi == 0):

        constr = [
            {'type': 'ineq', 'fun': lambda x: constraint1(x, A, b_vek)},
            {'type': 'ineq', 'fun': lambda x: constraint2(x, p_alle, K)},
            {'type': 'ineq', 'fun': lambda x: constraint3(x, p_alle, K)},
            {'type': 'ineq', 'fun': lambda x: constraint_j(x, p_alle, K, J_m, M)}]
            
    elif(sum(J_m) == M and simi == 0):
        constr = [
            {'type': 'ineq', 'fun': lambda x: constraint1(x, A, b_vek)},
            {'type': 'ineq', 'fun': lambda x: constraint2(x, p_alle, K)},
            {'type': 'ineq', 'fun': lambda x: constraint3(x, p_alle, K)},
            ]
                        
    return constr

def algorithm_max(cons, A, b_vek,p_alle,q_alle, K, J_m, simi, poss):
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
    M = len(p_alle[0])
    if(poss == "P1"):
        constr = constraints_all(A, b_vek, K, p_alle, q_alle, J_m, M, simi)
        x0 = [0]*len(x_bounds)
        result = minimize(objective, x0, args=(cons,), constraints=constr, bounds = x_bounds)
    elif(poss == "P2"):
        constr = constraints_all(A, b_vek, K, p_alle, q_alle, J_m, M, simi)
        x0 = [0]*len(x_bounds)
        result = minimize(entropy_max, x0, args=(q_alle,), constraints=constr, bounds = x_bounds)
    elif(poss == "P3"):
        constr = constraints_all(A, b_vek, K, p_alle, q_alle, J_m, M, simi)
        x0 = [0]*len(x_bounds)
        result = minimize(entropy, x0, args=(q_alle,), constraints=constr, bounds = x_bounds)
    
    return result.x

def algorithm_min(cons, A, b_vek,p_alle,q_alle, K, J_m, simi, poss):
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
    
    M = len(p_alle[0])
    constr = constraints_all(A, b_vek, K, p_alle, q_alle, J_m, M, simi)
    x0 = [0]*len(x_bounds)
    result = minimize(objective_min, x0, args=(cons,), constraints=constr, bounds = x_bounds)
    
    return result.x

# Finds the positions of every 0 in lst
def change_first_position(liste1):
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
def repeat_algo(q_vectors, p_alle, J_m, poss, simi):
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
    # K different opportunities for the maximization/minimization
    q_tausch = change_first_position(q_vectors)
    p_tausch = change_first_position(p_alle)
    for i in range(K):
        q_vec = q_tausch[i]
        p_vec = p_tausch[i]
        A = create_A(q_vec)
        conse = create_cons(A, poss, q_vectors)
        b_vek = create_b(q_vec)
        # Minimization
        if(poss == "P1"):
            res = algorithm_min(conse, A, b_vek, p_vec, q_vec, K, J_m, simi, poss)
            res_a.append(res)
        # Maximization
        res = algorithm_max(conse, A, b_vek, p_vec, q_vec,K, J_m,simi, poss)
        res_a.append(res)
        #print("res_a", res_a)
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
    matrix_inv = np.linalg.inv(matrix)    
    N = len(q_alle[0])
    M = len(p_alle[0])
    q_anders = np.array(rearrange_matrix(q_alle))
    p_anders = np.array(rearrange_matrix(p_alle))
    for i in range(N):
        q_temp = q_anders[i]@matrix
        res_q.append(q_temp)
    for m in range(M):
        p_temp = matrix_inv@p_anders[m]
        res_p.append(p_temp) 
    return res_p, res_q

# Transpose the matrix input_matrix
def transpose_matrix(input_matrix):
    return Matrix(input_matrix).transpose()

def save_values(temp, name, i):
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
    dateipfad = f"{name}_{i}.txt" 
    temp = transpose_matrix(temp)
    temp = temp.tolist()
    with open(dateipfad, 'w') as datei:
        for row in temp:
            zeilen_text = "\t".join(map(str, row)) + "\n"
            datei.write(zeilen_text)

def repeat_create_daten(q_alle, p_alle, K, parameters, poss, names):
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
    q_tausch = change_first_position(q_alle)
    p_tausch = change_first_position(p_alle)
    t = 0
    #for i in range(l):
    name_p = names[1]
    name_q = names[0]
    if(poss == "P1"):
        for i in range(l):
            if(i % 2 == 0):
                q_temp = q_tausch[t]
                p_temp = p_tausch[t]
                t += 1
            p,q = create_daten(q_temp, p_temp, K, parameters[i])
            p = transpose_matrix(p)
            q = transpose_matrix(q)
            save_values(q, name_q, i)
            save_values(p, name_p, i)
    else:
        q_temp = q_tausch[t]
        p_temp = p_tausch[t]
        p,q = create_daten(q_temp, p_temp, K, parameters[0])
        p = transpose_matrix(p)
        q = transpose_matrix(q)
        save_values(q, name_q, 0)
        save_values(p, name_p, 0)
    return

def algo_final(q_vectors, p_alle, K, J_m, poss, simi, names):
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
    K: Int
        number of populations
    J_m: List
        Number of alleles at marker m for all m = 1, \ldots, M
    poss: Int
        
    Returns
    -------
    None.

    '''
    temp = repeat_algo(q_vectors, p_alle, J_m, poss, simi)
    repeat_create_daten(q_alle, p_alle, K, temp, poss, names)
    return

q_alle, p_alle = correct_format(data_q, data_p)
# Final Result
K = len(data_p)
print(algo_final(q_alle, p_alle, 3, J_m, poss,simi, names))
