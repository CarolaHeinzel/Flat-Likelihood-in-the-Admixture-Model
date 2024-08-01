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
# Default (for only biallelic marker):
M = len(data_p)
J_m = [1]*M
poss = 1

# 2) Execute the functions
#%%
def umschreiben(array):
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
    elif(poss == "P2"):
        indices = get_indices_of_med_values_from_multiple_lists(umschreiben(q_vectors))
        conse = sum_selected_rows(A,K, indices)
    elif(poss == "P3"):
        indices = get_indices_of_max_values_from_multiple_lists(umschreiben(q_vectors))
        conse = sum_selected_rows(A,K, indices)
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
    # Initialisierung der Ergebnisliste
    result = []
    # Anzahl der Unterlisten (Anzahl der Zeilen in den Matrizen)
    num_rows = len(matrices[0])
    # Anzahl der Einträge in jeder Unterliste (Anzahl der Spalten in den Matrizen)
    num_cols = len(matrices[0][0])
    
    # Durch jede Zeile iterieren
    for row in range(num_rows):
        # Initialisierung der kombinierten Liste für die aktuelle Zeile
        combined_list = [max(matrix[row][col] for matrix in matrices) for col in range(num_cols)]
        # Den Index des maximalen Wertes in der kombinierten Liste finden und zur Ergebnisliste hinzufügen
        max_index = combined_list.index(max(combined_list))
        result.append(max_index)
    
    return result

matrix1 = [
    [0.402, 0.297, 0.3],
    [0.266, 0.352, 0.382]]

result = get_indices_of_med_values_from_multiple_lists(matrix1)
print("Die Indizes der maximalen Werte sind:", result)


def sum_selected_rows(matrix, K, indices):
    N = len(indices)
    result = np.zeros(len(matrix[0]))
    for i in range(N):
        row_index = 2*K * i + indices[i] * 2
        result += matrix[row_index]
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
    
    # Sort each individual permutation
    sorted_permuted_lists = [list(map(sorted, perm)) for perm in permuted_lists]
    
    return sorted_permuted_lists


def norm_similarity(q, q1):
    # Konvertieren der Listen in NumPy-Arrays für die Berechnung
    q = np.array(q)
    q1 = np.array(q1)
    temp = np.sum(np.sum((q - q1)**2, axis=0))
    #print(temp)
    #print(type(temp))
    sum_squared = sp.sqrt(temp)
    
    return sum_squared

def format_q(arrays):
    lists = [[arr.tolist() for arr in arrays]]
    return lists


def constraint4(x, q_alle, p_alle, K):
    #print("x", x)
    p_data, q_data = create_daten(q_alle, p_alle, K, x)
    q_hier = format_q(q_data)
    
    lists = q_hier[0]
    transposed_array = np.array(lists).T
    
    # Konvertieren des transponierten Arrays zurück in eine Liste von Listen
    q_data = transposed_array.tolist()
    
    q_permutiert = permute_and_sort_lists(q_alle)
    norm_comparisons = []
    #print(q_data)

    norm_target = norm_similarity(q_alle, q_data)
    for i in range(1,len(q_permutiert)):
        temp =  norm_similarity(q_permutiert[i], q_data)
        norm_comparisons.append(temp)
    constraints_similarity = min(norm_comparisons) - norm_target
    return constraints_similarity

#%%
def objective(x, conse):
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

#res_j = constraint_j(np.array([0,0,0,0,0,0]), p_alle, 3, [1]*20, 20)

#%%

def split_list(original_list, split_sizes):
    result = []
    current_index = 0

    for size in split_sizes:
        result.append(original_list[current_index:current_index + size])
        current_index += size

    return result

def split_all_lists(list_of_lists, split_sizes):
    return [split_list(sublist, split_sizes) for sublist in list_of_lists]

def sum_sublists(list_of_lists, split_sizes):
    split_result = split_all_lists(list_of_lists, split_sizes)
    list_of_lists = split_result
    return [[sum(sublist) for sublist in sublist_group] for sublist_group in list_of_lists]

# Gegebene Listen
list_of_lists = [[1, 2, 3, 4, 5, 6, 7, 8], [1, 2, 3, 4, 5, 6, 7, 8]]
split_sizes = [1, 2, 2, 3]

# Summieren der Elemente in den Unterlisten
sum_result = sum_sublists(list_of_lists, split_sizes)

# Ausgabe des Ergebnisses
print(sum_result)



#%%
def constraints_all(A, b_vek, K, p_alle, q_alle, J_m, M, simi):
    
    if(sum(J_m) > M & simi == 1):

        constr = [
            {'type': 'ineq', 'fun': lambda x: constraint1(x, A, b_vek)},
            {'type': 'ineq', 'fun': lambda x: constraint2(x, p_alle, K)},
            {'type': 'ineq', 'fun': lambda x: constraint3(x, p_alle, K)},
            {'type': 'ineq', 'fun': lambda x: constraint4(x, q_alle, p_alle, K)},
            {'type': 'ineq', 'fun': lambda x: constraint_j(x, p_alle, K, J_m, M)}]
            
    elif(sum(J_m) > M & simi == 1):

        constr = [
            {'type': 'ineq', 'fun': lambda x: constraint1(x, A, b_vek)},
            {'type': 'ineq', 'fun': lambda x: constraint2(x, p_alle, K)},
            {'type': 'ineq', 'fun': lambda x: constraint3(x, p_alle, K)},
            {'type': 'ineq', 'fun': lambda x: constraint4(x, q_alle, p_alle, K)}
            ]
    
    elif(sum(J_m) > M & simi == 0):

        constr = [
            {'type': 'ineq', 'fun': lambda x: constraint1(x, A, b_vek)},
            {'type': 'ineq', 'fun': lambda x: constraint2(x, p_alle, K)},
            {'type': 'ineq', 'fun': lambda x: constraint3(x, p_alle, K)},
            {'type': 'ineq', 'fun': lambda x: constraint_j(x, p_alle, K, J_m, M)}]
            
    elif(sum(J_m) > M & simi == 0):

        constr = [
            {'type': 'ineq', 'fun': lambda x: constraint1(x, A, b_vek)},
            {'type': 'ineq', 'fun': lambda x: constraint2(x, p_alle, K)},
            {'type': 'ineq', 'fun': lambda x: constraint4(x, q_alle, p_alle, K)}
            ]
                        
    return constr

def algorithm_max(cons, A, b_vek,p_alle,q_alle, K, J_m):
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
    constr = constraints_all(A, b_vek, K, p_alle, q_alle, J_m, M)
    x0 = [0]*len(x_bounds)
    result = minimize(objective, x0, args=(cons,), constraints=constr, bounds = x_bounds)
             
    return result.x

def algorithm_min(cons, A, b_vek,p_alle,q_alle, K, J_m):
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
    constr = constraints_all(A, b_vek, K, p_alle, q_alle, J_m, M)
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
def repeat_algo(q_vectors, p_alle, J_m, poss):
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
    q_tausch = tausche_erste_stelle(q_vectors)
    p_tausch = tausche_erste_stelle(p_alle)
    for i in range(K):
        q_vec = q_tausch[i]
        p_vec = p_tausch[i]
        A = create_A(q_vec)
        conse = create_cons(A, poss, q_vectors)
        b_vek = create_b(q_vec)
        # Minimization
        res = algorithm_min(conse, A, b_vek, p_vec, q_vec, K, J_m)
        res_a.append(res)
        # Maximization
        res = algorithm_max(conse, A, b_vek, p_vec, q_vec,K, J_m)
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
    print("N", N, param)
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
    #print(res_q)
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
        speichere_werte(q, "test_q1_K3_2", i)
        speichere_werte(p, "test_p1_K3_2", i)
    return

def algo_final(q_vectors, p_alle, K, J_m, poss):
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
    temp = repeat_algo(q_vectors, p_alle, J_m, poss)
    print("temp", temp)
    repeat_create_daten(q_alle, p_alle, K, temp)
    return

q_alle, p_alle = correct_format(data_q, data_p)
# Final Result
K = len(data_p)
print(algo_final(q_alle, p_alle, 3, J_m, poss))
