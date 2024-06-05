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
    [0.266, 0.352, 0.382],
    [0.365, 0.288, 0.347]
]
matrix2 = [
    [0.1, 0.5, 0.3],
    [0.4, 0.2, 0.35],
    [0.32, 0.4, 0.33]
]
result = get_indices_of_max_values_from_multiple_lists(matrix1, matrix2)
print("Die Indizes der maximalen Werte sind:", result)
#%%
import itertools
import numpy as np
import pandas as pd
from scipy.optimize import linprog
from sympy import symbols, Matrix
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
#file_path_p = 'p_values'
#file_path_q = 'q_values'
file_path_q = 'C:\\Users\\carol\\Downloads\\q_migtration7_mutation1'
file_path_p = 'C:\\Users\\carol\\Downloads\\p_migtration7_mutation1'

# 2) Execute the functions
#%%
data_q = pd.read_csv(file_path_q, sep=' ', header=None)
data_p = pd.read_csv(file_path_p, sep=' ', header=None)
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

# Evaluates whether every value in vector falls between 0 and 1
def check_vector_bounds(vector):
    for i in vector:
        if(i < 0 or i > 1):
            return 1
    return 0

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
    remaining_params = parameters
    for i in range(K):
        matrix[i, i] = 1 - sum(parameters[i * (K - 1): (i + 1) * (K - 1)])
    for i in range(K):
        for j in range(K):
            if j != i:
                matrix[i, j] = remaining_params.pop(0)
    return np.linalg.inv(matrix)

# Vektor with zeile_0 *x should be maximized or minimized
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
    value_range = [i * 0.01 for i in range(0,5)]
    zero = find_zero_positions(cons)
    rep = K*(K-1) - len(zero) -( K-1)
    combinations = list(itertools.product(value_range, repeat= rep))
    param = [0]*K*(K-1)
    combi_0 = add_zeros_to_tuples(combinations, zero,K)
    for combi in combi_0:
        #print("a")
        x_bounds = []
        for i in range(K*(K-1)):
            if(i < K-1):
                x_bounds.append((-2,0))
            else:
                x_bounds.append((-2,combi[i-K+1]))
        result = linprog(cons, A_ub=A, b_ub=b_vek, bounds=x_bounds, method='highs')
        resi = []
        for i in result.x:
            resi.append(i)
        test = p_korrekt(p_alle, K, resi) # = 1: Fehler, = 0 alles ok!
        if(result.fun < res_ok and test == 0):
            res_ok = result.fun
            param = result.x            
    return param# , res_ok

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
    #print(matrix)
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
# Finds the positions of every 0 in lst
def find_zero_positions(lst):
    zero_positions = [i for i, elem in enumerate(lst) if elem == 0]
    return zero_positions
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
def replace_first_k_minus_1_with_zero(lst, K):
    lst[:K-1] = [0] * (K-1)
    return lst

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

def sum_selected_rows(matrix, K, indices):
    N = len(indices)
    result = np.zeros(len(matrix[0]))
    for i in range(N):
        row_index = 2*K * i + indices[i] * 2
        result += matrix[row_index]
    return result

# Try every K IA of individual that should be maximized and minimized
def repeat_algo(q_vectors, p_alle, indices):
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
    A = create_A(q_vectors)
    print("A",A, indices)
    conse = sum_selected_rows(A,K, indices)
    b_vek = create_b(q_vectors)
    conse = np.multiply(conse, -1)
    res = algorithm_max(conse, A, b_vek, p_alle,K)
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

def speichere_werte(temp, name):
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
    dateipfad = f"{name}.txt"  # Füge i zum Dateinamen hinzu
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

    p,q = create_daten(q_alle, p_alle, K, parameters[0])
    p = transpose_matrix(p)
    q = transpose_matrix(q)
    speichere_werte(q, "population_q1_30_K3")
    speichere_werte(p, "population_p1_30_K3")
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
    indices = get_indices_of_max_values_from_multiple_lists(umschreiben(q_vectors))
    print(indices)
    temp = repeat_algo(q_vectors, p_alle, indices)
    print(temp)
    repeat_create_daten(q_alle, p_alle, K, temp)
    return

q_alle, p_alle = correct_format(data_q, data_p)
# Final Result
print(algo_final(q_alle, p_alle, 3))

