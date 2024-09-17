import numpy as np
import pandas as pd
from sympy import symbols, Matrix
from scipy.optimize import minimize
import itertools
import sympy as sp
import argparse
#%%
# Application
# 1) Load Data
#file_path_p = '/home/ch1158/Downloads/p_CEU_IBS_TSI_K2'
#file_path_q = '/home/ch1158/Downloads/q_CEU_IBS_TSI_K2.txt'
#file_path_pJ = '/home/ch1158/Downloads/p_CEU_IBS_TSI_K2_J'

#data_q = pd.read_csv(file_path_q, sep=' ', header=None)


def correct_formatJ(data_p):
    data_p = data_p.values.tolist()
    K = len(data_p[0]) # Number of populations
    p_alle = []
    for i in range(K):
        liste = data_p
        p1 = [float(subliste[i]) for subliste in liste]
        p_alle.append(p1)
    return p_alle

# poss names what we want to maximize/minimize
# Alternativly, you can also define the function that should be
# minimized in the function entropy() by yourself.
#poss = "P2"
simi = 0 # Do not take label switching into account (recommended)
k_specific = 0
n_trials = 10 # Number of different initial values for the minimization function
# Names of the output files

def main(data_q_input, data_p_input, pJ_input, data_q_output, data_p_output, poss):
    # Lese den Inhalt der Input-Dateien ein
    data_q = pd.read_csv(data_q_input, sep=' ', header=None)
    data_p = pd.read_csv(data_p_input, sep=' ', header=None)
    
    if pJ_input == "-1":
        pJ = None
    else:
        pJ = pd.read_csv(pJ_input, sep=' ', header=None)
    
    # Bereinigen von data_p: Entferne Zeilen, die nur 1 enthalten
    data_p = data_p[~(data_p == 1).all(axis=1)]
    
    # Speichere die Ergebnisse in den angegebenen Ausgabedateien
    data_q.to_csv(data_q_output, sep=' ', header=False, index=False)
    data_p.to_csv(data_p_output, sep=' ', header=False, index=False)
    
    # Erstelle die Liste der Dateinamen
    names = [data_q_output, data_p_output]
    
    return data_q, data_p, pJ, names, poss

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process input and output file names.")
    parser.add_argument("data_q_input", type=str, help="The input file name for data_q")
    parser.add_argument("data_p_input", type=str, help="The input file name for data_p")
    parser.add_argument("pJ_input", type=str, help="The input file name for pJ, or '-1' if not provided")
    parser.add_argument("data_q_output", type=str, help="The output file name for data_q")
    parser.add_argument("data_p_output", type=str, help="The output file name for data_p")
    parser.add_argument("poss", type=str, help="The output file name for data_p")

    parser.add_argument("--k_specific", type=str, help="Specific input for P4 or P5, required if poss is P4 or P5")

    args = parser.parse_args()

    # Check if k_specific is required
    if args.poss in ["P4", "P5"] and not args.k_specific:
        raise ValueError("k_specific is required when poss is P4 or P5")
    # Rufe die main-Funktion auf und übergebe die Input-Dateien
    data_q, data_p, pJ, names, poss = main(args.data_q_input, args.data_p_input, args.pJ_input, args.data_q_output, args.data_p_output, args.poss)
    
    # Optionale Verarbeitung für pJ
    if pJ is not None:
        pJ = correct_formatJ(pJ)
    data_p = data_p[~(data_p == 1).all(axis=1)]

    print("Input:", data_q, data_p, pJ)

#%%
# 2) Execute the functions
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

# Important to calculate the matrix A
def list_to_column_matrix(lst):
    return Matrix([[x] for x in lst])

# Important to calculate the matrix A
def extact_coeff(vector, K):
    sym = symbols('a0:%d' % (K * (K - 1) ))
    coefficients = [vector[0].coeff(symbol) for symbol in sym]
    return coefficients

# Calculates the matrix S_K for the symbols a
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

# Important to calculate the matrix A
def create_row(K,v, P_K, i): 
    prod = np.transpose(v) * P_K 
    temp = extact_coeff([prod[i]], K)
    return temp

# Matrix S_K for fixed parameter
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

# Invert matrix S_K for fixed values
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
# Only for P1 relevant
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

# Function that should be minimized with respect to x
def objective(x, conse):
    return - np.dot(conse, x)

# Function that should be minimized with respect to x
# Have an other minimization problem as above
def objective_min(x, conse):
    return np.dot(conse, x)

# Function that should be minimized with respect to x
def entropy(x, q_alle):
    q_hier = change_format(q_alle)
    S = create_S(len(q_alle), np.array(x))
    N = len(q_alle[0])
    ent = 0
    for i in range(N):
        p = np.dot(q_hier[i], S)
        b = all([x > 0.001 for x in p])
        if(b):
            ent += sum(p * np.log(p))
        else:
            ent -= 10**3
    return ent

# Function that should be minimized with respect to x
def entropy_pop_max(x, q_alle, k_specific):
    q_hier = change_format(q_alle)
    S = create_S(len(q_alle), np.array(x))
    N = len(q_alle[0])
    ent = 0
    for i in range(N):
        p = np.dot(q_hier[i], S)
        p = p[k_specific]
        if(p > 0):
            ent += - p 
    return ent

# Function that should be minimized with respect to x
def entropy_pop_min(x, q_alle, k_specific):
    q_hier = change_format(q_alle)
    S = create_S(len(q_alle), np.array(x))
    N = len(q_alle[0])
    ent = 0
    for i in range(N):
        p = np.dot(q_hier[i], S)
        p = p[k_specific]
        if(p > 0):
            ent += p 
    return ent

# Function that should be minimized with respect to x
def entropy_max(x, q_alle):
    q_hier = change_format(q_alle)
    
    S = create_S(len(q_alle), np.array(x))
    N = len(q_alle[0])
    ent = 0
    for i in range(N):
        p = np.dot(q_hier[i], S)
        b = all([x > 0.001 for x in p])
        if(b):
            ent -= sum(p * np.log(p))
        else:
            ent += 10**3
    return ent

# Define the constraints C1 - C4
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


def constraint_j(x, pJ, K):
    M_inv = create_matrix_p(K, x)    
    C = np.dot(M_inv, pJ)
    C = C.flatten()
    d = [1]*len(C)
    return d - C

# Define the constraint to avoid label switching
# Only use this for a small number of individuals and markers
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

# Constraints for the minimization problem
def constraints_all(A, b_vek, K, p_alle, q_alle, simi, pJ):
    if(pJ!= 0 and simi == 1):
        constr = [
            {'type': 'ineq', 'fun': lambda x: constraint1(x, A, b_vek)},
            {'type': 'ineq', 'fun': lambda x: constraint2(x, p_alle, K)},
            {'type': 'ineq', 'fun': lambda x: constraint3(x, p_alle, K)},
            {'type': 'ineq', 'fun': lambda x: constraint4(x, q_alle, p_alle, K)},
            {'type': 'ineq', 'fun': lambda x: constraint_j(x, pJ, K)}]
            
    elif(pJ == 0 and simi == 1):

        constr = [
            {'type': 'ineq', 'fun': lambda x: constraint1(x, A, b_vek)},
            {'type': 'ineq', 'fun': lambda x: constraint2(x, p_alle, K)},
            {'type': 'ineq', 'fun': lambda x: constraint3(x, p_alle, K)},
            {'type': 'ineq', 'fun': lambda x: constraint4(x, q_alle, p_alle, K)}]
            
    elif(pJ!= 0 and simi == 0):
        constr = [
            {'type': 'ineq', 'fun': lambda x: constraint1(x, A, b_vek)},
            {'type': 'ineq', 'fun': lambda x: constraint2(x, p_alle, K)},
            {'type': 'ineq', 'fun': lambda x: constraint3(x, p_alle, K)},
            {'type': 'ineq', 'fun': lambda x: constraint_j(x, pJ, K)}]
            
            
    elif(pJ == 0 and simi == 0):
        constr = [
            {'type': 'ineq', 'fun': lambda x: constraint1(x, A, b_vek)},
            {'type': 'ineq', 'fun': lambda x: constraint2(x, p_alle, K)},
            {'type': 'ineq', 'fun': lambda x: constraint3(x, p_alle, K)},
            ]
                        
    return constr


def initial(n_trials, cons, A, b_vek,p_alle,q_alle, K, simi, poss, k_specific, pJ, target_function):
    x_bounds = []
    for i in range(K*(K-1)):
        x_bounds.append((-2,2))
    constr = constraints_all(A, b_vek, K, p_alle, q_alle, simi, pJ)

    best = None
    for _ in range(n_trials):
        x0 = np.random.uniform(low=[bound[0] for bound in x_bounds], high=[bound[1] for bound in x_bounds])
        if(poss == "P1"):
            result = minimize(target_function, x0, args=(cons,), constraints=constr, bounds = x_bounds)
        elif(poss == "P2" or poss == "P3"):
            result = minimize(target_function, x0, args=(q_alle,), constraints=constr, bounds = x_bounds)
        else:
            result = minimize(target_function, x0, args=(q_alle, k_specific, ), constraints=constr, bounds = x_bounds)

        if result.success and (best is None or best.fun > result.fun):
            best = result
    return best

def algorithm_max(cons, A, b_vek,p_alle,q_alle, K, simi, poss, k_specific, pJ, n_trials):
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
    simi: Int
        0: Do not take label switching into account
        1: Take label switching into account
    poss: String
        Choice of the Target function
    k_specific: Names the population that should be considered (if any 
                specific population should be considered, i.e. P4 or P5).
    pJ: List
        Used to ensure that the sum over J equals 1.
    n_trails: Int
        Number of different initial values that are tried.    

    Returns
    -------
    List
        Parameters that maximize cons and make sure that the constraint 
        P_K^{-1}p \in [0,1] is fulfilled.

    '''

    if(poss == "P1"):
        result = initial(n_trials, cons, A, b_vek,p_alle,q_alle, K, simi, poss, k_specific, pJ, objective)
    elif(poss == "P2"):
        result = initial(n_trials, cons, A, b_vek,p_alle,q_alle, K, simi, poss, k_specific, pJ, entropy_max)
    elif(poss == "P3"):
        result = initial(n_trials, cons, A, b_vek,p_alle,q_alle, K, simi, poss, k_specific, pJ, entropy)
    elif(poss == "P4"):
        result = initial(n_trials, cons, A, b_vek,p_alle,q_alle, K, simi, poss, k_specific, pJ, entropy_pop_min)
    elif(poss == "P5"):
        result = initial(n_trials, cons, A, b_vek,p_alle,q_alle, K, simi, poss, k_specific, pJ, entropy_pop_max)
    return result.x

def algorithm_min(cons, A, b_vek,p_alle,q_alle, K, simi, poss, k_specific, pJ, n_trials):
    '''
    
    Calculates the optimal parameter for one case, e.g. find the parameters that
    minimize cons. This function is only used for poss == P1.
    
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
    simi: Int
        0: Do not take label switching into account
        1: Take label switching into account
    poss: String
        Choice of the Target function
    k_specific: Names the population that should be considered (if any 
                specific population should be considered, i.e. P4 or P5).
    pJ: List
        Used to ensure that the sum over J equals 1.
    n_trails: Int
        Number of different initial values that are tried.   

    Returns
    -------
    List
        Parameters that maximize cons and make sure that the constraint 
        P_K^{-1}p \in [0,1] is fulfilled.

    '''
    result = initial(n_trials, cons, A, b_vek,p_alle,q_alle, K, simi, poss, k_specific, pJ, objective_min)

    print(result.x)
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
def repeat_algo(q_vectors, p_alle, poss, simi, k_specific,pJ, n_trials):
    '''
    Parameters
    ----------
    q_vectors : List
        Estimated IAs.
    p_alle : List
        All allele frequencies in the STRUCTURE output. Format: 
        [[Allele Frequecies Population 1, Marker 1,...,M], ....,
         [Allele Frequecies Population 1, Marker 1,...,M]]
    poss: String
        Choice of the Target function
    simi: Int
        0: Do not take label switching into account
        1: Take label switching into account
    k_specific: Names the population that should be considered (if any 
                specific population should be considered, i.e. P4 or P5).
    pJ: List
        Used to ensure that the sum over J equals 1.
    n_trails: Int
        Number of different initial values that are tried.   

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
        # Minimization for P1
        if(poss == "P1"):
            res = algorithm_min(conse, A, b_vek,p_vec,q_vec, K, simi, poss, k_specific, pJ, n_trials)
            res_a.append(res)
        # Minimization for the other cases
        res = algorithm_max(conse, A, b_vek,p_vec,q_vec, K, simi, poss, k_specific, pJ, n_trials)
        res_a.append(res)
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
    poss: String
        Choice of the Target function
    names: List
        Names of the output lists
    Returns
    -------
    None.

    '''
    l = len(parameters)
    q_tausch = change_first_position(q_alle)
    p_tausch = change_first_position(p_alle)
    t = 0
    name_p = names[1]
    name_q = names[0]
    
    if(poss == "P1"):
        # Consider every population, i.e. 2*K possibilites
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

def algo_final(q_vectors, p_alle, K, poss, simi, names, k_specific,pJ, n_trials):
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
    poss: String
        Choice of the Target function
    names: List
        Names of the Output Files
    simi: Int
        0: Do not take label switching into account
        1: Take label switching into account
    k_specific: Names the population that should be considered (if any 
                specific population should be considered, i.e. P4 or P5).
    pJ: List
        Is used to make sure that the sum over J for every marker and every population
        equals 1.
    n_trials: Int
        Number of different initial values that the user should try
        
    Returns
    -------
    None.

    '''
    temp = repeat_algo(q_vectors, p_alle, poss, simi, k_specific,pJ, n_trials)
    repeat_create_daten(q_alle, p_alle, K, temp, poss, names)
    print("Successfully calculated the MLEs")
    return 

q_alle, p_alle = correct_format(data_q, data_p)
# Final Result
K = data_p.shape[1]
res_algo_final = algo_final(q_alle, p_alle, K, poss,simi, names, k_specific, pJ, n_trials)
