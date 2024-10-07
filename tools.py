import numpy as np
import pandas as pd
from scipy.stats import entropy
from scipy.optimize import minimize, LinearConstraint, NonlinearConstraint
import re


#################################################
## Importing data from a structure output file ## 
#################################################

def load_structure_file(path):
    with open(path, 'rb') as f:
        return f.readlines()

# Determine K from the structure File    
def get_K(structure_file):
    output_text = ''.join([line.decode('utf-8') for line in structure_file])
    match = re.search(r"(\d+)\s+populations assumed", output_text)
    if match:
        return int(match.group(1))  
    else:
        raise ValueError(f"K not found in {structure_file}")

# In a line of the form
# 0   (0.826) 0.542 0.999
# return [0.542, 0.999] : np.array
def get_data(line, K):
    # print(line)
    l = line.strip().split(" ")
    #print(l)
    return [float(x) for x in l[-K:]]
    
# get p from structure file (as lines)
# returns a dict key: value, where key is an identifier of a marker, and value is a K x J_m array
# with the allele frequencies in all populations at all allelels.
def get_p(lines):
    K = get_K(lines)
    lines_decoded = [line.decode('utf-8') for line in lines]  # Assuming 'utf-8' encoding
    output_text = ''.join([line.decode('utf-8') for line in lines])
    
    res = {}
    # The information for the next marker starts like this:
    #pattern = r'Locus\s+\d+\s*:\s*(\S+)'
    #pattern = r'Locus\s+\S+'
    for i, line in enumerate(lines_decoded):
        #match = re.search(pattern, line)
        if "Locus" in line: 
            l = line.strip().split(" ")
            marker = l[-1] if l[-1] != ":" else l[-2]
            j = next((k for k, line in enumerate(output_text.splitlines()[i:], start=i+1) if " missing data" in line), None)
            k = next((l for l, line in enumerate(output_text.splitlines()[j:], start=j) if line.strip() == ""), None)
            # allele frequencies are in lines j:k
            # print([get_data(line) for line in output_text.splitlines()[j:k]])
            res[marker] = np.array([get_data(line, K) for line in output_text.splitlines()[j:k]])
    return res

# get q from structure file (as lines)
# results in a dict { ind_name : q}
# can be transformed to a pd.DataFrame by usind
# pd.DataFrame.from_dict(get_q(lines), orient = 'index')
def get_q(lines):
    K = get_K(lines)
    # res stores the result
    res = {}
    # id stores if the identifiers were already used
    id = {}
    output_text = ''.join([line.decode('utf-8') for line in lines])
    start_index = next((i for i, line in enumerate(output_text.splitlines()) if "Inferred ancestry of individuals:" in line), None)+2
    end_index = next((i for i, line in enumerate(output_text.splitlines()) if "Estimated Allele Frequencies in each cluster" in line), None) - 1
    extracted_lines = output_text.splitlines()[start_index:end_index]
    # print(len(extracted_lines))
    for line in extracted_lines:
        if line != "":
            ind = line.split()[1]
            if ind in id.keys():
                id[ind] += 1
                # print(f"Identifier {ind} already taken. Using {ind}_{id[ind]} instead")
                ind = f"{ind}_{id[ind]}"
            else:
                id[ind] = 0
            res[ind] = np.array(line.split()[-K:], dtype = float)
    return res

# For a dict {key: value}, where values are np.arrays with the same number of columns, 
# this function returns a pd.DataFrame with all arrays stacked in top of each other.
# This function is applied to get_p(lines) and get_q(lines).
def to_df(res):
    return pd.DataFrame(np.vstack([value for _, value in res.items()]))

#####################################
## Importing data from other files ## 
#####################################

def load_file_notSTRUCTURE(path_q, path_p, path_pJ):
    data_p = pd.read_csv(path_p, delimiter=" ", header=None)
    data_q = pd.read_csv(path_q, delimiter=" ", header=None)
    if path_pJ:
        data_J = pd.read_csv(path_pJ, delimiter=" ", header=None)
    else:
         data_J = 0     
    return data_q, data_p, data_J

#################################################
## Optimizing routine for finding an optimal S ##
#################################################

# We start with a structure estimator 
# hatq = get_q(lines) (an NxK DataFrame), 
# hatp = get_p(lines) ((sum m*J_m)xK DataFrame) 
# from an output of structure
# We want to find a matrix S such that some target_function is minimized
# under the constraints 
# S.dot(ones(K))=1,
# 0 \leq hatq.dot S (<=1)
# 0\leq (S.linalg.inv).dot(p.T)\leq 1
# The matrix S is parametrized by the list x, where we will use
# S = np.reshape(np.array(x), (K,K))
np.array([[1,2],[3,4]]).dot(np.array([1,1]))
np.array([1,1]).dot(np.array([[1,2],[3,4]]))

# inds is a vector of 0/False and 1/True of length N. 
# This function is used in order to consider subsets of all N individuals
def subset_inds(hatq, inds):
    if hatq.shape[0] != len(inds):
        raise ValueError("In subset_inds, it must be hatq.shape[0] == len(inds)")
    # If inds = [1,1,0,1], A will be 
    # array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
    A = np.array([[1 if j == i else 0 for j in range(len(inds))] for i, x in enumerate(inds) if x == 1])
    return A.dot(hatq)

# The average entropy in q for individuals in inds   
def mean_entropy(x, hatq, inds = None):
    if inds:
        hatq = subset_inds(hatq, inds)
    K = hatq.shape[1]
    S = np.reshape(np.array(x), (K,K))
    res = entropy(hatq.dot(S), axis=1)
    return np.mean(np.maximum(res, 0))

# The average qS in population (column) pop along individuals in inds
# pop is a column index of hatq
def mean_size(x, hatq, pop, inds = None):
    if inds:
        hatq = subset_inds(hatq, inds)
    K = hatq.shape[1]
    S = np.reshape(np.array(x), (K,K))
    return np.mean((hatq.dot(S)).T[pop])

# We need the following two functions for the minimization
def neg_mean_entropy(x, hatq, inds = None):
    return - mean_entropy(x, hatq, inds)

def neg_mean_size(x, hatq, pop, inds = None):
    return - mean_size(x, hatq, pop, inds)

# This is a matrix A such that AS gives the constraint S1 = 1
def matrix_for_linear_constraint1(K):
    # Take the identity matrix, replace every i by [i]*K, and flatten the result
    # Example: K=2
    # Result is [[1, 1, 0, 0], [0, 0, 1, 1]]
    list_of_lists = [ [i]*K for i in np.identity(K).flatten()]
    y = [item for list in list_of_lists for item in list]
    # reshape in order to get K rows
    return np.reshape(y, (K, K*K))

# This is a matrix A such that AS gives the constraint hatq S \geq 0
def matrix_for_linear_constraint2(hatq):
    hatq = np.array(hatq)
    K = hatq.shape[1]
    # If hatq = [[0.5, 0.5], [0.8, 0.2]], the next line is
    # array([[0.5, 0. ], [0. , 0.5]]), 
    # array([[0.5, 0. ], [0. , 0.5]]), 
    # array([[0.8, 0. ], [0. , 0.8]]), 
    # array([[0.2, 0. ], [0. , 0.2]]
    arrays = [ np.hstack([np.identity(K) * qi for qi in q]) for q in hatq ]
    return np.vstack(arrays)

def function_for_nonlinear_constraint(x, hatp):
    K = hatp.shape[1]
    S = np.reshape(np.array(x), (K,K))
    T = np.linalg.inv(S)
    return T.dot(hatp.T).flatten()

# Find the best S for minimizing fun using all constraints
# For fun, we think of either of:
# fun = mean_size with args = (hatq, pop, inds)
# fun = neg_mean_size with args = (hatq, pop, inds)
# fun = mean_entropy with args = (hatq, inds)
# fun = neg_mean_entropy with args = (hatq, inds)
# n is the number of trials, the best is taken
def find_S(fun, hatq, hatp, n = 1, args = ()):
    K = hatq.shape[1]
    if K != hatp.shape[1]:
        raise ValueError("In find_S, hatq and hatp must have the same number of columns")
    # All entries in S cannot be smaller than -1 or larger than 1
    bounds = [(-1,1) for i in range(K*K)]    
    # The three constraints from above
    constraints = [
        LinearConstraint(matrix_for_linear_constraint1(K), 1, 1),
        LinearConstraint(matrix_for_linear_constraint2(hatq), 0, np.inf),
        NonlinearConstraint(lambda x: function_for_nonlinear_constraint(x, hatp), 0, np.inf)
        ]
    best = None
    for _ in range(n):
        # x0 is close to the identity, but some U[-.2,.2] away
        x0 = np.array([[np.random.uniform(-1/K, 1/K) for e in range(K)] for e in range(K)])
        x0 = np.identity(K) #+ x0
        # rows in x0 must be such that the corresponding S has row sum of 1
        x0 = (x0.dot(np.diag([1/z for z in x0.sum(axis=1)]))).flatten()
        res = minimize(fun = fun, args = args, x0 = x0, constraints = constraints, bounds = bounds)
        if res.success:
            try:
                if res.fun < best.fun:
                    best = res 
            except:
                best = res
    return best

def find_q_forKeq2(hatq, hatp, ind, pop, maximize = True):
    if hatq.shape[1] != 2 or hatp.shape[1] != 2 :
        raise ValueError("In find_q_forKeq2, hatq and hatp must have 2 columns.")
    if pop not in [0,1]:
        raise ValueError("In find_q_forKeq2, pop must be 0 or 1.")
    # We want to minimize/maximize column 0
    if pop == 1:
        hatq = hatq[:,[1,0]]
        hatp = hatp[:,[1,0]]

    U = max(hatp[:,1]/hatp[:,0]) # as u^\ast in the manuscript
    u = min(hatp[:,1]/hatp[:,0]) # as u_\ast in the manuscript
    V = max(hatq[:,0]/hatq[:,1]) # as v^\ast in the manuscript
    v = min(hatq[:,0]/hatq[:,1]) # as v_\ast in the manuscript

    if maximize:
        res = hatq[ind, 0] * (1 + V)/(u + V) + hatq[ind, 1] * (1 + V)/(1 + V / u)
    else:
        res = hatq[ind, 0] * (U - 1)/(U - v) + hatq[ind, 1] * (1 - U)/(1 + U / v)
    res = np.array([res, 1-res])

    if pop == 1:
        res = res[[1,0]]

    return res


# After finding the optimal S, we can also report the optimal q for minimizing fun
# We usually only report the optimal q for individuals in inds (which were used to minimize fun)
def find_q(fun, hatq, hatp, inds = None, n=1, args = ()):
    K = hatq.shape[1]
    # Does the formula apply?
    if K == 2 and sum(inds) == 1:        
        ind = inds.index(True)
        q_max = find_q_forKeq2(hatq, hatp, ind, 0, maximize = True)
        q_min = find_q_forKeq2(hatq, hatp, ind, 0, maximize = False)
        
        if fun.__name__ == mean_size:
            res = q_min
        elif fun.__name__ == neg_mean_size:
            res = q_max
        if fun.__name__ == mean_entropy:
            # entropy is minimal if individuals are as non-admixed as possible
            if min(q_min, 1-q_min) < min(q_max, 1-q_max): 
                res = q_min
            else:
                res = q_max
        if fun.__name__ == neg_mean_entropy:
            # Entropy is maximal if the solution closest to [0.5, 0.5]
            # If q_max[0] < 0.5, report q_max[0]
            # If q_min[0] > 0.5, report q_min[0]
            # If q_min[0] < 0.5 < q_max[0], report 0.5
            x = min(q_max[0], max(q_min[0],0.5))
            res = [[x, 1-x]]   

    else:
        res = find_S(fun, hatq, hatp, n, args)
        if res:
            S = np.reshape(res.x, (K,K))
        else:
            print("No optimum found. Proceeding with initial values.")
            S = np.identity(K)
        if inds:
            hatq = subset_inds(hatq, inds)
    
    return hatq.dot(S)

# Some tests
if __name__ == "__main__":
    np.random.seed(19)
#    default_STRUCTURE_path = 'Example_Input/output_structure_f'
    default_STRUCTURE_path = 'Example_Input/output_structure_K3_f'
#    default_STRUCTURE_path = 'Example_Input/output_structure_K3_not_biallelic_f'
    lines = load_structure_file(default_STRUCTURE_path)
    hatq_df = get_q(lines)
    hatq = to_df(hatq_df)
    # print(lines)
    hatp_df = get_p(lines)    
    # print(hatp_df)
    hatp = to_df(hatp_df)
    K = get_K(lines)
    N = hatq.shape[0]
    pop = 0
    first = 10
    inds = [1 if i < first else 0 for i in range(N)]
    print(f"Minimizing the contribution of population 0 in the first {first} individuals in {default_STRUCTURE_path}:")
    # fun = (lambda x : mean_size(x, hatq, pop, inds))
    print(find_S(mean_size, hatq, hatp, 10, (hatq, pop, inds)))
    print(find_q(mean_size, hatq, hatp, inds, 10, (hatq, pop, inds)))

    print(f"Minimizing the entropy in the first {first} individuals in {default_STRUCTURE_path}:")
    # fun = (lambda x : mean_entropy(x, hatq, inds))
    print(find_q(mean_entropy, hatq, hatp, inds, 10, (hatq, inds)))
    
    print("Generating some random data (with fixed seed) with K=3.")
    # np.random.seed(41)
    inds = None
    hatq1 = np.random.uniform(0.2, 0.8, size=(1000, 1)) 
    hatq2 = np.random.uniform(0.2, 0.8, size=(1000, 1)) 
    a = np.minimum(hatq1, hatq2)
    b = np.maximum(hatq1, hatq2)
    hatq = np.hstack([a, b-a, 1-b])
    print("hatq.T:")
    print(hatq.T)
    print("hatp.T:")
    hatp = array = np.random.uniform(0.3, 0.7, size=(200, 3)) 
    print(hatp.T)

    print("Maximizing the contribution of population 0 in all individuals:")
    # fun = (lambda x : -mean_size(x, hatq, pop))
    print("Initial IAs:")
    print(hatq[0:9])
    print("Optimized IAs:")
    print(find_q(neg_mean_size, hatq, hatp, inds, 10, (hatq, pop))[0:9])
    print("Maximizing the entropy in all individuals:")
    # fun = (lambda x : neg_mean_entropy(x, hatq, pop))
    print("Initial IAs:")
    print(hatq[0:9])
    print("Optimized IAs:")
    print(find_q(neg_mean_entropy, hatq, hatp, inds, 10, (hatq, pop))[0:9])
    

