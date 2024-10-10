import numpy as np
import pandas as pd
from scipy.stats import entropy
from scipy.optimize import minimize, LinearConstraint, NonlinearConstraint
import re, argparse, itertools

#################################################
## Importing data from a structure output file ## 
#################################################

def load_structure_file(path):
    with open(path, 'rb') as f:
        return f.readlines()

def load_structure_file(path):
    with open(path, 'rb') as f:
        return f.readlines()

# result is hatq_df
def load_q_file(path):
    return pd.read_csv(path, header = None, sep = " ")

# File must only contain bi-allelic marker and report the frequencies of a reference allele
# This is as in admixture output.
# result is hatp_dict
def load_p_file(path):
    df = pd.read_csv(path, header = None, sep = " ")
    res = { i : np.vstack([df.iloc[i,], 1-df.iloc[i,]]) for i in range(len(df[0]))}
    return res

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

## The following five functions are for the Jacobi matrices for the functions we want to optimize.

def entropy_jac(a):
    return -1 - np.log(np.maximum(a, 1e-12))

def mean_size_jac(x, hatq, pop, inds = None):
    K = hatq.shape[1]
    if inds:
        hatq = subset_inds(hatq, inds)
    I = hatq.shape[0]
    Apop = np.array([[1 if i == pop + j*K else 0 for i in range(K*K)] for j in range(K)])
    return np.ones(I).dot(hatq).dot(Apop) / I

def neg_mean_size_jac(x, hatq, pop, inds = None):
    return - mean_size_jac(x, hatq, pop, inds)

def mean_entropy_jac(x, hatq, inds = None):
    K = hatq.shape[1]
    if inds:
        hatq = subset_inds(hatq, inds)
    I = hatq.shape[0]
    A = { k: np.array([[1 if i == k + j*K else 0 for i in range(K*K)] for j in range(K)]) for k in range(K) }
    res = 0
    for i in range(I):
        for k in range(K):
            loc = hatq[i,:].dot(A[k])
            res = res + entropy_jac(loc.dot(x)) * loc
    return res / I

def neg_mean_entropy_jac(x, hatq, inds = None):
    return - mean_entropy_jac(x, hatq, inds)

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
    T = np.linalg.pinv(S) # more stable than .inv
    return T.dot(hatp.T).flatten()

def get_UuVv(hatp, hatq):
    U = max(hatp[:,1]/hatp[:,0]) # as u^\ast in the manuscript
    u = min(hatp[:,1]/hatp[:,0]) # as u_\ast in the manuscript
    V = max(hatq[:,0]/hatq[:,1]) # as v^\ast in the manuscript
    v = min(hatq[:,0]/hatq[:,1]) # as v_\ast in the manuscript
    U = min(U, 1e12)
    V = min(V, 1e12)
    u = max(u, 1e-12)
    v = max(v, 1e-12)
    return U, u, V, v

# Compute the distance of hatp and S^-1 hatp if usep and hatq S otherwise.
def dist(S, hat, minp = False):
    if minp:
        T = np.linalg.pinv(S) # more stable than .inv        
        return sum(abs((hat - T.dot(hat.T)).flatten()))
    else: 
        return sum(abs((hat - hat.dot(S)).flatten()))

# Returns S with switched labels, such that qS or S^-1 p are closest to hatq or hatp, respectively.
# hat is either hatq or hatp
def switch_labels_for_minimal_distance(S, hat, minp = False):
    K = hat.shape[1]
    best = S
    distbest = dist(S, hat, minp)
    for perm in itertools.permutations(range(K)):
        P = [[True if i == k else False for i in range(K)] for k in perm]
        res = S.dot(np.array(P))
        distres = dist(res, hat, minp)
        if distres < distbest:
            distbest = distres
            best = res
            #print("Labels are switched")
    return best

# Find the best S for minimizing fun using all constraints
# For fun, we think of either of:
# fun = mean_size with args = (hatq, pop, inds)
# fun = neg_mean_size with args = (hatq, pop, inds)
# fun = mean_entropy with args = (hatq, inds)
# fun = neg_mean_entropy with args = (hatq, inds)
# n is the number of trials, the best is taken
def find_S(fun, hatq, hatp, n = 1, args = (), jac = None, switch_labels = True, minp = False):
    K = hatq.shape[1]
    if K != hatp.shape[1]:
        raise ValueError("In find_S, hatq and hatp must have the same number of columns")
    # All entries in S cannot be smaller than -1 or larger than 1
    if K == 2 and (fun.__name__ == "mean_size" or fun.__name__ == "neg_mean_size"): 
        U, u, V, v = get_UuVv(hatp, hatq)
        # print(f"U = {U}, u = {u}, V = {V}, v = {v}")
        if fun.__name__ == "mean_size":
            a = (1 + v) / (U + v)
            b = (1 - U) / (1 + U / v)
        else:
            a = (u - 1) / (u + V)
            b = (1 + V) / (1 + V / u)
        best = np.array([[1-a, a], [b, 1-b]])
    else:
        bounds = [(-2,2) for i in range(K*K)]
        # The three constraints from above
        constraints = [
            LinearConstraint(matrix_for_linear_constraint1(K), 1, 1),
            LinearConstraint(matrix_for_linear_constraint2(hatq), 0, np.inf),
            NonlinearConstraint(lambda x: function_for_nonlinear_constraint(x, hatp), 0, np.inf)
            ]
        best = None
        for _ in range(n):
            # x0 is close to the identity, but some U[-.2,.2] away
            x0 = np.array([[np.random.uniform(-10/K, 10/K) for e in range(K)] for e in range(K)])
            x0 = np.identity(K) #+ x0
            # rows in x0 must be such that the corresponding S has row sum of 1
            x0 = (x0.dot(np.diag([1/z for z in x0.sum(axis=1)]))).flatten()
            res = minimize(fun = fun, args = args, x0 = x0, jac = jac, constraints = constraints, bounds = bounds)
            if res.success:
                try:
                    if res.fun < best.fun:
                        best = res 
                except:
                    best = res
            if best:
                best = np.reshape(best.x, (K,K))
        
        if switch_labels:
            best = switch_labels_for_minimal_distance(best, hatp, minp)

    return np.array(best)

# After finding the optimal S, we can also report the optimal q for minimizing fun
# n is the number of trials to find an optimum
def find_q(fun, hatq, hatp, n=1, args = (), jac = None):
    N = hatq.shape[0]
    K = hatq.shape[1]
    # Does the formula apply?
    if K == 2 and (fun.__name__ == "mean_size" or fun.__name__ == "neg_mean_size"):
        # We want to minimize/maximize column 0
        if args[1] == 1:
            hatq = hatq[:,[1,0]]
            hatp = hatp[:,[1,0]]
        U, u, V, v = get_UuVv(hatp, hatq)
        # print(f"U = {U}, u = {u}, V = {V}, v = {v}")
        if fun.__name__ == "neg_mean_size":
            res = hatq[:, 0] * (1 + V)/(u + V) + hatq[:, 1] * (1 + V)/(1 + V / u)
        else:
            res = hatq[:, 0] * (U - 1)/(U - v) + hatq[:, 1] * (1 - U)/(1 + U / v)
        # print(f"res vor b {res}")
        res = np.vstack([res, 1-res]).T
        if args[1] == 1:
            res = res[1,0]
        res = np.array(res).reshape(N,K)
    else:
        S = find_S(fun, hatq, hatp, n, args, jac = jac)
        if S is not None:
            S = np.reshape(S, (K,K))
        else:
            print("No optimum found. Proceeding with initial values.")
            S = np.identity(K)
        res = hatq.dot(S)    
    # print(f"res = {res}")
    return res

# After finding the optimal S, we can also report the optimal q for minimizing fun
# n is the number of trials to find an optimum
def find_p(fun, hatq, hatp, n=1, args = (), jac = None):
    M = hatp.shape[0]
    K = hatp.shape[1]
    # Does the formula apply?
    if K == 2 and (fun.__name__ == "mean_size" or fun.__name__ == "neg_mean_size"):
        # We want to minimize/maximize column 0
        if args[1] == 1:
            hatq = hatq[:,[1,0]]
            hatp = hatp[:,[1,0]]
        U, u, V, v = get_UuVv(hatp, hatq)
        # print(f"U = {U}, u = {u}, V = {V}, v = {v}")
        res1 = hatp[:, 0] * (U)/( U - 1 ) + hatp[:, 1] * ( -1 )/( U - 1 )
        res2 = hatp[:, 0] * ( -V )/( 1 + V ) + hatp[:, 1] * ( 1 )/( 1 + V )
        if fun.__name__ == "neg_mean_size":
            res = np.minimum(res1, res2)
        else:
            res = np.maximum(res1, res2)
        # print(f"res vor b {res}")
        res = np.vstack([res, 1-res]).T
        if args[1] == 1:
            res = res[1,0]
        res = np.array(res).reshape(M,K)
    else:
        S = find_S(fun, hatq, hatp, n, args, jac)
        T = np.linalg.pinv(S) # more stable than .inv
        res = T.dot(hatp.T)
    # print(f"res = {res}")
    return res

############################
## Functions for plotting ##
############################

def get_q_for_plot(q):
    N = q.shape[0]
    K = q.shape[1]    
    q_df = pd.DataFrame({
                    "ia": [x for qi in q for x in qi],
                    "ind": [x for i in range(N) for x in [i]*K],
                    "pop": [ j for i in range(N) for j in range(K)]
                })
    pivot_df = q_df.pivot(index='ind', columns='pop', values='ia')
    return q_df, pivot_df

def get_K_from_hatp_dict(hatp_dict):
    first_key = next(iter(hatp_dict))  # Erster Key
    first_value = hatp_dict[first_key]  # Erster Value
    return first_value.shape[1]

def get_p_for_plot(hatp_dict):
    K = get_K_from_hatp_dict(hatp_dict)
    M = len(hatp_dict.keys())
    # reformat marker names if they are ints
    try:
        hatp_dict = { int(key): value for key, value in hatp_dict.items()}
        m = max(hatp_dict.keys())
        hatp_dict = { f"{key:{len(str(m))}}": value for key, value in hatp_dict.items()}
    except:
        pass
    
    marker = []
    for k in range(K):
        for key, value in hatp_dict.items():
            J = value.shape[0]
            if J > 1:
                for j in range(J):
                    marker.append({
                        "id": key,
                        "p": hatp_dict[key][j,k],
                        "allele": j,
                        "pop": k,
                        "id.allele": f"{key}.{j}"
                    })


    p_df  = pd.DataFrame(marker)
    # st.write(p_df)
    pivot_df = p_df.pivot(index='id.allele', columns='pop', values='p')
    return p_df, pivot_df


# Command line tool
# Examples:
# python emalam.py --structure_filename Example_Input/CEU_IBS_TSI_enhanced_corr_K3_f --out out.Q out.P --fun entropy --min 
# python emalam.py --hatq_filename Example_Input/p_K3 --hatp_filename Example_Input/q_K3 --out out.Q out.P --fun entropy --min | less

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process input file names.")
    parser.add_argument("--hatq_filename", type=str, help="The input file name for hatq")
    parser.add_argument("--hatp_filename", type=str, help="The input file name for hatp")
    parser.add_argument("--structure_filename", type=str, help="The input file name for structure data")
    parser.add_argument("--out", nargs = 2, type=str, help="The output filenames for the optimized parameters. (Both are csv files.)")
    parser.add_argument("--fun", type=str, help="The target function to optimize, must be entropy or size", default = "entropy")
    parser.add_argument("--pop", type=str, help="If fun == size, the number of the population which is to be optimized.")
    parser.add_argument("--min", action='store_true', help="The target function is minimized.")
    parser.add_argument("--max", action='store_true', help="The target function is maximized.")
    parser.add_argument("--n", type=str, help="Number of iterations for the optimization.", default = 1)
    parser.add_argument("--inds", nargs = '+', type=str, help="The individuals which are used for the target function. If no names are given, a number starting with 0 is used. If missing, optimization is over all individuals.")
    parser.add_argument("--no_switch_labels", action='store_true', help="The optimum is given as is, and it is not checked if a relabeling of populations leads to a result which is closer to the input data. ")
    parser.add_argument("--minp", action='store_true', help="Only effective if no_swith_labels is set. The distance to the input data is computed with respect to the optimal p (rather than the optimal q). ")
    
    args = parser.parse_args()

    if not ((args.hatq_filename and args.hatp_filename) or args.structure_filename):
        raise ValueError("Either, structure_filename, or hatq_filename and hatp_filename must be given.")
    if args.structure_filename:
        lines = load_structure_file(args.structure_filename)
        hatq_dict = get_q(lines)
        hatq_df = to_df(hatq_dict)
        hatq = np.array(hatq_df)
        ind_ids = hatq_dict.keys()
        hatp_dict = get_p(lines)
    else:
        hatq = np.array(load_q_file(args.hatq_filename))
        N = hatq.shape[0]
        ind_ids = range(N) 
        hatp_dict = load_p_file(args.hatp_filename)

    hatp_df = to_df(hatp_dict)
    hatp = np.array(hatp_df)
    
    if args.min and args.max:
        raise ValueError("min and max must not be True simultaneously. ")
    if args.fun == "entropy":
        if args.min:
            f = {"fun": mean_entropy, "jac": mean_entropy_jac}
        elif args.max:
            f = {"fun": neg_mean_entropy, "jac": neg_mean_entropy_jac}
    elif args.fun == "size":
        if args.min:
            f = {"fun": mean_size, "jac": mean_size_jac}
        elif args.max:
            f = {"fun": neg_mean_size, "jac": neg_mean_size_jac}
    else:
        raise ValueError("fun must be 'entropy' or 'size'.")        
    
    if args.inds:
        wrong_inds = [id for id in args.names if id not in ind_ids]
        if wrong_inds:
            print(f"Warning: Individuals {wrong_inds} not in input file(s).")
        inds = [id in args.names for id in ind_ids]
    else:
        inds = None

    switch_labels = False if args.no_switch_labels else True
    
    S_opt = find_S(f["fun"], hatq, hatp, args.n, (hatq, inds), f["jac"], switch_labels, args.minp)
    q_opt = hatq.dot(S_opt)
    q_df, q_pivot = get_q_for_plot(q_opt)
    T_opt = np.linalg.pinv(S_opt)
    p_opt_dict = { key: value.dot(T_opt.T) for key, value in hatp_dict.items() } 
    p_df, p_pivot = get_p_for_plot(p_opt_dict)

    q_pivot.to_csv(args.out[0], header = False, sep = ' ')
    p_pivot.to_csv(args.out[1], header = False, sep = ' ')
    print(f"The result is in {args.out[0]} (for Q) and {args.out[1]} (for P).")    

