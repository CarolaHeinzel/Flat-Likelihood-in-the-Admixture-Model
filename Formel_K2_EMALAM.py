import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve
import pandas as pd
from sympy import symbols, Matrix
from scipy.optimize import minimize
import itertools
import sympy as sp
#%%
# Application
# 1) Load Data, i.e. insert here the correct names for q and p
file_path_p = 'C:\\Users\\carol\\Downloads\\p_example_strong.txt'
file_path_q = 'C:\\Users\\carol\\Downloads\\q_example_marker.txt'
data_q = pd.read_csv(file_path_q, sep=' ', header=None)
data_p = pd.read_csv(file_path_p, sep=' ', header=None)
#%%
def determine_p(data_p):
    result_list = data_p[0].tolist()
    result_list1 = data_p[1].tolist()
    M = len(result_list)
    s1 = 10**100
    sM = -10**100
    
    for m in range(M):
        temp_p = (1-result_list1[m])/(1-result_list[m])
        temp_p1 = result_list1[m]/result_list[m]
        temp_min = min(temp_p, temp_p1)
        temp_max = max(temp_p, temp_p1)
        if(s1 > temp_min):
            s1 = temp_min
        if(sM < temp_max):
            sM = temp_max
    return s1, sM
print(determine_p(data_p))

#%%
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


def calc_p_q(data_p, data_q, name_q, name_p):
    max_value = data_q[0].max()
    min_value = data_q[0].min()
    s1, sM = determine_p(data_p)

    a_max = (s1-1)/(max_value/(1-max_value) + s1)
    b_max = 1 + a_max * (max_value/(1-max_value))


    a_min = 1/(1-min_value)/(sM + min_value/(1-min_value))
    b_min =  min_value/(1-min_value)*(a_min - 1)
    
    N = len(data_q[0])
    for i in range(N):
        q_max = data_q[0]*(1-a_max) + (1-data_q[0])*b_max
        q_min = data_q[0]*(1-a_min) + (1-data_q[0])*b_min
        p_max_1 = (data_p[0] - b_max*data_p[0] - a_max*data_p[1])/(1 - a_max - b_max)
        p_min_1 = (data_p[0] - b_min*data_p[0] - a_min*data_p[1])/(1 - a_min - b_min)
        p_max_2 = (data_p[1] - a_max*data_p[1] - b_max*data_p[0])/(1 - a_max - b_max)
        p_min_2 = (data_p[1] - a_min*data_p[1] - b_min*data_p[0])/(1 - a_min - b_min)
    result_qmax = [q_max.tolist(), (1-q_max).tolist()]
    result_qmin = [q_min.tolist(), (1-q_min).tolist()]
    result_pmax = [p_max_1.tolist(), p_max_2.tolist()]
    result_pmin = [p_min_1.tolist(), p_min_2.tolist()]

    save_values(result_qmax, name_q, "max")
    save_values(result_pmax, name_p, "max")
    save_values(result_qmin, name_q, "min")
    save_values(result_pmin, name_p, "min")
    return 
print(calc_p_q(data_p, data_q, "formel_test_q_strong", "formel_test_p_strong"))

