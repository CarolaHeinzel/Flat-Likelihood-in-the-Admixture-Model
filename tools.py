import os
import re
import numpy as np
import pandas as pd
import io  

#################################################
## Importing data from a structure output file ## 
#################################################

def load_structure_file(path):
    with open(path, 'rb') as f:
        return f.readlines()

# Determine K from the structure File    
def get_population_count(structure_file):
    output_text = ''.join([line.decode('utf-8') for line in structure_file])
    match = re.search(r"(\d+)\s+populations assumed", output_text)
    if match:
        return int(match.group(1))  
    else:
        raise ValueError(f"K not found in {structure_file}")
        return None 

# In a line of the form
# 0   (0.826) 0.542 0.999
# return [0.542, 0.999] : np.array
def get_data(line):
    modified_line = re.sub(r"\s+[0-9]+\s*\(0\.[0-9]+\)", "", line)
    values = [val for val in modified_line.split() if val != ""]        
    return np.array(values, dtype=float)          

# get p from structure file (as lines)
# returns a dict key: value, where key is an identifier of a marker, and value is a K x J_m array
# with the allele frequencies in all populations at all allelels.
def get_p(lines):
    lines_decoded = [line.decode('utf-8') for line in lines]  # Assuming 'utf-8' encoding
    output_text = ''.join([line.decode('utf-8') for line in lines])
    
    res = {}
    # The information for the next marker starts like this:
    pattern = r'Locus\s+\d+\s*:\s*(\S+)'
    for i, line in enumerate(lines_decoded):
        match = re.search(pattern, line)
        if match:
            marker = match.group(1)
            j = next((k for k, line in enumerate(output_text.splitlines()[i:], start=i) if " missing data" in line), None) 
            k = next((l for l, line in enumerate(output_text.splitlines()[j:], start=j) if line.strip() == ""), None)-1
            # allele frequencies are in lines j:k
            print([get_data(line) for line in output_text.splitlines()[j:k]])
            res[marker] = np.array([get_data(line) for line in output_text.splitlines()[j:k]])
    return res

# get q from structure file (as lines)
# results in a dict { ind_name : q}
# can be transformed to a pd.DataFrame by usind
# df.fromDict(, orient = 'index)
def get_q(lines):
    K = get_population_count(lines)
    res = {}
    output_text = ''.join([line.decode('utf-8') for line in lines])
    start_index = next((i for i, line in enumerate(output_text.splitlines()) if "Inferred ancestry of individuals:" in line), None)+2
    end_index = next((i for i, line in enumerate(output_text.splitlines()) if "Estimated Allele Frequencies in each cluster" in line), None) - 1
    extracted_lines = output_text.splitlines()[start_index:end_index]
    for line in extracted_lines:
        if line != "":            
            print(f"line: {line}")
            ind = line.split()[1]
            q = np.array([val for val in line.split() if val != ""])
    res[ind] = q
    return res










default_STRUCTURE_path = 'Example_Input/output_structure_f'
lines = load_structure_file(default_STRUCTURE_path)
K = get_population_count(lines)
res = get_p(lines)
print(res)  
res = get_q(lines)
print(res)  
# data_pJ, data_p, marker_names, p_all = read_table_data(uploaded_file, K)
#print(p_all)
#data_q, individual_names = extract_q(uploaded_file, K)









