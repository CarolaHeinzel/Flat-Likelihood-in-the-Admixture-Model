import os
import re
import numpy as np
import pandas as pd
import io  

default_structure_path = 'Example_Input/output_structure_f'

def load_structure_file(path = default_structure_path):
    with open(path, 'rb') as f:
        return f.readlines()

# Extract p all
def split_array(arr, K):
    # print(arr)
    l = len(arr)
    l = int(l/K)
    temp = np.reshape(arr, (l, K))
    # print(temp, "t")
    return temp.tolist()

def format_and_save(data, K):
    # Flatten the nested list
    flattened_data = [item for sublist in data for row in sublist for item in row]
 
    # Split the flattened list into chunks of size K
    reshaped_data = [flattened_data[i:i+K] for i in range(0, len(flattened_data), K)]
    return reshaped_data

# get p from structure file (as lines)
def get_p(lines):
    K = get_population_count(lines)
    output_text = ''.join([line.decode('utf-8') for line in lines])
    start_indices = [i+1 for i, line in enumerate(output_text.splitlines()) if " missing data" in line]
    res_all = res = np.zeros((0, K))  
    marker = []
    pattern = r'Locus\s+\d+\s*:\s*(\S+)'
    lines_decoded = [line.decode('utf-8') for line in lines]  # Assuming 'utf-8' encoding
    for line in lines_decoded:
        match = re.search(pattern, line)
        if match:
            marker.append(match.group(1))
    for j, start_index in enumerate(start_indices):
        end_index = start_index + next((i for i, line in enumerate(output_text.splitlines()[start_index:], start=1) if line.strip() == ""), None)



        text = "\n".join(output_text.splitlines()[start_index:end_index])       
        if (end_index - start_index) > 2:  
            modified_text = re.sub(r"\s+[0-9]+\s*\(0\.[0-9]+\)", "", text)
            values = [val for val in modified_text.split() if val != ""]        
            numeric_values = np.array(values, dtype=float)          
            if len(numeric_values) >= 2 * K:
                res1 = numeric_values[:K] + numeric_values[K:2*K]
                if np.sum(res1) < K:  
                    res = np.vstack([res, res1])
                    res_all = np.vstack([res_all, numeric_values[K:2*K]])
                res_all = np.vstack([res_all, numeric_values[:K]])
    res_allp = []
    for j, start_index in enumerate(start_indices):
        # Find end index (first empty line after start index)
        end_index = start_index + next((i for i, line in enumerate(output_text.splitlines()[start_index:], start=1) if line.strip() == ""), None)
        # Combine the relevant lines into a single string
        text = "\n".join(output_text.splitlines()[start_index:end_index])
        if (start_index - end_index) < -2:
            # Remove the numbers followed by (0.x) as in R gsub() 
            modified_text = re.sub(r"\s+[0-9]+\s*\(0\.[0-9]+\)", "", text)
            # Split the modified text by spaces and filter out empty strings
            values = [val for val in modified_text.split() if val != ""]
            numeric_values = np.array(values, dtype=float)
            # Sum corresponding values for the two K sections
            res_allp.append(split_array(numeric_values, K))
   
    res_allp = format_and_save(res_allp, K)
    if(len(res) == 0):
        res = None
    return res, res_all, marker, res_allp

# get q from structure file (as lines)
def get_q(lines):
    K = get_population_count(lines)
    output_text = ''.join([line.decode('utf-8') for line in lines])
    start_index = next((i for i, line in enumerate(output_text.splitlines()) if "Inferred ancestry of individuals:" in line), None)+2
    end_index = next((i for i, line in enumerate(output_text.splitlines()) if "Estimated Allele Frequencies in each cluster" in line), None) - 1
    extracted_lines = output_text.splitlines()[start_index:end_index]
    df = pd.read_csv(io.StringIO("\n".join(extracted_lines)), delim_whitespace=True, header=None)
    df.set_index(1, inplace=True)
    df = df.iloc[:, -K:]  
    column_names = list(df.columns)
    df.rename({a: column_names.index(a) for a in column_names}, axis=1, inplace=True)
    return df

# Determine K for the structure File    
def get_population_count(structure_file):
    output_text = ''.join([line.decode('utf-8') for line in structure_file])
    match = re.search(r"(\d+)\s+populations assumed", output_text)
    if match:
        return int(match.group(1))  
    else:
        return None 

default_STRUCTURE_path = 'Example_Input/output_structure_f'
lines = load_structure_file()
K = get_population_count(lines)
res = get_p(lines)
print(res)  
# data_pJ, data_p, marker_names, p_all = read_table_data(uploaded_file, K)
#print(p_all)
#data_q, individual_names = extract_q(uploaded_file, K)
