import os
import re
import numpy as np
import pandas as pd
import io  

# Extract p
def read_table_data(lines, K):

    output_text = ''.join([line.decode('utf-8') for line in lines])
    

    start_indices = [i for i, line in enumerate(output_text.splitlines()) if " missing data" in line]
    
    res = np.zeros((0, K))  
    res_all = np.zeros((0, K))  
    locus_matches = []
    pattern = r'Locus\s+\d+\s*:\s*(\S+)'
    lines_decoded = [line.decode('utf-8') for line in lines]  # Assuming 'utf-8' encoding
    for line in lines_decoded:
        match = re.search(pattern, line)
        if match:
            locus_matches.append(match.group(1))
    for j, start_index in enumerate(start_indices):
        start_index += 1  


        end_index = start_index + next((i for i, line in enumerate(output_text.splitlines()[start_index:], start=1) if line.strip() == ""), None)


        text = "\n".join(output_text.splitlines()[start_index:end_index])
        print(text)
        
        if (end_index - start_index) > 2:  
            modified_text = re.sub(r"\s+[0-9]+\s*\(0\.[0-9]+\)", "", text)
            print(modified_text)
            

            values = [val for val in modified_text.split() if val != ""]
            print(values)
            

            numeric_values = np.array(values, dtype=float)
            

            if len(numeric_values) >= 2 * K:
                res1 = numeric_values[:K] + numeric_values[K:2*K]
                if np.sum(res1) < 2:  # Bedingung für die Summe
                    # Füge die Ergebnisse zum res-Array hinzu
                    res = np.vstack([res, res1])
                    res_all = np.vstack([res_all, numeric_values[K:2*K]])
                res_all = np.vstack([res_all, numeric_values[:K]])
    if(len(res) == 0):
    	res = None
    return res, res_all, locus_matches

# Extract q
def extract_q(lines, K):
    output_text = ''.join([line.decode('utf-8') for line in lines])

    
    start_index = [i for i, line in enumerate(output_text.splitlines()) if "Inferred ancestry of individuals:"  in line]
    end_index = [i for i, line in enumerate(output_text.splitlines()) if "Estimated Allele Frequencies in each cluster"  in line]
    output_text = ''.join([line.decode('utf-8') for line in lines])


    start_index = next((i for i, line in enumerate(output_text.splitlines()) if "Inferred ancestry of individuals:" in line), None)
    end_index = next((i for i, line in enumerate(output_text.splitlines()) if "Estimated Allele Frequencies in each cluster" in line), None)

    extracted_lines = output_text.splitlines()[start_index + 2:end_index - 1]


    res = pd.read_csv(io.StringIO("\n".join(extracted_lines)), delim_whitespace=True, header=None)
    individuals_names = res.iloc[:,1]

    res_final = res.iloc[:, -K:]  

    return res_final, individuals_names



