import os
import re
import numpy as np
import pandas as pd
import io  

file_pathK = "C:\\Users\\carol\\Downloads\\CEU_IBS_TSI_enhanced_corr_f" # Filepath
K = 2  # Number of populations
output_directory = "C:\\Users\\carol\\Downloads\\"  # Adjust this to your desired output directory

file_name_q = "q_STRUCTURE"
file_name_p = "p_CEU_IBS_TSI_K2"
file_name_pJ = "p_CEU_IBS_TSI_K2_J"
#%%
# Function to read table data
def read_table_data(file_path, K):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    # Find indices where "0.0% missing data" appears
    start_indices = [i for i, line in enumerate(lines) if "0.0% missing data" in line]
    
    res = np.zeros((0, K))  # Empty array to collect results
    res_all = np.zeros((0, K))
    for j, start_index in enumerate(start_indices):
        start_index += 1  # Move to the next line after "0.0% missing data"

        # Find end index (first empty line after start index)
        end_index = start_index + next((i for i, line in enumerate(lines[start_index:], start=1) if line.strip() == ""), None)
        
        # Combine the relevant lines into a single string
        text = "\n".join(lines[start_index:end_index])
        if (start_index - end_index) < -2:
            # Remove the numbers followed by (0.x) as in R gsub() 
            modified_text = re.sub(r"\s+[0-9]+\s*\(0\.[0-9]+\)", "", text)
            # Split the modified text by spaces and filter out empty strings
            values = [val for val in modified_text.split() if val != ""]
            numeric_values = np.array(values, dtype=float)
            # Sum corresponding values for the two K sections

            res1 = numeric_values[:K] + numeric_values[K:2*K]
            if(sum(res1) < 2):
            # Add the results to the res array
                res = np.vstack([res, res1])
                res_all = np.vstack([res_all, numeric_values[K:2*K]])
            res_all = np.vstack([res_all, numeric_values[:K]])

    return res, res_all

# Call the function to read data
p = read_table_data(file_pathK, K)

output_filepath = os.path.join(output_directory, file_name_pJ)
np.savetxt(output_filepath, p[0], delimiter=" ", fmt="%f")

output_filepath = os.path.join(output_directory, file_name_p)
np.savetxt(output_filepath, p[1], delimiter=" ", fmt="%f")

# Function to extract q matrix from the file
def extract_q(filepath, K):
    with open(filepath, 'r') as file:
        lines = file.readlines()

    # Find indices for "Inferred ancestry of individuals:" and "Estimated Allele Frequencies in each cluster"
    start_index = next(i for i, line in enumerate(lines) if "Inferred ancestry of individuals:" in line)
    end_index = next(i for i, line in enumerate(lines) if "Estimated Allele Frequencies in each cluster" in line)

    # Extract lines between these two markers
    extracted_lines = lines[start_index + 2:end_index - 1]

    # Read the extracted lines as a DataFrame
    res = pd.read_csv(io.StringIO("\n".join(extracted_lines)), delim_whitespace=True, header=None)

    C = res.shape[1]
    res_final = res.iloc[:, (C-K):C]

    return res_final

# Estimated IAs
q = extract_q(file_pathK, K)
# Save q matrix to file
output_filepath_q = os.path.join(output_directory, file_name_q)
q.to_csv(output_filepath_q, sep=" ", index=False, header=False)


