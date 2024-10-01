import os
import re
import numpy as np
import pandas as pd
import io  

file_pathK = "C:\\Users\\carol\\Downloads\\structure_output.txt_f" # Filepath
K = 2  # Number of populations
output_directory = "C:\\Users\\carol\\Downloads\\"  # Adjust this to your desired output directory

file_name_q = "q_STRUCTURE"
file_name_p = "p_CEU_IBS_TSI_K2"
file_name_pJ = "p_CEU_IBS_TSI_K2_J"


def split_array(arr, K):
    print(arr)
    l = len(arr)
    l = int(l/K)
    temp = np.reshape(arr, (l, K))
    print(temp, "t")
    return temp.tolist()

# Function to read table data
def read_table_data_all(file_path, K):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    # Find indices where "0.0% missing data" appears
    start_indices = [i for i, line in enumerate(lines) if " missing data" in line]
    
    res = []

    # Regex-Muster für das Extrahieren des Wertes nach "Locus"
    pattern = r'Locus\s+\d+\s*:\s*(\S+)'
# Liste, um die gefundenen Locus-Werte zu speichern
    locus_matches = []

# Durch jede Zeile in der Liste gehen
    for line in lines:
    # Suche nach "Locus" in der aktuellen Zeile
        match = re.search(pattern, line)
        if match:
        # Wenn ein Match gefunden wird, füge den Locus-Wert zur Liste hinzu
            locus_matches.append(match.group(1))
    for j, start_index in enumerate(start_indices):
        start_index += 1  # Move to the next line after "0.0% missing data"

        # Find end index (first empty line after start index)
        end_index = start_index + next((i for i, line in enumerate(lines[start_index:], start=1) if line.strip() == ""), None)
        # Combine the relevant lines into a single string
        text = "\n".join(lines[start_index:end_index])
        if (start_index - end_index) < -2:
            # Remove the numbers followed by (0.x) as in R gsub() 
            modified_text = re.sub(r"\s+[0-9]+\s*\(0\.[0-9]+\)", "", text)
            print(modified_text)
            # Split the modified text by spaces and filter out empty strings
            values = [val for val in modified_text.split() if val != ""]
            numeric_values = np.array(values, dtype=float)
            # Sum corresponding values for the two K sections
            print("n", numeric_values)
            res.append(split_array(numeric_values, K))

    return res

# Call the function to read data
p = read_table_data_all(file_pathK, K)

def format_and_save(data, K):
    # Flatten the nested list
    flattened_data = [item for sublist in data for row in sublist for item in row]

    # Split the flattened list into chunks of size K
    reshaped_data = [flattened_data[i:i+K] for i in range(0, len(flattened_data), K)]
    return reshaped_data
    # Save the reshaped data to a file

# Beispiel Daten
data = [
    [[0.13, 0.659]],
    [[0.735, 0.998], [0.265, 0.002]],
    [[0.708, 0.824], [0.292, 0.176]]
]

# Anzahl der Elemente pro Liste (K)
K = 2

# Funktion aufrufen und Datei speichern
res = format_and_save(p, K)

#%%

np.savetxt("C:\\Users\\carol\\Downloads\\all_p", res, delimiter=" ", fmt="%f")


