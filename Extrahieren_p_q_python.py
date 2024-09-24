import os
import re
import numpy as np
import pandas as pd
import io  


K = 2  # Number of populations

#%%
def read_table_data(lines, K):
    # Konvertiere die Byte-Strings in eine lesbare Form
    output_text = ''.join([line.decode('utf-8') for line in lines])
    print(output_text)
    
    # Finde die Startindices, die "0.0% missing data" enthalten
    start_indices = [i for i, line in enumerate(output_text.splitlines()) if "0.0% missing data" in line]
    print(start_indices)
    
    res = np.zeros((0, K))  # Leeres Array für die Ergebnisse
    res_all = np.zeros((0, K))  # Leeres Array für alle Ergebnisse

    for j, start_index in enumerate(start_indices):
        start_index += 1  # Gehe zur nächsten Zeile nach "0.0% missing data"

        # Finde das End-Index (erste leere Zeile nach dem Start-Index)
        end_index = start_index + next((i for i, line in enumerate(output_text.splitlines()[start_index:], start=1) if line.strip() == ""), None)

        # Kombiniere die relevanten Zeilen in einen einzelnen Text
        text = "\n".join(output_text.splitlines()[start_index:end_index])
        print(text)
        
        if (end_index - start_index) > 2:  # Überprüfe, ob genug Zeilen vorhanden sind
            # Entferne die Zahlen, die von (0.x) gefolgt werden
            modified_text = re.sub(r"\s+[0-9]+\s*\(0\.[0-9]+\)", "", text)
            print(modified_text)
            
            # Teile den modifizierten Text in Werte und filtere leere Strings
            values = [val for val in modified_text.split() if val != ""]
            print(values)
            
            # Konvertiere in ein numpy-Array vom Typ float
            numeric_values = np.array(values, dtype=float)
            
            # Summiere die entsprechenden Werte für die zwei K-Sektionen
            if len(numeric_values) >= 2 * K:
                res1 = numeric_values[:K] + numeric_values[K:2*K]
                if np.sum(res1) < 2:  # Bedingung für die Summe
                    # Füge die Ergebnisse zum res-Array hinzu
                    res = np.vstack([res, res1])
                    res_all = np.vstack([res_all, numeric_values[K:2*K]])
                res_all = np.vstack([res_all, numeric_values[:K]])
    if(len(res) == 0):
    	res = None
    return res, res_all



def extract_q(lines, K):
    output_text = ''.join([line.decode('utf-8') for line in lines])

    
    # Finde die Startindices, die "0.0% missing data" enthalten
    start_index = [i for i, line in enumerate(output_text.splitlines()) if "Inferred ancestry of individuals:"  in line]
    end_index = [i for i, line in enumerate(output_text.splitlines()) if "Estimated Allele Frequencies in each cluster"  in line]
    output_text = ''.join([line.decode('utf-8') for line in lines])

    # Finde den Start- und Endindex für die relevanten Abschnitte
    start_index = next((i for i, line in enumerate(output_text.splitlines()) if "Inferred ancestry of individuals:" in line), None)
    end_index = next((i for i, line in enumerate(output_text.splitlines()) if "Estimated Allele Frequencies in each cluster" in line), None)

    # Extrahiere die Zeilen zwischen diesen zwei Markierungen
    extracted_lines = output_text.splitlines()[start_index + 2:end_index - 1]


    # Lese die extrahierten Zeilen als DataFrame
    res = pd.read_csv(io.StringIO("\n".join(extracted_lines)), delim_whitespace=True, header=None)


    # Extrahiere die letzten K Spalten für das endgültige Ergebnis
    res_final = res.iloc[:, -K:]  # Negative Indizes, um die letzten K Spalten zu wählen

    return res_final



