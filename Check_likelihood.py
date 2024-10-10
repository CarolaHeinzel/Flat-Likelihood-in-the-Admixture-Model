import pandas as pd
import numpy as np

file1 = "C:\\Users\\carol\\Downloads\\q_minimal.csv"  
file2 = "C:\\Users\\carol\\Downloads\\p_minimal.csv"  
data1_min = pd.read_csv(file1).iloc[:, 1:4] 
data2_min = pd.read_csv(file2).iloc[:, 1:4]

#%%

results = []

for i, row1 in data1_min.iterrows():  
    row_results = []
    for j, row2 in data2_min.iterrows():  
        scalar_product = np.dot(row1.values, row2.values)
        row_results.append(scalar_product)
    results.append(row_results)

results_array = np.array(results)
print("Skalarprodukte:")
print(results_array)


#%%
file1 = "C:\\Users\\carol\\Downloads\\q_initial.csv"  
file2 = "C:\\Users\\carol\\Downloads\\p_initial.csv"  
data1 = pd.read_csv(file1).iloc[:, 1:4] 
data2 = pd.read_csv(file2).iloc[:, 1:4]


results_initial = []

for i, row1 in data1.iterrows():  
    row_results = []
    for j, row2 in data2.iterrows():  
        scalar_product = np.dot(row1.values, row2.values)
        row_results.append(scalar_product)
    results_initial.append(row_results)

results_array = np.array(results_initial)
print("Skalarprodukte:")
print(results_array)
