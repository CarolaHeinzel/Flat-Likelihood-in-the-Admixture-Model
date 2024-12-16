import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from itertools import permutations
#%%

def optimal_assignment(files):
    """
    Determines optimal assignment, i.e. minimal euclidian norm
    between the columns of files[i], i= 1,2,... compard to the one in files[0].
    
    """
    dataframes = [pd.read_csv(file, sep='\s+', header=None) for file in files]    
    # Baseline
    baseline_df = dataframes[0]
    K = baseline_df.shape[1]  
    reordered_data = []
    

    for df_idx, df in enumerate(dataframes[1:], start=1):
        # Calculate all permuations
        column_permutations = list(permutations(range(K)))
        
        min_distance = float('inf')
        best_permutation = None
        
        
        for perm in column_permutations:

            permuted_columns = df.iloc[:, list(perm)].to_numpy()
            
            # Calculate 2 norm
            total_distance = np.sum(np.linalg.norm(baseline_df.to_numpy() - permuted_columns, axis=0))
            
            if total_distance < min_distance:
                min_distance = total_distance
                best_permutation = perm
        
        # Save best permutation in reordered_df
        reordered_df = df.iloc[:, list(best_permutation)]
        reordered_data.append(reordered_df)
    
    return reordered_data

# Example
files = ["Q_output.txt", "Q_output_change.txt"]

results = optimal_assignment(files)

# Save results
for idx, result in enumerate(results):
    print(f"Ergebnisse für Datei {idx + 2}:")
    print(result)
    result.to_csv(f"data_alignment{idx + 2}.txt", index=False, header = False, sep = " ")

#%%
# Plots
def correct_format(data_q):
    data_q = data_q.values.tolist()
    K = len(data_q[0])
    q_alle = []
    for i in range(K):
        q1_vec = [float(subliste[i]) for subliste in data_q]
        q_alle.append(q1_vec)
    return q_alle

# Function to plot the results
def plot_q_values(data_files, colors, vertical_lines, output_file="plot.pdf"):
 
    data_list = []
    for file in data_files:
        data = pd.read_csv(file, sep="\s+", header=None)  
        q_values = correct_format(data)
        data_list.append(np.transpose(q_values))

    num_plots = len(data_list)
    num_arrays = len(data_list[0])  # N
    num_values = len(data_list[0][0])  # K

    # Set up figure and subplots
    fig, axes = plt.subplots(nrows=num_plots, ncols=1, figsize=(12, 5 * num_plots))

    if num_plots == 1:
        axes = [axes]

    for idx, (ax, q_alle) in enumerate(zip(axes, data_list)):
        for i in range(num_arrays):
            bottom = 0
            for j in range(num_values):
                ax.bar(i, q_alle[i][j], bottom=bottom, color=colors[j])
                bottom += q_alle[i][j]

        for vline in vertical_lines:
            ax.axvline(x=vline, color='black', linestyle='-', linewidth=4)

        ax.tick_params(axis='y', labelsize=20)
        ax.set_xlabel('Individuals', fontsize=20)
        ax.set_ylabel('Estimated Ancestry', fontsize=20)
        ax.set_xticklabels([])  
        
    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()

data_files = ["Q_output.txt", "data_alignment2.txt"]  
colors = ['red', 'blue', 'green', 'orange', "purple"]  
vertical_lines = [1.5, 2.5]  

plot_q_values(data_files, colors, vertical_lines, output_file="output_plot.pdf")

