import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#%%

file_path_p = 'C:\\Users\\carol\\p_K3_P1_0.txt'
file_path_q = 'C:\\Users\\carol\\q_K3_P5_0.txt'

#file_path_q = 'C:\\Users\\carol\\Downloads\\q_STRUCTURE_CEU_IBS_TSI'
data_q = pd.read_csv(file_path_q, sep='\t', header=None)
data_p = pd.read_csv(file_path_p, sep='\t', header=None)

file_path_q1 = 'C:\\Users\\carol\\q_K2_P1_2.txt'
data_q1 = pd.read_csv(file_path_q1, sep='\t', header=None)
data_q  = data_q[[0,1, 2]]

#%%
num_arrays = 313
num_values = 3

q_alle, p_alle = correct_format(data_q, data_p)

x = np.arange(num_arrays)

q_alle2, p_alle = correct_format(data_q1, data_p)

res_q_alle2 = sort_list(q_alle2, sorted_indices)

#%%
res_q_alle1 = sort_list(q_alle,  sorted_indices)
#%%
q_alle_2 = np.transpose(res_q_alle1)
q_alle_1 =  np.transpose(res_q_alle2)



value_labels = ['CEU', 'IBS/TSI']
colors = ["blue", "red"]

fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(8, 10))

# Plot 1
num_values = 2
for i in range(num_arrays):
    bottom = 0
    for j in range(num_values):
        ax1.bar(i, q_alle_1[i][j], bottom=bottom, color=colors[j], label=value_labels[j] if i == 0 else "")
        bottom += q_alle_1[i][j]
vertical_lines = [99.5, 209.5] 
for vline in vertical_lines:
    ax1.axvline(x=vline, color='black', linestyle='-', linewidth=4)
ax1.tick_params(axis='y', labelsize=20)   
ax1.set_xlabel('CEU                    IBS                  TSI\nIndividuals', fontsize=20)

ax1.set_ylabel('Estimated Ancestry', fontsize=20)
ax1.legend()
ax1.legend(fontsize = 20)

ax1.set_xticklabels([])
# Plot 2
value_labels = ['CEU', 'TSI', 'IBS']
colors = ["blue", "orange", "red"]

num_values = 3
for i in range(num_arrays):
    bottom = 0
    for j in range(num_values):
        ax2.bar(i, q_alle_2[i][j], bottom=bottom, color=colors[j], label=value_labels[j] if i == 0 else "")
        bottom += q_alle_2[i][j]
for vline in vertical_lines:
    ax2.axvline(x=vline, color='black', linestyle='-', linewidth=4)
ax2.set_xlabel('CEU                    IBS                  TSI\nIndividuals', fontsize=20)
ax2.set_ylabel('Estimated Ancestry', fontsize=20)
ax2.legend(fontsize = 20)
ax2.set_xticklabels([])
ax2.tick_params(axis='y', labelsize=20)    

plt.tight_layout()

plt.savefig('ancestry_estimation_mischung_K2_K3.pdf')

plt.show()
