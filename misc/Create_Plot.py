import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#%%

file_path_q = 'q_K2_P1_0.txt'

data_q = pd.read_csv(file_path_q, sep=' ', header=None)
data_q = data_q[[1,2,0]] # to change the order of the columns

file_path_q1 = 'q_K2_P1_0.txt' # Insert an other data file
data_q1 = pd.read_csv(file_path_q1, sep=' ', header=None)

value_labels = ['CEU', 'IBS', 'TSI']
colors = ["blue", "red", "orange"]

# Positions of the vertical lines
vertical_lines = [99.5, 209.5]
#%%

# num_arrays and num_values are defined based on the shape of data_q1
num_arrays, num_values = data_q1.shape 

def correct_format(data_q):

    data_q = data_q.values.tolist()
    K = len(data_q[0])  # Number of populations
    q_alle = []
    for i in range(K):
        q1_vec = [float(subliste[i]) for subliste in data_q]
        q_alle.append(q1_vec)
    return q_alle

# Convert data_q and data_q1 to the correct format
q_alle = correct_format(data_q)
q_alle_2 = np.transpose(q_alle)
q_alle2  = correct_format(data_q1)
q_alle_1 = np.transpose(q_alle2)

# Create Figure, with two rows and one column for the subplots
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(12, 10))

# Plot 1: Bar plots for the first dataset
for i in range(num_arrays):
    bottom = 0
    for j in range(num_values):
        ax1.bar(i, q_alle_1[i][j], bottom=bottom, color=colors[j], label=value_labels[j] if i == 0 else "")
        bottom += q_alle_1[i][j]

# Add vertical lines at specified positions
for vline in vertical_lines:
    ax1.axvline(x=vline, color='black', linestyle='-', linewidth=4)

# Customize the appearance of the first plot
ax1.tick_params(axis='y', labelsize=20)  # Increase font size for y-axis labels
ax1.set_xlabel('Individuals', fontsize=20)  # Set x-axis label
ax1.set_ylabel('Estimated Ancestry', fontsize=20)  # Set y-axis label
ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=20)  # Position legend to the right of the plot
ax1.set_xticklabels([])  # Remove x-axis tick labels

# Plot 2: Bar plots for the second dataset
for i in range(num_arrays):
    bottom = 0
    for j in range(num_values):
        ax2.bar(i, q_alle_2[i][j], bottom=bottom, color=colors[j], label=value_labels[j] if i == 0 else "")
        bottom += q_alle_2[i][j]

# Add vertical lines for the second plot
for vline in vertical_lines:
    ax2.axvline(x=vline, color='black', linestyle='-', linewidth=4)

# Customize the appearance of the second plot
ax2.set_xlabel('Individuals', fontsize=20)  # Set x-axis label
ax2.set_ylabel('Estimated Ancestry', fontsize=20)  # Set y-axis label
ax2.set_xticklabels([])  # Remove x-axis tick labels
ax2.tick_params(axis='y', labelsize=20)  # Increase font size for y-axis labels

# Save the figure as a PDF file
plt.savefig('plot.pdf')

# Display the figure
plt.show()



