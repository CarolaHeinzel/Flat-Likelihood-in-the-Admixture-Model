import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

num_arrays = len(data[1:5])
num_values = len(data[0])

x = np.arange(num_arrays)

fig, ax = plt.subplots()

colors = plt.cm.viridis(np.linspace(0, 1, num_values))

value_labels = ['AFR', 'AMR', 'EAS', 'EUR', 'MEA', 'OCE', 'SAS']
colors = ["blue", "red", "yellow", "black", "green", "orange", "brown"]

for i in range(num_arrays):
    bottom = 0
    for j in range(num_values):
        ax.bar(i, data[i][j], bottom=bottom, color=colors[j], label=value_labels[j] if i == 0 else "")
        bottom += data[i][j]
        
# Achsenbeschriftungen und Titel
ax.set_xlabel('Individuals')
ax.set_ylabel('Estimated Ancestry')
ax.legend()

plt.show()
