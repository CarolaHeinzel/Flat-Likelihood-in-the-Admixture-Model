#def calc_p_q(data_p, data_q, name_q, name_p):
    
import pandas as pd



data_q =  "C:\\Users\\carol\\Downloads\\random_values.txt" # Filepath
data = {
    0: [0.9, 0.8, 0.2],
    1: [0.82, 0.65, 0.3]
}

data_p = pd.DataFrame(data)


data_q = pd.read_csv(data_q, sep = '\t', header  = None)


def determine_p(data_p):
    result_list = data_p[0].tolist()
    result_list1 = data_p[1].tolist()
    M = len(result_list)
    s1 = 10**100
    sM = -10**100
    
    for m in range(M):
        temp_p = (1-result_list1[m])/(1-result_list[m])
        temp_p1 = result_list1[m]/result_list[m]
        temp_min = min(temp_p, temp_p1)
        temp_max = max(temp_p, temp_p1)
        if(s1 > temp_min):
            s1 = temp_min
        if(sM < temp_max):
            sM = temp_max
    return s1, sM

def calc_p_q(data_p, data_q):

    max_value = data_q[0].max()
    min_value = data_q[0].min()
    s1, sM = determine_p(data_p)

    a_max = (s1-1)/(max_value/(1-max_value) + s1)
    b_max = 1 + a_max * (max_value/(1-max_value))

    a_min = 1/(1-min_value)/(sM + min_value/(1-min_value))
    b_min =  min_value/(1-min_value)*(a_min - 1)


    q_max = data_q[0]*(1-a_max) + (1-data_q[0])*b_max
    q_min = data_q[0]*(1-a_min) + (1-data_q[0])*b_min
    p_max_1 = (data_p[0] - b_max*data_p[0] - a_max*data_p[1])/(1 - a_max - b_max)
    p_min_1 = (data_p[0] - b_min*data_p[0] - a_min*data_p[1])/(1 - a_min - b_min)
    p_max_2 = (data_p[1] - a_max*data_p[1] - b_max*data_p[0])/(1 - a_max - b_max)
    p_min_2 = (data_p[1] - a_min*data_p[1] - b_min*data_p[0])/(1 - a_min - b_min)
    result_qmax = pd.DataFrame([q_max.tolist(), (1-q_max).tolist()]).transpose()
    result_qmin = pd.DataFrame([q_min.tolist(), (1-q_min).tolist()]).transpose()
    result_pmax = pd.DataFrame([p_max_1.tolist(), p_max_2.tolist()]).transpose()
    result_pmin = pd.DataFrame([p_min_1.tolist(), p_min_2.tolist()]).transpose()

    
    return result_qmax, result_qmin, result_pmax, result_pmin

res = calc_p_q(data_p, data_q)


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Einlesen der Dateien
matrix_data_r1 = data_q


# Daten für m1 und m2
m1 = matrix_data_r1.iloc[0:50, 0]
m2 = matrix_data_r1.iloc[50:100, 0]

# Ordnung der Indizes
order_indices = np.argsort(m1)
order_indices1 = 50 + np.argsort(m2)

# Berechnung von q1 und Tabelle tbl1
q1 = matrix_data_r1.iloc[:, 0].values
res_max = np.zeros(len(q1))
res_min = np.zeros(len(q1))

tbl1 = np.column_stack((q1, 1 - q1))



t = np.sort(q1[0:50])
t1 = np.sort(q1[50:101])
t2 = np.concatenate((t, t1))

tbl1[:, 0] = t2
tbl1[:, 1] = 1 - tbl1[:, 0]



# Ordnung der Daten
list1_sorted = np.concatenate([matrix_data_r1.iloc[order_indices1, 0].values, matrix_data_r1.iloc[order_indices, 0].values])
res_max_sorted = np.concatenate([res[0].iloc[order_indices1,0].values, res[0].iloc[order_indices, 0].values])
res_min_sorted = np.concatenate([res[1].iloc[order_indices1, 0].values, res[1].iloc[order_indices, 0].values])
# Plot erstellen
x_values = np.arange(1, len(matrix_data_r1) + 1)
plt.figure(figsize=(10, 5))

plt.plot(x_values, res_max_sorted, label='Res Max', color='darkred')
plt.plot(x_values, res_min_sorted, label='Res Min', color='red')


plt.plot(x_values, list1_sorted, label='List 1', color='darkblue')

# Polygon zwischen res_max und res_min
plt.fill_between(x_values, res_max_sorted, res_min_sorted, color='red', alpha=0.3)

# Linien hinzufügen
plt.axvline(x=50.5, color='black', linewidth=7)

# Achsen und Labels anpassen
plt.xlim(1, 100)
plt.ylim(0, 1.01)
plt.xlabel('Individuals', fontsize = 14, labelpad=26)
plt.ylabel('Estimated Ancestry', fontsize = 14)
plt.xticks([])
plt.tick_params(axis='y', labelsize=14)  # 'labelsize' passt die Schriftgröße der Y-Achse an
# Mtext-Äquivalent in Python
plt.text(25, -0.05, "Population 1", ha='center', fontsize = 14)
plt.text(75, -0.05, "Population 2", ha='center', fontsize = 14)



# Plot anzeigen
plt.show()
#%%
#def calc_p_q(data_p, data_q, name_q, name_p):
    
import pandas as pd



data_q =  "C:\\Users\\carol\\Downloads\\random_values.txt" # Filepath
data = {
    0: [0.9, 0.8, 0.2, 0.0001, 0.999],
    1: [0.82, 0.65, 0.3, 0.999, 0.001]
}

data_p = pd.DataFrame(data)


data_q = pd.read_csv(data_q, sep = '\t', header  = None)


res = calc_p_q(data_p, data_q)

# Einlesen der Dateien
matrix_data_r1 = data_q


# Daten für m1 und m2
m1 = matrix_data_r1.iloc[0:50, 0]
m2 = matrix_data_r1.iloc[50:100, 0]

# Ordnung der Indizes
order_indices = np.argsort(m1)
order_indices1 = 50 + np.argsort(m2)

# Berechnung von q1 und Tabelle tbl1
q1 = matrix_data_r1.iloc[:, 0].values
res_max = np.zeros(len(q1))
res_min = np.zeros(len(q1))

tbl1 = np.column_stack((q1, 1 - q1))



t = np.sort(q1[0:50])
t1 = np.sort(q1[50:101])
t2 = np.concatenate((t, t1))

tbl1[:, 0] = t2
tbl1[:, 1] = 1 - tbl1[:, 0]



# Ordnung der Daten
list1_sorted = np.concatenate([matrix_data_r1.iloc[order_indices1, 0].values, matrix_data_r1.iloc[order_indices, 0].values])
res_max_sorted = np.concatenate([res[0].iloc[order_indices1,0].values, res[0].iloc[order_indices, 0].values])
res_min_sorted = np.concatenate([res[1].iloc[order_indices1, 0].values, res[1].iloc[order_indices, 0].values])
# Plot erstellen
x_values = np.arange(1, len(matrix_data_r1) + 1)
plt.figure(figsize=(10, 5))

plt.plot(x_values, res_max_sorted, label='Res Max', color='darkred')
plt.plot(x_values, res_min_sorted, label='Res Min', color='red')


plt.plot(x_values, list1_sorted, label='List 1', color='darkblue')

# Polygon zwischen res_max und res_min
plt.fill_between(x_values, res_max_sorted, res_min_sorted, color='red', alpha=0.3)

# Linien hinzufügen
plt.axvline(x=50.5, color='black', linewidth=7)
plt.tick_params(axis='y', labelsize=14)  # 'labelsize' passt die Schriftgröße der Y-Achse an
# Achsen und Labels anpassen
plt.xlim(1, 100)
plt.ylim(0, 1)
plt.xlabel('Individuals', fontsize = 14, labelpad=26)
plt.ylabel('Estimated Ancestry', fontsize = 14)
plt.xticks([])

# Mtext-Äquivalent in Python
plt.text(25, -0.05, "Population 1", ha='center', fontsize = 14)
plt.text(75, -0.05, "Population 2", ha='center', fontsize = 14)



# Plot anzeigen
plt.show()
