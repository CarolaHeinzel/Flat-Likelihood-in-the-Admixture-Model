import numpy as np
import pandas as pd

N = 50  #

first_p = np.random.uniform(0.2, 0.7, N)
second_p = np.random.uniform(0.3, 0.8, N)
first_p =  np.append(first_p, 0)
second_p = np.append(first_p, 1)
first_c = 1 - first_p
second_c = 1-second_p
M = 1000

p1 = np.random.uniform(0.25, 0.75, M)
p2 = p1
for i in range(M):
    e =  np.random.uniform(0, 0.15, 1)
    if(p1[i] + e < 0.8):
        p2[i] = p1[i] + e
    elif(p1[i] - e > 0.2):
        p2[i] = p1[i] - e
    else:
        p2[i] = 0.5
        
p1 = np.append(p1, 0.001)
p1 = np.append(p1, 0.999)
p2 = np.append(p2, 0.999)
p2 = np.append(p2, 0.001)


X = np.zeros((2*(N+1),M+2))

for i in range(N+1):
    for m in range(M+2):
        theta = first_p[i] * p1[m] + first_c[i] * p2[m]
        index = 2*i
        index1 = 2*i + 1
        x = np.random.binomial(2, theta)
        X[i,m] = x



for i in range(N+1):
    for m in range(M+2):
        theta = second_p[i] * p1[m] + second_c[i] * p2[m]
        index = 2*i
        index1 = 2*i + 1
        x = np.random.binomial(2, theta)
        X[(N+1)+i,m] = x

def convert_genotypes(array):
    output = np.empty(array.shape, dtype=object)
    output[array == 0] = 'GG'
    output[array == 1] = 'AG'
    output[array == 2] = 'AA'
    return output

converted_data = convert_genotypes(X)



ref_row = np.full((1, converted_data.shape[1]), 'A')  
alt_row = np.full((1, converted_data.shape[1]), 'G') 

new_array = np.vstack((ref_row, alt_row, converted_data))

df = pd.DataFrame(new_array)

output_filename = 'output_data.xlsx' 
df.to_excel(output_filename, index=False, header=False)
