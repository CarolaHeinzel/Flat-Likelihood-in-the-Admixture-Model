import pandas as pd
# This code convets a .vcf file to a file that has the correct input for STRUCTURE
#%%
vcf_file = 'combined_output_IBS_TSI_CEU.vcf'
        
# Read the Data
vcf_data = []

with open(vcf_file, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue  
        columns = line.strip().split("\t")  
        vcf_data.append(columns)

res_all = []
for i in range(len(vcf_data)):
    vcf_row = vcf_data[i]
    header = vcf_row[:9]
    genotypes = vcf_row[9:]

    genotype_sums = [sum(int(allele) for allele in gt.split('|')) for gt in genotypes]
    res_all.append(genotype_sums)
    

transposed = list(map(list, zip(*res_all))) 

def map_genotype(value):
    if value == 0:
        return [0, 0] 
    elif value == 1:
        return [1, 0]
    else:
        return [1, 1]  

mapped_genotypes = [map_genotype(value) for value in transposed[0]]

df = pd.DataFrame(mapped_genotypes, columns=[f"Marker_{i+1}" for i in range(len(mapped_genotypes[0]))])


for i in range(1, len(transposed)):
    mapped_genotypes = [map_genotype(value) for value in transposed[i]]
    temp_df = pd.DataFrame(mapped_genotypes, columns=[f"Marker_{i+1}_A", f"Marker_{i+1}_B"])
    df = pd.concat([df, temp_df], axis=1)



df_transposed = df.T

df_transposed.insert(0, "Sample_ID", [0]*len(df_transposed))

df_transposed.to_csv("genotype_data_IBS_TSI_CEU.txt", sep="\t", index=False, header=False)

