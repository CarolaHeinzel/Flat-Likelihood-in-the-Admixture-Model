# Use
# plink2 --vcf ET_trainingdata.vcf --vcf-half-call m --sort-vars --make-pgen --out ET_trainingdata
# plink2 --pfile ET_trainingdata --min-alleles 2 --max-alleles 2 --make-bed --out ET_trainingdata
# in order to obtain a bed/bim/fam file

import pandas as pd
import vcfpy 

# transform_to_biallelic = True
include = ["EUROPE", "MIDDLE"]
filename = 'ET_trainingdata_EUROPE_MIDDLEEAST.vcf'

def inc(pop):
    res = False
    for i in include:
        if i in pop:
            res = True 
    return res

df = pd.read_excel("/home/ch1158/Downloads/test_100.xlsx", header = 5, usecols="A:AAN")
samples = df.iloc[5:, 3].to_list()
pops = df.iloc[5:, 1].to_list()
#%%
samples_filtered = []
for i in range(102):
    #if inc(pops[i]):
    samples_filtered.append(samples[i])

rs = df.columns[5:].to_list()
chrom = df.iloc[0, 5:].to_list()
pos = df.iloc[2, 5:].to_list()
ref = df.iloc[3, 5:].to_list()
alt = df.iloc[4, 5:].to_list()

data = df.iloc[5:, 5:]

samples_info = df.iloc[5:, 0:4]
samples_info.to_csv('sample_infocmation.csv', header=False, index = False)

# s ist string, zB TC, ref ist Referenzallel, zB T, alt ist anderes Allel, zB C
def genotype(s, ref, alt):
    # print(s + " " + ref + " " + alt)
    refs = f"{ref}{alt.replace(',','')}"
    try:
        first = refs.index(s[0])
    except:
        first = "."
    try:
        second = refs.index(s[1])
    except:
        second = "."
    return f"{first}/{second}"


with open(filename, "w") as f:
    f.writelines('##fileformat=VCFv4.2\n')
    f.writelines('##FILTER=<ID=PASS,Description="All filters passed">\n')
    f.writelines('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.writelines("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + samples_filtered) + "\n")
    for i in range(len(rs)):
        print(i)    
        data = df[rs[i]][5:].to_list()
        genotypes = []
        for j in range(len(data)):
            if inc(pops[j]):
                # print(samples[j] + " " + rs[i])
                genotypes.append(genotype(data[j], ref[i], alt[i]))
        f.writelines("\t".join([str(chrom[i]), str(pos[i]), rs[i], ref[i], alt[i], "100", "PASS", "DP=100", "GT"] + genotypes) + "\n")
        

reader = vcfpy.Reader.from_path(filename)